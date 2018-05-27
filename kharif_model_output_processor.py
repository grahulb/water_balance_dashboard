import csv
import time
import itertools
from collections import OrderedDict
import numpy as np
from constants_dicts_lookups import *
from kharif_model_calculator import *
from qgis.core import QgsVectorLayer, QgsFeature, QgsField, QgsMapLayerRegistry, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter, QgsGeometry, QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer
from PyQt4.QtGui import QColor
import gdal, ogr, osr

class KharifModelOutputProcessor:

	def output_point_results_to_csv(self, output_grid_points, pointwise_output_csv_filepath, crops, end_date_indices):
		parameters = ['PET-AET', 'Soil Moisture', 'Infiltration', 'Runoff', 'GW Recharge']
		csvwrite = open(pointwise_output_csv_filepath,'w+b')
		writer = csv.writer(csvwrite)
		writer.writerow(['X', 'Y'] + [crops[i]+'-'+parameter+'-'+str(duration) for i in range(len(crops)) for parameter in parameters for duration in end_date_indices] + ['Vegetation-'+parameter+'-'+str(duration) for parameter in parameters for duration in end_date_indices])
		for point in output_grid_points:
			#~ if not point.zone_polygon:	continue
			if point.lulc_type in ['agriculture', 'fallow land']:
				writer.writerow([point.qgsPoint.x(), point.qgsPoint.y()] + 
								list(itertools.chain(*[
									list(itertools.chain(*[
										[point.budget.PET_minus_AET_till_date[end_date_index][i]  for end_date_index in end_date_indices],
										[point.budget.sm_on_date[end_date_index][i]   for end_date_index in end_date_indices],
										[point.budget.infil_till_date[end_date_index][i]  for end_date_index in end_date_indices],
										[point.budget.runoff_till_date[end_date_index][i] for end_date_index in end_date_indices],
										[point.budget.GW_rech_till_date[end_date_index][i]   for end_date_index in end_date_indices]
									]))
										for i in range(len(crops))
								])) +
								['']*5
							)
			else:
				writer.writerow([point.qgsPoint.x(), point.qgsPoint.y()] + 
								['']*(5*len(crops)*len(point.budget.sm_on_date)) +
								list(itertools.chain(*[
									[point.budget.PET_minus_AET_till_date[end_date_index][i]  for end_date_index in end_date_indices],
									[point.budget.sm_on_date[end_date_index][i] for end_date_index in end_date_indices],
									[point.budget.infil_till_date[end_date_index][i]    for end_date_index in end_date_indices],
									[point.budget.runoff_till_date[end_date_index][i]   for end_date_index in end_date_indices],
									[point.budget.GW_rech_till_date[end_date_index][i] for end_date_index in end_date_indices]
								]))
							)
		csvwrite.close()
	
	def render_and_save_pointwise_output_layer(self, points_values_dict, end_date_index, graduated_rendering_interval_points, shapefile_path=''):
		deficit_layer = QgsVectorLayer('Point?crs=epsg:32643', 'Total deficit by day ' + str(end_date_index), 'memory')
		deficit_layer.dataProvider().addAttributes([QgsField('Deficit', QVariant.Double)])
		deficit_layer.updateFields()
		fields = deficit_layer.fields()
		deficit_layer.startEditing()
		for point, deficits in points_values_dict.items():
			f = QgsFeature()
			f.setFields(fields)
			f.setGeometry(QgsGeometry.fromPoint(point.qgsPoint))

			f.setAttributes([float(deficits[end_date_index][0])])
		# for feature in deficit_layer.getFeatures():
		# 	feature['Deficit'] = points_values_dict[feature][end_date_index]
			deficit_layer.dataProvider().addFeatures([f])
			deficit_layer.updateFeature(f)
		deficit_layer.commitChanges()

		ET_D_max = max([v[end_date_index][0]   for v in points_values_dict.values()])
		graduated_symbol_renderer_range_list = []
		opacity = 1
		intervals_count = len(graduated_rendering_interval_points)
		for i in range(intervals_count):
			interval_min = 0 if graduated_rendering_interval_points[i] == 0 else (graduated_rendering_interval_points[i]*ET_D_max/100.0 + 0.01)
			interval_max = (graduated_rendering_interval_points[i+1] * ET_D_max / 100.0) if i+1 < intervals_count else ET_D_max
			label = "{0:.2f} - {1:.2f}".format(interval_min, interval_max)
			colour = QColor(int(255*(1-(i+1.0)/(intervals_count+1.0))), 0, 0)	# +1 done to tackle boundary cases
			symbol = QgsSymbolV2.defaultSymbol(deficit_layer.geometryType())
			symbol.setColor(colour)
			symbol.setAlpha(opacity)
			interval_range = QgsRendererRangeV2(interval_min, interval_max, symbol, label)
			graduated_symbol_renderer_range_list.append(interval_range)
		renderer = QgsGraduatedSymbolRendererV2('', graduated_symbol_renderer_range_list)
		renderer.setMode(QgsGraduatedSymbolRendererV2.EqualInterval)
		renderer.setClassAttribute('Deficit')
		deficit_layer.setRendererV2(renderer)
		QgsMapLayerRegistry.instance().addMapLayer(deficit_layer)

		if shapefile_path != '':
			QgsVectorFileWriter.writeAsVectorFormat(deficit_layer, shapefile_path, "utf-8", None, "ESRI Shapefile")
		
		return deficit_layer

	def save_pointwise_deficit_as_raster(self, iface, points_grid, end_date_index, save_file_path, layer_name):
		cols = len(points_grid[0]); rows = len(points_grid)
		originX = points_grid[0][0].qgsPoint.x();   originY = points_grid[0][0].qgsPoint.y()   #top-left

		driver = gdal.GetDriverByName('GTiff')
		outRaster = driver.Create(save_file_path, cols, rows, 1, gdal.GDT_Byte)
		outRaster.SetGeoTransform((originX, STEP, 0, originY, 0, STEP))
		outband = outRaster.GetRasterBand(1)

		deficit_values_array = np.array([[p.budget.PET_minus_AET_till_date[end_date_index] if not p.is_no_evaluation_point else -9999 for p in p_row] for p_row in points_grid])
		max_deficit = np.max(deficit_values_array)

		outband.WriteArray(deficit_values_array)
		outRasterSRS = osr.SpatialReference()
		outRasterSRS.ImportFromEPSG(32643)
		outRaster.SetProjection(outRasterSRS.ExportToWkt())
		outband.FlushCache()

		layer = iface.addRasterLayer(save_file_path, layer_name+'_max_deficit_'+str(int(max_deficit)))
		layer.dataProvider().setNoDataValue(1, -9999)
		fcn = QgsColorRampShader()
		fcn.setColorRampType(QgsColorRampShader.INTERPOLATED)

		lst = [QgsColorRampShader.ColorRampItem(0, QColor(240, 240, 240)), QgsColorRampShader.ColorRampItem(max_deficit, QColor(255, 0, 0))]
		fcn.setColorRampItemList(lst)
		shader = QgsRasterShader()
		shader.setRasterShaderFunction(fcn)
		renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)
		layer.setRenderer(renderer)
		layer.triggerRepaint()
