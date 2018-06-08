import numpy as np
from kharif_model_calculator import *
from qgis.core import QgsVectorLayer, QgsFeature, QgsField, QgsMapLayerRegistry, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter, QgsGeometry, QgsColorRampShader, QgsRasterShader, QgsSingleBandPseudoColorRenderer
from PyQt4.QtGui import QColor
import gdal, ogr, osr

class KharifModelOutputProcessor:

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
