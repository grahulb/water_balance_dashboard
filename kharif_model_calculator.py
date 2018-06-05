from __future__ import division
from qgis.gui import QgsMapToolEmitPoint
import qgis.gui
from qgis.core import QgsSpatialIndex, QgsPoint, QgsRectangle, QgsRaster, QgsVectorLayer, QgsFeature
from qgis.analysis import  QgsGeometryAnalyzer
from PyQt4.QtGui import *
from PyQt4.QtCore import QVariant
from PyQt4.QtCore import QFileInfo
import csv
import os
import time
import processing
import sys
import shutil
import numpy as np
from collections import OrderedDict
from configuration import *
from constants_dicts_lookups import *
from copy import deepcopy
import itertools

BOUNDARY_LABEL = 'Zones'
SOIL_LABEL = 'Soil'
LULC_LABEL = 'Land-Use-Land-Cover'
SLOPE_LABEL = 'Slope'
CADASTRAL_LABEL = 'Cadastral'

class Budget:
	
	def __init__(self):
		self.sm, self.runoff, self.infil, self.AET, self.GW_rech, self.sec_run_off = [],[],[],[],[], []
		
	def summarize(self, start_date_index, end_date_indices, PET):
		self.runoff = np.array(self.runoff)
		self.infil = np.array(self.infil)
		self.AET = np.array(self.AET)
		self.GW_rech = np.array(self.GW_rech)
		self.sm = np.array(self.sm)
		self.sec_run_off = np.array(self.sec_run_off)
		self.sm_on_date, self.runoff_till_date, self.infil_till_date, self.AET_till_date, self.GW_rech_till_date = {}, {}, {}, {}, {}
		self.PET_minus_AET_till_date = {}
		for end_date_index in end_date_indices:
			self.sm_on_date[end_date_index] = np.array([self.sm[end_date_index][0]])
			self.runoff_till_date[end_date_index] = [np.sum(self.runoff[start_date_index:end_date_index + 1,0]) + np.sum(self.sec_run_off[start_date_index:end_date_index+1,0])]
			self.infil_till_date[end_date_index] = [np.sum(self.infil[start_date_index:end_date_index+1,0]) - np.sum(self.sec_run_off[start_date_index:end_date_index+1,0])]
			self.AET_till_date[end_date_index] = [np.sum(self.AET[start_date_index:end_date_index+1,0])	]
			self.GW_rech_till_date[end_date_index] = [np.sum(self.GW_rech[start_date_index:end_date_index+1,0])	]
			self.PET_minus_AET_till_date[end_date_index] = np.array([np.sum(PET[start_date_index:end_date_index+1]) - self.AET_till_date[end_date_index][0]])
		# print self.GW_rech_monsoon_end
		# self.gwl=[float(k[0]) for k in self.GW_rech.tolist()]
		# self.ro =[float(k[0]) for k in self.runoff.tolist()]
		# self.ssm = [float(k[0]) for k in self.sm.tolist()]
		# self.infi = [float(k[0]) for k in self.infil.tolist()]
		# self.aet = [float(k[0]) for k in self.AET.tolist()]
		# a = [float(k[0]) for k in self.AET.tolist()]
		# for j in range (0,150):
		# 	print (ssm[j],ro[j],infi[j],a[j],gwl[j])
		# print ssm
		# print monsoon_end_date_index, end_date_index
		# print self.sm_crop_end, self.sm_monsoon_end
		# print "infil"
		# print self.infi
		# print "AET"
		# print self.aet
		# print "gw"
		# print self.gwl
		# print "sm"
		# print self.ssm
		# print "ro"
		# print self.ro


class Crop:
	def __init__(self, name):
		self.name = name
		self.PET = {};	self.end_date_index =dict(); self.PET_till_date = {}
	@property
	def root_depth(self):	return dict_crop[self.name][2] if self.name in dict_crop.keys() else dict_crop_current_fallow[self.name][2] if self.name in dict_crop_current_fallow.keys() else dict_LULC_pseudo_crop[self.name][2]
	@property
	def KC(self):	return dict_crop[self.name][0] if self.name in dict_crop.keys() else dict_crop_current_fallow[self.name][0]	if self.name in dict_crop_current_fallow.keys()	else dict_LULC_pseudo_crop[self.name][0]
	@property
	def depletion_factor(self):	return dict_crop[self.name][1] if self.name in dict_crop.keys() else dict_crop_current_fallow[self.name][1] if self.name in dict_crop_current_fallow.keys() else dict_LULC_pseudo_crop[self.name][1]


class Point:
	
	def __init__(self, qgsPoint, data=None):
		self.qgsPoint = qgsPoint
		self.data = data
		self.is_no_evaluation_point = False
		self.container_polygons = {}
		self.slope = None
		self._crop = None
		self.budget = Budget()
	
	@property
	def zone_polygon(self):	return self.container_polygons[BOUNDARY_LABEL]
	@property
	def soil_polygon(self):	return self.container_polygons[SOIL_LABEL]
	@property
	def lulc_polygon(self):	return self.container_polygons[LULC_LABEL]
	@property
	def cadastral_polygon(self):	return self.container_polygons[CADASTRAL_LABEL]
	@property
	def texture(self):
		if READ_FROM_POINTS_DATA_FILE:	return self.data['Soil Texture']
		else:                           return self.soil_polygon[TEX].lower()
	@property
	def district(self):
		if READ_FROM_POINTS_DATA_FILE:	return self.data['District']
		else:                           return self.container_polygons[BOUNDARY_LABEL]['District']
	@property
	def taluka(self):
		if READ_FROM_POINTS_DATA_FILE:  return self.data['Taluka']
		else:                           return self.container_polygons[BOUNDARY_LABEL]['Taluka']
	@property
	def rain_circle(self):
		if READ_FROM_POINTS_DATA_FILE:  return self.data['Rainfall Circle']
		else:                   		return self.container_polygons[BOUNDARY_LABEL][Circle]
	@property
	def crop(self):
		if self._crop is not None:  return self._crop
		if self.lulc_type in ['agriculture', 'fallow land']:
			self._crop = Crop(district_taluka_major_crops_dict[(self.district, self.taluka)])
		else:
			self._crop = Crop(self.lulc_type)
		return self._crop
	@property
	def depth_value(self):
		if READ_FROM_POINTS_DATA_FILE:  return float(self.data['Soil Depth']) if self.data['Soil Depth'] != '' else ''
		else:                           return dict_SoilDep[self.soil_polygon[Depth].lower()]
	@property
	def Ksat(self): return dict_SoilProperties[self.texture][7]
	@property
	def Sat(self):	return dict_SoilProperties[self.texture][6]
	@property
	def WP(self):   return dict_SoilProperties[self.texture][4]
	@property
	def FC(self):   return dict_SoilProperties[self.texture][5]
	@property
	def lulc_type(self):
		if READ_FROM_POINTS_DATA_FILE:  return self.data['LULC Type']
		else:                           return dict_lulc[self.lulc_polygon[Desc].lower()]
	@property
	def HSG(self):  return dict_SoilProperties[self.texture][0]
	
	@zone_polygon.setter
	def zone_polygon(self, polygon):	self.container_polygons[BOUNDARY_LABEL] = polygon
	@soil_polygon.setter
	def soil_polygon(self, polygon):	self.container_polygons[SOIL_LABEL] = polygon
	@lulc_polygon.setter
	def lulc_polygon(self, polygon):	self.container_polygons[LULC_LABEL] = polygon
	@cadastral_polygon.setter
	def cadastral_polygon(self, polygon):	self.container_polygons[CADASTRAL_LABEL] = polygon

	def set_PET_and_end_date_index_of_crop(self, rain, et0, sowing_threshold):
		def compute_sowing_index(daily_rain):
			rain_sum = 0
			for i in range(0, len(daily_rain)):
				if (rain_sum < sowing_threshold):	rain_sum += daily_rain[i]
				else:           					break
			return i
		pre_sowing_kc = [0] * compute_sowing_index(rain)
		kc = ([] if self.crop.name in dict_LULC_pseudo_crop else pre_sowing_kc) + self.crop.KC
		if (len(kc) < 365):	self.crop.end_date_index = len(kc) - 1
		else:   			self.crop.end_date_index = 364
		kc = kc + [0] * (365 - len(kc))
		self.crop.PET = np.array(et0[0:len(kc)]) * np.array(kc)

	def run_model(self, rain, crops, start_date_index, end_date_indices, lulc):
		self.setup_for_daily_computations(crops,lulc)
		self.SM1_fraction = self.layer2_moisture = self.WP

		for day in range (start_date_index, max(end_date_indices)+1):
			self.primary_runoff(day, rain)
			self.aet(day)
			self.percolation_below_root_zone(day)
			self.secondary_runoff(day)
			self.percolation_to_GW(day)

		self.budget.summarize(start_date_index, end_date_indices, self.crop.PET)
		for end_date_index in end_date_indices: self.budget.sm_on_date[end_date_index] -= self.WP_depth	# requirement expressed by users

	def setup_for_daily_computations(self, crops,lulc):
		"""
		"""
		self.lulc_type = lulc
		Sat_depth = self.Sat * self.depth_value * 1000
		self.WP_depth = self.WP * self.depth_value * 1000
		FC_depth = self.FC * self.depth_value * 1000
		
		root_depths = np.array([self.crop.root_depth])
		self.SM1 = np.where(self.depth_value <= root_depths, self.depth_value - 0.05, root_depths)
		self.SM2 = np.where(self.depth_value <= root_depths, 0.05, self.depth_value - root_depths)

		self.cn_val = dict_RO[lulc][self.HSG]
		
		cn_s = cn_val = self.cn_val
		cn3 = cn_s *np.exp(0.00673*(100-cn_s))
		if (self.slope > 5.0):
			cn_s = (((cn3-cn_val)/3)*(1-2*np.exp(-13.86*self.slope * 0.01))) + cn_val
		cn1_s = cn_s - 20*(100-cn_s)/(100-cn_s+np.exp(2.533-0.0636*(100-cn_s)))
		cn3_s = cn_s *np.exp(0.00673*(100-cn_s))
		
		self.Smax = 25.4 * (1000/cn1_s - 10)
		S3 = 25.4 * (1000/cn3_s - 10)
		self.W2 = (np.log((FC_depth- self.WP_depth)/(1-S3/self.Smax) - (FC_depth - self.WP_depth )) - np.log ((Sat_depth - self.WP_depth)/(1-2.54/self.Smax) - (Sat_depth - self.WP_depth)))/((Sat_depth- self.WP_depth) - (FC_depth - self.WP_depth))
		self.W1 = np.log((FC_depth- self.WP_depth)/(1- S3/self.Smax) - (FC_depth - self.WP_depth)) + self.W2 * (FC_depth -self.WP_depth)
		
		TT_perc = (Sat_depth- FC_depth)/self.Ksat	#SWAT equation 2:3.2.4
		self.daily_perc_factor = 1 - np.exp(-24 / TT_perc)	#SWAT equation 2:3.2.3
		
	def primary_runoff(self, day, rain):
		"""
		Retention parameter 'S_swat' using SWAT equation 2:1.1.6
		Curve Number for the day 'Cn_swat' using SWAT equation 2:1.1.11
		Initial abstractions (surface storage,interception and infiltration prior to runoff)
			'Ia_swat' derived approximately as recommended by SWAT
		Primary Runoff 'Swat_RO' using SWAT equation 2:1.1.1
		"""
		self.budget.sm.append((self.SM1_fraction * self.SM1 + self.layer2_moisture * self.SM2) * 1000)
		self.SW = self.budget.sm[-1] - self.WP_depth
		S_swat = self.Smax*(1 - self.SW/(self.SW + np.exp(self.W1 - self.W2 * self.SW)))
		Cn_swat = 25400/(S_swat+254)
		Ia_swat = 0.2 * S_swat
		self.budget.runoff.append(	np.where(rain[day] > Ia_swat,
											((rain[day]-Ia_swat)**2)/(rain[day] + 0.8*S_swat),
											0
									)	)
		self.budget.infil.append(rain[day] - self.budget.runoff[day])

	def aet(self, day):
		"""
		Water Stress Coefficient 'KS' using FAO Irrigation and Drainage Paper 56, page 167 and
			page 169 equation 84
		Actual Evapotranspiration 'AET' using FAO Irrigation and Drainage Paper 56, page 6 and 
			page 161 equation 81
		"""
		depletion_factors = np.array([self.crop.depletion_factor])
		KS = np.where(self.SM1_fraction < self.WP, 0,
						np.where(self.SM1_fraction > (self.FC *(1- depletion_factors) + depletion_factors * self.WP), 1,
							(self.SM1_fraction - self.WP)/(self.FC - self.WP) /(1- depletion_factors)
							)
						)
		PETs = np.array([self.crop.PET[day]	if day <= self.crop.end_date_index	else 0])
		self.budget.AET.append( KS * PETs )
	
	def percolation_below_root_zone(self, day):
		"""
		Calculate soil moisture (fraction) 'SM1_before' as the one after infiltration and (then) AET occur,
		but before percolation starts below root-zone. Percolation below root-zone starts only if
		'SM1_before' is more than field capacity and the soil below root-zone is not saturated,i.e.
		'layer2_moisture' is less than saturation. When precolation occurs it is derived as
		the minimum of the maximum possible percolation (using SWAT equation 2:3.2.3) and
		the amount available in the root-zone for percolation.
		"""
		self.SM1_before = (self.SM1_fraction*self.SM1 +((self.budget.infil[day]-self.budget.AET[day])/1000))/self.SM1
		#~ print np.logical_or(self.SM1_before < self.FC, self.layer2_moisture < self.Sat)
		#~ print np.minimum((self.Sat - self.layer2_moisture) * self.SM2 * 1000,
										 #~ (self.SM1_before - self.FC) * self.SM1 * 1000 * self.daily_perc_factor)
		self.R_to_second_layer = np.where(self.SM1_before < self.FC, 0,
									np.where(self.layer2_moisture < self.Sat,
									np.minimum((self.Sat - self.layer2_moisture) * self.SM2 * 1000,
										 (self.SM1_before - self.FC) * self.SM1 * 1000 * self.daily_perc_factor),
										0
									 ))
		self.SM2_before = (self.layer2_moisture*self.SM2*1000 + self.R_to_second_layer)/self.SM2/1000

	def secondary_runoff(self, day):
		"""
		
		"""
		self.budget.sec_run_off.append(np.where(
							((self.SM1_before*self.SM1 - self.R_to_second_layer/1000)/self.SM1) > self.Sat,
							(((self.SM1_before*self.SM1 - self.R_to_second_layer/1000)/self.SM1) - self.Sat) * self.SM1 * 1000,
							0
							))
		self.SM1_fraction = np.minimum((self.SM1_before*self.SM1*1000 - self.R_to_second_layer)/self.SM1/1000,self.Sat)
	
	def percolation_to_GW(self, day):
		"""
		
		"""
		self.budget.GW_rech.append(np.maximum((self.SM2_before - self.FC)*self.SM2*self.daily_perc_factor*1000,0))
		self.layer2_moisture = np.minimum(((self.SM2_before*self.SM2*1000- self.budget.GW_rech[day])/self.SM2/1000),self.Sat)
	

class VectorLayer:
	
	def __init__(self, qgsLayer, name=''):
		self.qgsLayer = qgsLayer
		self.name = name
		self.feature_dict = {f.id(): f for f in qgsLayer.getFeatures()}
		self.index = QgsSpatialIndex(qgsLayer.getFeatures())
	
	def get_polygon_containing_point(self, point):
		intersector_ids = self.index.intersects( QgsRectangle( point.qgsPoint, point.qgsPoint ) )
		for intersector_id in intersector_ids:
			polygon = self.feature_dict[intersector_id]
			if (polygon.geometry().contains(point.qgsPoint)):
				return polygon
		return None


class KharifModelCalculator:
	"""
	The actual algorithm for calculating results of the Kharif Model
	"""
	
	def __init__(self, path, et0, points_data=None, zones_layer=None, soil_layer=None, lulc_layer=None, slope_layer=None):

		self.path = path

		if points_data is not None:
			self.points_data = points_data
		else:
			self.zones_layer = VectorLayer(zones_layer, BOUNDARY_LABEL)
			self.soil_layer = VectorLayer(soil_layer, SOIL_LABEL)
			self.lulc_layer = VectorLayer(lulc_layer, LULC_LABEL)
			self.slope_layer = slope_layer
		
		self.et0 = et0
		assert et0 is not None
	
	@property
	def soil_types(self):	sts = dict_SoilProperties.keys();	sts.remove('soil type');	return sts
	@property
	def lulc_types(self):	return dict_RO.keys()
		
	def generate_output_points_grid(self):
		if READ_FROM_POINTS_DATA_FILE:
			grid_row_number = 1;    output_points_grid = [];    points_grid_row = []
			for points_data_row in self.points_data:
				if int(points_data_row['Grid Row Number']) > grid_row_number:
					grid_row_number += 1
					output_points_grid.append(points_grid_row)
					points_grid_row = []
				points_grid_row.append(Point(QgsPoint(float(points_data_row['X']), float(points_data_row['Y'])), points_data_row))
			output_points_grid.append(points_grid_row)
		else:
			xmin_on_aligned_grid = 1000.0 * int(self.zones_layer.qgsLayer.extent().xMinimum() / 1000)
			xmax_on_aligned_grid = 1000.0 * (int(self.zones_layer.qgsLayer.extent().xMaximum() / 1000) + 1)
			ymin_on_aligned_grid = 1000.0 * int(self.zones_layer.qgsLayer.extent().yMinimum() / 1000)
			ymax_on_aligned_grid = 1000.0 * (int(self.zones_layer.qgsLayer.extent().yMaximum() / 1000) + 1)
			print 'boundary min, max : ' , xmin_on_aligned_grid, xmax_on_aligned_grid, ymin_on_aligned_grid, ymax_on_aligned_grid
			def frange(start,end,step):
				i = start
				while i<=end :
					yield i
					i = i+step

			x_List = [x for x in frange(xmin_on_aligned_grid,xmax_on_aligned_grid,STEP)]
			y_List = [x for x in frange(ymin_on_aligned_grid,ymax_on_aligned_grid,STEP)]
			print len(x_List), len (y_List)
			output_points_grid = [[Point(QgsPoint(x,y))	for x in x_List]	for y in y_List]
		return output_points_grid
	
	def filter_out_points_outside_boundary(self):
		if READ_FROM_POINTS_DATA_FILE:
			for points_row in self.points_grid:
				for point in points_row:
					point.is_no_evaluation_point = (point.data['District'] == '')
		else:
			for point in self.output_grid_points:
				polygon = self.zones_layer.get_polygon_containing_point(point)
				if polygon is not None:	point.container_polygons[BOUNDARY_LABEL] = polygon
				else:   				point.is_no_evaluation_point = True


	def set_container_polygon_of_points_for_layers(self, points, polygon_vector_layers):
		for layer in polygon_vector_layers:
			print layer.name
			for p in points:
				p.container_polygons[layer.name] = layer.get_polygon_containing_point(p)
	
	def set_slope_at_points(self, points):
		if READ_FROM_POINTS_DATA_FILE:
			for point in points:
				point.slope = point.data['Slope']
		else:
			for point in points:
				point.slope = self.slope_layer.dataProvider().identify(
					point.qgsPoint, QgsRaster.IdentifyFormatValue).results()[1]

	def filter_out_points_with_incomplete_data(self, points):
		log_file = open(os.path.join(self.path, 'log'), 'a')
		log_file.write('\n' + time.ctime(time.time()) + '\n')
		if READ_FROM_POINTS_DATA_FILE:
			for points_row in points:
				for point in points_row:
					if '' in [point.texture, point.depth_value, point.lulc_type, point.rain_circle, point.slope]:
						point.is_no_evaluation_point = True
					if (point.district, point.taluka, point.rain_circle, YEAR) not in self.rain:
						point.is_no_evaluation_point = True
			return
		for point in points:
			if (None in [
					point.container_polygons[SOIL_LABEL],
					point.container_polygons[LULC_LABEL],
					point.container_polygons[BOUNDARY_LABEL],
					point.slope]
				):
				point.is_no_evaluation_point = True
				if point.container_polygons[SOIL_LABEL] is None:
					log_file.write('Soil polygon could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
				if point.container_polygons[LULC_LABEL] is None:
					log_file.write('LULC polygon could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
				if point.slope is None:
					log_file.write('Slope could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
				if point.container_polygons[BOUNDARY_LABEL] is None:
					log_file.write('Zone polygon could not be obtained for point at: x = '
								   + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
			if point.container_polygons[BOUNDARY_LABEL] is not None and point.container_polygons[BOUNDARY_LABEL][Circle] is None:
				point.is_no_evaluation_point = True
				log_file.write('Rainfall Circle could not be obtained for point at: x = '
				               + str(point.qgsPoint.x()) + ', y = ' + str(point.qgsPoint.y()))
		log_file.close()

	def calculate(self,
					rain,
					sowing_threshold,
					start_date_index=START_DATE_INDEX,
                    end_date_indices=[END_DATE_INDEX]
				):
		
		self.rain = rain

		self.points_grid = self.generate_output_points_grid()
		self.output_grid_points = list(itertools.chain(*self.points_grid))
		print 'Number of grid points to process : ', len([p for p in self.output_grid_points if not p.is_no_evaluation_point])
		self.filter_out_points_outside_boundary()
		print 'Number of grid points to process : ', len([p for p in self.output_grid_points if not p.is_no_evaluation_point])
		if not READ_FROM_POINTS_DATA_FILE:
			self.set_container_polygon_of_points_for_layers(self.output_grid_points, [self.soil_layer, self.lulc_layer, self.zones_layer])
			print 'Setting container polygons'
			self.set_slope_at_points(self.output_grid_points)
			print 'Setting slope'
		if READ_FROM_POINTS_DATA_FILE:
			self.filter_out_points_with_incomplete_data(self.points_grid)
		else:
			self.filter_out_points_with_incomplete_data(self.output_grid_points)
		print 'Number of grid points to process : ', len([p for p in self.output_grid_points if not p.is_no_evaluation_point])
		if READ_FROM_POINTS_DATA_FILE:
			for points_row in self.points_grid:
				for point in points_row:
					if not point.is_no_evaluation_point and point.lulc_type in ['water', 'habitation']:
						point.is_no_evaluation_point = True
		else:
			for point in self.output_grid_points:
				if not point.is_no_evaluation_point and point.lulc_type in ['water','habitation']:   point.is_no_evaluation_point = True
		total_points = len([p for p in self.output_grid_points if not p.is_no_evaluation_point])
		print 'Number of grid points to process : ', total_points
		count = 0
		if id(self.output_grid_points[0]) != id(self.points_grid[0][0]):    print "What!!"; return
		if READ_FROM_POINTS_DATA_FILE:
			# print "In READ_FROM_POINTS_DATA_FILE"
			# print len(self.points_grid[0])
			i = 0
			for points_row in self.points_grid:
				i += 1
				j = 0
				# summarized = []
				for point in points_row:
					j += 1
					# print point.is_no_evaluation_point
					# if count > 680:
					# 	# print j, ' before condition', len(summarized)
					# 	if j == 17:
					# 		print vars(point)
					if point.is_no_evaluation_point:
						# summarized.append(True)
						# if count > 680: print j, ' within condition', len(summarized)
						continue
					count += 1
					if count % 100 == 0:
						print count, '/', total_points
					try:
						rain_at_point = rain[(point.district, point.taluka, point.rain_circle, YEAR)]
					except:
						# point.is_no_evaluation_point = True
						# print "Helllllllo !!"
						continue
					point.set_PET_and_end_date_index_of_crop(rain_at_point, self.et0[point.district], sowing_threshold)
					point.run_model(rain_at_point, [point.crop], start_date_index, end_date_indices, point.lulc_type)
					# summarized.append((point.is_no_evaluation_point) or ('PET_minus_AET_till_date' in vars(point.budget).keys()))
					# if count > 680: print j, ' after running', len(summarized)
					# if 'PET_minus_AET_till_date' not in vars(point.budget).keys():
					# 	print "why1?"
				# if count > 680: print summarized
				# if not all([(point.is_no_evaluation_point) or ('PET_minus_AET_till_date' in vars(point.budget).keys())   for point in points_row]):
				# 	print [(point.is_no_evaluation_point) or ('PET_minus_AET_till_date' in vars(point.budget).keys())   for point in points_row]
				# 	print "why2?"
				# 	return
		else:
			for point in self.output_grid_points:
				if point.is_no_evaluation_point:    continue
				count += 1
				if count % 100 == 0:
					print count, '/', total_points
				try:
					rain_at_point = rain[(point.district, point.taluka, point.rain_circle, YEAR)]
				except:
					continue
				point.set_PET_and_end_date_index_of_crop(rain_at_point, self.et0[point.district], sowing_threshold)
				point.run_model(rain_at_point, [point.crop], start_date_index, end_date_indices, point.lulc_type)

		print "in kharif_model_calculator"
		for points_row in self.points_grid:
			for p in points_row:
				if (not p.is_no_evaluation_point) and ('PET_minus_AET_till_date' not in vars(p.budget).keys()):
					print p.is_no_evaluation_point, vars(p.budget).keys()
					return

