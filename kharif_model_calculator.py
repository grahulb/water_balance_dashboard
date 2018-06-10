from __future__ import division
from qgis.core import QgsSpatialIndex, QgsPoint, QgsRectangle, QgsRaster, QgsVectorLayer, QgsFeature
from PyQt4.QtGui import *
from math import *
from configuration import *
from constants_dicts_lookups import *
import itertools

class Budget:
	
	def __init__(self):
		self.sm, self.runoff, self.infil, self.AET, self.GW_rech, self.sec_run_off = [],[],[],[],[],[]
		
	def summarize(self, start_date_index, end_date_indices, PET):
		self.sm_on_date, self.runoff_till_date, self.infil_till_date, self.AET_till_date, self.GW_rech_till_date = {}, {}, {}, {}, {}
		self.PET_minus_AET_till_date = {}
		for end_date_index in end_date_indices:
			# self.sm_on_date[end_date_index] = self.sm[end_date_index]
			# self.sm_on_date[end_date_index] -= self.WP_depth  # requirement expressed by users
			# accumulated_sec_runoff = sum(self.sec_run_off[start_date_index:end_date_index + 1])
			# self.runoff_till_date[end_date_index] = sum(self.runoff[start_date_index:end_date_index + 1]) + accumulated_sec_runoff
			# self.infil_till_date[end_date_index] = sum(self.infil[start_date_index:end_date_index+1]) - accumulated_sec_runoff
			self.AET_till_date[end_date_index] = sum(self.AET[start_date_index:end_date_index+1])
			# self.GW_rech_till_date[end_date_index] = sum(self.GW_rech[start_date_index:end_date_index+1])
			self.PET_minus_AET_till_date[end_date_index] = sum(PET[start_date_index:end_date_index+1]) - self.AET_till_date[end_date_index]

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
	
	def __init__(self, qgsPoint, crop, data_object=None):
		self.qgsPoint = qgsPoint
		if INPUT_FROM_GRID_POINTS_ROSTER:   self.feature = data_object
		elif READ_FROM_POINTS_DATA_FILE:    self.data = data_object
		self.is_no_evaluation_point = False
		self.slope = None
		self.crop = crop
		self.budget = Budget()
	
	@property
	def texture(self):
		if INPUT_FROM_GRID_POINTS_ROSTER:   return str(self.feature['Soil Textu'])
		elif READ_FROM_POINTS_DATA_FILE:	return self.data['Soil Texture']
	@property
	def district(self):
		if INPUT_FROM_GRID_POINTS_ROSTER:   return str(self.feature['District']).lower()
		elif READ_FROM_POINTS_DATA_FILE:	return self.data['District'].lower()
	@property
	def taluka(self):
		# return 'Bhatkuli'
		if INPUT_FROM_GRID_POINTS_ROSTER:   taluka_code = str(self.feature['Taluka Cod']);   return taluka_code_to_name_dict[taluka_code] if taluka_code != 'NULL' else 'NULL'
		elif READ_FROM_POINTS_DATA_FILE:    return self.data['Taluka'].lower()
	@property
	def rain_circle(self):
		# return 'Bhatkuli'
		if INPUT_FROM_GRID_POINTS_ROSTER:   return str(self.feature['Rainfall C']).lower()
		elif READ_FROM_POINTS_DATA_FILE:    return self.data['Rainfall Circle'].lower()
	@property
	def depth_value(self):
		if INPUT_FROM_GRID_POINTS_ROSTER:   soil_depth_str = str(self.feature['Soil Depth']);   return float(soil_depth_str) if soil_depth_str != 'NULL' else 'NULL'
		elif READ_FROM_POINTS_DATA_FILE:    return float(self.data['Soil Depth']) if self.data['Soil Depth'] != '' else ''
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
		if INPUT_FROM_GRID_POINTS_ROSTER:   return str(self.feature['LULC Type'])
		elif READ_FROM_POINTS_DATA_FILE:    return self.data['LULC Type']
	@property
	def HSG(self):  return dict_SoilProperties[self.texture][0]
	
	def set_PET_and_end_date_index_of_crop(self, rain, et0, sowing_threshold):
		def compute_sowing_index(daily_rain):
			rain_sum = 0
			for i in range(0, len(daily_rain)):
				if (rain_sum < sowing_threshold):	rain_sum += daily_rain[i]
				else:           					break
			return i
		pre_sowing_kc = [0] * compute_sowing_index(rain)
		kc = ([] if self.crop.name in dict_LULC_pseudo_crop else pre_sowing_kc) + self.crop.KC
		if (len(kc) < 365): self.crop.end_date_index = len(kc) - 1
		else:   			self.crop.end_date_index = 364;     kc = kc[:365]
		kc = kc + [0] * (365 - len(kc))
		self.crop.PET = [et0[i]*kc[i]   for i in range(len(kc))]

	def run_model_numpy(self, rain, start_date_index, end_date_indices):
		self.setup_for_daily_computations()
		self.SM1_fraction = self.layer2_moisture = self.WP

		for day in range(start_date_index, max(end_date_indices) + 1):
			self.primary_runoff(day, rain)
			self.aet(day)
			self.percolation_below_root_zone(day)
			self.secondary_runoff()
			self.percolation_to_GW(day)

		self.budget.summarize(start_date_index, end_date_indices, self.crop.PET)
		for end_date_index in end_date_indices:
			self.budget.sm_on_date[end_date_index] -= self.WP_depth  # requirement expressed by users

	def run_model(self, rain, start_date_index, end_date_indices):
		self.setup_for_daily_computations()
		self.SM1_fraction = self.layer2_moisture = self.WP

		for day in range(start_date_index, max(end_date_indices) + 1):
			self.primary_runoff(day, rain)
			self.aet(day)
			self.percolation_below_root_zone(day)
			self.secondary_runoff()
			self.percolation_to_GW(day)
		self.budget.summarize(start_date_index, end_date_indices, self.crop.PET)

	def setup_for_daily_computations_numpy(self):
		"""
		"""
		Sat_depth = self.Sat * self.depth_value * 1000
		self.WP_depth = self.WP * self.depth_value * 1000
		FC_depth = self.FC * self.depth_value * 1000

		root_depths = np.array([self.crop.root_depth])
		self.SM1 = np.where(self.depth_value <= root_depths, self.depth_value - 0.05, root_depths)
		self.SM2 = np.where(self.depth_value <= root_depths, 0.05, self.depth_value - root_depths)

		self.cn_val = dict_RO[self.lulc_type][self.HSG]

		cn_s = cn_val = self.cn_val
		cn3 = cn_s * np.exp(0.00673 * (100 - cn_s))
		if (self.slope > 5.0):
			cn_s = (((cn3 - cn_val) / 3) * (1 - 2 * np.exp(-13.86 * self.slope * 0.01))) + cn_val
		cn1_s = cn_s - 20 * (100 - cn_s) / (100 - cn_s + np.exp(2.533 - 0.0636 * (100 - cn_s)))
		cn3_s = cn_s * np.exp(0.00673 * (100 - cn_s))

		self.Smax = 25.4 * (1000 / cn1_s - 10)
		S3 = 25.4 * (1000 / cn3_s - 10)
		self.W2 = (np.log((FC_depth - self.WP_depth) / (1 - S3 / self.Smax) - (FC_depth - self.WP_depth)) - np.log(
			(Sat_depth - self.WP_depth) / (1 - 2.54 / self.Smax) - (Sat_depth - self.WP_depth))) / (
				          (Sat_depth - self.WP_depth) - (FC_depth - self.WP_depth))
		self.W1 = np.log((FC_depth - self.WP_depth) / (1 - S3 / self.Smax) - (FC_depth - self.WP_depth)) + self.W2 * (
				FC_depth - self.WP_depth)

		TT_perc = (Sat_depth - FC_depth) / self.Ksat  # SWAT equation 2:3.2.4
		self.daily_perc_factor = 1 - np.exp(-24 / TT_perc)  # SWAT equation 2:3.2.3

	def setup_for_daily_computations(self):
		"""
		"""
		Sat_depth = self.Sat * self.depth_value * 1000
		self.WP_depth = self.WP * self.depth_value * 1000
		FC_depth = self.FC * self.depth_value * 1000

		if (self.depth_value <= self.crop.root_depth):  # thin soil layer
			self.SM1 = self.depth_value - 0.01;
			self.SM2 = 0.01
		else:
			self.SM1 = self.crop.root_depth;
			self.SM2 = self.depth_value - self.crop.root_depth

		self.cn_val = dict_RO[self.lulc_type][self.HSG]
		cn_s = self.cn_val
		cn3 = cn_s * exp(0.00673 * (100 - cn_s))
		if (self.slope > 5.0):
			cn_s = (((cn3 - self.cn_val) / float(3)) * (1 - 2 * exp(-13.86 * self.slope * 0.01))) + self.cn_val
		cn1_s = cn_s - 20 * (100 - cn_s) / float(100 - cn_s + exp(2.533 - 0.0636 * (100 - cn_s)))
		cn3_s = cn_s * exp(0.00673 * (100 - cn_s))

		self.Smax = 25.4 * (1000 / float(cn1_s) - 10)
		S3 = 25.4 * (1000 / float(cn3_s) - 10)
		self.W2 = (log((FC_depth - self.WP_depth) / (1 - float(S3 / self.Smax)) - (FC_depth - self.WP_depth)) - log(
			(Sat_depth - self.WP_depth) / (1 - 2.54 / self.Smax) - (Sat_depth - self.WP_depth))) / (
				          (Sat_depth - self.WP_depth) - (FC_depth - self.WP_depth))
		self.W1 = log((FC_depth - self.WP_depth) / (1 - S3 / self.Smax) - (FC_depth - self.WP_depth)) + self.W2 * (
				FC_depth - self.WP_depth)

		TT_perc = (Sat_depth - FC_depth) / self.Ksat  # SWAT equation 2:3.2.4
		self.daily_perc_factor = 1 - exp(-24 / TT_perc)  # SWAT equation 2:3.2.3

	def primary_runoff_numpy(self, day, rain):
		"""
		Retention parameter 'S_swat' using SWAT equation 2:1.1.6
		Curve Number for the day 'Cn_swat' using SWAT equation 2:1.1.11
		Initial abstractions (surface storage,interception and infiltration prior to runoff)
			'Ia_swat' derived approximately as recommended by SWAT
		Primary Runoff 'Swat_RO' using SWAT equation 2:1.1.1
		"""
		self.budget.sm.append((self.SM1_fraction * self.SM1 + self.layer2_moisture * self.SM2) * 1000)
		self.SW = self.budget.sm[-1] - self.WP_depth
		S_swat = self.Smax * (1 - self.SW / (self.SW + np.exp(self.W1 - self.W2 * self.SW)))
		Cn_swat = 25400 / (S_swat + 254)
		Ia_swat = 0.2 * S_swat
		self.budget.runoff.append(np.where(rain[day] > Ia_swat,
		                                   ((rain[day] - Ia_swat) ** 2) / (rain[day] + 0.8 * S_swat),
		                                   0
		                                   ))
		self.budget.infil.append(rain[day] - self.budget.runoff[day])

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
		S_swat = self.Smax * (1 - self.SW / (self.SW + exp(self.W1 - self.W2 * self.SW)))
		Cn_swat = 25400 / float(S_swat + 254)
		Ia_swat = 0.2 * S_swat
		self.budget.runoff.append(((rain[day] - Ia_swat) ** 2) / (rain[day] + 0.8 * S_swat) if (rain[day] > Ia_swat) else 0)
		self.budget.infil.append(rain[day] - self.budget.runoff[day])

	def aet_numpy(self, day):
		"""
		Water Stress Coefficient 'KS' using FAO Irrigation and Drainage Paper 56, page 167 and
			page 169 equation 84
		Actual Evapotranspiration 'AET' using FAO Irrigation and Drainage Paper 56, page 6 and
			page 161 equation 81
		"""
		depletion_factors = np.array([self.crop.depletion_factor])
		KS = np.where(self.SM1_fraction < self.WP, 0,
		              np.where(self.SM1_fraction > (self.FC * (1 - depletion_factors) + depletion_factors * self.WP), 1,
		                       (self.SM1_fraction - self.WP) / (self.FC - self.WP) / (1 - depletion_factors)
		                       )
		              )
		PETs = np.array([self.crop.PET[day] if day <= self.crop.end_date_index else 0])
		self.budget.AET.append(KS * PETs)

	def aet(self, day):
		"""
		Water Stress Coefficient 'KS' using FAO Irrigation and Drainage Paper 56, page 167 and
			page 169 equation 84
		Actual Evapotranspiration 'AET' using FAO Irrigation and Drainage Paper 56, page 6 and
			page 161 equation 81
		"""
		if (self.SM1_fraction < self.WP):
			KS = 0
		elif (self.SM1_fraction > (self.FC * (1 - self.crop.depletion_factor) + self.crop.depletion_factor * self.WP)):
			KS = 1
		else:
			KS = (self.SM1_fraction - self.WP) / (self.FC - self.WP) / (1 - self.crop.depletion_factor)
		self.budget.AET.append(KS * self.crop.PET[day])

	def percolation_below_root_zone_numpy(self, day):
		"""
		Calculate soil moisture (fraction) 'SM1_before' as the one after infiltration and (then) AET occur,
		but before percolation starts below root-zone. Percolation below root-zone starts only if
		'SM1_before' is more than field capacity and the soil below root-zone is not saturated,i.e.
		'layer2_moisture' is less than saturation. When precolation occurs it is derived as
		the minimum of the maximum possible percolation (using SWAT equation 2:3.2.3) and
		the amount available in the root-zone for percolation.
		"""
		self.SM1_before = (self.SM1_fraction * self.SM1 + (
				(self.budget.infil[day] - self.budget.AET[day]) / 1000)) / self.SM1
		self.R_to_second_layer = np.where(self.SM1_before < self.FC, 0,
		                                  np.where(self.layer2_moisture < self.Sat,
		                                           np.minimum((self.Sat - self.layer2_moisture) * self.SM2 * 1000,
		                                                      (
				                                                      self.SM1_before - self.FC) * self.SM1 * 1000 * self.daily_perc_factor),
		                                           0
		                                           ))
		self.SM2_before = (self.layer2_moisture * self.SM2 * 1000 + self.R_to_second_layer) / self.SM2 / 1000

	def percolation_below_root_zone(self, day):
		"""
		Calculate soil moisture (fraction) 'SM1_before' as the one after infiltration and (then) AET occur,
		but before percolation starts below root-zone. Percolation below root-zone starts only if
		'SM1_before' is more than field capacity and the soil below root-zone is not saturated,i.e.
		'layer2_moisture' is less than saturation. When precolation occurs it is derived as
		the minimum of the maximum possible percolation (using SWAT equation 2:3.2.3) and
		the amount available in the root-zone for percolation.
		"""
		self.SM1_before = (self.SM1_fraction * self.SM1 + (
				(self.budget.infil[day] - self.budget.AET[day]) / float(1000))) / self.SM1
		if (self.SM1_before < self.FC):
			self.R_to_second_layer = 0
		elif (self.layer2_moisture < self.Sat):
			self.R_to_second_layer = min((self.Sat - self.layer2_moisture) * self.SM2 * 1000,
			                             (self.SM1_before - self.FC) * self.SM1 * 1000 * self.daily_perc_factor)
		else:
			self.R_to_second_layer = 0
		self.SM2_before = (self.layer2_moisture * self.SM2 * 1000 + self.R_to_second_layer) / self.SM2 / 1000

	def secondary_runoff_numpy(self):
		"""

		"""
		self.budget.sec_run_off.append(np.where(
			((self.SM1_before * self.SM1 - self.R_to_second_layer / 1000) / self.SM1) > self.Sat,
			(((self.SM1_before * self.SM1 - self.R_to_second_layer / 1000) / self.SM1) - self.Sat) * self.SM1 * 1000,
			0
		))
		self.SM1_fraction = np.minimum((self.SM1_before * self.SM1 * 1000 - self.R_to_second_layer) / self.SM1 / 1000,
		                               self.Sat)

	def secondary_runoff(self):
		"""

		"""
		if (((self.SM1_before * self.SM1 - self.R_to_second_layer / 1000) / self.SM1) > self.Sat):
			self.budget.sec_run_off.append((((self.SM1_before * self.SM1 - self.R_to_second_layer / 1000) / self.SM1) - self.Sat) * self.SM1 * 1000)
		else:
			self.budget.sec_run_off.append(0)
		self.SM1_fraction = min((self.SM1_before * self.SM1 * 1000 - self.R_to_second_layer) / self.SM1 / 1000, self.Sat)

	def percolation_to_GW_numpy(self, day):
		"""

		"""
		self.budget.GW_rech.append(
			np.maximum((self.SM2_before - self.FC) * self.SM2 * self.daily_perc_factor * 1000, 0))
		self.layer2_moisture = np.minimum(
			((self.SM2_before * self.SM2 * 1000 - self.budget.GW_rech[day]) / self.SM2 / 1000), self.Sat)

	def percolation_to_GW(self, day):
		"""

		"""
		self.budget.GW_rech.append(max((self.SM2_before - self.FC) * self.SM2 * self.daily_perc_factor * 1000, 0))
		self.layer2_moisture = min(((self.SM2_before * self.SM2 * 1000 - self.budget.GW_rech[day]) / self.SM2 / 1000),
		                           self.Sat)


class KharifModelCalculator:

	def __init__(self, et0, crop_name, points_data=None, grid_points_roster_layer=None, zones_layer=None, soil_layer=None, lulc_layer=None, slope_layer=None):

		self.crop = Crop(crop_name)
		self.et0 = et0

		self.points_data = points_data  #only useful when READ_FROM_POINTS_DATA_FILE is True
		self.grid_points_roster_layer = grid_points_roster_layer    #only useful when INPUT_FROM_GRID_POINTS_ROSTER is True

	# @property
	# def soil_types(self):	sts = dict_SoilProperties.keys();	sts.remove('soil type');	return sts
	# @property
	# def lulc_types(self):	return dict_RO.keys()
		
	def generate_output_points_grid(self, progress_dlg):
		print 'Initializing output grid-points'
		if INPUT_FROM_GRID_POINTS_ROSTER:
			iterator = self.grid_points_roster_layer.getFeatures()
			grid_row_number_label = 'Grid Row N'
		elif READ_FROM_POINTS_DATA_FILE:
			iterator = self.points_data
			grid_row_number_label = 'Grid Row Number'
		grid_row_number = 1;    output_points_grid = [];	points_grid_row = []
		for it in iterator:
			if int(it[grid_row_number_label]) > grid_row_number:
				print progress_dlg.aborted
				if progress_dlg.aborted:    return None
				grid_row_number += 1
				output_points_grid.append(points_grid_row)
				points_grid_row = []
			points_grid_row.append(Point(QgsPoint(float(it['X']), float(it['Y'])), self.crop, it))
		output_points_grid.append(points_grid_row)
		return output_points_grid
	
	def filter_out_points_with_incomplete_data(self, points, year):
		self.rain_not_found_for = []
		null_token = 'NULL' if INPUT_FROM_GRID_POINTS_ROSTER else ''
		for points_row in points:
			for point in points_row:
				if null_token in [point.texture, point.depth_value, point.lulc_type, point.rain_circle, point.slope]:
					point.is_no_evaluation_point = True
				if (point.district, point.taluka, point.rain_circle, year) not in self.rain:
					point.is_no_evaluation_point = True
					if (point.district, point.taluka, point.rain_circle, year) not in self.rain_not_found_for:
						self.rain_not_found_for.append((point.district, point.taluka, point.rain_circle, year))

	def calculate(self,
					rain,
					sowing_threshold,
					year,
                    end_date_indices=[END_DATE_INDEX],
	                progress_dlg=None
				):
		
		self.rain = rain
		self.rain_year=year

		progress_dlg.progress_text.setText('Initializing grid-points data...')
		progress_dlg.progress_bar.setValue(0)
		self.points_grid = self.generate_output_points_grid(progress_dlg)
		if progress_dlg.aborted:    return
		self.output_grid_points = list(itertools.chain(*self.points_grid))
		self.filter_out_points_with_incomplete_data(self.points_grid,self.rain_year)
		for points_row in self.points_grid:
			for point in points_row:
				if point.lulc_type in ['water', 'habitation']:
					point.is_no_evaluation_point = True
		total_points = len([p for p in self.output_grid_points if not p.is_no_evaluation_point])
		count = 0
		for points_row in self.points_grid:
			for point in points_row:
				if point.is_no_evaluation_point:	continue
				if progress_dlg.aborted:    return
				count += 1
				if count % 100 == 0:
					progress_dlg.progress_text.setText('Running the model for grid-points : ' + str(count) + ' out of ' + str(total_points) + ' points completed')
					progress_dlg.progress_bar.setValue(int(100.0*count/total_points))
					print 'Processing output grid-points : ', count, '/', total_points, ' completed'
				rain_at_point = rain[(point.district, point.taluka, point.rain_circle, year)]
				point.set_PET_and_end_date_index_of_crop(rain_at_point, self.et0[point.district], sowing_threshold)
				point.run_model(rain_at_point, START_DATE_INDEX, end_date_indices)
