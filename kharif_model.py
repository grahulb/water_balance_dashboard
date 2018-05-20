# -*- coding: utf-8 -*-
"""
/***************************************************************************
 KharifModel
                                 A QGIS plugin
 Generates kharif season vulnerability map
                              -------------------
        begin                : 2017-11-18
        git sha              : $Format:%H$
        copyright            : (C) 2017 by IITB
        email                : sohoni@cse.iitb.ac.in
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
"""
from PyQt4.QtCore import QSettings, QTranslator, qVersion, QCoreApplication, QDate, QTimer
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QColor
# Initialize Qt resources from file resources.py
import resources
import time
# Import the code for the dialog
from kharif_model_dialog import KharifModelDialog
from kharif_model_output_processor import KharifModelOutputProcessor
import os.path, csv
# Import code for the calculation
from kharif_model_calculator import KharifModelCalculator
from qgis.core import QgsMapLayerRegistry, QgsVectorLayer, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter
from constants_dicts_lookups import *
from configuration import *
from collections import OrderedDict

class KharifModel:
	"""QGIS Plugin Implementation."""

	def __init__(self, iface):
		"""Constructor.

		:param iface: An interface instance that will be passed to this class
			which provides the hook by which you can manipulate the QGIS
			application at run time.
		:type iface: QgsInterface
		"""
		# Save reference to the QGIS interface
		self.iface = iface
		# initialize plugin directory
		self.plugin_dir = os.path.dirname(__file__)
		# initialize locale
		locale = QSettings().value('locale/userLocale')[0:2]
		locale_path = os.path.join(
			self.plugin_dir,
			'i18n',
			'KharifModel_{}.qm'.format(locale))

		if os.path.exists(locale_path):
			self.translator = QTranslator()
			self.translator.load(locale_path)

			if qVersion() > '4.3.3':
				QCoreApplication.installTranslator(self.translator)


		# Declare instance attributes
		self.actions = []
		self.menu = self.tr(u'&Water Balance Dashboard')
		# TODO: We are going to let the user set this up in a future iteration
		self.toolbar = self.iface.addToolBar(u'WaterBalanceDashboard')
		self.toolbar.setObjectName(u'WaterBalanceDashboard')
		

	# noinspection PyMethodMayBeStatic
	def tr(self, message):
		"""Get the translation for a string using Qt translation API.

		We implement this ourselves since we do not inherit QObject.

		:param message: String for translation.
		:type message: str, QString

		:returns: Translated version of message.
		:rtype: QString
		"""
		# noinspection PyTypeChecker,PyArgumentList,PyCallByClass
		return QCoreApplication.translate('WaterBalanceDashboard', message)


	def add_action(
		self,
		icon_path,
		text,
		callback,
		enabled_flag=True,
		add_to_menu=True,
		add_to_toolbar=True,
		status_tip=None,
		whats_this=None,
		parent=None):
		"""Add a toolbar icon to the toolbar.

		:param icon_path: Path to the icon for this action. Can be a resource
			path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
		:type icon_path: str

		:param text: Text that should be shown in menu items for this action.
		:type text: str

		:param callback: Function to be called when the action is triggered.
		:type callback: function

		:param enabled_flag: A flag indicating if the action should be enabled
			by default. Defaults to True.
		:type enabled_flag: bool

		:param add_to_menu: Flag indicating whether the action should also
			be added to the menu. Defaults to True.
		:type add_to_menu: bool

		:param add_to_toolbar: Flag indicating whether the action should also
			be added to the toolbar. Defaults to True.
		:type add_to_toolbar: bool

		:param status_tip: Optional text to show in a popup when mouse pointer
			hovers over the action.
		:type status_tip: str

		:param parent: Parent widget for the new action. Defaults None.
		:type parent: QWidget

		:param whats_this: Optional text to show in the status bar when the
			mouse pointer hovers over the action.

		:returns: The action that was created. Note that the action is also
			added to self.actions list.
		:rtype: QAction
		"""

		# Create the dialog (after translation) and keep reference
		self.dlg = KharifModelDialog(crops=dict_crop.keys())

		icon = QIcon(icon_path)
		action = QAction(icon, text, parent)
		action.triggered.connect(callback)
		action.setEnabled(enabled_flag)

		if status_tip is not None:
			action.setStatusTip(status_tip)

		if whats_this is not None:
			action.setWhatsThis(whats_this)

		if add_to_toolbar:
			self.toolbar.addAction(action)

		if add_to_menu:
			self.iface.addPluginToMenu(
				self.menu,
				action)

		self.actions.append(action)

		return action

	def initGui(self):
		"""Create the menu entries and toolbar icons inside the QGIS GUI."""

		icon_path = ':/plugins/water_balance_dashboard/icon.png'
		self.add_action(
			icon_path,
			text=self.tr(u'Water Balance Dashboard'),
			callback=self.run,
			parent=self.iface.mainWindow())


	def unload(self):
		"""Removes the plugin menu item and icon from QGIS GUI."""
		for action in self.actions:
			self.iface.removePluginMenu(
				self.tr(u'&Water Balance Dashboard'),
				action)
			self.iface.removeToolBarIcon(action)
		# remove the toolbar
		del self.toolbar


	def run(self):
		"""Run method that performs all the real work"""

		start_time = time.time()
		
		if PLUGIN_MODE == 'DEBUG':
			if not os.path.exists(DEBUG_BASE_FOLDER_PATH):	raise Exception('Set DEBUG_BASE_FOLDER_PATH for the debug dataset')
			paths = [DEBUG_BASE_FOLDER_PATH]
		elif PLUGIN_MODE == 'REAL':
			paths = ['']
		else:
			if not os.path.exists(TEST_SUITE_BASE_FOLDER_PATH):	raise Exception('Set TEST_SUITE_BASE_FOLDER_PATH for the debug dataset')
			paths = [os.path.join(TEST_SUITE_BASE_FOLDER_PATH, base_path)	for base_path in os.listdir(TEST_SUITE_BASE_FOLDER_PATH)	if os.path.isdir(os.path.join(TEST_SUITE_BASE_FOLDER_PATH, base_path))]

		points_values_dict = {};    i = 0
		for path in paths:
			i += 1
			print 'Cluster ', i, '/', len(paths), ' : ', path
			if self.fetch_inputs(path) is False:	return
			self.modelCalculator = KharifModelCalculator(self.path, self.et0, **self.input_layers)
			self.modelCalculator.calculate(self.rain, self.sowing_threshold, end_date_indices=END_DATE_INDICES)

			for point in self.modelCalculator.output_grid_points:
				points_values_dict[point] = point.budget.PET_minus_AET_till_date
			print 'Total points processed till now : ', len(points_values_dict)

			pointwise_output_csv_filepath = os.path.join(self.base_path, POINTWISE_OUTPUT_CSV_FILENAME)

			op = KharifModelOutputProcessor()
			op.output_point_results_to_csv	(
 				self.modelCalculator.output_grid_points,
				pointwise_output_csv_filepath,
				[crop.name for crop in self.modelCalculator.crops],
				END_DATE_INDICES
			)

			self.remove_layers()

		self.iface.addVectorLayer(ALL_CUSTERS_SHAPEFILE, 'Project Clusters', 'ogr')
		for end_date_index in END_DATE_INDICES:
			op_layer = op.render_and_save_pointwise_output_layer(points_values_dict, end_date_index, DEBUG_OR_TEST_GRADUATED_RENDERING_INTERVAL_POINTS, 'Deficit_by_day_'+str(end_date_index))
		self.iface.mapCanvas().setExtent(op_layer.extent())


			# zonewise_budgets = op.compute_zonewise_budget	(
			# 	self.modelCalculator.zone_points_dict ,
			# 	self.modelCalculator.zone_points_dict_ag_missing,
			# 	self.modelCalculator.zone_points_dict_current_fallow,
			# 	self.modelCalculator.zone_points_dict_non_ag_missing_LU,
			# 	self.modelCalculator.zones_layer
			# )
			# op.output_zonewise_budget_to_csv	(
			# 	zonewise_budgets,
			# 	self.modelCalculator.crops,
			# 	self.rabi_crop_names,
			# 	self.modelCalculator.currnet_fallow,
			# 	self.modelCalculator.LULC_pseudo_crops.values(),
			# 	os.path.join(self.base_path, ZONEWISE_BUDGET_CSV_FILENAME),
			# 	self.rain_sum_monsoon
			# )

			# op.compute_and_output_cadastral_vulnerability_to_csv(
			# 	self.crop_names,
			# 	self.modelCalculator.output_cadastral_points,
			# 	os.path.join(self.base_path, CADESTRAL_VULNERABILITY_CSV_FILENAME)
			# )
			# kharif_model_crop_end_output_layer = \
			# 	op.render_and_save_pointwise_output_layer(
			# 		pointwise_output_csv_filepath,
			# 		'Kharif Model Crop End Output',
			# 		'Crop duration PET-AET',
			# 		self.output_configuration['graduated_rendering_interval_points'],
			# 		shapefile_path=os.path.join(self.base_path, 'kharif_crop_duration_et_deficit.shp')
			# 	)
			# if(crop in long_kharif_crops):
			# 	kharif_model_monsoon_end_output_layer = \
			# 		op.render_and_save_pointwise_output_layer(
			# 			pointwise_output_csv_filepath,
			# 			'Kharif Model Monsoon End Output',
			# 			'Monsoon PET-AET',
			# 			self.output_configuration['graduated_rendering_interval_points'],
			# 			shapefile_path=os.path.join(self.base_path, 'kharif_monsoon_et_deficit.shp')
			# 		)
			# for i in range(len(self.crop_names)):
			# 	op.compute_and_display_cadastral_vulnerability(
			# 		self.modelCalculator.cadastral_layer,
			# 		self.modelCalculator.output_grid_points,
			# 		self.modelCalculator.output_cadastral_points,
			# 		i,
			# 		self.crop_names[i],
			# 		self.path
			# 	)
		print ("KM--- %s seconds ---" % (time.time() - start_time))
		print self.plugin_dir

	# self.iface.actionHideAllLayers().trigger()
	# 	self.iface.legendInterface().setLayerVisible(self.input_layers['zones_layer'], True)
	# 	if 'drainage_layer' in locals():	self.iface.legendInterface().setLayerVisible(self.input_layers['drainage_layer'], True)
	# 	if (crop in long_kharif_crops):		self.iface.legendInterface().setLayerVisible(kharif_model_monsoon_end_output_layer, True)
	# 	self.iface.legendInterface().setLayerVisible(kharif_model_crop_end_output_layer, True)
	# 	self.iface.mapCanvas().setExtent(self.input_layers['zones_layer'].extent())
	# 	self.iface.mapCanvas().mapRenderer().setDestinationCrs(self.input_layers['zones_layer'].crs())

		#~ if self.dlg.save_image_group_box.isChecked():
			#~ QTimer.singleShot(1000, lambda :	self.iface.mapCanvas().saveAsImage(self.dlg.save_image_filename.text()))
	
	def fetch_inputs(self, path):
		def set_et0_from_et0_file_data(et0_file_data):
			et0 = []
			for i in range (0,len(et0_file_data)):
				if (i in [0,3,5,10]):	et0.extend([et0_file_data[i]]*30)
				elif i == 8:		et0.extend([et0_file_data[i]]*28)
				else:					et0.extend([et0_file_data[i]]*31)
			return et0
		self.rain=OrderedDict()
		if path != '':
			self.base_path = self.path = path
			self.input_layers = {}
			self.input_layers['zones_layer'] = self.iface.addVectorLayer(os.path.join(path, 'Villages.shp'), 'Zones', 'ogr')
			self.input_layers['soil_layer'] = self.iface.addVectorLayer(os.path.join(path, 'Soil.shp'), 'Soil Cover', 'ogr')
			self.input_layers['lulc_layer'] = self.iface.addVectorLayer(os.path.join(path, 'LULC.shp'), 'Land-Use-Land-Cover', 'ogr')
			self.input_layers['cadastral_layer'] = self.iface.addVectorLayer(os.path.join(path, 'Cadastral.shp'), 'Cadastral Map', 'ogr')
			self.input_layers['slope_layer'] = self.iface.addRasterLayer(os.path.join(path, 'Slope.tif'), 'Slope')
			#~ self.input_layers['drainage_layer'] = self.iface.addRasterLayer(os.path.join(path, 'Drainage.shp'), 'Drainage', 'ogr')

			# data_dir = os.path.join(self.plugin_dir,'Data')

			# self.input_layers['soil_layer'] = self.iface.addVectorLayer(os.path.join(data_dir, 'soil utm.shp'), 'Soil Cover', 'ogr')
			# self.input_layers['lulc_layer'] = self.iface.addVectorLayer(os.path.join(data_dir, 'lulc utm.shp'), 'Land-Use-Land-Cover', 'ogr')
			csvreader=csv.reader(open(os.path.join(path, RAINFALL_CSV_FILENAME)))
			next(csvreader)
			self.rain = OrderedDict((row[0].lower(),{'year': row[1],'daily_rain':[float(val) for val in row[2:]]}) for row in csvreader)
			self.rain['0'] = self.rain[next(iter(self.rain.keys()))]

			et0_file_data =  [float(row["ET0"]) for row in csv.DictReader(open(os.path.join(path, ET0_CSV_FILENAME)))]
			self.et0 = set_et0_from_et0_file_data(et0_file_data)
			self.sowing_threshold = DEFAULT_SOWING_THRESHOLD
			self.monsoon_end_date_index = MONSOON_END_DATE_INDEX
			self.rain_sum_monsoon={key:{'year':self.rain[key]['year'],'sum':sum(self.rain[key]['daily_rain'][START_DATE_INDEX : self.monsoon_end_date_index+1])} for key in self.rain.keys()}

			if not OVERRIDE_FILECROPS_BY_DEBUG_OR_TEST_CROPS and os.path.exists(os.path.join(path, CROPS_FILENAME)):
				self.crop_names = open(os.path.join(path, CROPS_FILENAME), 'r').read().split(',')
				print (self.crop_names)
				if len(self.crop_names) == 0 :	raise Exception('No crop selected')
			else:
				self.crop_names = DEBUG_OR_TEST_CROPS
			self.rabi_crop_names = DEBUG_OR_TEST_RABI_CROPS
			self.output_configuration = {}
			self.output_configuration['graduated_rendering_interval_points'] = DEBUG_OR_TEST_GRADUATED_RENDERING_INTERVAL_POINTS
			
		else:
			self.dlg.show()
			if self.dlg.exec_() == QFileDialog.Rejected:	return False
			
			path = self.path = self.base_path = self.dlg.folder_path.text()
			self.input_layers = {}
			self.input_layers['zones_layer'] = self.iface.addVectorLayer(self.dlg.zones_layer_filename.text(), 'Zones', 'ogr')
			self.input_layers['soil_layer'] = self.iface.addVectorLayer(self.dlg.soil_layer_filename.text(), 'Soil Cover', 'ogr')
			self.input_layers['lulc_layer'] = self.iface.addVectorLayer(self.dlg.lulc_layer_filename.text(), 'Land-Use-Land-Cover', 'ogr')
			self.input_layers['cadastral_layer'] = self.iface.addVectorLayer(self.dlg.cadastral_layer_filename.text(), 'Cadastral Map', 'ogr')
			self.input_layers['slope_layer'] = self.iface.addRasterLayer(self.dlg.slope_layer_filename.text(), 'Slope')
			if self.dlg.drainage_layer_filename.text() != '':
				self.drainage_layer = self.iface.addVectorLayer(self.dlg.drainage_layer_filename.text(), 'Drainage', 'ogr')


			csvreader=csv.reader(open(str(self.dlg.rainfall_csv_filename.text())))
			next(csvreader)
			self.rain = OrderedDict((row[0].lower(),{'year': row[1],'daily_rain':[float(val) for val in row[2:]]}) for row in csvreader)
			self.rain['0'] = self.rain[next(iter(self.rain.keys()))]

			et0_file_data = [float(row["ET0"]) for row in csv.DictReader(open(os.path.join(path, ET0_CSV_FILENAME)))]
			self.et0 = set_et0_from_et0_file_data(et0_file_data)
			self.sowing_threshold = self.dlg.sowing_threshold.value()
			self.monsoon_end_date_index = self.dlg.monsoon_end.value()+122
			self.rain_sum_monsoon={key:{'year':self.rain[key]['year'],'sum':sum(self.rain[key]['daily_rain'][START_DATE_INDEX : self.monsoon_end_date_index+1])} for key in self.rain.keys()}

			self.crop_names = self.dlg.crops
			self.rabi_crop_names = self.dlg.rabi_crops
			if len(self.crop_names) == 0:    raise Exception('No crop selected')
			self.output_configuration = {}
			self.output_configuration['graduated_rendering_interval_points'] = [
				int(self.dlg.colour_code_intervals_list_widget.item(i).text().split('-')[0])
					for i in range(1,self.dlg.colour_code_intervals_list_widget.count())
			]

	def remove_layers(self):
		registry = QgsMapLayerRegistry.instance()
		registry.removeMapLayer(self.input_layers['zones_layer'].id())
		registry.removeMapLayer(self.input_layers['soil_layer'].id())
		registry.removeMapLayer(self.input_layers['lulc_layer'].id())
		registry.removeMapLayer(self.input_layers['cadastral_layer'].id())
		registry.removeMapLayer(self.input_layers['slope_layer'].id())
