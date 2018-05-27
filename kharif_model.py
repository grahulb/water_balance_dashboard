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
import itertools
# Import the code for the dialog
from kharif_model_dialog import KharifModelDialog
from kharif_model_output_processor import KharifModelOutputProcessor
import os.path, csv
# Import code for the calculation
from kharif_model_calculator import KharifModelCalculator
from qgis.core import QgsMapLayerRegistry, QgsVectorLayer, QgsSymbolV2, QgsRendererRangeV2, QgsGraduatedSymbolRendererV2, QgsVectorFileWriter
from constants_dicts_lookups import *
from configuration import *
import numpy

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

		if READ_FROM_POINTS_DATA_FILE:
			reader = csv.DictReader(open(REGION_POINTS_DATA_FILE_PATH, 'r'))
			points_data = list(reader)
		self.fetch_inputs(REGION_DATA_FILES_PATH)
		self.modelCalculator = KharifModelCalculator(self.path, self.et0, points_data=points_data, **self.input_layers)
		self.modelCalculator.calculate(self.rain, self.sowing_threshold, end_date_indices=END_DATE_INDICES)

		# points_values_dict = {point: point.budget.PET_minus_AET_till_date   for point in self.modelCalculator.output_grid_points if not point.is_no_evaluation_point}

		# pointwise_output_csv_filepath = os.path.join(self.base_path, POINTWISE_OUTPUT_CSV_FILENAME)
		#
		op = KharifModelOutputProcessor()
		# op.output_point_results_to_csv	(
         #    self.modelCalculator.output_grid_points,
		# 	pointwise_output_csv_filepath,
		# 	[crop.name for crop in self.modelCalculator.crops],
		# 	END_DATE_INDICES
		# )
		#
		# self.remove_layers()

		for end_date_index in END_DATE_INDICES:
			layer_name = 'Deficit_on_day_' + str(end_date_index)
			save_file_path = os.path.join(REGION_DATA_FILES_PATH, layer_name + '.tif')
			op.save_pointwise_deficit_as_raster(self.iface, self.modelCalculator.points_grid, end_date_index, save_file_path, layer_name)

			# op_layer = op.render_and_save_pointwise_output_layer(points_values_dict, end_date_index, DEBUG_OR_TEST_GRADUATED_RENDERING_INTERVAL_POINTS, 'Deficit_by_day_'+str(end_date_index))
		# self.iface.mapCanvas().setExtent(op_layer.extent())

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
		self.base_path = self.path = path
		self.input_layers = {}

		self.iface.addVectorLayer(os.path.join(path, REGION_SHAPEFILE), 'Zones', 'ogr')
		if not READ_FROM_POINTS_DATA_FILE:
			print 'Loading Zones Layer'
			self.input_layers['zones_layer'] = self.iface.addVectorLayer(os.path.join(path, 'Villages.shp'), 'Zones', 'ogr')
			print 'Loading Soil Layer'
			self.input_layers['soil_layer'] = self.iface.addVectorLayer(os.path.join(path, 'Soil.shp'), 'Soil Cover', 'ogr')
			print 'Loading LULC Layer'
			self.input_layers['lulc_layer'] = self.iface.addVectorLayer(os.path.join(path, 'LULC.shp'), 'Land-Use-Land-Cover', 'ogr')
			print 'Loading Slope Layer'
			self.input_layers['slope_layer'] = self.iface.addRasterLayer(os.path.join(path, 'Slope.tif'), 'Slope')

		print 'Fetching rainfall data'
		with open(os.path.join(ET0_AND_RAINFALL_MASTER_FILES_PATH, 'Rainfall.csv'), 'r') as f:
			reader = csv.reader(f)
			reader.next()
			self.rain = {(row[0], row[1], row[2], row[3]): numpy.array([float(v)    for v in row[4:]])   for row in reader}

		print 'Fetching ET0 data'
		with open(os.path.join(ET0_AND_RAINFALL_MASTER_FILES_PATH, 'ET0.csv'), 'r') as f:
			reader = csv.DictReader(f)
			self.et0 = {
				row['District']: list(itertools.chain(*[
					[float(row['Jun'])] * 30, [float(row['Jul'])] * 31, [float(row['Aug'])] * 31,
					[float(row['Sep'])] * 30, [float(row['Oct'])] * 31, [float(row['Nov'])] * 30,
					[float(row['Dec'])] * 31, [float(row['Jan'])] * 31,	[float(row['Feb'])] * 28,
					[float(row['Mar'])] * 31, [float(row['Apr'])] * 30, [float(row['May'])] * 31
				]))	for row in reader
			}

		self.sowing_threshold = DEFAULT_SOWING_THRESHOLD
		self.monsoon_end_date_index = MONSOON_END_DATE_INDEX

		self.output_configuration = {}
		self.output_configuration['graduated_rendering_interval_points'] = DEBUG_OR_TEST_GRADUATED_RENDERING_INTERVAL_POINTS
