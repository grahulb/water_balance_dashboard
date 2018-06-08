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
from kharif_model_dialog import KharifModelDialog, KharifModelProgressDialog
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
		self.dlg = KharifModelDialog()

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
		# show the dialog and get user response by executing its event-loop
		self.dlg.show()
		if self.dlg.exec_() == QFileDialog.Rejected:    return False

		start_time = time.time()

		self.progress_dlg = KharifModelProgressDialog()
		self.progress_dlg.show()
		self.progress_dlg.progress_text.setText('Reading roster of grid-points...')
		self.progress_dlg.progress_bar.setValue(0)
		points_data = None
		if READ_FROM_POINTS_DATA_FILE:
			reader = csv.DictReader(open(REGION_POINTS_DATA_FILE_PATH, 'r'))
			points_data = list(reader)

		self.fetch_inputs()
		self.modelCalculator = KharifModelCalculator(self.et0, self.crop_name, points_data=points_data, **self.input_layers)
		self.modelCalculator.calculate(self.rain, self.sowing_threshold, self.rain_year, self.end_date_indices, self.progress_dlg)

		if self.progress_dlg.aborted:   self.progress_dlg.close(); return
		self.progress_dlg.close();

		op = KharifModelOutputProcessor()

		for end_date_index in self.end_date_indices:
			layer_name = 'Deficit_on_day_' + str(end_date_index)
			save_file_path = os.path.join(REGION_DATA_FILES_PATH, layer_name + '.tif')
			op.save_pointwise_deficit_as_raster(self.iface, self.modelCalculator.points_grid, end_date_index, save_file_path, layer_name)

		print ("KM--- %s seconds ---" % (time.time() - start_time))

	def fetch_inputs(self):
		if INPUT_FROM_GRID_POINTS_ROSTER:
			district = self.dlg.district.currentText()
			self.input_layers = {
				'grid_points_roster_layer':
					self.iface.addVectorLayer(ROSTER_SHAPEFILES_PATH + '/' + district + '/Grid_points.shp', 'Grid-points', 'ogr')
			}
			self.iface.addVectorLayer(ROSTER_SHAPEFILES_PATH+'/'+district+'/District.shp', district+' District', 'ogr')

		if self.progress_dlg.aborted:   return

		self.progress_dlg.progress_text.setText('Fetching rainfall data...')
		self.progress_dlg.progress_bar.setValue(0)
		print 'Fetching rainfall data'
		with open(os.path.join(ET0_AND_RAINFALL_MASTER_FILES_PATH, 'Rainfall.csv'), 'r') as f:
			reader = csv.reader(f)
			reader.next()
			self.rain = {(row[0], row[1], row[2], row[3]): numpy.array([float(v)    for v in row[4:]])   for row in reader}

		self.progress_dlg.progress_text.setText('Fetching ET0 data...')
		self.progress_dlg.progress_bar.setValue(0)
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

		self.sowing_threshold = self.dlg.sowing_threshold.value()
		self.monsoon_end_date_index = self.dlg.monsoon_end.value()+122
		self.crop_name = self.dlg.selected_crop.currentText()
		self.rain_year=self.dlg.rainfall_year.currentText()
		self.end_date_indices = self.dlg.end_date_indices
		return