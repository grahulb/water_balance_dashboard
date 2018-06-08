
"""
/***************************************************************************
 KharifModelDialog
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

import os

from PyQt4 import QtGui, uic
from PyQt4.QtGui import QFileDialog
from configuration import *

from constants_dicts_lookups import dict_crop, district_list, rain_years

FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'kharif_model_dialog_base.ui'))

PROGRESS_FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'progress_dialog.ui'))

class KharifModelDialog(QtGui.QDialog, FORM_CLASS):
	def __init__(self, parent=None):
		"""Constructor."""
		super(KharifModelDialog, self).__init__(parent)
		# Set up the user interface from Designer.
		# After setupUI you can access any designer object by doing
		# self.<objectname>, and you can use autoconnect slots - see
		# http://qt-project.org/doc/qt-4.8/designer-using-a-ui-file.html
		# #widgets-and-dialogs-with-auto-connect
		self.setupUi(self)
		
		self.sowing_threshold.setValue(DEFAULT_SOWING_THRESHOLD)
		self.monsoon_end.setValue(MONSOON_END_DATE_INDEX-122)
		
		self.district.clear()
		self.district.addItems(district_list)

		self.rainfall_year.clear()
		self.rainfall_year.addItems(rain_years)

		self.selected_crop.clear()
		self.selected_crop.addItems(dict_crop.keys())

		self.end_date_indices = []

		self.output_deficit_add_button.clicked.connect(self.on_add)
		self.output_deficit_remove_button.clicked.connect(self.on_remove)

	def on_add(self):
		accumulated_deficit_at = self.output_deficit_date.selectedDate()
		self.deficit_dates_list_widget.addItem(accumulated_deficit_at.toString())
		self.end_date_indices.append(accumulated_deficit_at.dayOfYear()-152)

	def on_remove(self):
		selection = self.deficit_dates_list_widget.currentRow()
		self.deficit_dates_list_widget.takeItem(selection)
		del self.end_date_indices[selection]


class KharifModelProgressDialog(QtGui.QDialog, PROGRESS_FORM_CLASS):
	def __init__(self, parent=None):
		super(KharifModelProgressDialog, self).__init__(parent)
		self.setupUi(self)
		self.progress_text.setText('')
		self.progress_bar.setValue(0)
		self.aborted = False
		def set_aborted():
			self.aborted = True
		self.abort_button.clicked.connect(set_aborted)