# -*- coding: utf-8 -*-

"""
******************************************************************
                    LI4MOHID QGIS Plugin
******************************************************************

**li4mohid_dockwidget.py**
Generated by Plugin Builder: http://g-sherman.github.io/Qgis-Plugin-Builder/

* *Project:* li4mohid QGIS plugin
* *author:*
  + Carlos F. Balseiro (4Gotas, cfbalseiro@4gotas.com)
  + Pedro Montero (INTECMAR, pmontero@intecmar.gal)
* *license:* Copyright (c) 2020 INTECMAR 2020. Lincesed under MIT
* *funding:* MYCOAST  Interreg Atlantic Programme, Project nr. EAPA 285/2016
             http://www.mycoas-project.org
* *version:* 0.0.1

* *Purpose:* A QGIS plugin to manage MOHID Lagrangian model

"""


import os

from PyQt5 import QtGui, QtWidgets, uic
from PyQt5.QtCore import pyqtSignal
from PyQt5.QtWidgets import QFileDialog

# My imports
from PyQt5.QtCore import QDate, QTime, QDateTime, Qt, QVariant
from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsField, QgsMessageLog, Qgis

from .utils import THREDDS_parser, outputReader, ModelGrid, Application
from .li4mohid_dockwidget_base import Ui_li4mohidDockWidgetBase

'''
FORM_CLASS, _ = uic.loadUiType(os.path.join(
    os.path.dirname(__file__), 'li4mohid_dockwidget_base.ui'))
'''


# class li4mohidDockWidget(QtWidgets.QDockWidget, FORM_CLASS):
class li4mohidDockWidget(QtWidgets.QDockWidget, Ui_li4mohidDockWidgetBase):

    closingPlugin = pyqtSignal()

    def __init__(self, iface, parent=None):
        """Constructor."""
        super(li4mohidDockWidget, self).__init__(parent)
        # Set up the user interface from Designer.
        # After setupUI you can access any designer object by doing
        # self.<objectname>, and you can use autoconnect slots - see
        # http://doc.qt.io/qt-5/designer-using-a-ui-file.html
        # #widgets-and-dialogs-with-auto-connect
        self.setupUi(self)

        application_path = QFileDialog.getExistingDirectory(self, "Select working directory")

        if application_path:
            self.application_path = application_path
        else:
            pass  # Needs to handle closing plugin in case of no dir selected

        self.exe = None

        self.comboBoxHydro.addItems(['artabro', 
                                     'arousa',
                                     'vigo',
                                     'noia',
                                     'iberia',
                                     'tamar',
                                     'portugal'])

        self.comboBoxWind.addItems(['wrf04km', 
                                    'wrf12km'])

        self.wind_state = False
        self.checkBoxWind.stateChanged.connect(self.checkBoxState)

        self.comboBoxHydro.currentIndexChanged.connect(self.enable_calendar)
        self.comboBoxWind.currentIndexChanged.connect (self.enable_calendar)

        self.applyButton.clicked.connect(self.apply)
        self.runButton.clicked.connect(self.run)

        # Access to iface from plugin:
        self.iface = iface

        # Enable initial dates in calendar:
        self.enable_calendar()

    def closeEvent(self, event):

        self.closingPlugin.emit()
        event.accept()

    def apply(self):

        hydro = self.comboBoxHydro.currentText()
        self.app = Application(self.application_path, hydro, self.iface)
        # Auxiliary layers:
        self.app.hydro.get_vectorLayer()
        self.app.define_input_layer()

    def run(self):

        if self.exe is None:
            self.exe, _ = QFileDialog.getOpenFileName(self, "Choose MOHID executable", "",
                                                      "All Files (*);;Exe Files (*.exe)")

        if self.exe is None:
            pass  # Needs to handle closing plugin in case of no executable is selected
        else:
            QgsMessageLog.logMessage('MOHID executable: %s' % self.exe, 'li4mohid', level=Qgis.Info)

        # Get simulation info:
        output = 3600  # Needs input from user
        start = QDateTime(self.dateTimeEditStart.dateTime()).toPyDateTime()
        end = QDateTime(self.dateTimeEditEnd.dateTime()).toPyDateTime()

        # Add wind forcing to application if enabled:
        if self.wind_state:
            meteo = self.comboBoxWind.currentText()
        else:
            meteo = None
        self.app.setDates(start, end, output, meteo)
        # Get defined sources by user:
        self.app.getSources()
        # Write configuration file:
        self.app.write()
        # Write auxiliary data:
        self.app.aux_data()
        # Forcing XML:
        self.app.build_hydro_xml()
        # Add wind forcing if enabled
        if self.wind_state:
            self.app.build_meteo_xml()
        hydro = self.comboBoxHydro.currentText()
        QgsMessageLog.logMessage('%s -i %s/%s.xml -o %s' % (self.exe, self.application_path, hydro, self.application_path), 'li4mohid', level=Qgis.Info)
        os.chdir(self.application_path)
        os.system('%s -i %s/%s.xml -o %s' % (self.exe, self.application_path, hydro, self.application_path))

        # Model results:
        reader = outputReader(self.app.application_path, hydro)
        reader.get_layer()

    def checkBoxState(self):

        self.wind_state = not self.wind_state
        self.comboBoxWind.setEnabled(self.wind_state)

        # Change dates in calendars:
        self.enable_calendar()

    def enable_calendar(self, **kwargs):

        grid = ModelGrid(self.comboBoxHydro.currentText())
        hydro_dates = grid.get_dates()
        print('as hydro_dates = ', hydro_dates)
        timespan = grid.time_span

        if self.wind_state:

            grid = ModelGrid(self.comboBoxWind.currentText())
            wind_dates = grid.get_dates()

            dates = [date for date in hydro_dates if date in wind_dates]
            timespan = min(timespan, grid.time_span)

        else:

            dates = hydro_dates

        self.start = dates[-1]
        self.end  = dates[0]
        print(self.start, self.end)

        # Available date output files:
        self.calendarWidget.setDateRange(QDate(self.start), QDate(self.end))

        self.end += timespan

        QgsMessageLog.logMessage('---> Start: %s End: %s' % (self.start.isoformat(), self.end.isoformat()),
                                 'li4mohid', level=Qgis.Info)

        self.dateTimeEditStart.setDateRange(QDate(self.start), QDate(self.end))
        self.dateTimeEditEnd.setDateRange(QDate(self.start), QDate(self.end))

        self.dateTimeEditStart.setDateTime(QDateTime(self.end))
        self.dateTimeEditEnd.setDateTime(QDateTime(self.end))


