# -*- coding: utf-8 -*-
"""
******************************************************************
                    LI4MOHID QGIS Plugin
******************************************************************

**procesa_input.py**


* *Project:* li4mohid QGIS plugin
* *author:*
  + Carlos F. Balseiro (4Gotas, cfbalseiro@4gotas.com)
  + Pedro Montero (INTECMAR, pmontero@intecmar.gal)
* *license:* Copyright (c) 2020 INTECMAR 2020. Lincesed under MIT
* *funding: MYCOAST European Project
* *version:* 0.0.1

* *Purpose:* All classes to process input from THREDDS

"""
import os
import re

from collections import OrderedDict
from datetime import datetime, timedelta
from glob import glob
from urllib import request

from xml.etree import ElementTree
from xml.etree.ElementTree import Element, SubElement, Comment, tostring
from xml.dom import minidom

import numpy as np
import ogr
import vtk

from netCDF4 import Dataset, num2date

from PyQt5.QtCore import QDate, QTime, QDateTime, Qt, QVariant
from PyQt5.QtGui import QColor

from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsField, QgsGeometry, QgsMessageLog, Qgis

from scipy.spatial import cKDTree
from vtk.util.numpy_support import vtk_to_numpy

# import argparse

PLUGIN_NAME = 'li4mohid'


class THREDDS_parser:
    """Class for parse thredds catalog"""

    URL_XML = {
        'artabro': 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/artabro/catalog.xml',
        'arousa': 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/arousa/catalog.xml',
        'vigo': 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/vigo/catalog.xml',
        'noia': 'http://193.144.35.143/thredds/catalog/MyCoast/MOHID/noia/catalog.xml',
        'iberia': 'http://193.144.35.143/thredds/catalog/MyCoast/ROMS/iberia/catalog.xml',
        'tamar': 'https://data.ecosystem-modelling.pml.ac.uk/thredds/catalog/mycoast-all-files/Model/FVCOM/tamar/catalog.xml',
        'portugal': 'http://thredds.maretec.org/thredds/catalog/portugal/catalog.xml',
        'wrf12km': 'http://193.144.35.143/thredds/catalog/MyCoast/WRF/iberia/catalog.xml',
        'wrf04km': 'http://193.144.35.143/thredds/catalog/MyCoast/WRF/galicia/catalog.xml',
        }

    def __init__(self, model):

        self.URL = self.URL_XML[model]

    def parse_dates(self):

        request.urlopen(self.URL)
        contenido = ''.join([linea.decode("utf-8") for linea in request.urlopen(self.URL).readlines()])
        XML = ElementTree.fromstring(contenido)
        filtered = XML.findall('{http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset/{' +
                               'http://www.unidata.ucar.edu/namespaces/thredds/InvCatalog/v1.0}dataset')
        dates = [datetime.strptime(re.findall(r'\d{10}', element.attrib['name'])[0],
                                   '%Y%m%d%H') for element in filtered if '.nc' in element.attrib['name']]

        if dates[0]<dates[-1]:
            dates.reverse()
        return dates


class modelGrid:
    """An abstraction of different model grids that will be used"""

    url_templates = {
        'artabro': 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/artabro/MyCOAST_V1_MeteoGalicia_MOHID_artabro_01hr_%Y%m%d00_PR.ncml',
        'arousa': 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/arousa/MyCOAST_V1_MeteoGalicia_MOHID_arousa_01hr_%Y%m%d00_PR.ncml',
        'vigo': 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/vigo/MyCOAST_V1_MeteoGalicia_MOHID_vigo_01hr_%Y%m%d00_PR.ncml',
        'noia': 'http://193.144.35.143/thredds/dodsC/MyCoast/MOHID/noia/MyCOAST_V1_MeteoGalicia_MOHID_noia_01hr_%Y%m%d00_PR.ncml',
        'iberia': 'http://193.144.35.143/thredds/dodsC/MyCoast/ROMS/iberia/MyCOAST_V1_MeteoGalicia_ROMS_iberia_01hr_%Y%m%d00_PR.ncml',
        'tamar': 'https://data.ecosystem-modelling.pml.ac.uk/thredds/dodsC/mycoast-all-files/Model/FVCOM/tamar/MyCOAST_V0_PML_FVCOM_tamar_01hr_%Y%m%d00_ANPR.ncml',
        'portugal': 'http://thredds.maretec.org/thredds/dodsC/portugal/MyCOAST_V1_IST_MOHID_portugal_03hr_%Y%m%d00_PR.ncml',
        'wrf12km': 'http://193.144.35.143/thredds/dodsC/MyCoast/WRF/iberia/MyCOAST_V1_MeteoGalicia_WRF_iberia_01hr_%Y%m%d00_PR.ncml',
        'wrf04km': 'http://193.144.35.143/thredds/dodsC/MyCoast/WRF/galicia/MyCOAST_V1_MeteoGalicia_WRF_galicia_01hr_%Y%m%d00_PR.ncml',
    }

    def __init__(self, model):

        self.gridName = model
        self.THREDDS_parser = THREDDS_parser(model)

        if model not in self.url_templates.keys():
            QgsMessageLog.logMessage('No template for %s' % model, PLUGIN_NAME, level=Qgis.Critical)
            exit()
        else:
            self.template = self.url_templates[model]

        # Last element of available dates in THREDDS server for grid:
        origen = Dataset(self.THREDDS_parser.parse_dates()[-1].strftime(self.template))

        # This translation is based on standard_name attribute (based on THREDDS data standardisation):
        standard_names_to_var = {}
        for key in origen.variables.keys():
            try:
                standard_names_to_var[origen.variables[key].standard_name] = key
            except:
                pass

        # Search by standard_name attribute:
        self.lon = origen.variables[standard_names_to_var['longitude']][:].astype('double')
        self.lat = origen.variables[standard_names_to_var['latitude']][:].astype('double')
        # Time span of grid file predictions:
        time = origen.variables[standard_names_to_var['time']]
        time = num2date(time[:], time.units)
        self.timespan = time[-1]-time[0]

        if len(self.lon.shape) == 1:
            self.lon, self.lat = np.meshgrid(self.lon, self.lat)
        self.Xmin, self.Ymin = self.lon.min(), self.lat.min()
        self.Xmax, self.Ymax = self.lon.max(), self.lat.max()

    def get_vectorLayer(self):
        """Get model grid contour as a layer for QGIS"""
        vectorlayer = QgsVectorLayer("Linestring?crs=EPSG:4326", "Bounding box", "memory")
        segmento = ogr.Geometry(ogr.wkbLineString)

        for X, Y in zip(self.lon[0, :], self.lat[0, :]):
            segmento.AddPoint(X, Y)
        for X, Y in zip(self.lon[:, -1], self.lat[:, -1]):
            segmento.AddPoint(X, Y)
        for X, Y in zip(self.lon[-1, ::-1], self.lat[-1, ::-1]):
            segmento.AddPoint(X,Y)
        for X, Y in zip(self.lon[::-1, 0], self.lat[::-1, 0]):
            segmento.AddPoint(X, Y)

        geom = QgsGeometry.fromWkt(segmento.ExportToWkt())
        feature = QgsFeature()
        feature.setGeometry(geom)

        pr = vectorlayer.dataProvider()
        pr.addAttributes([QgsField("id", QVariant.Int)])

        vectorlayer.updateFields()

        feature.setAttributes([int(0)])
        pr.addFeature(feature)

        vectorlayer.renderer().symbol().setWidth(0.7)
        vectorlayer.renderer().symbol().setColor(QColor.fromRgb(0, 137, 0))

        proyecto = QgsProject.instance()
        proyecto.addMapLayer(vectorlayer)

    def get_boundingBox(self):
        """ Get aproximate bounding box"""

        return self.Xmin, self.Ymin, self.Xmax, self.Ymax

    def get_dates(self):

        return list(self.THREDDS_parser.parse_dates())


class outputReader:

    def __init__(self, path, model):

        self.path = path
        self.model = model
        self.xml_file = '%s.xml' % model

        root = ElementTree.parse('%s/%s' % (self.path, self.xml_file)).getroot()

        for parameter in root.findall('execution/parameters/parameter'):

            if parameter.get('key') == 'Start':
                start_time = datetime.strptime(parameter.get('value'), '%Y %m %d %H %M %S')
            if parameter.get('key') == 'End':
                end_time = datetime.strptime(parameter.get('value'), '%Y %m %d %H %M %S')
            if parameter.get('key') == 'OutputWriteTime':
                dt = np.float(parameter.get('value'))

        self.ficheros = glob('%s/%s_?????.vtu' % (self.path, self.model))

        fechas = [start_time + timedelta(seconds=int(re.findall(r'\d{5}', fichero)[0])*dt) for fichero in self.ficheros]
        self.fechas = [fecha.strftime('%Y/%m/%d %H:%M') for fecha in fechas]

    def get_layer(self):

        # Feature store:
        features = []

        for fichero, fecha in dict(zip(self.ficheros, self.fechas)).items():

            reader = vtk.vtkXMLUnstructuredGridReader()
            reader.SetFileName(fichero)
            reader.Update()
            nvars = reader.GetOutput().GetPointData().GetNumberOfArrays()

            vars = []

            for i in range( nvars ):
                var = reader.GetOutput().GetPointData().GetArrayName(i)
                # print(var)
                vars.append(var)

            if len(vars) == 1:
                QgsMessageLog.logMessage("No data in file: %s date: %s " %
                                         (fichero, fecha), PLUGIN_NAME, level=Qgis.Info)
                continue
            else:
                QgsMessageLog.logMessage("Processing file: %s date: %s " %
                                         (fichero, fecha), PLUGIN_NAME, level=Qgis.Info)

            # Las coordenadas se cogen así:
            coordenadas = vtk_to_numpy(reader.GetOutput().GetPoints().GetData())
            x = coordenadas[:, 0]
            y = coordenadas[:, 1]

            arrays = {}
            for var in vars:
                arrays[var] = vtk_to_numpy(reader.GetOutput().GetPointData().GetArray(var))

            i = 0
            for X, Y in zip(x, y):
                punto = ogr.Geometry(ogr.wkbPoint)
                punto.AddPoint(X, Y)

                geom = QgsGeometry.fromWkt(punto.ExportToWkt())
                feature = QgsFeature()
                feature.setGeometry(geom)

                feature.setAttributes([int(arrays['id'][i]),
                                       fecha,
                                       int(arrays['source'][i]),
                                       # float(arrays['velocity'] [i]),
                                       int(arrays['state'][i]),
                                       float(arrays['age'][i]),
                                       ])
                features.append(feature)
                i += 1

        vectorlayer = QgsVectorLayer("Point?crs=epsg:4326", "temporary_points", "memory")
        # vectorlayer = QgsVectorLayer("%s/output.shp&Point?crs=epsg:4326" % self.path, "temporary_points", "ogr")

        pr = vectorlayer.dataProvider()
        pr.addAttributes([QgsField("id", QVariant.Int),
                          QgsField("time", QVariant.String),
                          QgsField("source", QVariant.Int),
                          # QgsField("velocity", QVariant.Double),
                          QgsField("state", QVariant.Int),
                          QgsField("age", QVariant.Double),
                          ])
        vectorlayer.updateFields()
        pr.addFeatures(features)
        proyecto = QgsProject.instance()
        proyecto.addMapLayer(vectorlayer)


class application:

    # The template string is the same for all applications::
    input_string = '''
<?xml version="1.0" encoding="UTF-8" ?>
<case>
    <execution>
        <parameters>
            <!--Space for general set up of application-->
        </parameters>
        <outputFields>
            <file name="data/outputFields.xml"/>
        </outputFields>
        <variableNaming>
            <file name="data/NamesLibrary.xml"/>
        </variableNaming>
    </execution>
    <caseDefinitions>
        <inputData>
            <inputDataDir name="nc_fields/hydro/" type="hydrodynamic"/>
            <inputDataDir name="nc_fields/meteo/" type="meteorology"/>
        </inputData>
        <simulation>
            <!--Space for grid and simulation time step-->
        </simulation>
        <sourceDefinitions>
            <!--Space for point sources definition-->
        </sourceDefinitions>
        <constants>
            <BeachingLevel value="-3.0" comment="Level above which beaching can occur. Default = -3.0" units_comment="m" />
            <BeachingStopProb value="80" comment="Probablity of beaching stopping a tracer. Default = 50%" units_comment="%" />
            <DiffusionCoeff value="0.75" comment="Horizontal diffusion coefficient. Default = 1.0" units_comment="m2/s" />
        </constants>
    </caseDefinitions>
</case>
'''

    def __init__(self, APPLICATION_PATH, hydro_in_use, iface):

        self.APPLICATION_PATH = APPLICATION_PATH  # Set working path for application
        self.iface = iface  # Access to QGIS interface from this class
        self.hydro = modelGrid(hydro_in_use)

    def setDates(self, start, end, output, meteo_in_use):

        # Checks whether wind forcing is available
        if meteo_in_use is not None:
            self.meteo = modelGrid(meteo_in_use)

        self.start_time, self.end_time, self.dt = start, end, output

        contenido = re.sub(r"[\n\t]*", "", self.input_string) # Get rid of tabs and new lines
        self.XML  = ElementTree.fromstring(contenido)
        self.XML_INPUTS = Element('file_collection')

        # Execution parameters:
        parameters = self.XML.findall('execution/parameters')[0] # Only one group per file

        parameter = SubElement(parameters, 'parameter', {'key': "Start", 'value': start.strftime('%Y %m %d %H %M %S'), 'comment':"Date of initial instant", 'units_comment':"space delimited ISO 8601 format up to seconds"})
        parameter = SubElement(parameters, 'parameter', {'key': "End", 'value': end.strftime('%Y %m %d %H %M %S')  , 'comment':"Date of final instant"  , 'units_comment':"ISO format"})
        parameter = SubElement(parameters, 'parameter', {'key': "Integrator", 'value': "3"                                , 'comment':"Integration Algorithm 1:Euler, 2:Multi-Step Euler, 3:RK4 (default=1)"})
        parameter = SubElement(parameters, 'parameter', {'key': "Threads", 'value': "4"                                , 'comment':"Computation threads for shared memory computation (default=auto)"})
        parameter = SubElement(parameters, 'parameter', {'key': "OutputWriteTime", 'value':"%d" % output                      , 'comment':"Time out data (1/Hz)"    , 'units_comment':"seconds"})

        # Simulation parameters:
        simulation     = self.XML.findall('caseDefinitions/simulation')[0] # Only one group per file

        resolution     = SubElement(simulation, 'resolution'     ,{'dp':"50"    , 'units_comment':"metres (m)"})
        timestep       = SubElement(simulation, 'timestep'       ,{'dt':"1200.0", 'units_comment':"seconds (s)"})

        # At first, only hydro limits the geographical span of sims:
        Xmin, Ymin, Xmax, Ymax = self.hydro.get_boundingBox()
        BoundingBoxMin = SubElement(simulation, 'BoundingBoxMin' ,{'x':"%f" % Xmin   , 'y':"%f" % Ymin, 'z':"-1", 'units_comment':"(deg,deg,m)"})
        BoundingBoxMax = SubElement(simulation, 'BoundingBoxMax' ,{'x':"%f" % Xmax   , 'y':"%f" % Ymax, 'z': "1", 'units_comment':"(deg,deg,m)"})

    def getSources(self):

        # Incoming data from input layer:
        features = self.iface.activeLayer().getFeatures()

        points = []

        # Accedemos a la linea para obtener los vertices:
        for current, feature in enumerate(features):
            
            listado = {}
            
            listado['id'] = feature.attributes()[feature.fieldNameIndex('id')]
            listado['name'] = feature.attributes()[feature.fieldNameIndex('name')]
            listado['rate'] = feature.attributes()[feature.fieldNameIndex('rate')]
            listado['start'] = feature.attributes()[feature.fieldNameIndex('start')]
            listado['end'] = feature.attributes()[feature.fieldNameIndex('end')]
            
            for point in feature.geometry().vertices():
                listado['geometry'] = (point.x(), point.y())
            points.append(listado)

        # Source definition:
        sourceDefinitions = self.XML.findall('caseDefinitions/sourceDefinitions')[0]  # Only one group per file

        # Remove existing child nodes from XML if any:
        for child in list(sourceDefinitions):
            sourceDefinitions.remove(child) 

        for point in points:

            source = SubElement(sourceDefinitions, 'source')
            setsource = SubElement(source, 'setsource', {'id':'%d' % point['id'], 'name': point['name']})
            rate       = SubElement(source, 'rate', {'value':'%f' % point['rate'], 'comment':'emission rate (Hz)'})
            active     = SubElement(source, 'active', {'start':'%f' % point['start'], 'end':'%f' % point['end'], 'comment':"example: start='12.7' end='end'; start='0.0' end='95' ", 'units_comment':'seconds (s)'})
            point      = SubElement(source, 'point', {'x':'%f' % point['geometry'][0], 'y':'%f' % point['geometry'][1], 'z':'0', 'units_comment':'(deg,deg,m)'})

    # XML prettifier:
    @staticmethod
    def prettify(elem):

        """Return a pretty-printed XML string for the Element."""

        rough_string = ElementTree.tostring(elem, 'utf-8')
        reparsed = minidom.parseString(rough_string)
        return reparsed.toprettyxml(indent="  ")

    def write(self):

        # Writing configuration XML:
        if not os.path.exists(self.APPLICATION_PATH):
            os.makedirs(self.APPLICATION_PATH)

        f = open('%s/%s.xml' % (self.APPLICATION_PATH, self.hydro.gridName),'w')
        f.write(self.prettify(self.XML))
        f.close()

        if DEBUG:
            print(self.prettify(self.XML))


    def aux_data(self):

        outputFields = '''<?xml version="1.0" encoding="UTF-8" ?>
<!-- Basic output fields (always printed) are
-id
-source
-position
-velocity -->
<!--Optional output fields can be listed here
-"yes" - field is written, even if it doesn't exist
-"no" - field is not written -->
<output>
    <field name="age" output="yes" />
    <field name="condition" output="no" />
</output>'''

        NamesLibrary = '''<?xml version="1.0" encoding="UTF-8" ?>
<naming>
    <variables>
        <eastward_wind name="u10">
            <variant name="u10" comment="used in ECWMF"/>
        </eastward_wind>
        <northward_wind name="v10">
            <variant name="v10" comment="used in ECWMF"/>
        </northward_wind>
        <eastward_sea_water_velocity name="u">
            <variant name="u" comment="used in MOHID" />
            <variant name="uu" />
            <variant name="U" />
            <variant name="uo" comment="used in CMEMS" />
        </eastward_sea_water_velocity>
        <northward_sea_water_velocity name="v">
            <variant name="v" comment="used in MOHID" />
            <variant name="vv" />
            <variant name="V" />
            <variant name="vo" comment="used in CMEMS" />
        </northward_sea_water_velocity>
        <upward_sea_water_velocity name="w">
            <variant name="w" comment="used in MOHID and CMEMS" />
            <variant name="W" />
        </upward_sea_water_velocity>
        <sea_water_temperature name="temp">
            <variant name="temp" comment="used in MOHID and CMEMS" />
            <variant name="Temp" />
            <variant name="temperature" />
            <variant name="Temperature" />
        </sea_water_temperature>
        <sea_water_salinity name="salt">
            <variant name="salt" comment="used in MOHID and CMEMS" />
            <variant name="salinity" />
            <variant name="Salt" />
            <variant name="Salinity" />
        </sea_water_salinity>
        <emission_rate name="rate">
            <variant name="rate" />
            <variant name="Rate" />
            <variant name="RATE" />
        </emission_rate>
    </variables>
    <dimensions>
        <longitude name="lon">
            <variant name="lon" />
            <variant name="Lon" />
            <variant name="LON" />
            <variant name="longitude" />
            <variant name="Longitude" />
            <variant name="LONGITUDE" />
        </longitude>
        <latitude name="lat">
            <variant name="lat" />
            <variant name="Lat" />
            <variant name="LAT" />
            <variant name="latitude" />
            <variant name="Latitude" />
            <variant name="LATITUDE" />
        </latitude>
        <vertical name="level">
            <variant name="depth" />
            <variant name="Depth" />
            <variant name="DEPTH" />
            <variant name="level" />
            <variant name="Level" />
        </vertical>
        <time name="time">
            <variant name="time" />
            <variant name="Time" />
            <variant name="TIME" />
        </time>
    </dimensions>
</naming>'''

        if not os.path.exists('%s/data' % self.APPLICATION_PATH):
            os.makedirs('%s/data' % self.APPLICATION_PATH)

        f = open('%s/data/outputFields.xml' % self.APPLICATION_PATH,'w')
        f.write(outputFields)
        f.close()

        f = open('%s/data/NamesLibrary.xml' % self.APPLICATION_PATH,'w')
        f.write(NamesLibrary)
        f.close()

    @staticmethod

    def descarga(f_origen, f_destino, full_flag):

        nt = 24
        if full_flag:
            nt = None
        
        while True:

            try:

                origen = Dataset(filename=f_origen, mode='r', set_auto_mask=False)
                break

            except Exception as e:

                QgsMessageLog.logMessage('Error en el procesamiento de la URL: %s' % f_origen, PLUGIN_NAME, level=Qgis.Critical)
                QgsMessageLog.logMessage('Tipo de la excepcion es: %s' % type(e), PLUGIN_NAME, level=Qgis.Critical)


        variables_origen = [u'time' ,
                  u'longitude'  ,
                  u'latitude'  ,
                  u'uo'    ,
                  u'vo'     ]

        # Optional translation of varnames:
        variables_destino = variables_origen

        if 'iberia' in f_origen:

            lon   = origen.variables['longitude'][0,:]
            lat   = origen.variables['latitude'][:,0]

        else:
            lon   = origen.variables['longitude'][:]
            lat   = origen.variables['latitude'][:]

        
        variables = OrderedDict(zip(variables_destino, variables_origen))

        destino = Dataset(filename=f_destino, mode='w', format='NETCDF4', clobber=True)
        
        destino.createDimension('time', None)
        destino.createDimension('longitude', len(lon))
        destino.createDimension('latitude', len(lat))
       
        for local, remoto in variables.items():

            QgsMessageLog.logMessage('---> Storing variable %s --> %s' % (remoto, local), PLUGIN_NAME, level=Qgis.Info)

            # Variable de origen:       
            variable_origen  = origen.variables[remoto]

            dimensiones = len(variable_origen.dimensions)

            ## Desactivamos el rescalado automatico:
            # variable_origen.set_auto_scale(False)

            # Variable destino:
            print(local)
            if local=='time':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('time',))
                times = num2date(variable_origen[0:nt], variable_origen.units)
                variable_destino[:] = variable_origen[0:nt]

            elif local=='longitude':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('longitude',))
                variable_destino[:] = lon[:]

            elif local == 'latitude':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('latitude',))
                variable_destino[:] = lat[:]

            else:
                variable_destino = destino.createVariable(local, variable_origen.dtype,('time', 'latitude', 'longitude',), fill_value='-9999.')
                # variable_destino = destino.createVariable(local,variable_origen.dtype,('time','latitude','longitude',), fill_value=variable_origen._FillValue)

                if 'iberia' in f_origen:
                    variable_destino[:] = variable_origen[0:nt, -1, :, :]
                elif 'tamar' in f_origen:
                    variable_destino[:] = variable_origen[0:nt, -1, :, :]
                elif 'portugal' in f_origen:
                    variable_destino[:] = variable_origen[0:nt, -1, :, :]
                else:
                    variable_destino[:] = variable_origen[0:nt,:,:]

            atributos = variable_origen.ncattrs()

            for atributo in atributos:

                if atributo not in ['scale_factor','add_offset','_FillValue','_ChunkSize']:

                    variable_destino.setncattr(atributo, variable_origen.getncattr(atributo))


        destino.close()
        origen.close()

        return times[0], times[-1]

    @staticmethod
    def descarga_wrf_alt(f_origen, f_destino, Lon, Lat, full_flag):
        
        origen = Dataset(filename=f_origen, mode='r', set_auto_mask=False)

        variables_origen = [u'time' ,
                            u'longitude'  ,
                            u'latitude'  ,
                            u'u'    ,
                            u'v'     ]

        # Optional translation of varnames:
        variables_destino = [u'time' ,
                            u'longitude'  ,
                            u'latitude'  ,
                            u'u10'  ,
                            u'v10'   ]

        lon = origen.variables['longitude'][:]
        lat = origen.variables['latitude'][:]

        nt = 24

        if full_flag:
            nt = len(origen.dimensions['time'])

        original_points     = np.column_stack((lon.flatten(), lat.flatten()))
        destination_points  = np.column_stack((Lon.flatten(), Lat.flatten()))
        kd     = cKDTree(original_points)
        distancia, indice = kd.query(destination_points)

        variables = OrderedDict(zip(variables_destino, variables_origen))

        destino = Dataset(filename=f_destino, mode='w', format='NETCDF4', clobber=True)
        
        destino.createDimension('time', None)
        destino.createDimension('longitude', len(Lon[0,:]))
        destino.createDimension('latitude', len(Lat[:,0]))
       
        for local, remoto in variables.items():

            QgsMessageLog.logMessage('---> Storing variable %s --> %s' % (remoto, local), PLUGIN_NAME, level=Qgis.Info)

            # Variable de origen:       
            variable_origen  = origen.variables[remoto]

            if local=='time':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('time',))
                times = num2date(variable_origen[0:nt], variable_origen.units)
                variable_destino[:] = variable_origen[0:nt]

            elif local=='longitude':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('longitude',))
                variable_destino[:] = Lon[0,:]

            elif local=='latitude':

                variable_destino = destino.createVariable(local,variable_origen.dtype,('latitude',))
                variable_destino[:] = Lat[:,0]

            else:

                variable_destino = destino.createVariable(local,variable_origen.dtype,('time','latitude','longitude',))

                destination_tmp = np.empty_like(variable_destino[:])
                origin_tmp      = variable_origen[0:nt,:]

                for i in range(nt):

                    destination_tmp[i,:] =  origin_tmp[i,:].flatten()[indice].reshape(Lon.shape)

                variable_destino[:] = destination_tmp[:]

            atributos = variable_origen.ncattrs()

            for atributo in atributos:

                if atributo not in ['scale_factor','add_offset','_FillValue','_ChunkSizes']:

                    variable_destino.setncattr(atributo, variable_origen.getncattr(atributo))


        destino.close()
        origen.close()

        return times[0], times[-1]

    def build_hydro_xml(self):

        dates = [self.start_time + timedelta(days=i) for i in range((self.end_time-self.start_time).days+1)]

        today = datetime.today()
        today = datetime(today.year, today.month, today.day)

        # Correct dates:
        dates = [min(date, today) for date in dates]
        print('dates = ', dates)
        # XML generation:
        file_collection = self.XML_INPUTS
        hydrodynamic    = SubElement(file_collection, 'hydrodynamic')
        print('hydrodynamic = ', hydrodynamic)
        # Remove existing child nodes from XML if any:
        for child in list(hydrodynamic):
            hydrodynamic.remove(child) 

        if not os.path.exists('%s/nc_fields/hydro' % self.APPLICATION_PATH):
            os.makedirs('%s/nc_fields/hydro' % self.APPLICATION_PATH)

        # Loop to add files:
        full_flag = False

        for date in dates:

            if date == today:
                full_flag = True

            QgsMessageLog.logMessage('Downloading date: %s' % date.date().isoformat(), PLUGIN_NAME, level=Qgis.Info)

            fichero_in   = date.strftime(self.hydro.template)
            print('fichero_in ', fichero_in)

            fichero_out  = '%s.nc' % fichero_in.split('/')[-1].split('.')[0]
            print('antes de descarga', date)
            start_hydro, end_hydro = self.descarga(fichero_in, '%s/nc_fields/hydro/%s' % (self.APPLICATION_PATH, fichero_out), full_flag)
            print('despois de descarga', start_hydro)
            dt_inicio = ( start_hydro - self.start_time).total_seconds()
            dt_fin    = ( end_hydro   - self.start_time).total_seconds()

            file            = SubElement(hydrodynamic, 'file')
            name            = SubElement(file, 'name', {'value': 'nc_fields/hydro/%s' % fichero_out})
            startTime       = SubElement(file, 'startTime', {'value': '%10.1f' % dt_inicio})
            endTime         = SubElement(file, 'endTime',   {'value': '%10.1f' % dt_fin})

        # XML resultante:
        if DEBUG:
            pass
        print('adfafa ', self.prettify(file_collection))

        # Stores xml for inputs for further processing:
        self.XML_INPUTS = file_collection

        # Stores xml for inputs that eventually will be overwritted:
        f = open('%s/%s_inputs.xml' % (self.APPLICATION_PATH, self.hydro.gridName), 'w')
        f.write(self.prettify(file_collection))
        f.close()


    def build_meteo_xml(self):

        dates = [self.start_time + timedelta(days=i) for i in range((self.end_time-self.start_time).days+1)]

        today = datetime.today()
        today = datetime(today.year, today.month, today.day)

        # Correct dates:
        dates = [min(date,today) for date in dates]

        # XML generation:
        file_collection = self.XML_INPUTS
        meteorology     = SubElement(file_collection, 'meteorology')

        # Remove existing child nodes from XML if any:
        for child in list(meteorology):
            meteorology.remove(child) 

        if not os.path.exists('%s/nc_fields/meteo' % self.APPLICATION_PATH):
            os.makedirs('%s/nc_fields/meteo' % self.APPLICATION_PATH)

        # Loop to add files:
        full_flag = False

        for date in dates:
            if date == today:
                full_flag = True

            QgsMessageLog.logMessage('Downloading date: %s' % date.date().isoformat(), PLUGIN_NAME, level=Qgis.Info)

            fichero_in = date.strftime(self.meteo.template)
            fichero_out = '%s.nc' % fichero_in.split('/')[-1].split('.')[0]

            start_meteo, end_meteo = self.descarga_wrf_alt(fichero_in, '%s/nc_fields/meteo/%s' % (self.APPLICATION_PATH, fichero_out), self.hydro.lon, self.hydro.lat, full_flag)

            dt_inicio = (start_meteo - self.start_time).total_seconds()
            dt_fin = (end_meteo - self.start_time).total_seconds()

            file = SubElement(meteorology, 'file')
            name = SubElement(file, 'name', {'value': 'nc_fields/meteo/%s' % fichero_out})
            startTime = SubElement(file, 'startTime', {'value': '%10.1f' % dt_inicio})
            endTime = SubElement(file, 'endTime',   {'value': '%10.1f' % dt_fin})

        # XML resultante:
        if DEBUG:
            print(prettify(file_collection))

        # Stores xml for inputs for further processing:
        self.XML_INPUTS = file_collection

        # Stores xml for inputs that eventually will be overwritted:
        f = open('%s/%s_inputs.xml' % (self.APPLICATION_PATH, self.hydro.gridName),'w')
        f.write(self.prettify(file_collection))
        f.close()

    @staticmethod
    def defineInputLayer():

        # Input vector point layer
        vectorlayer = QgsVectorLayer("Point?crs=epsg:4326", "Input points", "memory")

        pr = vectorlayer.dataProvider()

        # Creamos aquí los campos que sea necesario introducir cuando se define un origen:
        pr.addAttributes([QgsField("id"    , QVariant.Int),
                          QgsField("name"  , QVariant.String),
                          QgsField("rate"  , QVariant.Double),
                          QgsField("start" , QVariant.Double),
                          QgsField("end"   , QVariant.Double),
                          ])

        vectorlayer.updateFields()

        proyecto = QgsProject.instance()
        proyecto.addMapLayer(vectorlayer)


DEBUG = False

'''
# Model results:
reader = outputReader(app.APPLICATION_PATH, app.hydro.gridName)
reader.get_layer()
'''
