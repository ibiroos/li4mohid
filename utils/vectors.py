import netCDF4
import h5py
import numpy as np
import datetime
from math import isnan, pi

from qgis.core import QgsProject, QgsVectorLayer, QgsFeature, QgsField, QgsGeometry, QgsMessageLog, Qgis, QgsPointXY
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QMessageBox


class hydroDataset(self):
    def __init__(self, file, type = 'nc'):
        if file==None:
            quit('A netcdf or hdf5 file is needed')
        else:
            self.file
            self.type = type
            if self.type != 'nc' and self.type !='hdf5':
                quit('type must be nc or hdf5')

    def read(self):
        if self.type=='nc':
            self.ds = netCDF4.Dataset(self.file)
            times_in = getvar_standardname(self.ds, ['time'])
            self.time = netCDF4.num2date(times_in[:], units=times_in.units)[0]
            self.lat = getvar_standardname(self.ds, ['latitude'])[:]
            self.lon = getvar_standardname(self.ds, ['longitude'])[:]
            self.u = getvar_standardname(self.ds, ['surface_eastward_sea_water_velocity',
                                           'eastward_sea_water_velocity'])[:]
            self.v = getvar_standardname(self.ds, ['surface_northward_sea_water_velocity',
                                           'northward_sea_water_velocity'])[:]
        else:
            self.ds = h5py.File(self.file, 'r')
            self.time = self.ds.get('time')
            self.lat = self.ds.get('latitude')
            self.lon = self.ds.get('longitude')
            self.u = self.ds.get('uo')
            self.v = self.ds.get('vo')

    def uv2md(self):
        self.mod = pow((pow(self.u, 2) + pow(self.v, 2)), .5)
        self.dir = (180 * np.arctan2(self.u, self.v)) / pi


def add_current_layer(h_ds):

    project = QgsProject.instance()

    layer = QgsVectorLayer("Point?crs=epsg:4326", "currents", "memory")
    layer.dataProvider().addAttributes([QgsField("water_u", QVariant.Double),
                                        QgsField("water_v", QVariant.Double),
                                        QgsField("mod", QVariant.Double),
                                        QgsField("dir", QVariant.Double)])
    layer.updateFields()
    features = []
    for i, lon in enumerate(h_ds.lon):
        for j, lat in enumerate(h_ds.lat):
            if not isnan(h_ds.u[j][i]):
                feature = QgsFeature()
                feature.setFields(layer.fields())
                pt = QgsPointXY(lon, lat)
                geom = QgsGeometry.fromPointXY(pt)
                feature.setGeometry(geom)
                feature.setAttributes(
                        [float(h_ds.u[j][i]), float(h_ds.v[j][i]), float(h_ds.mod[j][i]), float(h_ds.dir[j][i])])
                features.append(feature)
    layer.dataProvider().addFeatures(features)
    project.addMapLayer(layer)
    return layer





def unix_time(dt):
    """ Seconds since 01_01_1970."""
    epoch = datetime.datetime.utcfromtimestamp(0)
    delta = dt - epoch
    return delta.total_seconds()


def getvar_standardname(f, nome_standards):
    """Return values using the CF standard name of a variable in a netCDF file."""
    for var in f.variables:
        for atributo in (f.variables[var].ncattrs()):
            if atributo == 'standard_name':
                nome_atributo = (getattr(f.variables[var], 'standard_name'))
                for nome_standar in nome_standards:
                    if nome_atributo == nome_standar:
                        return f.variables[var]
    print('standard_name = {0} not found'.format(nome_standar))


def getvar_longname(f, nome_longs):
    """Return values using the CF long name of a variable in a netCDF file."""
    for var in f.variables:
        for atributo in (f.variables[var].ncattrs()):
            if atributo == 'long_name':
                nome_atributo = (getattr(f.variables[var], 'long_name'))
                for nome_long in nome_longs:
                    if nome_atributo == nome_long:
                        return f.variables[var]
    print('long_name = {0} not found'.format(nome_long))









def hf():
    file_in = r'http://150.145.136.27:8080/thredds/dodsC/Ibiza_NRT212/2020/2020_02/2020_02_12/HFR-Ibiza-Total_2020_02_12_1700.nc'
    print('vou a ler {0}'.format(file_in))

    f = netCDF4.Dataset(file_in)

    nc_attrs = f.ncattrs()

    detailed_text = 'NetCDF Global Attributes:\n\n'
    for nc_attr in nc_attrs:
        value = '%s' % repr(f.getncattr(nc_attr), )
        spam = f'- {nc_attr}: {value}; \n'
        detailed_text += spam

    # Radial or Total file
    total = False
    # if f.getncattr('data_type') == 'HF radar total data':
    total = True

    # Extension

    lat_min = float(f.getncattr('geospatial_lat_min'))
    lat_max = float(f.getncattr('geospatial_lat_max'))
    lon_min = float(f.getncattr('geospatial_lon_min'))
    lon_max = float(f.getncattr('geospatial_lon_max'))
    extension = [lat_min, lat_max, lon_min, lon_max]

    # Variables with time

    times_in = getvar_standardname(f, ['time'])
    tempos = netCDF4.num2date(times_in[:], units=times_in.units)
    lat_in = getvar_standardname(f, ['latitude'])[:]
    lon_in = getvar_standardname(f, ['longitude'])[:]

    u_in = getvar_standardname(f, ['surface_eastward_sea_water_velocity',
                                   'eastward_sea_water_velocity'])[:]
    # if u_in is None:
    # u_in = getvar_standardname(f, 'eastward_sea_water_velocity')[:]
    v_in = getvar_standardname(f, ['surface_northward_sea_water_velocity',
                                   'northward_sea_water_velocity'])[:]
    # if v_in is None:
    # v_in = getvar_standardname(f, 'northward_sea_water_velocity')[:]
    print(v_in.shape)

    tempo = netCDF4.num2date(times_in[:], units=times_in.units)[0]
    times = unix_time(tempo)
    print(tempo)

    water_u = u_in[:]
    water_v = v_in[:]
    water_u = water_u[0][0]
    water_v = water_v[0][0]

    mod = pow((pow(water_u, 2) + pow(water_v, 2)), .5)
    dir = (180 * np.arctan2(water_u, water_v)) / pi

    f.close()

    # msg = QMessageBox()
    # msg.setText(f'file: {file_in}')
    # msg.setInformativeText("fe")
    # msg.setWindowTitle("Open NetCDF")
    # msg.setDetailedText(detailed_text)
    # msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
    # msg.exec_()

    proyecto = QgsProject.instance()
    print(proyecto)

    layer = QgsVectorLayer("Point?crs=epsg:4326", "capa_radar", "memory")
    layer.dataProvider().addAttributes([QgsField("water_u", QVariant.Double),
                                        QgsField("water_v", QVariant.Double),
                                        QgsField("mod", QVariant.Double),
                                        QgsField("dir", QVariant.Double)])
    layer.updateFields()
    features = []

    if total:

        for i, lon in enumerate(lon_in):
            for j, lat in enumerate(lat_in):
                if not isnan(water_u[j][i]):
                    feature = QgsFeature()
                    feature.setFields(layer.fields())
                    pt = QgsPointXY(lon, lat)
                    geom = QgsGeometry.fromPointXY(pt)
                    feature.setGeometry(geom)
                    feature.setAttributes(
                        [float(water_u[j][i]), float(water_v[j][i]), float(mod[j][i]), float(dir[j][i])])
                    features.append(feature)
    layer.dataProvider().addFeatures(features)
    proyecto.addMapLayer(layer)

    return layer


def sinrampa(layer):
    # You can alter a single property...
    symbol = QgsMarkerSymbol.createSimple({'name': 'square', 'color': 'red'})
    renderer = QgsSingleSymbolRenderer(symbol)
    layer.setRenderer(renderer)

    layer.renderer().symbol().symbolLayer(0).setSize(3)
    # ... but not all properties are accessible from methods,
    # you can also replace the symbol completely:
    props = layer.renderer().symbol().symbolLayer(0).properties()
    props['color'] = 'red'
    props['name'] = 'arrow'
    layer.renderer().setSymbol(QgsMarkerSymbol.createSimple(props))
    # show the changes
    # layer.triggerRepaint()

    layer.triggerRepaint()


def rampa(layer):
    vals = []
    fld = 'mod'
    for f in layer.getFeatures():
        vals.append(f[fld])
    # If you don't like these colors, change them out for ones you do, using hexcodes,
    # RGB codes etc. as long as the items in this list are valid strings you
    # can pass to a QColor constructor
    colors = ['#0011FF', '#0061FF', '#00D4FF', '#00FF66', '#00FF00', '#E5FF32', '#FCFC0C', '#FF9F00', '#FF3F00',
              '#FF0000']
    lower = sorted(vals)[0]
    upper = sorted(vals)[-1]
    step = (upper - lower) / len(colors)
    range_list = []
    for c in colors:
        cat = [lower, lower + step, c]
        # sym = QgsSymbol.defaultSymbol(layer.geometryType())
        sym = QgsMarkerSymbol.createSimple({'name': 'arrow'})

        sym.setColor(QColor(cat[2]))
        rng = QgsRendererRange(cat[0], cat[1], sym, '{0:.1f}-{1:.1f}'.format(cat[0], cat[1]))
        range_list.append(rng)
        lower = (lower + step) + 0.1
        # sym.setSize(3)
        # sym.setName('arrow')
        # sym.symbolLayer(0)
        print(QgsProperty.fromField(fld).asExpression())
        sym.symbolLayer(0).setDataDefinedProperty(QgsSymbolLayer.PropertyAngle, QgsProperty.fromField('dira'))
    renderer = QgsGraduatedSymbolRenderer(fld, range_list)

    renderer.setSymbolSizes(0, 10)
    layer.setRenderer(renderer)
    layer.triggerRepaint()

    # proyecto.addMapLayer(layer)


def rampa2(layer):
    vals = []
    fld = 'mod'
    for f in layer.getFeatures():
        vals.append(f[fld])
    print(sorted(vals)[-1])
    # If you don't like these colors, change them out for ones you do, using hexcodes,
    # RGB codes etc. as long as the items in this list are valid strings you
    # can pass to a QColor constructor
    colors = ['#0011FF', '#0061FF', '#00D4FF', '#00FF66', '#00FF00', '#E5FF32', '#FCFC0C', '#FF9F00', '#FF3F00',
              '#FF0000']
    lower = sorted(vals)[0]
    upper = sorted(vals)[-1]
    step = (upper - lower) / len(colors)
    print(lower, upper, step)
    range_list = []
    for c in colors:
        cat = [lower, lower + step, c]
        # sym = QgsSymbol.defaultSymbol(layer.geometryType())
        sym = QgsMarkerSymbol.createSimple({'name': 'arrow'})

        sym.setColor(QColor(cat[2]))
        rng = QgsRendererRange(cat[0], cat[1], sym, '{0:.1f}-{1:.1f}'.format(cat[0], cat[1]))
        print(cat[0], cat[1], cat[2])
        range_list.append(rng)
        lower = (lower + step) + 0.001
        # sym.setSize(3)
        # sym.setName('arrow')
        # sym.symbolLayer(0)
        sym.symbolLayer(0).setDataDefinedProperty(QgsSymbolLayer.PropertyAngle, QgsProperty.fromField('dir'))
    renderer = QgsGraduatedSymbolRenderer(fld, range_list)

    renderer.setSymbolSizes(1, 5)
    layer.setRenderer(renderer)
    layer.triggerRepaint()

    # proyecto.addMapLayer(layer)
