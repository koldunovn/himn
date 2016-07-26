from netCDF4 import Dataset
import geojson 
from geojson import Feature, Point, FeatureCollection, Polygon
import matplotlib.pylab as plt
import numpy as np
import PyRemo.OoPlot.Grid as gr
import geojson 
import pandas as pd
import sh
import os

def _minmax(v):
    return np.min(v), np.max(v)

def get_lonlat(ifile, lonname, latname):
    fl = Dataset(ifile)
    lon = fl.variables[lonname][:]
    lat = fl.variables[latname][:]
    return lon, lat

def shift(lon, lat, lon_resolution=None, lat_resolution=None):
    ''' resulution - in degrees, half of the resolution will be used.
                     assume regular regional (no global) grid
    '''
    if not lon_resolution:
        lon_resolution = (lon[1:]-lon[:-1]).mean()
        print('lon_resolution not provided, guessing: {}'.format(str(lon_resolution)))
    
    if not lat_resolution:
        lat_resolution = (lat[1:]-lat[:-1]).mean()
        print('lat_resolution not provided, guessing: {}'.format(str(lat_resolution)))
        
    grid_lat = lat-(lat_resolution/2.)
    grid_lat = np.append(grid_lat, lat[-1]+(lat_resolution/2.))
    
    grid_lon = lon-(lon_resolution/2.)
    grid_lon = np.append(grid_lon, lon[-1]+(lon_resolution/2.))
    
    return grid_lon, grid_lat

def make_2d(lon, lat):
    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    return lon, lat

def rotated2geo(lon, lat, pole_latitude, pole_longitude):
    grid_lon, grid_lat = gr.rotated_grid_transform(lon, lat, pole_latitude,pole_longitude)
    return grid_lon, grid_lat

def grid2gis(ifile, lonname='rlon', latname='rlat', lon_resolution=None, lat_resolution=None,\
             pole_latitude=None, pole_longitude=None):
    
    lon, lat = get_lonlat(ifile, lonname, latname)
    lon, lat = shift(lon, lat, lon_resolution, lat_resolution)
    lon, lat = make_2d(lon, lat)
    if pole_latitude and pole_longitude:
        glon, glat = rotated2geo(lon, lat, pole_longitude, pole_latitude)
    else:
        print('Assume non rotated grid (pole_latitude and pole_longitude not provided)')
        glon, glat = lon, lat
    return glon, glat

def cut_indexes(lon, lat, bbox):
    #bbox = [70,100,25,40]

    inregion = np.logical_and(np.logical_and(lon > bbox[0],
                                             lon < bbox[1]),
                              np.logical_and(lat > bbox[2],
                                             lat < bbox[3]))
    region_inds = np.where(inregion)
    imin, imax = _minmax(region_inds[0])
    jmin, jmax = _minmax(region_inds[1])
    return imin, imax, jmin, jmax

def grid2geojson(lon, lat, bbox=None,  \
       ofile = 'grid2json_output.geojson', field=None, writecsv=False):
    
    if ofile.split('.')[-1] != 'geojson':
        ofile = ofile+'.geojson'
    
    if bbox:
        imin, imax, jmin, jmax = cut_indexes(lon, lat, bbox)
    else:
        imin = 0
        imax=lon.shape[0]-1
        jmin=0
        jmax = lon.shape[1]-1
    
    
    iid = 0
    id_list = []
    temp_list = []

    if field is not None:
        tvar = field
    else:
        tvar = np.zeros((lon.shape[0]-1,lon.shape[1]-1))
      

    b = FeatureCollection([])
    for x in range(imin,imax):
        for y in range(jmin,jmax):
            
            

            inc = 0.0000001

            pp = Polygon([[(float(lon[x,y]), float(lat[x,y])), \
                       (float(lon[x,y+1]), float(lat[x,y+1])), \
                       (float(lon[x+1,y+1]), float(lat[x+1,y+1])), \
                       (float(lon[x+1,y]), float(lat[x+1,y])), \
                       (float(lon[x,y]), float(lat[x,y])),]])

            my_feature = Feature(geometry=pp, id =iid,  properties={'value':float(tvar[x,y]),'gridid':float(iid)})
            b['features'].append(my_feature)
            id_list.append(iid)
            temp_list.append(tvar[x,y])
            iid = iid+1

    with open(ofile, 'w') as outfile:
        geojson.dump(b, outfile)
        
    if writecsv:
        csvofile=os.path.splitext(ofile)[0]+'.csv'
        df = pd.DataFrame({'temp':temp_list}, index=id_list)
        df.index.name = 'id'
        df.to_csv(csvofile)

def geojson2shp(gfile, sfile, zipit=False):
    ogr2ogr = sh.ogr2ogr
    ogr2ogr('-nlt', 'POLYGON', '-skipfailures', sfile, gfile, 'OGRGeoJSON')
    if zipit:
        zipfile=os.path.splitext(sfile)[0]+'.zip'
        zipdir=os.path.splitext(sfile)[0]
        os.system('mkdir '+zipdir)
        os.system('cp {}* {}'.format(zipdir, zipdir))
        os.system('zip -r {} {}'.format(zipfile, zipdir))




