from netCDF4 import Dataset
import geojson 
from geojson import Feature, Point, FeatureCollection, Polygon
import matplotlib.pylab as plt
import numpy as np
import Grid as gr
import geojson 
import pandas as pd
import sh
import os

def _minmax(v):
    return np.min(v), np.max(v)

def get_lonlat(ifile, lonname, latname):
    '''
    Load lon lat data from netCDF file. 
    '''
    fl = Dataset(ifile)
    lon = fl.variables[lonname][:]
    lat = fl.variables[latname][:]
    return lon, lat

def shift(lon, lat, lon_resolution=None, lat_resolution=None):
    '''
    Shifts lons and lats by half a gridbox to get lon.lat of verticies. 
    Assumes that lon and lat are 1d and the grid is regular and not global.

    lon - longitude, 1d
    lat - latitude, 1d
    lon_resolution - grid resolution in lon direction (in degrees)
    lat_resolution - grid resolution in lat direction (in degrees)

    '''
    # If resilition is not provided, guessing it from the mean of differences
    if not lon_resolution:
        lon_resolution = (lon[1:]-lon[:-1]).mean()
        print('lon_resolution not provided, guessing: {}'.format(str(lon_resolution)))
    
    if not lat_resolution:
        lat_resolution = (lat[1:]-lat[:-1]).mean()
        print('lat_resolution not provided, guessing: {}'.format(str(lat_resolution)))
    
    # Shift by half a gridbox
    grid_lat = lat-(lat_resolution/2.)
    grid_lat = np.append(grid_lat, lat[-1]+(lat_resolution/2.))
    
    grid_lon = lon-(lon_resolution/2.)
    grid_lon = np.append(grid_lon, lon[-1]+(lon_resolution/2.))
    
    return grid_lon, grid_lat

def make_2d(lon, lat):
    '''
    Create 2d mesh from 1d coordinates
    '''
    if lon.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    return lon, lat

def rotated2geo(lon, lat, pole_latitude, pole_longitude):
    '''
    Use function from PyRemo to convert from rotated to geographical coordinates.

    pole_latitude  - latitude of rotated north pole
    pole_longitude - longitude of rotated north pole

    '''
    grid_lon, grid_lat = gr.rotated_grid_transform(lon, lat, pole_latitude,pole_longitude)
    return grid_lon, grid_lat

def grid2gis(ifile, lonname='rlon', latname='rlat', lon_resolution=None, lat_resolution=None,\
             pole_latitude=None, pole_longitude=None):
    '''
    Combune several functions to convert from lon/lat of grid centers to lon/lat of verticies
    in one go. Assumes that lon and lat are 1d and the grid is regular and not global.

    ifile          - input netCDF file
    lonname        - name of the lon coordinate in the netCDF4 file
    latname        - name of the lat coordinate in the netCDF4 file

    lon_resolution - grid resolution in lon direction (in degrees, will be guessed if not set)
    lat_resolution - grid resolution in lat direction (in degrees, will be guessed if not set)

    pole_latitude  - latitude of rotated north pole (if provided, convertion from rotated to 
                                                     geographical coordinates will be prformed)
    pole_longitude - longitude of rotated north pole (if provided, convertion from rotated to 
                                                     geographical coordinates will be prformed)

    '''
    
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
    '''
    Return indecies for 2d lon lat arrays that surround
    bounding box.

    lon - 2d
    lat - 2d
    bbox - [miLon, maLon, miLat, maLat], where miLon - minimum longitude
                                               maLon - maximum longitude 
                                               miLat - minimum latitude
                                               maLat - maximum latitude
    '''
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
    '''
    Convert lon/lat of verticies to geojson and
    assighn optional values to the 'value' property.

    lon - 2d lons of verticies
    lat - 2d lats of verticies

    bbox - optional bounding box. [miLon, maLon, miLat, maLat],
                                         where miLon - minimum longitude
                                               maLon - maximum longitude 
                                               miLat - minimum latitude
                                               maLat - maximum latitude

    ofile - name of the output file.
    field - field of the shape (lon.shape[0]-1, lon.shape[1]-1),
            that will be assighned to the 'value' property of each grid Polygon.
    writecsv - if True write additional csv file with grid ids and values (from the fiels),
               nessesary for some applications, e.g. plotting with leaflet/folium
    '''
    #fix file extention
    if ofile.split('.')[-1] != 'geojson':
        ofile = ofile+'.geojson'
    
    #Process bbox
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

    # Generate fake 'field' filled with zeros if the 'field' is not provided
    if field is not None:
        tvar = field
    else:
        tvar = np.zeros((lon.shape[0]-1,lon.shape[1]-1))
      
    
    #Create and fill geojson structures
    b = FeatureCollection([])
    
    for x in range(imin,imax):
        for y in range(jmin,jmax):
             
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
    '''
    Just a wraper around ogr2ogr, converts from geojson to shapefile.
    gfile - geojson file name
    sfile - shapefile name
    zipit - optianally zip the result

    '''
    ogr2ogr = sh.ogr2ogr
    ogr2ogr('-nlt', 'POLYGON', '-skipfailures', sfile, gfile, 'OGRGeoJSON')
    if zipit:
        zipfile=os.path.splitext(sfile)[0]+'.zip'
        zipdir=os.path.splitext(sfile)[0]
        os.system('mkdir '+zipdir)
        os.system('cp {}* {}'.format(zipdir, zipdir))
        os.system('zip -r {} {}'.format(zipfile, zipdir))




