# coding: utf-8
#
# This file is part of PyRemo. PyRemo is a toolbox to facilitate
# treatment and plotting of REMO or other rotated and non-rotated
# data.
#
# Copyright (C) 2010-2014 REMO Group
# See COPYING file for copying and redistribution conditions.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; version 2 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

"""
Defining a grid for rotated and non-rotated coordinates.
"""
import numpy as np
import math


class Grid(object):
    """This class contains gridded geographic coordinates.

    Variables ending with _geo are describing the points in the geographical
    non-rotated coordinates system you find on a globe.

    **Attributes:**
        *lon_arr_geo:*
            longitudes in geographical system (2d-array)
        *lat_arr_geo:*
            latitudes in geographical system (2d-array)
        *cyclic:*
            True if the data is on a cyclic (global) grid
    """
    lon_arr_geo = None  # (array of float)
    lat_arr_geo = None  # (array of float)
    cyclic      = None  # (boolean)


    # Methods
    def __init__(self, lon_arr_geo, lat_arr_geo):
        """Setting lon/lat-array

        **Arguments:**
            *lon_arr_geo:*
                longitudes in geographical coordinate system (2d-array)
            *lat_arr_geo:*
                latitudes in geographical coordinate system (2d-array)
        """
           
        self.lon_arr_geo, self.lat_arr_geo = self.init_lon_lat_arr(
                lon_arr_geo, lat_arr_geo)
        assert(self.lon_arr_geo.shape == self.lat_arr_geo.shape)
        self._check_cyclic(self.lon_arr_geo)


    def __eq__(self, other):
        """Check for equality.

        Two objects are equal if the grids are equal.
        """
        if self.get_dimensions() == other.get_dimensions():
            is_equal = (np.allclose(self.lon_arr_geo, other.lon_arr_geo) and
                        np.allclose(self.lat_arr_geo, other.lat_arr_geo))
        else:
            is_equal = False
        return is_equal


    def __str__(self):
        """Details of the Grid to a string.
        """
        out_tmplt = (
            "lon_arr_geo:\n{lon_arr_geo}\n"
            "lat_arr_geo:\n{lat_arr_geo}\n"
            )
        dic = {'lon_arr_geo': self.lon_arr_geo,
               'lat_arr_geo': self.lat_arr_geo
               }
        return out_tmplt.format(**dic)


    def _check_cyclic(self, lon_arr):
        """Sets the *cyclic* attribute.
        
        Checks if a grid is cylic, which means that the longitude span
        is 360 degree. If it is cyclic the *cyclic* attribute will be set
        to *True*.

        **Arguments:**
            *lon_arr:*
                The longitude array of the grid.

        """
        lon_step = abs(lon_arr[0, 0] - lon_arr[0, 1])
        xdim = lon_arr.shape[1]
        if xdim * lon_step < 360:
            self.cyclic = False
        else:
            self.cyclic = True


    def is_cyclic(self):
        """Returns *True* if the grid is cyclic.

        **Returns:**
            *cyclic:*
                *True* if the grid is cyclic. Otherwise *False*.
        """
        return self.cyclic


    @staticmethod
    def init_lon_lat_arr(lon_arr, lat_arr):
        """Creates two 2d-arrays of two 1d-lon/lat-arrays

        Two 2d-arrays are created if the dimension of both input-arrays is one,
        otherwise the arrays are just returned.
        Values are filled up to the size of the other array to create 2d-arrays
        with dimensions:
        (len(lon_arr), len(lat_arr))

        **Arguments:**
            *lon_arr:*
                array of longitudes (1d or 2d)
            *lat_arr:*
                array of latitudes (1d or 2d)

        **Returns:**
            *lon_arr:*
                array of longitudes (1d or 2d)
            *lat_arr:*
                array of latitudes (1d or 2d)
        """
        tmp_lon = np.array(lon_arr).squeeze()
        tmp_lat = np.array(lat_arr).squeeze()
        if np.ndim(tmp_lon) == np.ndim(tmp_lat) == 1:
            my_lon_arr = np.vstack(len(tmp_lat)*(tmp_lon,))
            my_lat_arr = np.hstack(len(tmp_lon)*(tmp_lat[:, np.newaxis],))
        else:
            my_lon_arr = tmp_lon
            my_lat_arr = tmp_lat
        return my_lon_arr, my_lat_arr

 
    def get_coordinates_geo(self):
        """Returns non-rotated geographical coordinates

        **Returns:**
            *lon_arr:*
                longitudes in geographical coordinates
            *lat_arr:*
                latitudes in geographical coordinates
        """
        return self.lon_arr_geo, self.lat_arr_geo


    def get_dimensions(self):
        """Returns the dimensions of the grid.

        **Returns:**
            *dimensions:*
                dimensions of the grid
        """
        return self.lon_arr_geo.shape


    def get_boundary_as_polygon(self, do_geo=True):
        """Returns lon-lat information at the boundary.

        **Returns:**
            *lon_square:*
                list of longitudes describing the boundary-polygon
            *lat_square:*
                list of latitudes describing the boundary-polygon
        """
        if do_geo:
            xhor, yhor = self.get_coordinates_geo()
        else:
            xhor, yhor = self.get_coordinates_rot()
        dimensions = xhor.shape
        xbottom = xhor[0, :]
        xright = xhor[:, dimensions[1]-1]
        xtop = xhor[dimensions[0]-1, :][::-1]
        xleft = xhor[:, 0][::-1]

        ybottom = yhor[0, :]
        yright = yhor[:, dimensions[1]-1]
        ytop = yhor[dimensions[0]-1, :][::-1]
        yleft = yhor[:, 0][::-1]

        lon_square = np.concatenate((xbottom, xright, xtop, xleft))
        lat_square = np.concatenate((ybottom, yright, ytop, yleft))

        return lon_square, lat_square


    def get_grid_box(self, lon, lat):
        """Returns the grid box containing the given point.

        """
        lon_1d = self.lon_arr_geo[1,:]
        lat_1d = self.lat_arr_geo[:,1]
        box_number_x = self._ret_box_position(lon_1d, lon)
        box_number_y = self._ret_box_position(lat_1d, lat)

        return (box_number_x, box_number_y)


    @staticmethod
    def _ret_box_position(coord_arr, coord):
        """
        """
        min_distance = 360.0
        box_number = 0
        for position, coordinate in enumerate(coord_arr, start=1):
            distance = abs(coordinate - coord)
            if distance < min_distance:
                min_distance = distance
                box_number = position
            elif distance > min_distance:
                break

        return box_number


class RotGrid(Grid):
    """This class contains rotated or non-rotated grid-information.

    We are looking at the same points in different coordinate-systems:

    * Variables ending with _geo are describing the points in the
      geographical non-rotated coordinates system you find on a globe.

    * Variables ending with _rot are describing the same points in the
      rotated coordinates system.

    **Attributes:**
        *lon_arr_rot:*
            longitudes in rotated coordinate system (2d-array)
        *lat_arr_rot:*
            latitudes in rotated coordinate system (2d-array)
        *pol_lon_geo:*
            longitude of rotated North Pole
        *pol_lat_geo:*
            latitude of rotated North Pole
        *lon_arr_geo:*
            longitudes in geographical system (2d-array)
        *lat_arr_geo:*
            latitudes in geographical system (2d-array)
    """
    # Attributes:
    lon_arr_rot = None  # (array of float)
    lat_arr_rot = None  # (array of float)
    pol_lon_geo = None  # (float)
    pol_lat_geo = None  # (float)


    # Methods
    def __init__(self, lon_arr_rot, lat_arr_rot, pol_lon_geo, pol_lat_geo):
        """Setting lon/lat-array and rotated North Pole
        
        **Arguments:**
            *lon_arr_rot:*
                longitudes in rotated coordinate system (2d-array)
            *lat_arr_rot:*
                latitudes in rotated coordinate system (2d-array)
            *pol_lon_geo:*
                longitude of rotated North Pole
            *pol_lat_geo:*
                latitude of rotated North Pole
        """
        self.lon_arr_rot, self.lat_arr_rot = \
                self.init_lon_lat_arr(lon_arr_rot, lat_arr_rot)
        self.pol_lon_geo = pol_lon_geo
        self.pol_lat_geo = pol_lat_geo
        (lon_arr_geo, lat_arr_geo) = self.__calc_coordinates_geo()
        Grid.__init__(self, lon_arr_geo, lat_arr_geo)
        self._check_cyclic(self.lon_arr_rot)


    def __eq__(self, other):
        """Check if two RotGrid objects are equal.

        Two objects are equal if the grids and the rotation are equal
        """
        if Grid.__eq__(self, other):
            is_equal = (np.allclose(self.lon_arr_rot, other.lon_arr_rot) and
                        np.allclose(self.lat_arr_rot, other.lat_arr_rot))
        else:
            is_equal = False
        return is_equal

 
    def __str__(self):
        """Details of the Grid to a string.
        """
        out_tmplt = (
                "Rotated pole (lon/lat): {pollon}/{pollat}\n"
                "lon_arr_rot:\n{lon_arr_rot}\n"
                "lat_arr_rot:\n{lat_arr_rot}\n"
                "lon_arr_geo:\n{lon_arr_geo}\n"
                "lat_arr_geo:\n{lat_arr_geo}\n"
                )
        dic = {'pollon': self.pol_lon_geo,
               'pollat': self.pol_lat_geo,
               'lon_arr_rot': self.lon_arr_rot,
               'lat_arr_rot': self.lat_arr_rot,
               'lon_arr_geo': self.lon_arr_geo,
               'lat_arr_geo': self.lat_arr_geo
               }
        return out_tmplt.format(**dic)


    @staticmethod
    def get_real_coord(phis, rlas, polphi, pollam):
        '''Returns the regular lat/lon coordinate of a rotated coordinate.

        This definition was taken from the REMO model and translated
        into python code.

        **Arguments:**
            *phis:*
                Latitude coordinate of the rotated grid.
            *rlas:*
                Longitude coordinate of the rotated grid.
            *polphi:*
                Latitude coordinate of the rotated pole.
            *pollam:*
                Longitude coordinate of the rotated pole.

        **Returns:**
            *(phstoph, rlstorl):*
                Tuple of the regular coordinates (lat, lon)

        Written by Kevin Sieck

        Last changes 31.10.2010
        '''
        rpi18 = 57.2957795
        pir18 = 0.0174532925

        sinpolp = math.sin(pir18*polphi)
        cospolp = math.cos(pir18*polphi)

        sinpoll = math.sin(pir18*pollam)
        cospoll = math.cos(pir18*pollam)

        sinphis = math.sin(pir18*phis)
        cosphis = math.cos(pir18*phis)

        if rlas > 180.0:
            rlas = rlas - 360.0

        sinrlas = math.sin(pir18*rlas)
        cosrlas = math.cos(pir18*rlas)

        # compute latitude coordinate
        arg = cospolp*cosphis*cosrlas + sinpolp*sinphis
        phstoph = rpi18*math.asin(arg)

        # compute longitude coordinate
        arg1 = sinpoll*(- sinpolp*cosrlas*cosphis +
                          cospolp*sinphis) - cospoll*sinrlas*cosphis
        arg2 = cospoll*(- sinpolp*cosrlas*cosphis +
                          cospolp*sinphis) + sinpoll*sinrlas*cosphis

        if arg2 == 0.0:
            arg2 = pow(10, -20)

        rlstorl = rpi18*math.atan2(arg1, arg2)

        return(phstoph, rlstorl)

    
    def __calc_coordinates_geo(self):
        """Calculates the geographical coordinates.
        
        Geographical (non-rotated) coordinates are calculated.

        **Returns:**
            *lon_arr_geo:*
                longitudes on geographical non-rotated grid.
            *lat_arr_geo:*
                latitudes on geographical non-rotated grid.
        """
        dimensions = self.lon_arr_rot.shape
        lon_arr_geo = np.zeros(dimensions)
        lat_arr_geo = np.zeros(dimensions)
        for i in range(dimensions[1]):       # x-direction
            for j in range(dimensions[0]):   # y-direction
                lon = self.lon_arr_rot[j, i]
                lat = self.lat_arr_rot[j, i]
                #my_coord = Coordinates(lam=lon, phi=lat,
                #                       lam_pol=self.pol_lon_geo,
                #                       phi_pol=self.pol_lat_geo)
                #lon_rot, lat_rot = my_coord.getRotCoordinates()
                lat_rot, lon_rot = self.get_real_coord(
                    lat, lon, self.pol_lat_geo, self.pol_lon_geo)
                lon_arr_geo[j, i] = lon_rot
                lat_arr_geo[j, i] = lat_rot
        return lon_arr_geo, lat_arr_geo


    def get_coordinates_rot(self):
        """Returns rotated coordinates
        
        **Returns:**
            *lon_rot_arr:*
                longitudes in rotated coordinates
            *lat_rot_arr:*
                latitudes in rotated coordinates
        """
        return self.lon_arr_rot, self.lat_arr_rot


    def get_rot_pole(self):
        """Returns the coordinates of the rotated North Pole

        **Returns:**
            *pol_lon_geo:*
                longitude of the rotated North Pole
            *pol_lat_geo:*
                latitude of the rotated North Pole
        """
        return (self.pol_lon_geo, self.pol_lat_geo)


    def get_grid_box(self, lon, lat):
        """Returns the grid box containing the given point.

        """
        lon_1d = self.lon_arr_rot[1,:]
        lat_1d = self.lat_arr_rot[:,1]
        lon_rot, lat_rot = rotated_coord_transform(lon, lat,
                                                   self.pol_lon_geo,
                                                   self.pol_lat_geo,
                                                   direction='geo2rot')
        box_number_x = self._ret_box_position(lon_1d, lon_rot)
        box_number_y = self._ret_box_position(lat_1d, lat_rot)

        return (box_number_x, box_number_y)


def rotated_coord_transform(lon, lat, np_lon, np_lat,
                            direction='rot2geo'):
    """Transforms a coordinate into a rotated grid coordinate and vice versa.

    The coordinates have to given in degree and will be returned in degree.

    **Arguments:**
        *lon:*
            Longitude coordinate.
        *lat:*
            Latitude coordinate.
        *np_lon:*
            Longitude coordinate of the rotated pole.
        *np_lat:*
            Latitude coordinate of the rotated pole.
        *direction:*
            Direction of the rotation.
            Options are: 'rot2geo' (default) for a transformation to regular
            coordinates from rotated. 'geo2rot' transforms regular coordinates
            to rotated.

    **Returns:**
        *lon_new:*
            New longitude coordinate.
        *lat_new:*
            New latitude coordinate.

    Written by Kevin Sieck
    """

    # Convert degrees to radians
    lon = (lon * math.pi) / 180.
    lat = (lat * math.pi) / 180.

#    SP_lon = SP_coor(1)
#    SP_lat = SP_coor(2)

    theta = 90. - np_lat # Rotation around y-axis
    phi = np_lon + 180.  # Rotation around z-axis

    # Convert degrees to radians
    phi = (phi * math.pi) / 180.
    theta = (theta * math.pi) / 180.

    # Convert from spherical to cartesian coordinates
    x = math.cos(lon) * math.cos(lat)
    y = math.sin(lon) * math.cos(lat)
    z = math.sin(lat)

    # Regular -> Rotated
    if direction == 'geo2rot':

        x_new = (math.cos(theta) * math.cos(phi) * x +
                 math.cos(theta) * math.sin(phi) * y +
                 math.sin(theta) * z)
        y_new = (- math.sin(phi) * x +
                   math.cos(phi) * y)
        z_new = (- math.sin(theta) * math.cos(phi) * x -
                   math.sin(theta) * math.sin(phi) * y +
                   math.cos(theta) * z)

    # Rotated -> Regular
    elif direction == 'rot2geo':
    
        phi = - phi
        theta = - theta
    
        x_new = (math.cos(theta) * math.cos(phi) * x +
                 math.sin(phi) * y +
                 math.sin(theta) * math.cos(phi) * z)
        y_new = (- math.cos(theta) * math.sin(phi) * x +
                   math.cos(phi) * y -
                   math.sin(theta) * math.sin(phi) * z)
        z_new = (- math.sin(theta) * x +
                   math.cos(theta) * z)

    # Convert cartesian back to spherical coordinates
    lon_new = math.atan2(y_new, x_new)
    lat_new = math.asin(z_new)

    # Convert radians back to degrees
    lon_new = (lon_new * 180.) / math.pi
    lat_new = (lat_new * 180.) / math.pi;

    return (lon_new, lat_new)


def rotated_grid_transform(lon_arr, lat_arr, np_lon, np_lat,
                           direction='rot2geo'):
    """Transforms a grid into a rotated grid and vice versa.

    The grid coordinates have to given in degree and will be returned in degree.

    **Arguments:**
        *lon_arr:*
            Array with longitude coordinates (at least 2D).
        *lat_arr:*
            Array with latitude coordinates (at least 2D).
        *np_lon:*
            Longitude coordinate of the rotated pole.
        *np_lat:*
            Latitude coordinate of the rotated pole.
        *direction:*
            Direction of the rotation.
            Options are: 'rot2geo' (default) for a transformation to regular
            coordinates from rotated. 'geo2rot' transforms regular coordinates
            to rotated.

    **Returns:**
        *lon_arr_new:*
            New longitude coordinates.
        *lat_arr_new:*
            New latitude coordinates.

    Written by Kevin Sieck
    """
    dimensions = lon_arr.shape
    lon_arr_new = np.zeros(dimensions)
    lat_arr_new = np.zeros(dimensions)
    for i in range(dimensions[1]):       # x-direction
        for j in range(dimensions[0]):   # y-direction
            lon = lon_arr[j, i]
            lat = lat_arr[j, i]
            lon_new, lat_new = rotated_coord_transform(
                lon, lat, np_lon, np_lat, direction=direction)
            lon_arr_new[j, i] = lon_new
            lat_arr_new[j, i] = lat_new

    return (lon_arr_new, lat_arr_new)


def from_cdo_griddes(griddes):
    """Returns a Grid instance from reading a cdo griddes file.

    **Arguments:**
        *griddes:*
            Textfile with a grid description (cdo style).

    **Returns:**
        *grid:*
            Grid instance.

    Written by Kevin Sieck
    """
    
    with open(griddes) as grid_file:
        grid_file_lines = grid_file.readlines()

    grid_dic = {}

    for line in grid_file_lines:
        words = line.split()
        if words[0] == '#':
            continue
        else:
            length = len(words)
            if length == 3:
                grid_dic[words[0]] = words[2]
            else:
                value_string = ' '.join(words[2:length-1])
                grid_dic[words[0]] = value_string

    if grid_dic['gridtype'] != 'lonlat':
        print 'Gridtype {0} not supported'.format(grid_dic['gridtype'])
        return ''

    lon = np.zeros(int(grid_dic['xsize']))
    lat = np.zeros(int(grid_dic['ysize']))

    for i in range(len(lon)):
        lon[i] = float(grid_dic['xfirst']) + i * float(grid_dic['xinc'])
    for j in range(len(lat)):
        lat[j] = float(grid_dic['yfirst']) + j * float(grid_dic['yinc'])

    if grid_dic['xname'] == 'rlon':
        pol_lon = float(grid_dic['xnpole'])
        pol_lat = float(grid_dic['ynpole'])
        grid = RotGrid(lon, lat, pol_lon, pol_lat)
    else:
        grid = Grid(lon, lat)

    return grid
