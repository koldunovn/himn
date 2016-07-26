# himn
Convert lon/lat grid to geojson

# Basic usage
```python
import himn

# Read lon lat coordinates, shift them to get verticies and convert to geographical coordinates if needed 
lon, lat = himn.grid2gis('./input.nc', pole_latitude=79.95, pole_longitude=-123.34)

# Convert to geojson and cut specific region if needed
himn.grid2geojson(lon, lat, [70,100,25,38], ofile = 'remo_grid.geojson')

# Convert from geojson to shapefile. Have to have working ogr2org. Alternativelly can use web service http://ogre.adc4gis.com/
himn.geojson2shp('remo_grid.geojson', 'remo_grid.shp', zipit=True)
```

# WRF grid case
```python
import himn
from netCDF4 import Dataset
import numpy as np

flw = Dataset('./data/geo_em.d01.nc')

lonwU = flw.variables['XLONG_U'][0,:,:]
latwV = flw.variables['XLAT_V'][0,:,:]

#expand grid
ddlon = lonwU[-1,:] + (lonwU[-2,:]-lonwU[-1,:])
ddlat = latwV[:,-1] + (latwV[:,-2]-latwV[:,-1])

newlon = np.row_stack((lonwU, ddlon[None,:]))
newlat = np.column_stack((latwV, ddlat[:, None]))

himn.grid2geojson(newlon, newlat, [70,100,25,38], ofile='wrg_glacindia25km_grid.geojson')
himn.geojson2shp('wrg_glacindia25km_grid.geojson', 'wrg_glacindia25km_grid.shp', zipit=True)
```

# Requirements

- netCDF4
- geojson 
- numpy 
- PyRemo
- pandas
- sh (https://github.com/amoffat/sh)
