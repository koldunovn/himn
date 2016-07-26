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
himn.geojson2shp('remo_cordex_grid_RCA4.geojson', 'remo_cordex_grid_RCA4.shp', zipit=True)
```
# Requirements

- netCDF4
- geojson 
- numpy 
- PyRemo
- pandas
- sh (https://github.com/amoffat/sh)
