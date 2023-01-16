import pytest
import geopandas as gpd
import shapely 
from shapely import wkt
from percentage_cover import bens_donuts
def test_constructor():
    assert bens_donuts.PercentageInDonut()
#self intersecting poly:
#gon
# test data
test_trees = gpd.read_file("test_trees.geojson")
test_trees = test_trees.reset_index(drop = True).geometry 
test_buildings = gpd.read_file("test_buildings.geojson")
test_buildings = test_buildings.reset_index(drop = True).geometry
