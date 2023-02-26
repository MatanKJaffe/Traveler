
from osgeo import gdal, osr
import os
import geopandas as gpd
import rasterio
from shapely.wkt import loads
from pandas import to_numeric
from geopandas import GeoDataFrame, GeoSeries
from random import choices
from string import ascii_lowercase, digits
from math import atan, degrees
from itertools import combinations
from shapely.geometry import Point, LineString
from matplotlib.pyplot import subplots, show
from warnings import filterwarnings
from haversine import haversine, Unit
from rasterio.sample import sample_gen
filterwarnings("ignore")


def get_all_tifs(path, extension='.tif'):
    output = []
    for root, dirs, files in os.walk(path):
        for file in files:
            if file.endswith(extension):
                output.append(os.path.join(root, file))
    return output

def build_vrt_file(path, vrt_path = "out.vrt",out_vrt = None, extension='.tif', epsg = 3857):
    filepaths = get_all_tifs(path, extension)# list of paths to raster files
    # vrt_options = gdal.BuildVRTOptions(resampleAlg='cubic', addAlpha=True)
    from_dataset = gdal.BuildVRT(vrt_path, filepaths)#, options=vrt_options)
    
    from_dataset = reproject_dataset(from_dataset , epsg, out_vrt)
    from_dataset = None

def reproject_dataset(from_dataset, epsg = 3857, out_vrt="tiles.vrt"):
    """
    Returns the input dataset in the expected "destination" SRS.
    If the dataset is already in the correct SRS, returns it unmodified
    """
    from_srs = from_dataset.GetSpatialRef()
    to_srs = osr.SpatialReference()
    to_srs.ImportFromEPSG(epsg)
    
    if not from_srs or not to_srs:
        raise ValueError("from and to SRS must be defined to reproject the dataset")

    if (from_srs.ExportToProj4() != to_srs.ExportToProj4()) or (from_dataset.GetGCPCount() != 0):
        to_dataset = gdal.AutoCreateWarpedVRT(from_dataset,
                                              from_srs.ExportToWkt(), to_srs.ExportToWkt())

        # if options and options.verbose:
        print(f"Warping of the raster by AutoCreateWarpedVRT (result saved into '{out_vrt}')")
        to_dataset.GetDriver().CreateCopy(out_vrt, to_dataset)
        return to_dataset
    else:
        return from_dataset 
    
if __name__ == "__main__":
    path= "/Users/osn/Code_Library/Traveler/data/italy_dem"
    vrt_path = "/Users/osn/Code_Library/Traveler/data/italy_dem/italy_dem.vrt"
    epsg = 3857
    out_vrt = f"/Users/osn/Code_Library/Traveler/data/italy_dem/italy_dem_{epsg}.vrt" 
    extension='.tif'
    build_vrt_file(path, vrt_path, out_vrt, extension, epsg)
    