import rasterio as rio
import rasterio.mask
import geopandas as gpd
import shapely
from Anomaly_Index import AnomalyIndex
import numpy as np
import fiona
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
from numba import jit


@jit(forceobj=True)
def _mask_to_polygons_layer(mask, out_transform, boundary):
    """AI is creating summary for mask_to_polygons_layer

    Args:
        mask ([type]): [description]
        src ([type]): [description]

    Returns:
        [type]: [description]
    """
    all_polygons = []
    for shape, value in rio.features.shapes(mask.astype(np.uint8), mask=(mask >0), transform = out_transform):
        #return shapely.geometry.shape(shape)
        all_polygons.append(shapely.geometry.shape(shape))
    all_polygons = shapely.geometry.MultiPolygon(all_polygons)
    if not all_polygons.is_valid:
        all_polygons = all_polygons.buffer(0)
        # Sometimes buffer() converts a simple Multipolygon to just a Polygon,
        # need to keep it a Multi throughout
        if all_polygons.type == 'Polygon':
            all_polygons = shapely.geometry.MultiPolygon([all_polygons])
    gdf = gpd.GeoDataFrame({'geometry':all_polygons})
    return gpd.clip(gdf, shapely.geometry.shape(boundary))

def generate_anomaly_contours(in_tiff:str, roof_db_shapefile:str,out_shapefile:str):
    """_summary_

    Args:
        in_tiff (str): _description_
        roof_db_shapefile (str): _description_
        out_shapefile (str): _description_
    """    
    Radius_px_lst = [[5, 50],[2,100]]
    bright_dict = [gpd.GeoDataFrame() for _ in range(len(Radius_px_lst))]
    dark_dict  = [gpd.GeoDataFrame() for _ in range(len(Radius_px_lst))]
    with fiona.Env(CPL_DEBUG=True):
        # This ensures that all drivers are registered in the global
        # context. Within this block *only* GDAL's debugging messages
        # are turned on and the raster block cache size is set to 1 gb.
        with rio.Env(CPL_DEBUG=True, GDAL_CACHEMAX=1073741824):
            with fiona.open(roof_db_shapefile, "r") as shapefile:
                shapes = [feature["geometry"] for feature in shapefile]
                #with alive_bar(len(shapes)) as bar:
                for shape in tqdm(shapes):
                    with rio.open(in_tiff,"r",format='GTiff',dtype = np.uint8) as src:
                        out_image, out_transform = rasterio.mask.mask(src, [shape], crop=True)
                        #sensitivity = 0.25
                        #sensitivity = 0.5
                        #sensitivity = 1
                        #[2,100]
                        
                        for idx, radius_lst in enumerate(Radius_px_lst):
                            anomi_array1  = AnomalyIndex().execute(out_image[0],radius_lst).astype(np.uint8)
                            anomi_array2  = AnomalyIndex().execute(out_image[1],radius_lst).astype(np.uint8)
                            anomi_array3  = AnomalyIndex().execute(out_image[2],radius_lst).astype(np.uint8)
                            anomi_array = anomi_array1 + anomi_array2 + anomi_array3
                            
                            #bright stains 22
                            mask = anomi_array >= 22
                            mask[mask== 1] = 255
                            mask = mask.astype(np.uint8)
                            bright_dict[idx] = bright_dict[idx].append(_mask_to_polygons_layer(mask, out_transform, shape)).reset_index(drop = True)
                            bright_dict[idx].set_crs(src.crs)
                            bright_dict[idx].to_file(driver = 'ESRI Shapefile', filename=out_shapefile+f"bright_{radius_lst}.shp")
                            
                            #dark stains 18
                            mask = anomi_array <= 18
                            mask[mask== 1] = 255
                            mask = mask.astype(np.uint8)
                            dark_dict[idx] = dark_dict[idx].append(_mask_to_polygons_layer(mask, out_transform, shape)).reset_index(drop = True)
                            dark_dict[idx].set_crs(src.crs)
                            dark_dict[idx].to_file(driver = 'ESRI Shapefile', filename=out_shapefile+f"dark_{radius_lst}.shp")