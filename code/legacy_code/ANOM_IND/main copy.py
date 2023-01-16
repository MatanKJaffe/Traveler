import rasterio as rio
import rasterio.mask
import geopandas as gpd
from itertools import repeat
import shapely
from Anomaly_Index import AnomalyIndex
import numpy as np
import fiona
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
from numba import jit
from multiprocess import Pool
import multiprocessing
from parallelbar import progress_map
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

def run_calc(ziped_args):
    shapes, Radius_px_lst, in_tiff = ziped_args
    #Radius_px_lst = [[5, 50],[2,100]]
    with rio.open(in_tiff,"r",format='GTiff',dtype = np.uint8) as src:
        bright_lst = [gpd.GeoDataFrame({'geometry':[]}) for _ in range(len(Radius_px_lst))]
        dark_lst  = [gpd.GeoDataFrame({'geometry':[]}) for _ in range(len(Radius_px_lst))]
        #for shape in tqdm(shapes):
        for shape in shapes:
            subsample_image, subsample_transform = rasterio.mask.mask(src, [shape], crop=True)
            #sensitivity = 0.25
            #sensitivity = 0.5
            #sensitivity = 1
            #[2,100]
            
            for idx, radius_lst in enumerate(Radius_px_lst):
                anomi_array1  = AnomalyIndex().execute(subsample_image[0],radius_lst).astype(np.uint8)
                anomi_array2  = AnomalyIndex().execute(subsample_image[1],radius_lst).astype(np.uint8)
                anomi_array3  = AnomalyIndex().execute(subsample_image[2],radius_lst).astype(np.uint8)
                anomi_array = anomi_array1 + anomi_array2 + anomi_array3
                
                #bright stains 22
                mask = anomi_array >= 22
                mask[mask== 1] = 255
                mask = mask.astype(np.uint8)
                bright_lst[idx] = bright_lst[idx].append(_mask_to_polygons_layer(mask, subsample_transform, shape)).reset_index(drop = True)
                
                #dark stains 18
                mask = anomi_array <= 18
                mask[mask== 1] = 255
                mask = mask.astype(np.uint8)
                dark_lst[idx] = dark_lst[idx].append(_mask_to_polygons_layer(mask, subsample_transform, shape)).reset_index(drop = True)

            return [bright_lst,dark_lst]
            
def generate_anomaly_contours(in_tiff:str, roof_db_shapefile:str,out_shapefile:str):
    """_summary_

    Args:
        in_tiff (str): _description_
        roof_db_shapefile (str): _description_
        out_shapefile (str): _description_
    """    
    Radius_px_lst = [[5, 50],[2,100]]
    with fiona.Env(CPL_DEBUG=True):
        # This ensures that all drivers are registered in the global
        # context. Within this block *only* GDAL's debugging messages
        # are turned on and the raster block cache size is set to 1 gb.
        with rio.Env(CPL_DEBUG=True, GDAL_CACHEMAX=1073741824):
            with fiona.open(roof_db_shapefile, "r") as shapefile:
                shapes = [feature["geometry"] for feature in shapefile]
                #with alive_bar(len(shapes)) as bar:
                split_shapes = np.array_split(shapes, 16)
                with rio.open(in_tiff,"r",format='GTiff',dtype = np.uint8) as src:
                    ziped_args = list(zip(split_shapes, repeat(Radius_px_lst), repeat(in_tiff)))
                    
                    #result_list_tqdm = progress_map(func = run_calc, tasks = ziped_args, n_cpu= multiprocessing.cpu_count() , chunk_size= multiprocessing.cpu_count(), core_progress=True, #context='spawn', total=len(split_shapes[0]), bar_step=1,
                    #disable=False)
                    with Pool(multiprocessing.cpu_count()) as pool:
                        result_list_tqdm = []
                        #for position in range(multiprocessing.cpu_count()+1):
                        for result in tqdm(pool.imap(func=run_calc, iterable=ziped_args)):
                            result_list_tqdm.append(result)
                        pool.close()
                        pool.join()

                    #with rio.open(in_tiff,"r",format='GTiff',dtype = np.uint8) as src:
                    for idx in tqdm(range(len(Radius_px_lst))):
                        bright_gdf = gpd.GeoDataFrame({'geometry':[]},crs = src.crs)
                        dark_gdf = gpd.GeoDataFrame({'geometry':[]},crs = src.crs)
                        
                        print(f"--------------------{Radius_px_lst(idx)} Radius Discoloration shapefile export started--------------------")
                        for index in tqdm(range(len(result_list_tqdm))):
                            bright_gdf = bright_gdf.concat(result_list_tqdm[index][0]).reset_index(drop = True)
                            dark_gdf = dark_gdf.concat(result_list_tqdm[index][1]).reset_index(drop = True)
                        
                        print(f"--------------------{Radius_px_lst(idx)} Radius Discoloration function mapping finished--------------------")
                        bright_gdf.to_file(driver = 'ESRI Shapefile', filename=out_shapefile+f"bright_{Radius_px_lst(idx)}.shp")
                        dark_gdf.to_file(driver = 'ESRI Shapefile', filename=out_shapefile+f"dark_{Radius_px_lst(idx)}.shp")
                        
                        print(f"--------------------{Radius_px_lst(idx)} Radius Discoloration shapefile export finished--------------------")

def run_imap_multiprocessing(func, argument_list):
    with Pool(multiprocessing.cpu_count()) as pool:
        result_list_tqdm = []
        for result in tqdm(pool.imap(func=func, iterable=argument_list), total=len(argument_list)):
            result_list_tqdm.append(result)
    return result_list_tqdm
                           
@jit(forceobj=True)
def main():
    in_tiff = "Q:/Projects/roof_material_20220118/Miami_Part9/Ortho_Miami_Part1_ecw.ecw.tiff"
    roof_db_shapefile = "Z:/Projects/AI/Model - JM/Deliveries/Annada/Roof Material Data Set/Roof Material/Miami Part 9/Miami_Part_9_Ortho_Miami_Part_1_Buildings_Roof_Material/Miami_Part_9_Ortho_Miami_Part_1_Buildings_Roof_Material.shp"
    out_shapefile = "Z:/Product/Data Preparation/United States/Production for AI/Roof Condition Segmentation/multiprocess test/stains_"
    print ('--------------------Discoloration Detector Started---------------')
    generate_anomaly_contours(in_tiff, roof_db_shapefile,out_shapefile)
    print('---------------------All Done---------------------')
if __name__ == "__main__":
    main()