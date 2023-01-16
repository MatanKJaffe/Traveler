from geopandas import read_file, GeoDataFrame, GeoSeries
import geopandas as gpd
from collections import defaultdict
from warnings import filterwarnings
from numpy import ndarray
from multiprocessing import cpu_count
import inspect

filterwarnings("ignore")

def aux_retrieve_name(var):
    callers_local_vars = inspect.currentframe().f_back.f_back.f_locals.items()
    return [var_name for var_name, var_val in callers_local_vars if var_val is var]

def _zip_components(first_results:ndarray, second_results:ndarray) -> defaultdict:
    """
    zips the contents of two np arrays and return a dictionary
    Args:
        first_results (ndarray): ndarray of elemnts to be matched
        second_results (ndarray): ndarray of elemnts to be matched

    Returns:
        defaultdict: a defaultdict of the shape {first_results element as key:second_results elements as values}
    
    Example: 
            >>> df1 = gpd.GeoDataFrame(geometry = [Point(0, 0), Point(0, 0)])
            >>> df2 = gpd.GeoDataFrame(geometry = [Point(1, 0), Point(0, 0)])
            >>> df2_results, df1_results = df1.sindex.query_bulk(df2.geometry,predicate = 'intersects')
            >>> _zip_components(df1_results, df2_results)
            >>>  #Will result in a dictionary of 0: [1] 1: [1]
    """
    #define output defaultdict
    collections_dict = defaultdict(list)
    #loop to match value from first_results to second_results index
    for v, j in zip(first_results, second_results):
        collections_dict[v].extend([j])
    
    #return matched defaultdict
    return collections_dict

def orginal_percentage_in_donut(building_gdf:GeoDataFrame, object_gdf:GeoDataFrame, radius_lst:list = []) ->GeoDataFrame:
    """
    calculates the percentage of overlap between one geodataframe and another by value in a specified range. if only an outer range is required than set one of the radius_lst values as 1

    Args:
        building_gdf (GeoDataFrame): a geodataframe containing footprints 
        object_gdf (GeoDataFrame):  a geodataframe containing objects to be searched (right now its trees) 
        radius_lst (list, optional): list of radius to be calculated. Defaults to empty list.

    Returns:
        GeoDataFrame: added intersect_percent_at_radius_lst_m column to building_gdf
    """

    assert len(radius_lst) == 2, f"please enter two values into the search radius {radius_lst} argument" 
    assert building_gdf.crs.equals(object_gdf.crs), "make sure the given geodataframes have the same crs"
    assert (building_gdf.crs.axis_info[0].unit_name == "metre") & (object_gdf.crs.axis_info[0].unit_name == "metre"), "make sure the given geodataframes crs is in meteres"
    #assert any(building_gdf.is_empty), f"building_gdf has an empty geometry"
    #assert any(object_gdf.is_empty), f"object_gdf has an empty geometry " 

    #set our querry polygons
    #diff_envelope =  building_gdf.geometry.envelope.map(lambda  geom : geom.buffer(max(radius_lst) , join_style = 2, cap_style = 3).difference(geom.buffer(min(radius_lst), join_style = 2, cap_style = 3)))
    
    small_envelope =  building_gdf.geometry.envelope.buffer(min(radius_lst), join_style = 2, cap_style = 3)
    
    big_envelope =  building_gdf.geometry.envelope.buffer(max(radius_lst), join_style = 2, cap_style = 3)
    #querry polygons intersection with our object_gdf
    objects_results, R_results = big_envelope.sindex.query_bulk(object_gdf.geometry, predicate='intersects')

    # zip the index of each querry polygon to the object_gdf polygons that intersect with it 
    indx_dict = _zip_components(R_results, objects_results)

    #initialize a column in the building_gdf to contain the functions output
    building_gdf[f"intersect_percent_at_{radius_lst}_m"] = None

    #iterate over each polygon and calculate its overlap percentage with polygons from the object_gdf
    for i in range(len(building_gdf)):
        #define the search polygon for this iteration
        #donut =  building_gdf.geometry.envelope.buffer(min(radius_lst), join_style = 2, cap_style = 3)[i].difference(small_envelope[i])
        donut = big_envelope[i].difference(small_envelope[i])
        donut_area = donut.area

        #get the intersection polygons
        intersecty = object_gdf.geometry[indx_dict[i]].intersection(donut)
        intersecty = intersecty[~intersecty.is_empty]

        #if any intersections were found than add the raltional area to the output
        print()
        if intersecty is not None:
            building_gdf[f"intersect_percent_at_{radius_lst}_m"][i] = int((sum(intersecty.area)/donut_area)*100)
    
    #return building_gdf with relevant data
    return building_gdf

def simplified_percentage_in_donut(building_gdf:GeoDataFrame, object_gdf:GeoDataFrame, radius_lst:list = []) ->GeoDataFrame:
    """
    calculates the percentage of overlap between one geodataframe and another by value in a specified range. if only an outer range is required than set one of the radius_lst values as 1

    Args:
        building_gdf (GeoDataFrame): a geodataframe containing footprints 
        object_gdf (GeoDataFrame):  a geodataframe containing objects to be searched (right now its trees) 
        radius_lst (list, optional): list of radius to be calculated. Defaults to empty list.

    Returns:
        GeoDataFrame: added intersect_percent_at_radius_lst_m column to building_gdf
    """

    assert len(radius_lst) == 2, f"please enter two values into the search radius {radius_lst} argument" 
    assert building_gdf.crs.equals(object_gdf.crs), "make sure the given geodataframes have the same crs"
    assert (building_gdf.crs.axis_info[0].unit_name == "metre") & (object_gdf.crs.axis_info[0].unit_name == "metre"), "make sure the given geodataframes crs is in meteres"
    #assert any(building_gdf.is_empty), f"building_gdf has an empty geometry"
    #assert any(object_gdf.is_empty), f"object_gdf has an empty geometry " 

    #set our querry polygons
    diff_envelope =  building_gdf.geometry.envelope.map(lambda  geom : geom.buffer(max(radius_lst) , join_style = 2, cap_style = 3).difference(geom.buffer(min(radius_lst), join_style = 2, cap_style = 3)))
    
    #querry polygons intersection with our object_gdf
    objects_results, R_results = diff_envelope.sindex.query_bulk(object_gdf.geometry, predicate='intersects')

    # zip the index of each querry polygon to the object_gdf polygons that intersect with it 
    indx_dict = _zip_components(R_results, objects_results)
    diff_envelope = diff_envelope.rename("geometry").pipe(GeoDataFrame)
    #initialize a column in the building_gdf to contain the functions output
    intersecty = diff_envelope.apply(lambda row: int(sum(object_gdf.geometry[indx_dict[row.name]].intersection(row["geometry"]).area) / row["geometry"].area *100) , axis= 1)
    building_gdf[f"intersect_percent_at_{radius_lst}_m"] = intersecty
    
    #return building_gdf with relevant data
    return building_gdf

def parallel_percentage_in_donut(building_gdf:GeoDataFrame, object_gdf:GeoDataFrame, radius_lst:list = []) ->GeoDataFrame:
    from dask_geopandas import from_geopandas
    """
    calculates the percentage of overlap between one geodataframe and another by value in a specified range. if only an outer range is required than set one of the radius_lst values as 1

    Args:
        building_gdf (GeoDataFrame): a geodataframe containing footprints 
        object_gdf (GeoDataFrame):  a geodataframe containing objects to be searched (right now its trees) 
        radius_lst (list, optional): list of radius to be calculated. Defaults to empty list.

    Returns:
        GeoDataFrame: added intersect_percent_at_radius_lst_m column to building_gdf
    """

    assert len(radius_lst) == 2, f"please enter two values into the search radius {radius_lst} argument" 
    assert building_gdf.crs.equals(object_gdf.crs), "make sure the given geodataframes have the same crs"
    assert (building_gdf.crs.axis_info[0].unit_name == "metre") & (object_gdf.crs.axis_info[0].unit_name == "metre"), "make sure the given geodataframes crs is in meteres"
    #assert any(building_gdf.is_empty), f"building_gdf has an empty geometry"
    #assert any(object_gdf.is_empty), f"object_gdf has an empty geometry " 
    building_gdf = from_geopandas(building_gdf, npartitions=cpu_count())
    #set our querry polygons
    diff_envelope =  building_gdf.geometry.envelope
    diff_envelope = diff_envelope.map_partitions(lambda  geom : geom.buffer(max(radius_lst) , join_style = 2, cap_style = 3).difference(geom.buffer(min(radius_lst), join_style = 2, cap_style = 3)))

    #querry polygons intersection with our object_gdf
    objects_results, R_results = diff_envelope.compute().sindex.query_bulk(object_gdf.geometry, predicate='intersects')

    # zip the index of each querry polygon to the object_gdf polygons that intersect with it 
    indx_dict = _zip_components(R_results, objects_results)
    diff_envelope = diff_envelope.pipe(GeoDataFrame).rename(columns= {0:"geometry"})

    #initialize a column in the building_gdf to contain the functions output
    building_gdf[f"intersect_percent_at_{radius_lst}_m"] = from_geopandas(diff_envelope, npartitions=cpu_count()).apply(lambda row: int(sum(object_gdf.geometry[indx_dict[row.name]].intersection(row["geometry"]).area) / row["geometry"].area *100) , axis= 1, meta =(None, int)).compute()

    #return building_gdf with relevant data
    return building_gdf.compute()
    
def percentage_in_buffer(building_gdf:gpd.GeoDataFrame, object_gdf:gpd.GeoDataFrame, radius_lst:list, column_name: str) ->gpd.GeoDataFrame:
    """
    calculates the percentage of overlap between one geodataframe and another by value in a specified range. if only an outer range is required than set one of the radius_lst values as 1

    Args:
        building_gdf (GeoDataFrame): a geodataframe containing footprints 
        object_gdf (GeoDataFrame):  a geodataframe containing objects to be searched (right now its trees) 
        radius_lst (list): list of radius to be calculated. Defaults to empty list.
        column_name (str): tag for the name of the column based on object type

    Returns:
        GeoDataFrame: added intersect_percent_at_radius_lst_m column to building_gdf
    """

    #assert len(radius_lst) == 2, f"please enter two values into the search radius {radius_lst} argument" 
    #assert building_gdf.crs.equals(object_gdf.crs), "make sure the given geodataframes have the same crs"
    #assert (building_gdf.crs.axis_info[0].unit_name == "metre") & (object_gdf.crs.axis_info[0].unit_name == "metre"), "make sure the given geodataframes crs is in meteres"
    #assert not any(building_gdf.geometry.is_empty), f"{aux_retrieve_name(building_gdf)} has an empty geometry"
    #assert not any(object_gdf.geometry.is_empty), f"{aux_retrieve_name(object_gdf)} has an empty geometry " 

    full_column_name = f"{column_name}_{'-'.join(map(str, radius_lst))}"
    #convert to meters
    #radius_lst = [ConvertToMeters(type = "LINEAR", value = r, unit = "FEET_US") for r in radius_lst]

    #set our query polygons
    diff_envelope =  building_gdf.geometry.envelope.map(lambda  geom : geom.buffer(max(radius_lst) , join_style = 2, cap_style = 3).difference(geom.buffer(min(radius_lst), join_style = 2, cap_style = 3)))
    
    #query polygons intersection with our object_gdf
    diff_envelope.reset_index(inplace = True, drop = True)
    object_gdf.reset_index(inplace = True, drop = True)

    objects_results, R_results = diff_envelope.sindex.query_bulk(object_gdf.geometry, predicate='intersects')

    # zip the index of each querry polygon to the object_gdf polygons that intersect with it 
    indx_dict = _zip_components(R_results, objects_results)
    diff_envelope = diff_envelope.rename("geometry").pipe(gpd.GeoDataFrame)

      #initialize a column in the building_gdf to contain the functions output
    #building_gdf[full_column_name] = diff_envelope.apply(lambda row: apply_sum_of_intersection(row, indx_dict) , axis= 1)

    building_gdf[full_column_name] = diff_envelope.apply(lambda row: int(object_gdf.geometry[indx_dict[row.name]].intersection(row["geometry"]).dropna().cascaded_union.area / row["geometry"].area *100) if len(indx_dict[row.name]) > 0 else 0 , axis= 1)

    #building_gdf[full_column_name] = diff_envelope.apply(lambda row: int(sumu(GeoSeries((row["geometry"]).difference(object_gdf.geometry[indx_dict[row.name]].intersection(row["geometry"]).dropna())).area) / row["geometry"].area *100) if len(indx_dict[row.name]) > 0 else 0 , axis= 1)
    #print(building_gdf.loc[building_gdf[full_column_name] > 100, full_column_name])
    #return building_gdf with relevant data
    return building_gdf
def percentage_in_buffer(building_gdf:gpd.GeoDataFrame, object_gdf:gpd.GeoDataFrame, radius_lst:list, column_name: str) ->gpd.GeoDataFrame:
    """
    calculates the percentage of overlap between one geodataframe and another by value in a specified range. if only an outer range is required than set one of the radius_lst values as 1

    Args:
        building_gdf (GeoDataFrame): a geodataframe containing footprints 
        object_gdf (GeoDataFrame):  a geodataframe containing objects to be searched (right now its trees) 
        radius_lst (list): list of radius to be calculated. Defaults to empty list.
        column_name (str): tag for the name of the column based on object type

    Returns:
        GeoDataFrame: added intersect_percent_at_radius_lst_m column to building_gdf
    """

    #assert len(radius_lst) == 2, f"please enter two values into the search radius {radius_lst} argument" 
    #assert building_gdf.crs.equals(object_gdf.crs), "make sure the given geodataframes have the same crs"
    #assert (building_gdf.crs.axis_info[0].unit_name == "metre") & (object_gdf.crs.axis_info[0].unit_name == "metre"), "make sure the given geodataframes crs is in meteres"
    #assert not any(building_gdf.geometry.is_empty), f"{aux_retrieve_name(building_gdf)} has an empty geometry"
    #assert not any(object_gdf.geometry.is_empty), f"{aux_retrieve_name(object_gdf)} has an empty geometry " 

    full_column_name = f"{column_name}_{'-'.join(map(str, radius_lst))}"
    #convert to meters
    radius_lst = [ConvertToMeters(type = "LINEAR", value = r, unit = "FEET_US") for r in radius_lst]

    #set our query polygons
    diff_envelope =  building_gdf.geometry.envelope.map(lambda  geom : geom.buffer(max(radius_lst) , join_style = 2, cap_style = 3).difference(geom.buffer(min(radius_lst), join_style = 2, cap_style = 3)))
    
    #query polygons intersection with our object_gdf
    diff_envelope.reset_index(inplace = True, drop = True)
    object_gdf.reset_index(inplace = True, drop = True)
    objects_results, R_results = diff_envelope.sindex.query_bulk(object_gdf.geometry, predicate='intersects')

    # zip the index of each querry polygon to the object_gdf polygons that intersect with it 
    indx_dict = _zip_components(R_results, objects_results)
    diff_envelope = diff_envelope.rename("geometry").pipe(gpd.GeoDataFrame)

    building_gdf[full_column_name] = diff_envelope.apply(lambda row: int(object_gdf.geometry[indx_dict[row.name]].intersection(row["geometry"]).dropna().cascaded_union.area / row["geometry"].area *100) if len(indx_dict[row.name]) > 0 else 0 , axis= 1)
    #return building_gdf with relevant data
    return building_gdf

#def main():
building_gdf = read_file("C:/Users/matan/OneDrive/Documents/percent_in_donut/alabama_0_footprint_df/alabama_0_footprint_df.shp")
object_gdf = read_file("C:/Users/matan/OneDrive/Documents/percent_in_donut/alabama_0_ai_df/alabama_0_ai_df.shp")
radius_lst = [5,20]
#%load_ext line_profiler
#%lprun -f percentage_in_buffer 
percentage_in_buffer(building_gdf, object_gdf, radius_lst, column_name = "veg").to_file("veg.shp")
#if __name__ == "__main__": main()