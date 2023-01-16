from geopandas import GeoDataFrame, sjoin
from shapely.geometry import MultiPolygon
from string import ascii_uppercase, digits
from random import choices

def geo_groupby(gdf , element_id , buffer_m = 0, return_simplified_dataframe = False ):
    """
    group intersecting goemetrys in geodataframe and give them a gorup ID

    Args:
        gdf (geodataframe): a geodataframe containing overlaping polygons in the geometry column
        element_id (string): (column name) a single single value given to each polygon row
        buffer_m (int): the buffer to give the elements in the array in meters
        return_simplified_dataframe (bool, optional): :sets the export value data type

    Returns:
        geodataframe: function returns a geodataframe , either simplified by group id or with added group id and group geometry values

    """
    #TODO: optimize unary_union and sjoin with sindex to save RAM 
    #notes: spatial indeX INTERSECTs
    #TODO: write out tests!!!!
    crs = gdf.crs
    if buffer_m > 0: 
        if crs is None:
            raise ValueError("Please set a CRS for the input geodataframe so that the function can reproject back to it after adding the buffer")
        
        if crs != 3857:
            gdf = gdf.to_crs("EPSG:3857")
            gdf  = gdf.geometry.buffer(buffer_m)
            gdf = gdf.to_crs(crs)
        else:
            gdf  = gdf.geometry.buffer(buffer_m)

    union = gdf.unary_union
    if union.type == "Polygon":
        data = {"Group_ID":"".join(choices(ascii_uppercase + digits ,k = 30)), "geometry":union}
        gdf_union = GeoDataFrame(data, geometry= "geometry", index = [0], crs = crs)

    else:
        data = {"Group_ID": ["".join(choices(ascii_uppercase + digits ,k = 30)) for _ in union] , "geometry":[x for x in union]}
        gdf_union = GeoDataFrame(data, geometry= "geometry" , crs = crs)
    
    gdf_union = sjoin(gdf_union, GeoDataFrame(gdf[[element_id, "geometry"]]), how= "left", op = "intersects").drop(columns = ["index_right"])
    group_gdf = gdf.merge(gdf_union, on = element_id)
    
    if return_simplified_dataframe is True:
        merged_df = group_gdf.groupby("Group_ID").agg(lambda x : list(x))
        
        for i in merged_df.index:
            merged_df["geometry"][i] = MultiPolygon(merged_df["geometry_y"][i])
        merged_gdf = GeoDataFrame(merged_df ,crs = crs )
        merged_gdf = merged_gdf.rename(columns = {"geometry_y": "Group_Geometry","geometry_x":"element_geometry"}).reset_index(drop = True)
        return merged_gdf
        
    group_gdf = group_gdf.rename(columns = {"geometry_y": "Group_Geometry","geometry_x":"element_geometry"})
    group_gdf = group_gdf.sort_values("Group_ID")
    return group_gdf

