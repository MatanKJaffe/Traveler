"#!/usr/bin/env python"
#-*- coding: utf-8- -*-
# from osgeo import gdal, osr
# import os
import geopandas as gpd
import rasterio
import numpy as np
from shapely.wkt import loads
from pandas import to_numeric
from geopandas import GeoDataFrame, GeoSeries
from random import choices
from string import ascii_lowercase, digits
from math import atan, degrees
from itertools import combinations
from shapely.geometry import Point, LineString
from warnings import filterwarnings
from haversine import haversine, Unit
from rasterio.sample import sample_gen

import rasterio

from rasterio.vrt import WarpedVRT
from rasterio.mask import mask
from shapely.geometry import LineString, mapping, Point, box

filterwarnings("ignore")

# group lines by relative angle
class Trail(object):
    """
    Class object for runing Straightest_Lines grouping algorithm
    input:
        lines ; list of linestrings as strings.
        segment_lines ; bool value for whether or not to split the linestrings into component
        straight lines.

    output:
        lines ; GeoSeries of shapely linestring geometrys
        nodes ; GeoSeries shapely point geometrys representing intersections between at least 2 linestrings
        gdf ; GeoDataFrame with geometry and unique label for generated straightest line groups.
    """

    def __init__(self, lines, raster_file, segment_lines=True, is_wkt = True):
        self.raster_file = raster_file
        # linestring as strings to shapely objects
        if segment_lines is True:
            if not is_wkt:
                self.lines = GeoSeries([segment for line in lines for segment in self.segment_linestring(loads(line))])

            self.lines = GeoSeries([segment for line in lines for segment in self.segment_linestring(line)])

        else:
            self.lines = GeoSeries([loads(line) for line in lines])
        # add lines to output gdf
        self.gdf = GeoDataFrame(geometry=self.lines)
        self.full_line = lines[0]
        self.get_line_data()
        
    def get_line_data(self):
        self.generate_nodes()
        self.calculate_angle_between_all_lines()
        self.line_length = self.gdf["geometry"].length
        self.get_line_lengths()
        self.calculate_angle_between_all_lines()
        self.get_raster_data()
        self.snap_raster_points_to_line()
        
    def segment_linestring(self, line):
        """
        Function for segmenting linestring  into individual straight lines
        """
        if len(line.coords) > 2:
            return list(map(LineString, zip(line.coords[:-1], line.coords[1:])))

        return [line]

    def slope(self, x1, y1, x2, y2):
        """
        Function for calculating Line slope given two points:
                 y2 - y1
        slope = -------
                 x2 - x1
        """
        if (x2 - x1) == 0:
            return 0
        
        return (y2 - y1) / (x2 - x1)

    def angle(self, slope1, slope2):
        """
        Function for calculating given two slopes:
               |  slope2 - slope1  |
        tan0 = |-------------------|
               |1 + (slope2*slope1)|
        *angle is in degrees
        """
        return degrees(atan(abs((slope2 - slope1) / (1 + (slope2 * slope1)))))

    def calculate_ang_between_2_lines(self, lineA, lineB):
        """
        Function: calculating the angle between two lines using the formula
               |  slope2 - slope1  |                                       y2 - y1
        tan0 = |-------------------| where slope is calculated as: slope = -------
               |1 + (slope2*slope1)|                                       x2 - x1
        *angle is in degrees
        """
        # set shapely linestrings (lineA,lineB) as lists of point coordinates from
        lineA = [i for i in lineA.coords]
        lineB = [i for i in lineB.coords]

        # calculate slope for each line
        slope1 = self.slope(lineA[0][0], lineA[0][1], lineA[1][0], lineA[1][1])
        slope2 = self.slope(lineB[0][0], lineB[0][1], lineB[1][0], lineB[1][1])

        # calculate and return angle between the two lines in degrees
        return self.angle(slope1, slope2)

    def generate_inter_pts(self, inter, inter_pts):
        """
        Function for finding intersection points between lines given the intersection geometry
        between those lines
        inter : intersection shapely geometry
        inter_pts : list of all intersection points between the lines
        """
        if "Point" == inter.type:
            inter_pts.append(inter)

        elif "MultiPoint" == inter.type:
            inter_pts.extend([pt for pt in inter])

        elif "MultiLineString" == inter.type:
            multiline = [line for line in inter]
            first_coords = multiline[0].coords[0]
            last_coords = multiline[len(multiline) - 1].coords[1]
            inter_pts.append(Point(first_coords[0], first_coords[1]))
            inter_pts.append(Point(last_coords[0], last_coords[1]))

        return inter_pts

    def generate_nodes(self):
        """
        Function generating intersection points between at least two geometrys (nodes) given
        a list of linestrings.
        """
        # get all end points for all lines
        endpts = [(Point(list(line.coords)[0]), Point(list(line.coords)[-1])) for line in self.lines]
        # flatten the resulting list to a simple list of points
        endpts = [pt for sublist in endpts for pt in sublist]
        # remove duplicates from endpoint list
        endpts = list(GeoSeries(endpts).unique())

        # get all generateable intersection points
        inter_pts = []
        for line1, line2 in combinations(self.lines, 2):
            if line1.intersects(line2):
                inter = line1.intersection(line2.buffer(2))
                inter_pts = self.generate_inter_pts(inter, inter_pts)
                if "GeometryCollection" == inter.type:
                    for geom in inter:
                        inter_pts = self.generate_inter_pts(geom, inter_pts)

        # add all points to the same list for iteration
        endpts.extend([pt for pt in inter_pts if pt not in endpts])

        # from list of intersection and endpoints find all nodes
        nodepts = []
        for node in endpts:
            if self.gdf.intersects(node).sum() > 1:
                nodepts.append(node)

        # save the nodes to the class object
        self.nodes = GeoSeries(GeoSeries(nodepts).unique())


    def calculate_angle_between_all_lines(self):
        """
        Function calculating straightest continuous lines in a list of  intersecting lines
        """
        angels = []
        
        # subsample geodataframe by intersection with each node
        for node in self.nodes:
            node_gdf = self.gdf[(self.gdf.geometry.intersects(node.buffer(5)))]
            node_gdf.reset_index(drop= True, inplace= True)
            # if the node only intersects two geometrys they are obviusly the straightest line
            if len(node_gdf) == 2:
                angle = self.calculate_ang_between_2_lines(node_gdf["geometry"][1], node_gdf["geometry"][0])
                angels.append(angle)
                
        self.angels = gpd.GeoDataFrame()
        self.angels["angle"] = np.ceil(np.array(angels))
        self.angels["angle"] = 180 - self.angels["angle"]


        
    def get_raster_data_of_nodes(self, vrt, points, out_var):
        """
        Get the elevation profile (distance vs. altitude) of a path segment from the list of coordinates.
        Args:
            dataset: The opened rasterio dataset for the topography data.

        Returns: The distance (in meters) and elevation (in meters) vectors.
        """
        # coordinates are [lon, lat], flip for rasterio
        coords = [[c.x, c.y] for c in points]
        # convert meters to feet and use rasterio.sample.sample_gen to query each point
        out_var = [e[0] for e in sample_gen(vrt, coords)]

    def get_line_lengths(self):
        coords = [[c.x, c.y] for c in self.nodes.set_crs(3857).to_crs(4326)]
        line_length = [0.0]
        line_length.extend([haversine(coords[i], coords[i+1], Unit.METERS) for i in range(len(coords) - 1)])
        self.line_length_haversine = line_length

    def pixel_to_coords(self, pixel_x, pixel_y, transform):
        #decompose geotransform
        xoffset, px_w, rot1, yoffset, px_h, rot2 = transform
        # supposing x and y are your pixel coordinate this 
        # is how to get the coordinate in space.
        posX = px_w * pixel_x + rot1 * pixel_y + xoffset
        posY = rot2 * pixel_x + px_h * pixel_y + yoffset

        # shift to the center of the pixel
        posX += px_w / 2.0
        posY += px_h / 2.0

        return Point(posX, posY)

    def get_raster_data(self):
        # Open raster file
        with rasterio.open(self.raster_file) as src:
            # Open VRT file as a WarpedVRT object
            with WarpedVRT(src) as vrt:
                # Get bounds of raser object
                self.raster_bounds = vrt.bounds
                self.get_raster_points_intersecting_linestring(self.full_line, vrt)
                
                # coordinates are [lon, lat], flip for rasterio
                self.nodes_elev = None
                self.raster_points_elev =  None
                self.nodes_elev = self.get_raster_data_of_nodes(vrt, self.nodes, )
                self.get_raster_data_of_nodes(vrt, self.raster_points, self.raster_points_elev)
                
    def get_raster_points_intersecting_linestring(self, linestring, dataset):
        # Get bounds of linestring
        bounds = linestring.bounds
        # Mask raster with linestring geometry
        out_image, out_transform = mask(dataset, [mapping(linestring)], crop=True, all_touched=True)
        #get coordinates of intersecting raster cells
        coords = np.argwhere(out_image != -9999)[:,[1,2]]
        self.raster_xy_coords = coords
        # Get coordinates of pixels with value > 0
        x_res, y_res = out_transform.a, out_transform.e
        x_min, y_min, x_max, y_max = bounds
        #set geotransform      
        geotransform = (x_min, x_res, 0, y_max, y_res, 0)
        #generate point from array coordinates
        point_coords = [self.pixel_to_coords(x[1], x[0], geotransform) for x in coords] 
        self.raster_points = GeoSeries(point_coords)

    def snap_raster_points_to_line(self):
        self.snaped_raster_points = GeoSeries(linestring.interpolate(linestring.project(point)) for point in line_map.raster_points).set_crs(3857)


if __name__ == "__main__":
    shapes = gpd.read_file("/Users/osn/Code_Library/Traveler/code/website_mockup/hiker_site/static/edges/aosta_valley_italy.geojson")
    shapes.to_crs(3857, inplace= True)
    
    raster_file = "/Users/osn/Code_Library/Traveler/data/italy_dem/italy_dem_3857.vrt"
    
    test_line = GeoSeries([shapes.loc[[shapes["length"].idxmax()]]["geometry"].values[0]])
    linestring = test_line[0]
    
    line_map = Trail(test_line, raster_file, segment_lines= True)
