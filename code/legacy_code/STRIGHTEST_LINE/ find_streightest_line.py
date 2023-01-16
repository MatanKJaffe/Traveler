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
filterwarnings("ignore")

lines = ['LINESTRING (-10797337.23321533 3873376.289489746, -10797338.03521729 3873369.30078125)',
'LINESTRING (-10797338.03521729 3873369.30078125, -10797344.22186279 3873376.174865723)',
'LINESTRING (-10797338.03521729 3873369.30078125, -10797349.60656738 3873365.63470459)',
'LINESTRING (-10797335.74387592 3873335.732715198, -10797342.27421039 3873343.523120471)',
'LINESTRING (-10797342.27421039 3873343.523120471, -10797353.84556049 3873339.857043811)',
'LINESTRING (-10797344.90928119 3873326.223437854, -10797342.27421039 3873343.523120471)',
'LINESTRING (-10797342.27421039 3873343.523120471, -10797348.4608559 3873350.397204944)',
'LINESTRING (-10797342.27423096 3873343.523132324, -10797338.03521729 3873369.30078125)',
'LINESTRING (-10797326.92220234 3873337.336473834, -10797324.28713154 3873354.636156451)',
'LINESTRING (-10797317.75679706 3873346.845751178, -10797324.28713154 3873354.636156451)',
'LINESTRING (-10797324.28713154 3873354.636156451, -10797335.85848163 3873350.970079791)',
'LINESTRING (-10797288.31299486 3873332.753778585, -10797285.67792406 3873350.053461202)',
'LINESTRING (-10797285.67792406 3873350.053461202, -10797297.24927416 3873346.387384542)',
'LINESTRING (-10797299.42600981 3873364.718086001, -10797310.99735991 3873361.05200934)',
'LINESTRING (-10797298.62400785 3873371.706794497, -10797299.42600981 3873364.718086001)',
'LINESTRING (-10797285.67794462 3873350.053473055, -10797281.43893095 3873375.831121981)',
'LINESTRING (-10797306.30007371 3873321.640742605, -10797303.66500291 3873338.940425222)',
'LINESTRING (-10797303.66500291 3873338.940425222, -10797309.85164842 3873345.814509694)',
'LINESTRING (-10797303.66500291 3873338.940425222, -10797315.23635301 3873335.274348562)',
'LINESTRING (-10797279.14758959 3873342.263055929, -10797285.67792406 3873350.053461202)',
'LINESTRING (-10797299.42600981 3873364.718086001, -10797305.61265532 3873371.592170473)',
'LINESTRING (-10797297.13466844 3873331.150019948, -10797303.66500291 3873338.940425222)',
'LINESTRING (-10797303.66502348 3873338.940437075, -10797301.58660621 3873351.579393888)',
'LINESTRING (-10797301.58660621 3873351.579393888, -10797299.42600981 3873364.718086001)',
'LINESTRING (-10797324.28710938 3873354.636230469, -10797332.53413939 3873361.277603304, -10797338.03521729 3873369.300781248)',
'LINESTRING (-10797324.2871521 3873354.636168303, -10797324.0795451 3873366.992879203, -10797320.4043048 3873379.295513599)',
'LINESTRING (-10797285.67790189 3873350.053535215, -10797291.82450636 3873358.430903191, -10797299.42600981 3873364.718086012)']


# group lines by relative angle
class StraightestLines:
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

    def __init__(self, lines, segment_lines=True):
        # linestring as strings to shapely objects
        if segment_lines is True:
            self.lines = GeoSeries([segment for line in lines for segment in self.segment_linestring(loads(line))])

        else:
            self.lines = GeoSeries([loads(line) for line in lines])
        # add lines to output gdf
        self.gdf = GeoDataFrame(geometry=self.lines)

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
        endptsd = list(GeoSeries(endpts).unique())

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

    def find_lowest_angle_line_index(self, line1, line1_indx, intersecting_lines):
        """
        Function to find the index of the line with the lowest angel relative to line1
        """
        for line2_indx, line2 in intersecting_lines.geometry.items():
            intersecting_lines["angle"][line2_indx] = self.calculate_ang_between_2_lines(line1, line2)
        # return the index of the lowest angle value (which isnt the angle between the line and itself)
        return to_numeric(intersecting_lines.loc[intersecting_lines.index != line1_indx, "angle"]).idxmin()

    def label_line_group(self, group_indexs):
        """Function give unique label to each straightest line group.
        group_indexs ; list of straightest line indexes
        """
        if any(self.gdf.loc[group_indexs, "labels"].notna()):
            uniques = list(self.gdf.loc[group_indexs, "labels"].dropna().unique())
            lst = self.gdf.loc[self.gdf["labels"].isin(uniques)].index.to_list()
            group_indexs.extend(lst)
            group_indexs = list(set(group_indexs))
        self.gdf.loc[group_indexs, "labels"] = "".join(choices(ascii_lowercase + digits, k=30))

    def calculate_straightest_lines(self):
        """
        Function calculating straightest continuous lines in a list of  intersecting lines
        """
        # instance labels column
        self.gdf["labels"] = None

        # subsample geodataframe by intersection with each node
        for node in self.nodes:
            node_gdf = self.gdf[(self.gdf.geometry.intersects(node.buffer(1)))]
            # if the node only intersects two geometrys they are obviusly the straightest line
            if len(node_gdf) == 2:
                group_indexs = node_gdf.index.to_list()
                self.label_line_group(group_indexs)
                continue

            # iterate over all lines in node_gdf
            for line1_indx, line1 in node_gdf.geometry.items():
                intersecting_lines = node_gdf[(node_gdf.geometry.intersects(line1.buffer(1)))]
                intersecting_lines["angle"] = None

                # get index of the line with the lowest angle relative to current line
                min_angle_idx = self.find_lowest_angle_line_index(line1, line1_indx, intersecting_lines)
                # get the index of the line with the lowest angle relative to line at min_angle_idx
                check_angle_idx = self.find_lowest_angle_line_index(self.gdf.geometry[min_angle_idx], min_angle_idx,
                                                                    intersecting_lines)

                # if the two are equal than they should be grouped into the same straightest line group
                if check_angle_idx == line1_indx:
                    group_indexs = [check_angle_idx, min_angle_idx]
                    self.label_line_group(group_indexs)

        # add label to all ungrouped lines
        for i in self.gdf[self.gdf["labels"].isna()].index:
            self.gdf["labels"][i] = "".join(choices(ascii_lowercase + digits, k=30))

    def fit(self):
        """
        Function running straightest line grouping algorithm
        """
        self.generate_nodes()
        self.calculate_straightest_lines()

    def visualize(self, with_nodes=True):
        """
        Function called to visualize output.
        with_nodes ; bool value whether to plot nodes
        """
        fig, ax = subplots()
        if with_nodes is True:
            self.nodes.plot(ax=ax, color="r", marker="o", markersize=30, zorder=2, edgecolor="k")

        self.gdf.plot(ax=ax, cmap="rainbow", column="labels", zorder=1)
        show()



def main():
    line_map = StraightestLines(lines, segment_lines= False)
    line_map.fit()
    line_map.visualize(with_nodes= False)

if __name__ == "__main__":
    main()