from geopandas import read_file
from warnings import filterwarnings
filterwarnings("ignore")
from percentage_in_buffers import percentage_in_donut
#def main():
building_gdf = read_file("C:/Users/matan/OneDrive/Documents/percent_in_donut/alabama_0_footprint_df/alabama_0_footprint_df.shp")
object_gdf = read_file("C:/Users/matan/OneDrive/Documents/percent_in_donut/alabama_0_ai_df/alabama_0_ai_df.shp")
radius_lst = [5,20]
percentage_in_donut(building_gdf, object_gdf, radius_lst)

#if __name__ == "__main__": main()
