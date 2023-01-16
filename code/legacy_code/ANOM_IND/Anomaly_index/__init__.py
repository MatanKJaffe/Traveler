
import warnings
warnings.filterwarnings('ignore')
from numba import jit
from poly_from_contours import generate_anomaly_contours
@jit(forceobj=True)
def main():
    in_tiff = "Q:/Projects/roof_material_20220118/Miami_Part9/Ortho_Miami_Part1_ecw.ecw.tiff"
    roof_db_shapefile = "Z:/Projects/AI/Model - JM/Deliveries/Annada/Roof Material Data Set/Roof Material/Miami Part 9/Miami_Part_9_Ortho_Miami_Part_1_Buildings_Roof_Material/Miami_Part_9_Ortho_Miami_Part_1_Buildings_Roof_Material.shp"
    out_shapefile = "Z:/Product/Data Preparation/United States/Production for AI/Roof Condition Segmentation/full_data_set/stains_"
    generate_anomaly_contours(in_tiff, roof_db_shapefile,out_shapefile)
if __name__ == "__main__":
    main()