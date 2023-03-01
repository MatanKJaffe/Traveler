from Anomaly_Index import AnomalyIndex
from Difference_of_Gaussians import DoG
from numpy import gradient, ndarray, uint8 ,pi ,arctan ,arctan2 , sin, cos, sqrt, uint16, degrees, array, float32 as nparray
from osgeo import gdal, gdalconst, gdal_array

gdal.SetCacheMax(1000000000)
import os
import numpy as np
from tqdm import tqdm
import numpy as np
from osgeo import gdal, gdal_array, ogr, osr
from tqdm import tqdm
gdal.AllRegister()
from math import atan
from scipy.ndimage import convolve
class VrtBuilder:
    def __init__(self, vrt_path = "out.vrt", extension='.tif', epsg = 3857):
        self.path_list = []
        self.vrt_path = vrt_path
        self.extension = extension
        self.epsg = epsg
        
    def _add_path(self,path):
        self.path_list.append(path)

    def set_filepaths(self):
        filepaths = []
        for path in self.path_list:
            filepaths.append(self.get_all_tifs(path))
            
        self.filepaths = filepaths

    def build_vrt_file(self):
        filepaths = self.get_all_tifs()# list of paths to raster files
        vrt_options = gdal.BuildVRTOptions(resampleAlg='cubic', addAlpha=True)
        vrt_dataset = gdal.BuildVRT(self.vrt_path, filepaths, options=vrt_options)
        vrt_dataset = self.normalize_projections(vrt_dataset)
        vrt_dataset.FlushCache()
        vrt_dataset = None
        
    def normalize_projections(self,from_dataset):
        """
        Returns the input dataset in the expected "destination" SRS.
        If the dataset is already in the correct SRS, returns it unmodified
        """
        from_srs = from_dataset.GetSpatialRef()
        to_srs = osr.SpatialReference()
        to_srs.ImportFromEPSG(self.epsg)
        out_vrt = f'{os.path.splitext(self.vrt_path)[0]}_{self.epsg}.tif'
        
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

    def get_all_tifs(self, path):
        output = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(self.extension):
                    output.append(os.path.join(root, file))
        return output

class RasterProcessor(VrtBuilder):
    def __init__(self, input_file, func, output_file = None, block_size = 256, description = "Elevation", name = "elevation"):
        self.input_file = input_file
        self.block_size = block_size
        self.func = func
        
        self.input_dataset = gdal.Open(self.input_file, gdal.GA_ReadOnly)
        self.num_bands = self.input_dataset.RasterCount
        
        self.xsize = self.input_dataset.RasterXSize
        self.ysize = self.input_dataset.RasterYSize
        
        self.resolution = self.input_dataset.GetGeoTransform()[1]
        
        driver = gdal.GetDriverByName('GTiff')
        if output_file is None:
            output_file = f'{os.path.splitext(self.input_file)[0]}_{name}.tif'
        self.output_dataset = driver.Create(output_file, self.xsize, self.ysize, 1, gdal.GDT_Float32)
        self.output_dataset.SetProjection(self.input_dataset.GetProjection())
        self.output_dataset.SetGeoTransform(self.input_dataset.GetGeoTransform())
        self.output_band = self.output_dataset.GetRasterBand(1)
        self.output_band.SetNoDataValue(-9999)
        self.output_band.SetDescription(description)
        
        
                
    def process_block(self, x, y, cols, rows):
        band = self.input_dataset.GetRasterBand(1)
        data = band.ReadAsArray(x, y, cols, rows)
        if data is not None:
            processed_data = self.func(data)
            output_band = self.output_band
            output_band.WriteArray(processed_data, x, y)
    
    def process_band(self):
        band = self.input_dataset.GetRasterBand(1)
        x_block_size, y_block_size = band.GetBlockSize()
        xsize = band.XSize
        ysize = band.YSize
        blocks = 0
        for y in range(0, ysize, y_block_size):
            #print blocks
            if y + y_block_size < ysize:
                rows = y_block_size
            else:
                rows = ysize - y
                
            for x in range(0, xsize, x_block_size):
                if x + x_block_size < xsize:
                    cols = x_block_size
                else:
                    cols = xsize - x
                self.process_block(x, y, cols, rows)
                #print array.shape
                blocks += 1
        
    def process(self):
        self.process_band()
        self.output_dataset.FlushCache()
        self.output_band.FlushCache()
        self.input_dataset.FlushCache()
        self.input_dataset = None
        self.output_dataset = None
        self.output_band = None

class GdalGeomorphology(RasterProcessor):
    def __init__(self):
        pass
    def gdal_slope(self, outfile,dsc):
        # Compute the slope using GDAL
        gdal.DEMProcessing(outfile, dsc, 'slope', computeEdges=True, slopeFormat='degree', scale=111120, alg='Horn')

    def gdal_aspect(self, outfile, dsc):
        # Compute the aspect using GDAL
        gdal.DEMProcessing(outfile, dsc, 'aspect', computeEdges=True, alg='ZevenbergenThorne')

    def gdal_tri_wilson(self, outfile, dsc):
        gdal.DEMProcessing(outfile, dsc, 'tri', computeEdges=True, alg='Wilson')

    def gdal_tri_riley(self, outfile, dsc):
        gdal.DEMProcessing(outfile, dsc, 'tri', computeEdges=True, alg='Riley')
    
    def gdal_tpi(self, outfile, dsc):
        # Compute the TPI using GDAL
        gdal.DEMProcessing(outfile, dsc, 'tpi', computeEdges=True)

    def gdal_roughness(self, outfile, dsc):
        # Compute the roughness using GDAL
        gdal.DEMProcessing(outfile, dsc, 'roughness', computeEdges=True)
        
#TODO: simplify the NpGeomorphology class
class NpGeomorphology(RasterProcessor):
    def __init__(self):
        pass
    
    def rmss(data):
        """Root Mean Square Slope (RMSS)"""
        dy, dx = np.gradient(data)
        slope = np.sqrt(pi/2. - arctan(sqrt(dx*dx + dy*dy)))
        return slope**2
    
    def msr(data):
        scales = [1, 2, 4, 8, 16, 32] # Define scales to compute roughness at
        roughness = np.zeros(data.shape)
        for s in scales:
            dx, dy = np.gradient(data, s, s)
            dzdx, dzdy = np.gradient(data, s, s)
            r = np.sqrt(dx**2 + dy**2 + dzdx**2 + dzdy**2)
            roughness += r
        return roughness/len(scales)

    def vrm(data):
        """Vector Ruggedness Measure (VRM)"""
        dy, dx = np.gradient(data)
        dxy, _ = np.gradient(dx)
        _, dyy = np.gradient(dy)
        return np.sqrt(dxy**2 + dyy**2).mean()
    
    def np_aspect(self, data):
        x, y = gradient(data)
        return np.array(degrees(arctan2(-x, y)), dtype = uint16)

    def np_slope(self, data):
        x, y = gradient(data)
        return  np.array(degrees(pi/2. - arctan(sqrt(x*x + y*y))), dtype = uint16)
    
    def aspect_cosine(self, data):
        dx, dy = np.gradient(data)
        return np.cos(np.arctan2(dy, -dx))
    
    def aspect_sine(self,data):
        dx, dy = np.gradient(data)
        return np.sin(np.arctan2(dy, -dx))

    def eastness(self, data):
        dx, dy = np.gradient(data)
        return np.arctan2(dy, dx)

    def northness(self,data):
        dx, dy = np.gradient(data)
        return np.arctan2(-dx, dy)
    
    def convergence(self, data):
        dx, dy = np.gradient(data)
        return np.arctan2(-dx, -dy)

#TODO: add more functions to the Hydromorphology class
#TODO: build soil comp class
#TODO: build land cover class 
#TODO: build expected meteorological conditions class
#TODO: build meteorology topography class

class Hydromorphology(RasterProcessor):
    def __init__(self,):
        pass
    
    def sca(self,data):
        data = data.astype(np.float32)
        dx, dy = np.gradient(data)
        slope = np.arctan(np.sqrt(dx * dx + dy * dy))
        aspect = np.arctan2(dy, -dx)

        # Compute the specific catchment area (SCA)
        sca = np.empty_like(data)
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                sca[y, x] = np.sum(np.where((slope > 0) & (aspect >= 0) & (aspect < np.pi / 2), np.sin(slope[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= np.pi / 2) & (aspect < np.pi), np.cos(np.pi - aspect[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= -np.pi) & (aspect < -np.pi / 2), np.sin(slope[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= -np.pi / 2) & (aspect < 0), np.cos(aspect[y, x]), 0))
        return sca.astype(np.float32)


    def twi(self,data):
        # Convert the elevation data to float32 and compute slope and aspect
        data = data.astype(np.float32)
        dx, dy = np.gradient(data)
        slope = np.arctan(np.sqrt(dx * dx + dy * dy))
        aspect = np.arctan2(dy, -dx)

        # Compute the specific catchment area (SCA)
        sca = np.empty_like(data)
        for y in range(data.shape[0]):
            for x in range(data.shape[1]):
                sca[y, x] = np.sum(np.where((slope > 0) & (aspect >= 0) & (aspect < np.pi / 2), np.sin(slope[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= np.pi / 2) & (aspect < np.pi), np.cos(np.pi - aspect[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= -np.pi) & (aspect < -np.pi / 2), np.sin(slope[y, x]), 0))
                sca[y, x] += np.sum(np.where((slope > 0) & (aspect >= -np.pi / 2) & (aspect < 0), np.cos(aspect[y, x]), 0))

        # Compute the topographic wetness index (TWI)
        twi = np.log10(sca / np.sin(slope))
        twi[np.isinf(twi)] = -9999.0

        return twi.astype(np.float32)
     

class LandFormDelineation(AnomalyIndex, DoG):
    def __init__(self):
        pass
    
class Geomorphology(NpGeomorphology, GdalGeomorphology):
    def __init__(self):
        pass
    
if __name__ == '__main__':
    path= "/Users/osn/Code_Library/Traveler/code/geomorpho/SRTM_GL3/SRTM_GL3_srtm/North/North_30_60"
    extension='.tif'
    
    def get_all_tifs(path,extension):
        output = []
        for root, dirs, files in os.walk(path):
            for file in files:
                if file.endswith(extension):
                    output.append(os.path.join(root, file))
        return output

    tiff_list = get_all_tifs(path,extension)
    processor = RasterProcessor(input_file=tiff_list[0],
                                block_size=256, func=  ,output_file= "out.tif", description = "Slope", name ="Slope").process()

    # processor = RasterProcessor(input_file=tif,
    #                             block_size=256, func=aspect, description = "Aspect", name ="aspect").process()
    # processor = RasterProcessor(input_file=tif,
    #                             block_size=256, func=AnomalyIndex().execute, description = "TPI_Landforms", name ="tpi_L").process()
        