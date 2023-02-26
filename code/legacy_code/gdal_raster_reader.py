from osgeo import gdal, gdal_array
import numpy as np
import rasterio as rio
import rasterio.mask
import geopandas as gpd
import shapely
import numpy as np
import fiona
from tqdm import tqdm
import warnings
warnings.filterwarnings('ignore')
from numba import jit
import cv2

class ReadRaster:
    #TODO: add contextmanager 
    def __init__(self,src_img , func=None, dst_img = None, obj = None, with_debug = True, **func_kwargs):
        self.src_img = src_img
        self.dst_img = dst_img
        self.func_kwargs = func_kwargs
        self.func = func
        self.obj = obj
        if with_debug is True:
            self._add_osgeo_debug()
        self._setup_gdal_env()

    def _add_osgeo_debug(self):
        global logging
        import logging
        logging.basicConfig(level=logging.DEBUG)
    
    def _setup_gdal_env(self):
        # tell gdal to use its inbuilt exceptions
        gdal.UseExceptions()
        #set the gdal block ram cache max to 2GB
        gdal.SetCacheMax(2^60)
        # register all of the drivers
        gdal.AllRegister()

    def _instance_dst_ds(self):
        driver = gdal.GetDriverByName('GTiff')
        self.dst_ds = driver.Create(self.dst_img , self.xsize , self.ysize, self.in_ds.RasterCount , gdal.GDT_Byte )
        self.dst_ds.SetGeoTransform(self.in_ds.GetGeoTransform())
        self.dst_ds.SetProjection(self.in_ds.GetProjection())

    def _write_to_dst_ds(self, out_data:np.ndarray, band_num:int, row:int, column:int):
        #dst_ds.GetRasterBand(band_num).WriteArray(out_data, row, column)
        gdal_array.BandWriteArray(self.dst_ds.GetRasterBand(band_num),out_data, row, column)

    def _open_in_ds(self):
        self.in_ds = gdal.Open(self.src_img, gdal.GA_ReadOnly)
        # Get total raster XY size

    def coords_to_pixel(self, lat:float or np.ndarray, lon:float or np.ndarray) ->tuple:
        """
        _summary_

        Parameters
        ----------
        gt : tuple
            _description_
        lat : float
            _description_
        lon : float
            _description_

        Returns
        -------
        _type_
            _description_
        """    
        transform = self.in_ds.GetGeoTransform()
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = -transform[5]
        if not isinstance(lat, np.ndarray):
            # get x,y by compute pixel offset
            xOffset = int((lat - xOrigin) / pixelWidth)
            yOffset = int((yOrigin - lon) / pixelHeight)
            return xOffset, yOffset

        # get x,y by compute pixel offset
        xOffset = np.array([*map(int,np.floor((lat - xOrigin) / pixelWidth))])
        yOffset = np.array([*map(int, np.floor((yOrigin - lon ) / pixelHeight))])
        return xOffset, yOffset 

    def pixel_to_coords(self, pixel_x, pixel_y):
        xoffset, px_w, rot1, yoffset, px_h, rot2 = self.in_ds.GetGeoTransform()
        # supposing x and y are your pixel coordinate this 
        # is how to get the coordinate in space.
        posX = px_w * pixel_x + rot1 * pixel_y + xoffset
        posY = rot2 * pixel_x + px_h * pixel_y + yoffset

        # shift to the center of the pixel
        posX += px_w / 2.0
        posY += px_h / 2.0

        return posX, posY

    def read_block_from_pixel(self, x:float, y:float, is_coord:bool = False, block_shape:list = [256, 256],overlap:int = 0):
        if is_coord is False:
            if overlap == 0 :
                return self.in_ds.ReadAsArray(x ,y ,block_shape[0] ,block_shape[1] )
            return self.in_ds.ReadAsArray(xOffset, yOffset, block_shape[0], block_shape[1], overlap, overlap)

        xOffset, yOffset = self.coords_to_pixel(x,y)
        if overlap == 0 :
                return self.in_ds.ReadAsArray(x ,y ,block_shape[0] ,block_shape[1] )
        return self.in_ds.ReadAsArray(xOffset, yOffset, block_shape[0], block_shape[1], overlap, overlap)
    
    def process_block_bands_descretely(self, x, y, cols, rows, overlap):
        for band_idx in range(1,self.in_ds.RasterCount + 1):
            array = self.in_ds.ReadAsArray(band_idx ,x, y, cols, rows, overlap, overlap)
            self.process(array)

    def process_block_bands_together(self,x, y, cols, rows, overlap):
        return self.process(self.in_ds.ReadAsArray(x, y, cols, rows, overlap, overlap))
    
    def process_whole_raster_bands_descretely(self):
        for band_idx in range(1,self.in_ds.RasterCount + 1):
            array = self.in_ds.ReadAsArray(band_idx)
            self.process(array)

    def process_whole_raster_bands_together(self):
        return self.process(self.in_ds.ReadAsArray())
    
    def process(self ,array ,func ,**kwargs ):
        return func(array,**kwargs)
        
    # Function to read the raster as arrays of size 256
    def gdal_block_reader(self,window_size:list = [256, 256] , overlap:int = 1, windowing_type:str or None = "windows", descreete:bool = True):

        if windowing_type == "windows":
            x_window_size, y_window_size = window_size[0], window_size[1]
        
        if (windowing_type == "blocks") | (windowing_type == "tiles"):
            block_sizes = self.in_ds.GetRasterBand(1)
            x_window_size, y_window_size = block_sizes[0], block_sizes[1]

        for y in range(0, self.ysize, y_window_size):
            if y + y_window_size < self.ysize:
                rows = y_window_size
            else:
                rows = self.ysize - y
            for x in range(0, self.ysize, x_window_size):
                if x + x_window_size < self.ysize:
                    cols = x_window_size
                else:
                    cols = self.ysize - x

                if descreete is True:
                    self.process_block_bands_descretely(x, y, cols, rows, overlap)
                else:
                    self.process_block_bands_together(x, y, cols, rows, overlap)
    
    def gdal_read_simple(self, descreete:bool = True):
        if descreete is True:
                    self.process_whole_raster_bands_descretely()
        else:
            self.process_whole_raster_bands_together()
    @jit
    def _generate_thresh(mask:np.ndarray) -> np.ndarray:
        """
        Function converting binary bask values  to 127 and 255 for easier contour detection

        Args:
            mask (ndarray): binary mask

        Returns:
            ndarray: threshold array from binary mask
        """        
        return cv2.threshold(mask ,127,255,cv2.THRESH_BINARY+cv2.THRESH_OTSU)[1]
    
    
    @jit(forceobj=True)
    def _greyscale(self):
        """
        Returns:
            grayscale of orginal image array
        """
        self.img  = cv2.cvtColor(self.img, cv2.COLOR_BGR2GRAY)

    @jit(forceobj=True)
    def mask_to_polygons_layer(self, mask, out_transform, boundary = None):
        """_summary_

        Args:
            mask (_type_): _description_
            out_transform (_type_): _description_
            boundary (_type_, optional): _description_. Defaults to None.

        Returns:
            _type_: _description_
        
        all_polygons = []
        for shape, value in rio.mask.astype(np.uint8):
             (mask >0)
             transform = out_transform:
            #return shapely.geometry.shape(shape)
            all_polygons.append(shapely.geometry.shape(shape))
        all_polygons = shapely.geometry.MultiPolygon(all_polygons)
        if not all_polygons.is_valid:
            all_polygons = all_polygons.buffer(0)
            # Sometimes buffer() converts a simple Multipolygon to just a Polygon,
            # need to keep it a Multi throughout
            if all_polygons.type == 'Polygon':
                all_polygons = shapely.geometry.MultiPolygon([all_polygons])
        if boundary is not None:
            gdf = gpd.GeoDataFrame({'geometry'(boundary))"""
        return gpd.GeoDataFrame({'geometry':all_polygons})

    def _get_polygons_from_shp(self, shapefile:str):
        with fiona.Env(CPL_DEBUG=True):
            with fiona.open(shapefile, "r") as shapefile:
                shapes = [feature["geometry"] for feature in shapefile]
        return shapes
    
    def _read_and_process(self, src, band_id:int or None = None, window = None):
            src_data = src.read(band_id, window = window)
            return self._run_func(src_data)# Do the Processing Here
      
    def _write(self,dst, dst_data, band_id = None, window = None):
        if band_id is None:
            if window is None:
                dst.write(dst_data)
            else:
                dst.write(dst_data, window = window)
        if window is None:
                dst.write(dst_data)
        else:
             dst.write(dst_data, band_id, window = window)

    def _process_subsample_and_export_as_shape(self,src , gdf, mask, out_transform, shape):
            mask = mask.astype(np.uint8)
            gdf.set_crs(src.crs)
            gdf = gdf.append(self.mask_to_polygons_layer(mask, out_transform, shape)).reset_index(drop = True)
            gdf.to_file(self.dst_img)

    def _get_polygons_from_shp(self, shapefile:str):
            with fiona.Env(CPL_DEBUG=True):
                with fiona.open(shapefile, "r") as shapefile:
                    shapes = [feature["geometry"] for feature in shapefile]
            return shapes

    def subsample_tiff_by_shape(self,src, shapes):
        for shape in tqdm(shapes):
                out_image, out_transform = rio.mask.mask(src, [shape], crop=True)
        return out_image, out_transform