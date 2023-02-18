from numpy import column_stack
from numpy import full
from math import sqrt
from skimage.draw import circle_perimeter
from skimage.draw import line
from shapely.geometry import Point
#the viewshed function will work out the endpoint for each line of sight (each cell on the perimiter of the circle)
#and than call the line_of_sight() function to work out the visibility of each point on a line from the origin to the location of interest
class Visibility:
    def __init__(self, obserever_point, radius_m, observer_height):
        self.xObserver = obserever_point.x
        self.yObserver = obserever_point.y
        self.radius = radius_m
        self.observer_height = observer_height
        self.target_height = target_height
    def viewshed(self,xObserver, yObserver, radius, observer_height , target_height, dem_data, dem_array, transform):
        #get base data from raster in order to preform conversion calculations
        xOrigin = transform[0]
        yOrigin = transform[3]
        pixelWidth = transform[1]
        pixelHeight = transform[5]
        
        #make sure the observer point coordinates are in the raster
        if ((xObserver - xOrigin) < 0 ) | ((yOrigin - yObserver) < 0 ):
            print(f"Sorry , {xOrigin ,yOrigin} is not within the elevation dataset.")
            return
        
        #convert radius to pixels
        radius_px = int(radius / pixelWidth)
        
        #convert observer point coordinates to row and column indexes
        row = int((yOrigin - y0) / pixelHeight)
        col = int((x0 - xOrigin) / pixelWidth)

        #get observer height (above sea level)
        height0 = dem_array[(row,col)] + observer_height
        
        #define output array
        output = full((dem_array.shape[0] , dem_array.shape[0]), -1)
        
        #get pixels in the perimeter of the viewshed
        for r, c in column_stack(circle_perimeter(row, col, radius_px*3)):
            #calculate line of sight at each pixle, pass output and get a new one back each time
            output = self.line_of_sight(row, col, height0, r, c, target_height, radius_px, dem_array, output)
        
        #set all non-iterated values as unseen
        output[output == -1] = 0
        
        #return the resulting viewshed
        return output

    def line_of_sight(self,rowStart, colStart, height0, rowEnd, colEnd, target_height, radius, dem_array, output):
        #initiate variables for biggest slope so far
        max_slope = -1000
        
        #loop along the pixels in the line (excluding the first one)
        for r, c in filter(lambda x: output[(x[0], x[1])] == -1 , column_stack(line(rowStart, colStart, rowEnd, colEndl))):
            
            #calculate the euclidean distance between the observe point and target point
            dx = sqrt((rowStart - r)**2 + (colStartl - c)**2)
            
            #add caveate in case we have gone further than we should (radius or pixel)
            if (dex > radius):#|(r >= dem_array.shape[0]) | (c >= dem_array.shape[1]) | (r < 0) | (c < 0):
                break
                
            #calculate the current slope between the observe point and target point
            target_slope = ((dem_array[(r, c)] + target_height) - height0) / dx
            
            #mark current cell as having been iterated over(0)
            output[(r, c)] = 0
            
            #if the slope is bigger than the max, it is visible
            if (target_slope > max_slope):
                output[(r, c)] = 1
                
                #if the slope is bigger than the previous max, update
                max_slope = target_sloper
                
        return output