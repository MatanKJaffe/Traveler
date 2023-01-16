"#!/usr/bin/env python"
#-*- coding: utf-8- -*-
from numpy import ndarray, array, pad, s_, ones, zeros, ndenumerate, isnan
from numba import jit

class AnomalyIndex:
    """
    Class prefroming AnomalyIndex calculation based on the work of Andrew D. Weiss as described in his 2001 paper and associated poster 
    (http://www.jennessent.com/downloads/TPI-poster-TNC_18x22.pdf)
    
    """
    def __init__(self):
        pass
    @jit(forceobj=True)
    def _padder(self):
        """
        Function to pad an array by linearly streatching the data at the adge of the array in the necessary axis 

        """           
        self.orig_shape = array(self.array.shape)
        new_shape = array([axis if axis >= max(self.Radius_px_lst) else axis + (max(self.Radius_px_lst)-axis) for axis in self.array.shape])
        #print(f"array axis has been padded to maximum radius e.g. {self.orig_shape} is now {new_shape}")
        pad_width = array([(0, j - i) for i, j in zip(self.orig_shape, new_shape)])
        # example:
        # orig_shape = [79, 123]
        # new_shape = [100, 123]
        # pad_width = [(0, 21), (0, 0)]
        self.array = pad(self.array, pad_width, mode = "edge")
    
    @jit(forceobj=True)
    def _array_shape_to_radius_validator(self):
        """
        Function checking the shape of the input array against the radius given as input, if any radius is larger than one of the array axis than either the arra
        is lengthend and the space filled with a linear regression from the edge or the radius is made to mach the size of the array
        
        """        
        if min(self.array.shape) < max(self.Radius_px_lst):
            if self.padding is True:
                self._padder()
            else:
                new_rad_lst = array([radius - (radius - min(self.array.shape)) if radius > min(self.array.shape) else radius  for radius in self.Radius_px_lst ])
                self.Radius_px_lst = new_rad_lst
        else:
            self.padding = False
    
    @jit(forceobj=True)
    def _sliding_window_view(self, offset_y:ndarray, offset_x:ndarray, step:int = 1):
        """
        Function returning two matching numpy views for moving window routines.

        Args:
            offset_y (ndarray):  refer to the shift in relation to the analysed (central) cell y axis
            offset_x (ndarray):  refer to the shift in relation to the analysed (central) cell x axis
            step (int, optional): size of center point. Defaults to 1.

        Returns:
                ndarray:  window view (in) and main view (out)
        """        
       
        size_y, size_x = self.array.shape
        x, y = abs(offset_x), abs(offset_y)

        x_in = slice(x , size_x, step) 
        x_out = slice(0, size_x - x, step)

        y_in = slice(y, size_y, step)
        y_out = slice(0, size_y - y, step)

        # the swapping trick    
        if offset_x < 0: x_in, x_out = x_out, x_in                                 
        if offset_y < 0: y_in, y_out = y_out, y_in

        # return window view (in) and main view (out)
        return s_[y_in, x_in], s_[y_out, x_out]
    
    @jit(forceobj=True)
    def _calculate_Anomaly_Index(self, radius_px:int)  -> ndarray:
        """
        Function calculating anomI value of dem array at input radius 
        
        Args:
            radius_m (int): radius in meters for the Anomaly index calculation
        
        Returns:
                ndarray: normalized AnomI(spot height – average neighbourhood height) as z score 
        """
   
        #radius in pixels
        #r = int(floor(radius_m/self.resolution))
        r = radius_px
        #instance window
        win = ones((2* r +1, 2* r +1))
        # win = array( [    [0, 1, 1, 1, 0]
        #                      [1, 1, 1, 1, 1],
        #                      [1, 1, 0, 1, 1],
        #                      [1, 1, 1, 1, 1],
        #                      [0, 1, 1, 1, 0]  ])

        # window radius is needed for the function,
        # deduce from window size (can be different for height and width…)
        r_y, r_x  = win.shape[0]//2, win.shape[1]//2
        win[r_y, r_x  ] = 0  # let's remove the central cell 

        #matrices for temporary data
        mx_temp = zeros(self.array.shape)
        mx_count = zeros(self.array.shape)

        # loop through window and accumulate values
        for (y,x), weight in ndenumerate(win):
            if weight == 0 : continue  #skip zero values !
            # determine views to extract data
            #'view_in' is the shifted view and 'view_out' is the position of central cells
            view_in, view_out = self._sliding_window_view(y - r_y, x - r_x)
            
            # using window weights (eg. for a Gaussian function)
            mx_temp[view_out] += self.array[view_in]  * weight

           # track the number of neighbours 
           # (this is used for weighted mean : Σ weights*val / Σ weights)
            mx_count[view_out] += weight

        # this is AnomI (spot value – average neighbourhood value)
        out = self.array - mx_temp / mx_count
        #normalize output as z score
        out = ((((out - out.mean()) / out.std()) * 100) + 0.5)
        
        return out
    
    @jit(forceobj=True)
    def _set_AnomI_classes(self, big_AnomI_array:ndarray ,small_AnomI_array:ndarray )  -> ndarray:
        """
        Use two differently radiused Anomaly Index arrays to generate an 8 classed array.
        
        *disclaimer: the project this code is based on was originally built to run on elevation data and classify landforms
        the label descriptions have not been change from this configuration as i think they describe the labels quite well.

        Args:
            big_AnomI_array (ndarray): larger  and smaller anomI arrays of the same shape 
            small_AnomI_array (ndarray): larger  and smaller anomI arrays of the same shape 
        
        Returns:
                AnomI_class_array (ndarray): classified 
        """
        

        #calculate standard deviation with sensitivity modifier
        std = big_AnomI_array.std()*self.sensitivity

        #instantiate zeros filled array of same size as input to hold tpi calculation valuse
        AnomI_class_array = zeros(big_AnomI_array.shape)

        # Anomaly_Label = 8 | mountain tops | High narrow ridges
        AnomI_class_array[(small_AnomI_array >= std) & (big_AnomI_array >= std)] = 8

        # Anomaly_Label = 7 | flat ridge tops | mesa tops
        AnomI_class_array[((small_AnomI_array > -std) & (small_AnomI_array < std)) & (big_AnomI_array >= std)] = 7

        # Anomaly_Label = 6 | Lateral midslope |drainage divides
        AnomI_class_array[(small_AnomI_array >= std) & ((big_AnomI_array > -std) & (big_AnomI_array < std))] = 6

        # Anomaly_Label = 5 |Local ridge/hilltops within | broad valleys
        AnomI_class_array[(small_AnomI_array >= std) & (big_AnomI_array <= -std)] = 5

        # Anomaly_Label = 4 | U-shape valleys
        AnomI_class_array[((small_AnomI_array > -std ) & (small_AnomI_array < std)) & (big_AnomI_array <= -std)] = 4

        # Anomaly_Label = 3 | Upland incised |drainages |Stream headwaters
        AnomI_class_array[ (small_AnomI_array <= -std) & (big_AnomI_array >= std)] = 3

        # Anomaly_Label = 2 | Lateral midslope | incised drainages
        AnomI_class_array[(small_AnomI_array <= -std) & ((big_AnomI_array > -std ) & (big_AnomI_array < std))] = 2

        # Anomaly_Label = 1 | V-shape river valleys | Deep narrow canyons
        AnomI_class_array[(small_AnomI_array <= -std) & (big_AnomI_array <= -std)] = 1

        return AnomI_class_array

    def execute(self, array:ndarray ,Radius_px_lst:list , padding:bool = True, sensitivity:float = 0.25) -> ndarray :
        """
        Funcution running the AnomalyIndex algorithm
        
        Args:
            array (ndarray): list of radius in px values to comput for feature recognition
            Radius_px_lst (list): array of values for anomaly detection
            padding (bool, optional): _description_. Defaults to True.
        
        Returns:
                ndarray: ndarray of classified pixels 
        """

        assert (type(Radius_px_lst) == list) | (len(Radius_px_lst) == 2)  , "Please enter a list of two integer values to Radius_px_lst"
        self.array = array
        self.Radius_px_lst = Radius_px_lst
        self.padding = padding
        self.sensitivity = sensitivity
        #validate shape and array size compatibility
        self._array_shape_to_radius_validator()

        #calculate large_tpi
        big_radius = self._calculate_Anomaly_Index(radius_px = max(self.Radius_px_lst))
        big_radius[isnan(big_radius)] = 0

        #calculate smaller_tpi
        small_radius = self._calculate_Anomaly_Index(radius_px = min(self.Radius_px_lst))
        small_radius[isnan(small_radius)] = 0
        
        #calculate landform classes
        anomi_array = self._set_AnomI_classes(big_radius , small_radius)
        #resize output as relating to padding
        if self.padding is True:
            return anomi_array[:self.orig_shape[0],:self.orig_shape[1]]

        return anomi_array