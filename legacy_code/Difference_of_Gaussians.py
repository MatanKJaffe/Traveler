import cv2
import numpy as np
from numba import jit

class DoG:
    """
    Class preforming difference between two Gausian shifts applied to an image 

    Returns:
        ndarray: ndarray of 
    """
    
    def __init__(self, img:str or np.ndarray, read_image:bool = False, grayscale:bool = False , k1:int = 5, s1:int = 2, k2:int = 3, s2:int = 1):
        """
        Base init function containing boolians which decide object configuration

        Args:
            img (str or np.ndarray): either a ndarray of an image or the str directory to that image
            read_image (bool, optional): wheather to read img as a str directory. Defaults to False.
            grayscale (bool, optional): weather to procces the input array into grayscael representation . Defaults to False.
            k1 (int, optional): kernal size. Defaults to 5.
            s1 (int, optional): sigma. Defaults to 2.
            k2 (int, optional): kernal size. Defaults to 3.
            s2 (int, optional): sigma. Defaults to 1.
        """        
        if read_image is True:
            self.img = cv2.imread(img)
        else:
            self.img = img

        if grayscale is True:
            self.img  = self._greyscale(self.img)
        self.k1, self.s1, self.k2, self.s2 = k1, s1, k2, s2
    
    @jit(forceobj=True)
    def _greyscale(self):
        """
        Returns:
            grayscale of orginal image array
        """
        self.img  = cv2.cvtColor(self.img, cv2.COLOR_BGR2GRAY)
    
    @jit(forceobj=True)
    def difference_of_Gaussians(self)  -> np.ndarray:
        """
        Function to run diffrence of Gaussians on a given image array

        Returns:
            ndarray: the diffrence between two gaussian blur on an image 
        """
        # Apply Gaussian blur
        b1 = cv2.GaussianBlur(self.img,(self.k1, self.k1), self.s1, borderType=cv2.BORDER_REPLICATE)
        b2 = cv2.GaussianBlur(self.img,(self.k2, self.k2), self.s2, borderType=cv2.BORDER_REPLICATE)

        # Calculate the DoG by subtracting
        return b1 - b2