from numba import jit
from numpy import gradient, ndarray, uint8 ,pi ,arctan ,arctan2 , sin, cos, sqrt, uint8, degrees, array as nparray

@jit
def hillshade(array: ndarray ,azimuth: int ,angle_altitude:int):
    azimuth = 360.0 - azimuth 
    x, y = gradient(array)
    slope = pi/2. - arctan(sqrt(x*x + y*y))
    aspect = arctan2(-x, y)
    azimuthrad = azimuth*pi/180.
    altituderad = angle_altitude*pi/180.
    shaded = sin(altituderad)*sin(slope) + cos(altituderad)*cos(slope)*cos((azimuthrad - pi/2.) - aspect)
    return nparray(255*(shaded + 1)/2 , dtype = uint8) 

@jit
def aspect(array: ndarray, deg: int = False) -> ndarray:  
    x, y = gradient(array)
    if deg is True:
        return nparray(degrees(arctan2(-x, y)), dtype = uint8)
    return  nparray(arctan2(-x, y), dtype = uint8)
    
@jit
def slope(array: ndarray, deg: int = False):
    x, y = gradient(array)
    if deg is True:
        return  nparray(degrees(pi/2. - arctan(sqrt(x*x + y*y))), dtype = uint8)
    return nparray(pi/2. - arctan(sqrt(x*x + y*y)), dtype = uint8)


