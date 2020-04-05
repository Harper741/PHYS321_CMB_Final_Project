import numpy as np
import math
import healpy as hp
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

def splitting(regions, data_map):
'function that takes in the splitting regions (this array must be in the format value,rel_range,abs_range) and data maps and returns an array of regions cut out from the initial map as well as the boolean map'
    value, rel_range, abs_range = regions
    
    cut_map = [] #cut out pieces of original data map
    bool_map = [] #entire map as bool values indicating whether a value is contained on a piece of a map or not
    size = [] #size of cut-out region
    
    for i in range(len(regions)):
        bool_map.append(ealpy.pixelfunc.mask_bad(data_map, badval=value[i], rtol=rel_range[i], atol=abs_range[i]))
        cut_map.append(data_map[bool_map[i] == True])
        size.append(len(cut_map))
    
    if np.sum(size =! data_map):
        print('Error: a pixel on a boundary was excluded, try changing the absolute tolerance')
        break 
    else:
        return cut_map, bool_map
    
    
    
def find_weight(regions,maps):
    'function that takes the splitting regions