import numpy as np
import math
import healpy as hp
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import functions_harper as hf


   
def splitting_step1(regions, data_map): #FOR JUNK MAP
    #function that takes in the splitting regions and data maps and returns an array of regions cut out from the initial map as well as the boolean map
    #assumes regions is a list of intervals that partition the real line ex: [[-np.inf,-10],[-10,10],[10,np.inf]]
       
    cut_map = [] #cut out pieces of original data map
    bool_map = [] #entire map as bool values indicating whether a value is contained on a piece of a map or not
    size = [] #size of cut-out region
    
    for i in range(len(regions)):
        
        temp_map = hf.mask_outside_of_interval(regions[i], data_map)
        bool_map.append(hp.pixelfunc.mask_bad(temp_map))
        cut_map.append(data_map[bool_map[i] == True])
        size.append(len(cut_map))
    
    if (np.sum(size) != len(data_map)):
        print('Error: a pixel on a boundary was excluded or included in two or more regions, try changing the absolute tolerance')
         
    else:
        return cut_map, bool_map
    
def splitting_step2(bool_maps, data_map): #for other maps, after junk map
    #function that takes in the splitting regions and data maps and returns an array of regions cut out from the initial map
       
    cut_map = [] #cut out pieces of original data map
    size = [] #size of cut-out region
    
    for i in range(len(bool_maps)):
        cut_map.append(data_map[bool_map[i] == True])
        
    if (len(np.flatten(cut_map)) != len(data_map)):
        print('Error: a pixel on a boundary was excluded or included in two or more regions, try changing the absolute tolerance')
       
    else:
        return cut_map #array of cut-out maps 
