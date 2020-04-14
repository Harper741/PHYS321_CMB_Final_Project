import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename


#checks that the boolean masks partition the sky
#takes a list masks = [[mask region 1],...,[mask region n]] where n is the number of regions
#returns true if the masks partition the sky, false if they don't
def test_masks(masks):
    

    #masks are by True for excluded points, and false for included points, so we need to reverse this
    #for this test to work
    negated_masks = []
    for i in masks:
        negated_masks.append(np.logical_not(i))

    #making an array of the proper length to hold boolean values
    vals = np.copy(negated_masks[0])

    #doing the logical XOR operation across each element of the masks so that for each pixel,
    #if there is not exactly one true element, there will be a false value for that pixel
    for i in range(len(masks)-1):
        vals = np.logical_xor(vals, negated_masks[i+1])

    #vals should be an array of all true values if the mask partitions the sky
    #checking that indeed all elements are true
    if(np.all(vals)):
        print('The masks partition the sky')
        return(np.all(vals))
    else:
        print('ERROR: Masks do not partition the sky')

    return(np.all(vals))


#function that combines a map and a mask, and makes healpy ignore the masked values
#takes a mask, list of boolean values, 1 for each sky pixel, if mask[i]=True, then we set map[i]=hp.UNSEEN
#returns a map with the masked values set to hp.UNSEEN, so that healpy will ignore those pixels.
def map_and_mask(mask,data_map):
    
    a = np.copy(data_map) #to avoid pointer errors

    #check each mask value and set to hp.UNSEEN if the value is masked
    for i in range(len(data_map)):     
        if(mask[i]):
            a[i] = hp.UNSEEN

    return(a)

#function that combines a map and a mask, and sets the masked values to zero
#takes a mask, list of boolean values, 1 for each sky pixel, if mask[i]=True, then we set map[i]=hp.UNSEEN
#returns a map with the masked values set to hp.UNSEEN, so that healpy will ignore those pixels.
def map_and_zero(mask,data_map):
    
    a = np.copy(data_map) #to avoid pointer errors
    #check each mask value and set to zero if the value is masked
    for i in range(len(data_map)):     
        if(mask[i]):
            a[i] = 0

    return(a)

#function that combines all of the seperate maps into one
#takes maps = [[region 1 temps],[region 2 temps],...,[region n temps]] where n is the number of regions
#takes masks = [[region 1 boolean mask],...,[region n boolean mask]] where n is the number of regions
#returns a sky map that is all of the regions put togather
def recombine_maps(maps, masks):
    if(not test_masks(masks)):
       return()

    #making an array of new masks to add together
    new_maps = np.array([map_and_zero(masks[i], maps[i]) for i in range(len(masks))])

    new_map = np.sum(new_maps, axis=0)
    return(new_map)

#function that applies the weights to a given frequency range
#takes weights = [weight 1, weight 2, ..., weight n] n is number of frequencies
#takes freqs = [[freq1 temps], [freq 2 temps], ... ,[freq n temps]] n is number of frequencies
#returns [temp 1 weighted, temp 2 weighted, ..., temp n weighted] n is number of sky pixels
def weight_freqs(weights, freqs):

    #checking that the weights and frequencies are the same length
    if(not (len(weights)==len(freqs))):
        print("ERROR:different number of weights and frequencies")
        return()

    #applying weights and returning the summed result
    weighted = [weights[i]*freqs[i] for i in range(len(weights))]
    return(np.sum(weighted,axis=0))

#function that takes edges of masked intervals and sets the masked
#values to hp.UNSEEN so masking is easy
def mask_outside_of_interval(interval, the_map):
    map2 = np.copy(the_map) #avoiding pointer errors

    #Masking everything outsied of the given interval
    for i in range(len(map2)):
        if( (interval[0] >= map2[i]) or (map2[i] > interval[1]) ):
            map2[i] = hp.UNSEEN

    return(map2)

#function to compute weights that minimize the variance of the sum of the channels

#takes a list of temperatures for a region, one for each frequency being used for weight computation
#takes [[freq 1 temps],[freq 2 temps],...,[freq n temps]]
#returns [weight 1, weight 2, ..., weight n]
def compute_weights(freq_temps):

    n = len(freq_temps)
    H = np.ones((n,n))
    
    #constructiong the symmetric H matrix
    for i in range(n):
        for j in range(i,n):
            H[i][j] = np.sum(freq_temps[i]*freq_temps[j])

            #since H is symmetric we use this to avoid computing matrix elements more than is necessary
            if(not (i==j)):
                H[j][i] = H[i][j]
            
    #computing the weights using the formula derived from minimizing the variance while
    #requiring the weights to sum to 1
    Hinv = np.linalg.inv(H)
    
    unit_vec = np.ones(n)

    #formula derived from lagrange multipliers
    weights = (Hinv @ unit_vec)/(np.transpose(unit_vec) @ Hinv @ unit_vec )
    return(weights)

#FOR JUNK MAP
    #function that takes in the splitting regions and data maps and returns an array of regions cut out from the initial map as well as the boolean map
    #assumes regions is a list of intervals that partition the range of available data ex: [[-20,-10],[-10,10],[10,20]] for data
    #that is valued from -20 to 20
def splitting_step1(regions, data_map): 
    
    cut_map = [] #cut out pieces of original data map
    bool_map = [] #entire map as bool values indicating whether a value is contained on a piece of a map or not
    size = [] #size of cut-out region
    
    #for each interval that defines a region, we construct a boolean mask with the healpy function mask_bad, which can be
    #used by other healpy functions. We also make a cut map with all masked valued deleted for computing region weights

    for i in range(len(regions)):
        temp_map = mask_outside_of_interval(regions[i], data_map)
        bool_map.append(hp.pixelfunc.mask_bad(temp_map))
        #hp.pixelfunc.mask_bad is a function that sets the hp.UNSEEN values to True
        cut_map.append(data_map[bool_map[i] == False])

        #so we can chack that the regions are appropriately sized
        size.append(len(cut_map[-1]))
        
    #checking that the cut regions' lengths add up to the total length of the sky map
    if (np.sum(size) != len(data_map)):
        print('Error: a pixel on a boundary was excluded or included in two or more regions, try changing the absolute tolerance')
        print('all cut out regions sum up to', np.sum(size))
        print('the pixel number:', len(data_map))
    else:
        return bool_map

def splitting_step2(bool_map, data_map): #for other maps, after junk map
    #function that takes in the boolean masks and data map and returns an array of regions cut out from that map
       
    cut_map = [] #cut out pieces of original data map
    len_map = 0
    
    #cutting the maps
    for i in range(len(bool_map)):
        cut_map.append(data_map[bool_map[i] == False])
        len_map += len(cut_map[-1])
    
    #checking pixel lengths
    if (len_map != len(data_map)):
        print('Error: a pixel on a boundary was excluded or included in two or more regions')
        print(len_map, len(data_map))
       
    else:
        return cut_map #array of cut-out maps 
    
    #takes in a list of bin edges and returns a list of pairs giving the upper and lower bounds for each region
    #for example [1,2,3] -> [[1,2],[2,3]]
def bins_to_regions(bin_edges):
    
    regions = []
    
    #extending the outermost regions a bit to make sure that there are no 
    #points that don't get included due to some kind of floating point rounding error etc.
    bin_edges[0] -= 1
    bin_edges[-1] += 1
    
    for i in range(len(bin_edges)-1):
        regions.append([bin_edges[i],bin_edges[i+1]])
    
    return(regions)


#takes in the power spectrum c and returns the y-axis suitable value to plot that is scaled by L(L+1)
def yaxis_pow_spec(c, l=1025):
    y = []
    for i in range(l):
        y.append(i*(i+1)*c[i]/(2*np.pi))
    return y
