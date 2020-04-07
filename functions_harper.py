import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

#checks that the boolean masks partition the sky
#takes a list masks = [[mask region 1],...,[mask region n]] where n is the number of regions
#returns true if the masks partition the sky, false if they don't
def test_masks(masks):
	#making an array of the proper length to hold boolean values
	vals = masks[0]

	#doing the logical XOR operation across each element of the masks to make sure
	#that only one of them is true for each index
	for i in range(len(masks)-1):
		vals = np.logical_xor(vals, masks[i+1])

	if(np.all(vals)):
		print('The masks partition the sky')
		return(np.all(vals))
	else:
		print('ERROR: Masks do not partition the sky')

	return(np.all(vals))


#function that combines a map and a mask, and sets the masked values to zero
#takes a map of temp values for the whole sky
#takes a mask, list of boolean values, 1 for each sky pixel, if mask[i]=True, then we set map[i]=0
#returns a map with the masked values set to zero
def map_and_mask(map,mask):

	#check each mask value and set to zero if the value is masked
	for i in range(len(map)):
		if(mask[i]):
			map[i] = 0

	return(map)

#function that combines all of the seperate maps into one
#takes maps = [[region 1 temps],[region 2 temps],...,[region n temps]] where n is the number of regions
#takes masks = [[region 1 boolean mask],...,[region n boolean mask]] where n is the number of regions
#returns a sky map that is all of the regions put togather
def recombine_maps(maps, masks):
	if(not test_masks(masks)):
		return()

	#making an array of new masks to add together
	new_maps = [map_and_mask(maps[i],masks[i]) for i in range(len(masks))]

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

	#applying weights and returning the result
	weighted = [weights[i]*freqs[i] for i in range(len(weights))]
	print(weighted)
	return(np.sum(weighted,axis=0))

#function that takes edges of intervals and sets the masked
#values to hp.UNSEEN so masking is easy
def mask_outside_of_interval(interval, map):
	print(map)
	for i in range(len(map)):
		if((not (interval[0] < map[i])) or (not (map[i] <= interval[1])) ):
			map[i] = hp.UNSEEN

	return(map)

#function to compute weights that minimize the variance of the sum of the channels

#takes a list of temperatures for a region, one for each frequency being used for weight computation
#takes [[freq 1 temps],[freq 2 temps],...,[freq n temps]]
#returns [weight 1, weight 2, ..., weight n]
def compute_weights(freq_temps):
    #making matrix
    n = len(freq_temps)
    
    H = np.ones((n,n))
    
    #constructiong the symmetric H matrix
    for i in range(n):
        for j in range(i,n-i):
            H[i][j] = np.sum(freq_temps[i]*freq_temps[j])
            if(not (i==j)):
                H[j][i] = H[i][j]
            
    #computing the weights using the formula derived from minimizing the variance while
    #requiring the weights to sum to 1
    Hinv = np.linalg.inv(H)
    print(H,"\n" ,Hinv)
    unit_vec = np.ones(n)
    weights = (Hinv @ unit_vec)/(np.transpose(unit_vec) @ Hinv @ unit_vec )

    return(weights)








