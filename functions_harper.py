import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

#checks that the boolean masks partition the sky
def test_masks(masks):
	#making an array of the proper length to hold boolean values
	vals = masks[0]

	#doing the logical XOR operation across each element of the masks to make sure
	#that only one of them is true for each index
	for i in range(len(masks)-1):
		vals = np.logical_xor(vals, masks[i+1])

	if(np.all(vals)):
		print('The masks partition the sky')
	else:
		print('ERROR: Masks do not partition the sky')

	return(np.all(vals))


#function that combines a map and a mask, and sets the masked values to zero
def map_and_mask(map,mask):

	#check each mask value and set to zero if the value is masked
	for i in range(len(map)):
		if(mask[i]):
			map[i] = 0

	return(map)

#function that combines all of the seperate maps into one
def recombine_maps(maps, masks):
	if(not test_masks(masks)):
		return()

	#making an array of new masks to add together
	new_maps = [map_and_mask(maps[i],masks[i]) for i in range(len(masks))]

	new_maps = np.sum(new_maps, axis=0)
	return(new_maps)

#function that applies the weights to a given frequency range
def weight_freqs(weights, freqs):

	#checking that the weights and frequencies are the same length
	if(not (len(weights)==len(freqs))):
		print("ERROR:different number of weights and frequencies")
		return()

	#applying weights and returning the result
	weighted = [weights[i]*freqs[i] for i in range(len(weights))]
	return(np.sum(weighted))

#function that takes edges of intervals and sets the masked
#values to hp.UNSEEN so masking is easy
def mask_outside_of_interval(interval, map):
	print(map)
	for i in range(len(map)):
		if((not (interval[0] < map[i])) or (not (map[i] <= interval[1])) ):
			map[i] = hp.UNSEEN

	return(map)







