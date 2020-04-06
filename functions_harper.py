import numpy as np
import matplotlib.pyplot as plt
import healpy as hp

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


def recombine_maps(maps, masks):
	if(not test_masks(masks)):
		return()

	#making an array of new masks to add together
	new_maps = [map_and_mask(maps[i],masks[i]) for i in range(len(masks))]

	new_maps = np.sum(new_maps, axis=0)
	return(new_maps)



