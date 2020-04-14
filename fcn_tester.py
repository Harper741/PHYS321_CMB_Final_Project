import unittest
import nose.tools as nt
import numpy as np
import functions as f
import matplotlib.pyplot as plt
import healpy as hp

class tester():
    
    def setUp(self):
        
        #test boolean masks used to test the boolean mask function
        #that checks that the regions nicely partition the sky
        n_map = 3
        pix_in_map = 100
        self.bool_lis = np.ones((pix_in_map,n_map))
        for i in range(pix_in_map):
            self.rand_temp = np.random.randint(0,3)
            for j in range(n_map):
                if(j == self.rand_temp):
                    self.bool_lis[i][j] = False
                else:
                    self.bool_lis[i][j] = True
                    
        self.bool_masks = np.transpose(self.bool_lis)
        
        #making a mask that shouldn't pass since the whole map is included in every mask
        self.bad_mask_1 = np.zeros((n_map, pix_in_map))
        
        #making a mask that shouldn't pass since the whole map is masked in every mask
        self.bad_mask_2 = np.ones((n_map, pix_in_map))
        
        #changing one value in the good mask list in order to see if the testing function notices
        
        
    def tearDown(self):
        pass
    
    def test_boolean_mask_tester1(self):
        #checking that a good mask passes the function
        nt.assert_true(f.test_masks(self.bool_masks))
        
    def test_boolean_mask_tester2(self):
        #checking that a bad mask doesn't pass the function
        nt.assert_false(f.test_masks(self.bad_mask_1))
        
    def test_boolean_mask_tester3(self):
        #checking that a bad mask doesn't pass the function
        nt.assert_false(f.test_masks(self.bad_mask_2))
        
    def test_boolean_mask_tester4(self):
        #checking that a bad mask doesn't pass the function
        self.bool_lis[0][0] = not self.bool_lis[0][0]
        nt.assert_false(f.test_masks(self.bool_masks))
        
    def test_boolean_mask_maker(self):
        #making a mask based on some regions
        self.loaded_map = hp.read_map("../CMB_maps/LFI_SkyMap_030_1024_R2.01_full.fits")
        hp.pixelfunc.ud_grade(self.loaded_map, 64, order_in='RING')
        self.regions = [[-np.inf,-100],[-100,100],[100,1000],[1000,np.inf]]
        self.made_bool_map = f.splitting_step1(self.regions, self.loaded_map)
        
        #using the test_masks function on this mask to see if it partitions the sky
        nt.assert_true(f.test_masks(self.made_bool_map))
        
    
        
    