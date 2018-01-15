import unittest
from density_of_inverse import BlockFunction, GridFunction
from scaling_parameters import LmdPdf, Qoi
import numpy as np
import numbers
import pickle

class TestScalingParameters(unittest.TestCase):
    """
    """
    def test01(self):
        """
        energies['computed'] = []
        energies['dim'] = 6
        energies['n']   = lmd['n']
        energies['nist_data'] = [4524640.0, 4585620.0, 4585680.0,
                                 4586230.0, 4588380.0, 4629201.0]
        # nist_energies = np.mean(energies, axis=0)                          
        energies['nist_ratings'] = ['AAA','AAA','AAA','AAA','AAA','AAA']  # made-up
        energies['tags'] = [2,3,4,5,6,7]
        energies['search_labels'] = ["1S1 2S1           (3)0( 1.0)",
                                     "1S1 2P1           (3)1( 0.0)",
                                     "1S1 2P1           (3)1( 1.0)",
                                     "1S1 2P1           (3)1( 2.0)",
                                     "1S1 2S1           (1)0( 0.0)",
                                     "1S1 2P1           (1)1( 1.0)"]
        """
        #
        # Quantities of interest
        # 
        e1 = Qoi(category='Energy', tag=2, 
                 search_label='1S1 2S1           (3)0( 1.0)', 
                 nist_value=4524640.0, nist_rating='AAA')
    
        e2 = Qoi(category='Energy', tag=3, 
                 search_label='1S1 2P1           (3)1( 0.0)', 
                 nist_value=4585620.0, nist_rating='AAA')
        
        output_qois = [e1, e2]
        with open('e2.pickle', 'wb') as f:
            pickle.dump(e2, f, pickle.HIGHEST_PROTOCOL)
        
        with open('e2.pickle', 'rb') as f:
            pickle.load(f)
            
        #
        # LmdPdf
        # 
        tags = ['1s', '2s', '2p']
        rng = np.array([[0.8, 1.2], [0.8, 1.2], [0.8, 1.2]])
        resolution = (2,2,2)
        path_to_input = '/home/hans-werner/Dropbox/work/projects'+\
                        '/atomic_data_uncertainty/code/icft/o_6/'
         
        lmd = LmdPdf(tags, rng, resolution, path_to_input, output_qois)
        
        #
        # Interpolants/gridfunctions
        #
        # Constuct
        lmd.construct_interpolants()
        
        # Evaluate interpolant
        points = np.array([[1.1,0.99, 1.01],[1.19, 0.85, 1.0]])
        y1 = lmd.interpolate(points)
        y2 = lmd.interpolate(points, qoi_index=0) 
        y3 = lmd.interpolate(points, qoi_index=[0,1])
        
        # Evaluate gridfunction
        # TODO: Grid Functions cannot be evaluated. 
        #print(lmd.grid_functions[0].eval(points))
        
if __name__ == "__main__":
    unittest.main()