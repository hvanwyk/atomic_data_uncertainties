import unittest
from density_of_inverse import BlockFunction, GridFunction
from scaling_parameters import LmdPdf, Qoi
import numpy as np
import numbers
import pickle
import os


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
        # Energies
        e2 = Qoi(category='Energy', tag=2, 
                 search_label='1S1 2S1           (3)0( 1.0)', 
                 nist_value=4524640.0, nist_rating='AAA')
    
        e3 = Qoi(category='Energy', tag=3, 
                 search_label='1S1 2P1           (3)1( 0.0)', 
                 nist_value=4585620.0, nist_rating='AAA')
        
        e4 = Qoi(category='Energy', tag=4, 
                 search_label='1S1 2P1           (3)1( 1.0)', 
                 nist_value=4585680.0, nist_rating='AAA')
        
        e5 = Qoi(category='Energy', tag=5, 
                 search_label='1S1 2P1           (3)1( 2.0)', 
                 nist_value=4586230.0, nist_rating='AAA')
        
        e6 = Qoi(category='Energy', tag=6, 
                 search_label='1S1 2S1           (1)0( 0.0)', 
                 nist_value=4588380.0, nist_rating='AAA')
        
        e7 = Qoi(category='Energy', tag=7, 
                 search_label='1S1 2P1           (1)1( 1.0)', 
                 nist_value=4629201.0, nist_rating='AAA')
        
        # A-values
        """
        avalues['nist_data'] = [1.04e3, None, 3.31e5, None, 3.309e12]       
        avalues['nist_ratings'] = ['AA', None, 'A', None, 'AA']
        avalues['tags'] = [2,4,5,6,7]
        avalues['search_labels'] = \
                           ["   2   1", 
                            "   4   1",
                            "   5   1",
                            "   6   1",
                            "   7   1"]
        """
        a2 = Qoi(category='A-value', tag=2, 
                 search_label='   2   1', 
                 nist_value=1.04e3, nist_rating='AA')
        
        a5 = Qoi(category='A-value', tag=5, 
                 search_label='   5   1', 
                 nist_value=3.31e5, nist_rating='AA')
        
        a7 = Qoi(category='A-value', tag=7, 
                 search_label='   7   1', 
                 nist_value=3.309e12, nist_rating='AA')
        #
        # Test pickle
        # 
        output_qois = [e1, e2]
        with open('e2.pickle', 'wb') as f:
            pickle.dump(e2, f, pickle.HIGHEST_PROTOCOL)
        
        with open('e2.pickle', 'rb') as f:
            pickle.load(f)
            
        #
        # LmdPdf
        #
        lmd_saved = True
        if not lmd_saved:
            #
            # Initialize
            #
            tags = ['1s', '2s', '2p']
            rng = np.array([[0.8, 1.2], [0.8, 1.2], [0.8, 1.2]])
            resolution = (20,20,20)
            path_to_input = '/home/hans-werner/Dropbox/work/projects'+\
                            '/atomic_data_uncertainty/code/icft/o_6/'
            lmd = LmdPdf(tags, rng, resolution, path_to_input, output_qois)
            
            #
            # Construct interpolants/grid_
            #
            lmd.construct_interpolants()
            
            #
            # Pickle
            # 
            with open('lmd.pickle', 'wb') as f:
                pickle.dump(lmd, f, pickle.HIGHEST_PROTOCOL)
        else:
            #
            # Load pickle
            # 
            assert 'lmd.pickle' in os.listdir(os.getcwd()), \
            'Cannot find file "lmd.pickle" in current working directory.'
            with open('lmd.pickle', 'rb') as f:
                lmd = pickle.load(f)
    
        #
        # Determine indices of qois
        # 
        self.assertEqual(lmd.get_qoi_index(category='Energy', tag=2), 0, \
                         'Qoi index incorrectly identified.')
        self.assertEqual(lmd.get_qoi_index(qoi=e1), 0, \
                         'Qoi index incorrectly identified.')
        
        
        # Evaluate interpolant
        points = np.array([[1.1,0.99, 1.01],[1.19, 0.85, 1.0]])
        #y1 = lmd.interpolate(points)
        #y2 = lmd.interpolate(points, qoi_index=0) 
        #y3 = lmd.interpolate(points, qoi_index=[0,1])
        
        # Evaluate gridfunction
        # TODO: Grid Functions cannot be evaluated. 
        #print(lmd.grid_functions[0].eval(points))
        
if __name__ == "__main__":
    unittest.main()