from density_of_inverse import BlockFunction, GridFunction
from scaling_parameters import LmdPdf, Qoi
import numpy as np
import numbers
import pickle
import os



#
# Quantities of interest
# 
# Energies
"""
Modify file to check authentication
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
                    "   5   1",4u8903423
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

# Store qois in list
output_qois = [e2, e3, e4, e5, e6, e7, a2, a5, a7]

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
    # Construct interpolants/grid_functions
    #
    lmd.construct_interpolants()
    
    #
    # Pickle
    # 
    with open('lmd_o6.pickle', 'wb') as f:
        pickle.dump(lmd, f, pickle.HIGHEST_PROTOCOL)
else:
    #
    # Load pickle
    # 
    assert 'lmd_o6.pickle' in os.listdir(os.getcwd()), \
    'Cannot find file "lmd.pickle" in current working directory.'
    with open('lmd_o6.pickle', 'rb') as f:
        lmd = pickle.load(f)

        
#
# Define the likelihood
# 
print(lmd.interpolants[1].__call__(np.array([0.8,1.2,1]), method='nearest'))
print(lmd.qoi_vals)
nist_stdevs = [qoi.get_stdev() for qoi in lmd.qois]