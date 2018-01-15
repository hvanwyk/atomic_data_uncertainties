"""
This module is used to estimate and sample from the joint density function
of the atomic orbital scaling parameters. 

The density function
 
- specify lambda range and resolution
- loop over lambda space
- change the input data in "input.dat"
- run ../adas803.testern.pl --proc=pp --inp --born input.dat 8
- go into the file born/adf04ic and extract energies and avalues
   
- store a-values and energies in variables. 
"""
# =============================================================================
# Import Libraries
# =============================================================================
# System tools
import os
import sys
import subprocess

# Scientific Comptuting
import numpy as np
from scipy import linalg as la
from scipy.interpolate import RegularGridInterpolator
from itertools import product
from scipy.stats import norm

# Density estimation
from density_of_inverse import GridFunction

import matplotlib.pyplot as plt


class Qoi(object):
    """
    Quantity of interest (such as an energy or A-value)


    Attributes

        category: 
        
        tag:
        
        search_label:
        
        nist_data:
        
        nist_rating:
        

    Methods
    
        compute_histogram: 
        
        
    """
    def __init__(self, category, tag, search_label, nist_data, nist_rating):
        """
        Constructor

        Inputs:

            category: str, 'A-value' or 'Energy'

            tag: str/int used to identify specific qoi for user

            search_label: str, used to search for qoi in adf04ic file

            nist_data: double, nominal NIST value (from database)

            nist_rating: str, rating for the estimated measurement error

               'AAA': 0.003
               'AA' : 0.01 
               'A+' : 0.02 
               'A'  : 0.03 
               'B+' : 0.07 
               'B'  : 0.1 
               'C+' : 0.18
               'C'  : 0.25
               'D+' : 0.4 
               'D'  : 0.5 
               'E'  : 0.5
        """
        self.category = category
        self.tag = tag
        self.search_label = search_label
        self.nist_data = nist_data
        self.nist_rating = nist_rating


    def set_histogram(self, density):
        """
        Compute the histogram
        """
        pass
    
    def sample(self, n):
        """
        Sample from density function associated with quantity
        """
        pass


class LmdPdf(object):
    """
    Class for storing the joint distribution of orbital scaling parameters


    Attributes:
    
        dim: int, number of orbital scaling parameters considered
        
        n_points: int, number of gridpoints used. 
        
        oned_vals: double, dim-length list of one dimensional grids, one 
            for each lambda. 
        
        path_to_input: str, path to input file used to call the perl script
            samples the adas code.

        qois: Qoi, list of quantities of interest (energies and A-values) 
            dependent on lambda. 
        
        range: double, (dim,2) array, the ith row of which contains the
            endpoints for the ith component of lambda
        
        resolution: int, dim-tuple specifying the number of cells in grid
            for each direction

        tags: str/int, list of dim names used to identify scaling parameters.
        
        vals: double, (n,dim) array of lambda values on the grid.
        
        
    Methods: 

        compute_output_qois: 

        construct_grid_functions:
        
        get_qoi_index: Determine the position of a qoi in the "qois" list, based on
            its category and tag. 
    
    """
    def __init__(self, tags, rng, resolution, 
                 path_to_input_file, output_qois):
        """
        Constructor: 


        Inputs:
            
            tags: str/int, list of dim names used to identify scaling parameters.

            rng: double, (dim,2) array the ith row of which specifies the range
                of the ith lambda parameter. 

            resolution: int, dim-tuple specifying the number of subintervals
                for each lambda parameter. 
            
            path_to_input_file: str, path to input file used by the adas code.

            output_qois: Qoi, list of quantities of interest 
            
        """
        self.tags = tags
        self.resolution = resolution
        self.dim = len(resolution)  
        self.n_points = np.prod(resolution)
        self.range = rng
        oned_grids = []
        for i in range(self.dim):
            oned_grids.append(np.linspace(self.range[i,0], self.range[i,1],\
                             self.resolution[i]))
        self.oned_vals = oned_grids
        L = np.meshgrid(*oned_grids)
        self.vals = np.array([L[i].ravel() for i in range(self.dim)]).T
        self.path_to_input = path_to_input_file
        self.qois = output_qois
    
        
    def construct_interpolants(self):
        """
        """
        #
        # Compute Qoi's at all grid points.
        # 
        n_qois = len(self.qois)
        qoi_vals = np.empty((self.n_points, n_qois))
        for i in range(self.n):
            point = self.vals[i,:]
            qoi_vals[i,:] = self.map_to_qois(point)
        # Store values
        self.qoi_vals = qoi_vals
        
        interpolants = []
        grid_functions = []
        for j in range(n_qois):
            #
            # Interpolate j'th output for quantity of interest
            #
            qoi_vals_grid = np.reshape(qoi_vals[:,j], self.resolution)
            f = RegularGridInterpolator(tuple(self.oned_grids), qoi_vals_grid, 
                                        bounds_error=False)
            interpolants.append(f)
            
            #
            # Estimate Grid Functions
            # 
            n_intervals = tuple([i-1 for i in self.resolution])
            bnd = self.range
            g = GridFunction(bnd, n_intervals, f)
            grid_functions.append(g)
            
        self.interpolants = interpolants
        self.grid_functions = grid_functions
           
        
    def map_to_qois(self, point):
        """
        Returns the 
        
        Input:
            
            point: double, (dim,) array specifying a single sample point
        
            qoi_list: Qoi, list of quantities of interest to be computed
            
            
        Output:
        
            qoi_vals: double, length-n list of quantities of interests
                computed.
        """
        tags = self.tags.copy()
        point = list(point)
        n_qois = len(self.qois)
        #
        # Modify lambda parameters in autostructure input file
        #     
        with open(self.path_to_input,'r+') as infile:
            lines = infile.readlines()
            infile.seek(0)
            infile.truncate()
            for line in lines:
                modified = False
                if len(tags)>0 and not modified:
                    for j in range(len(tags)):
                        if tags[j] + ' = ' in line:
                            tag, l = tags.pop(j), point.pop(j)
                            line = tag + ' = ' + str(l) + '\n'
                            modified = True
                            break
                infile.write(line)
        # 
        # Call autostructure function
        # 
        # TODO: Suppress autostructure stdout
        # TODO: Fix folder
        subprocess.call(['../adas803.testern.pl', '--proc=pp ', 
                         '--inp', '--born', 'input.dat', '8'])
    
        #
        # Search output files for energies and A-values
        #
        with open('born/adf04ic','r') as infile:
            # Initialize temporary storage
            qoi_vals = [None]*n_qois
            qoi_extracted = [False]*n_qois
            for line in infile:
                while not all(qoi_extracted): 
                    #
                    # Extract computed qoi
                    #
                    for i in range(n_qois):
                        if self.qois[i].search_label in line:
                            #
                            # Found Qoi in file
                            # 
                            qoi = self.qois[i]
                            words = line.split() 
                            if qoi.category == 'Energy':
                                #
                                # Energy
                                #
                                qoi_vals[i] = float(words[-1])
                                qoi_extracted[i] = True
                                break 
                            elif qoi.category == 'A-value':
                                #
                                # A-value
                                #
                                
                                # To convert data to floats, we use the fact 
                                # that A-values are stored to 2 decimals. This
                                # may have to be changed in the generic case. 
                                aval = words[2]
                                dec_pos = aval.find('.')
                                a = float(aval[:dec_pos+3] + 'e' + aval[dec_pos+3:])
                                qoi_vals[i] = a
                                qoi_extracted[i] = True
                                break        
    
    
    def interpolate(self, points):
        """
        Use the interpolants constructed via "construct_interpolants" to 
        interpolate the mapping from lambda parameters to qois.
        
        Inputs:
        
            points: double, (n, dim) array of points
        """
        pass
                     
        
    def get_qoi_index(self):
        """
        """
        pass
    
    


'''

# =============================================================================
# Import Libraries
# =============================================================================
# System tools
import os
import sys
#import stat
import subprocess
#from tempfile import mkstemp
#from shutil import move


# Scientific Comptuting
import numpy as np
from scipy import linalg as la
from scipy.interpolate import RegularGridInterpolator
from itertools import product

# Calibration
sys.path.insert(0, '/home/hans-werner/Dropbox/work/projects/'+
                    'atomic_data_uncertainty/code/src/')
from density_of_inverse import GridFunction

# Plotting tools
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

# =============================================================================
# Initialize
# =============================================================================

# .............................................................................
#  Scaling Parameters
# .............................................................................
#
# Initialize lambda dictionary
#
lmd = dict.fromkeys(['dim', 'lmd_1d', 'n', 'range', 'resolution', 
                     'tags', 'vals'])

#
# Specify sampling parameters
#
lmd['dim'] = 3  
lmd['resolution'] = (2,2,2) 
lmd['tags'] = ['1s', '2s', '2p'] 
lmd['range'] = np.array([[0.5, 3], [0.5, 3], [0.5, 3]])
lmd['n'] = np.prod(lmd['resolution'])
#
# Generate lambda values
#
lmd_samples = []
for i in range(lmd['dim']):
    lmd_samples.append(np.linspace(lmd['range'][i,0],\
                                   lmd['range'][i,1],\
                                   lmd['resolution'][i]))
lmd['lmd_1d'] = lmd_samples
L = (1,)*lmd['dim']  # initialize tuple
L = np.meshgrid(*lmd_samples)
lmd['vals'] = np.array([L[i].ravel() for i in range(lmd['dim'])]).T


# ------------------------------------------------------------------------------
# Energies
# ------------------------------------------------------------------------------
energies = dict.fromkeys(['computed', 'dim', 'grid_functions', 'interpolants', 
                          'n', 'nist_data', 'nist_ratings', 'tags', 
                          'search_labels'])
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

# ..............................................................................
# A-Values
# ..............................................................................
avalues = dict.fromkeys(['computed', 'dim', 'grid_functions', 'interpolants', 
                         'n', 'nist_data', 'nist_ratings', 'tags', 
                         'search_labels'])

avalues['computed'] = []
avalues['dim'] = 5
avalues['n']   = lmd['n']
avalues['nist_data'] = [1.04e3, None, 3.31e5, None, 3.309e12]       
avalues['nist_ratings'] = ['AA', None, 'A', None, 'AA']
avalues['tags'] = [2,4,5,6,7]
avalues['search_labels'] = ["   2   1", 
                            "   4   1",
                            "   5   1",
                            "   6   1",
                            "   7   1"]


# ==============================================================================
# Sample Energies and A-Values (Autostructure Code)
# ==============================================================================
path_to_input = 'input.dat'
for i in range(lmd['n']):
    # ..........................................................................
    # Loop over lambda samples
    # ..........................................................................
    tags = lmd['tags'].copy()    
    lmd_sample = list(lmd['vals'][i,:])
    
    #
    # Modify lambda parameters in autostructure input file
    #     
    with open(path_to_input,'r+') as infile:
        lines = infile.readlines()
        infile.seek(0)
        infile.truncate()
        for line in lines:
            modified = False
            if len(tags)>0 and not modified:
                for j in range(len(tags)):
                    if tags[j] + ' = ' in line:
                        tag, l = tags.pop(j), lmd_sample.pop(j)
                        line = tag + ' = ' + str(l) + '\n'
                        modified = True
                        break
            infile.write(line)
    # 
    # Call autostructure function
    # 
    # TODO: Suppress autostructure stdout
    subprocess.call(['../adas803.testern.pl', '--proc=pp ', 
                     '--inp', '--born', 'input.dat', '8'])

    #
    # Search output files for energies and A-values
    #
    with open('born/adf04ic','r') as infile:
        # Initialize temporary storage
        energies_data = [None]*energies['dim']
        avalues_data = [None]*energies['dim']

        # Keeps track of which values have been extracted
        energies_extracted = [False]*energies['dim']
        avalues_extracted = [False]*avalues['dim']
        for line in infile:
            line_processed = False
            if not all(energies_extracted):
                #
                # Extract computed energies
                #
                for i in range(energies['dim']):
                    if energies['search_labels'][i] in line:
                        words = line.split()
                        energies_data[i] = float(words[-1])
                        energies_extracted[i] = True
                        line_processed = True
                        break

            if not all(avalues_extracted) and not line_processed:
                #
                # Extract computed A-values
                #
                for i in range(avalues['dim']):
                    if avalues['search_labels'][i] in line:
                        words = line.split()
                        #
                        # To convert data to floats, we use the fact that 
                        # A-values are stored to 2 decimals. This may have to be
                        # changed in the generic case. 
                        #
                        aval = words[2]
                        dec_pos = aval.find('.')
                        a = float(aval[:dec_pos+3] + 'e' + aval[dec_pos+3:])
                        avalues_data[i] = a
                        avalues_extracted[i] = True
                        break

    energies['computed'].append(energies_data)
    avalues['computed'].append(avalues_data)

energies['computed'] = np.array(energies['computed'])
avalues['computed'] = np.array(avalues['computed']) 

#
# Store all data in Structure
# 
structure = {'lambda': lmd, 'energies': energies, 'avalues': avalues} 

# ==============================================================================
# Approximate the lambda to QoI (autostructure) mappings
# ==============================================================================
lmd_1d = structure['lambda']['lmd_1d']
for qoi in ['energies', 'avalues']:
    interpolants = []
    grid_functions = [] 
    for i in range(structure[qoi]['dim']):
        #
        # Interpolate i'th output for quantity of interest
        #
        output = structure[qoi]['computed'][:,i]
        output_grid = np.reshape(output, structure['lambda']['resolution'])
        f = RegularGridInterpolator(tuple(lmd_1d), output_grid, 
                                    bounds_error=False)
        interpolants.append(f)
        
        #
        # Estimate Grid Functions
        # 
        resolution = list(structure['lambda']['resolution'])
        n_intervals = tuple([i-1 for i in resolution])
        bnd = structure['lambda']['range']
        g = GridFunction(bnd, n_intervals, f)
        grid_functions.append(g)
        
    structure[qoi]['interpolants'] = interpolants   
    structure[qoi]['grid_functions'] = grid_functions


# ==============================================================================
# Check whether the computed quantities cover the reported NIST ranges
# ==============================================================================
conversion_table = {'AAA': 0.003, 
                    'AA' : 0.01, 
                    'A+' : 0.02, 
                    'A'  : 0.03, 
                    'B+' : 0.07, 
                    'B'  : 0.1, 
                    'C+' : 0.18, 
                    'C'  : 0.25,
                    'D+' : 0.4, 
                    'D'  : 0.5, 
                    'E'  : 0.5}


fig, ax = plt.subplots(2,1,figsize=(10,6))
iqoi = 0
for qoi in ['energies', 'avalues']:    
    count = 0
    n_plots = structure[qoi]['dim']
    for i in range(n_plots):
        #
        # Compute sample range
        # 
        qoi_min = np.amin(structure[qoi]['computed'][:,i])
        qoi_max = np.amax(structure[qoi]['computed'][:,i])

        if structure[qoi]['nist_data'][i] is not None:
            #
            # Compute NIST range
            # 
            nist_value  = structure[qoi]['nist_data'][i]
            nist_rating = structure[qoi]['nist_ratings'][i]
            nist_relerr = conversion_table[nist_rating] 
            nist_range  = np.array([(1-nist_relerr)*nist_value, \
                                    (1+nist_relerr)*nist_value])
            #
            # Plot overlap
            # 
            ax[iqoi].plot((count,count), tuple(nist_range), 'r',
                          (count,count), (qoi_min, qoi_max), 'y', 
                          linewidth=16, alpha=0.4)
        else:
            #
            # Plot sample range
            ax[iqoi].plot((count,count), (qoi_min, qoi_max), 'y', 
                          linewidth=16, alpha=0.4)
            
            print([qoi_min, qoi_max])

        count += 1  
    # --------------------------------------------------------------------------
    # Format Figure
    # --------------------------------------------------------------------------
    #
    # Xticks
    # 
    ax[iqoi].set_xticks([i for i in range(n_plots)])
    ax[iqoi].set_xticklabels(structure[qoi]['tags'])
    ax[iqoi].set_xlim([-0.5,n_plots-0.5])
    
    #
    # Adjust y limits
    #     
    if iqoi==0:
        ax[iqoi].set_ylim([4.5e6, 4.66e6])
        ax[iqoi].set_yscale('symlog')
        ax[iqoi].set_yticks(np.linspace(4.5e6,4.65e6,4))

        ax[iqoi].set_title('Energies')
    elif iqoi==1:
        ax[iqoi].set_ylim([1e1,1e14])                            
        ax[iqoi].set_yscale('symlog')
        ax[iqoi].set_yticks([1e1, 1e3, 1e5, 1e7, 1e9, 1e11, 1e13])
        ax[iqoi].set_title('A-Values')
    iqoi += 1
fig.savefig('overlap.pdf')


# =============================================================================
# Estimate input histograms for each qoi
# =============================================================================

# -----------------------------------------------------------------------------
# Uniform Histogram on NIST values
# -----------------------------------------------------------------------------    
n_partition = 10
f_prob = np.ones(n_partition)/n_partition
for i in range(n_avalues):
    x_min, g_min, b_min, x_max, g_max, b_max = g_avalue[i].global_extrema()
    nist_min, nist_max = nist_avalue_range[i,:]
    #print(g_min, g_max)
    #print(nist_min, nist_max)
    f_grid = np.linspace(nist_avalue_range[i,0], \
                         nist_avalue_range[i,1], n_partition+1)                         
    g_avalue[i].set_output_histogram(f_grid, f_prob)    
    g_avalue[i].compute_histogram()
    
for i in range(n_energies):
    f_grid = np.linspace(nist_energy_range[i,0], \
                         nist_energy_range[i,1], n_partition+1)                         
    g_energy[i].set_output_histogram(f_grid, f_prob)    
    g_energy[i].compute_histogram()
'''

