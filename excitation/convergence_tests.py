
import sys
if ".." not in sys.path:
    sys.path.append("..")

from state import State
import numpy as np
from lambdas import make_lambda_distribution, make_rates_distribution


atom = "o" # nucleus
seq = "he" # isoelectronic sequence
corex = "1-2" # Core excitation shells

ion = State(atom, seq, corex)

nmax = 3 # Highest n shell used in R-matrix input file
n_lambdas = 2 # Number of lambda parameters to vary, starting with 2s

# Min and max lambda parameter values to use in generating interpolators
lambda_min = 0.2
lambda_max = 1.9
x_bnd = []
for i in range(n_lambdas):
    x_bnd.append([lambda_min, lambda_max])
x_bnd = np.array(x_bnd)

for i in range(6):
    grid_size_structure = 2**i # Number of points per lambda used in autostructure runs to generate interpolators for level energies
    x_res_structure = np.array([grid_size_structure]*n_lambdas)
    print(x_res_structure)

    

