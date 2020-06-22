#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:21:10 2020

@author: kyle
"""

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


#---------------------------------------------------------------------------
"""
Generating distribution of lambda parameters.
"""

# Min and max lambda parameter values to use in generating interpolators
lambda_min = 0.8
lambda_max = 1.2
x_bnd = []
for i in range(n_lambdas):
    x_bnd.append([lambda_min, lambda_max])
x_bnd = np.array(x_bnd)


grid_size_structure = 5 # Number of points per lambda used in autostructure runs to generate interpolators for level energies
x_res_structure = np.array([grid_size_structure]*n_lambdas)

nist_cutoff = 0.05 # Tolerance for computed energies relative to NIST values
potential_type = 1 # Central field potential used by autostructure. +1 for Thomas-Fermi potential, -1 for Slater-type

grid_size_rates = 5 # Number of points per lambda used in full R-matrix runs to generate rate interpolators
x_res_rates = np.array([grid_size_rates]*n_lambdas)


lambda_samples_tf = make_lambda_distribution(ion=ion, x_bnd=x_bnd, x_res=x_res_structure, 
                                                 n_lambdas=n_lambdas, nmax=nmax, nist_cutoff=nist_cutoff, potential_type=1)

T, rate_samples_tf, df_base = make_rates_distribution(ion=ion, lambda_samples=lambda_samples_tf, 
                                                          x_bnd=x_bnd, x_res=x_res_rates, n_lambdas=n_lambdas, potential_type=1)

# Computes level energies over lambda grid, generates lambda distribution, and then samples from distribution using MCMC
lambda_samples_slater = make_lambda_distribution(ion=ion, x_bnd=x_bnd, x_res=x_res_structure, 
                                                 n_lambdas=n_lambdas, nmax=nmax, nist_cutoff=nist_cutoff, potential_type=-1)


# Compute excitation rates for each sample set of lambdas using interpolators
T, rate_samples_slater, df_base = make_rates_distribution(ion=ion, lambda_samples=lambda_samples_slater, 
                                                          x_bnd=x_bnd, x_res=x_res_rates, n_lambdas=n_lambdas, potential_type=-1)


lambda_samples = np.array([lambda_samples_tf, lambda_samples_slater])
np.save(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/lambda_samples.npy", arr=lambda_samples)

np.save(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rate_samples.npy", np.array([T, rate_samples_tf, rate_samples_slater]))