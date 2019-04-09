#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:03:40 2018

@author: kyle
"""
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner 
from recombination_methods import State, structure, structure_dr, postprocessing_rates
from bayesian_methods import log_posterior, lambdas_grid, energy_grid, rates_grid, interpolators
import time
from graphing import graph_rates_from_file

"""
This code snippet illustrates how to wrap a Bayesian sampler around
the mapping from lambda's to the quantities to be used for calibration. 
We want to use the EMCEE package for calibration, which is based on Monte Carlo
sampling. To accommodate a possibly large sample size, we first build a 
surrogate model by creating a linear interpolant over a regular lambda grid.   
In particular, we want to achieve the following.

    1. Interpolate a multivariate mapping f from R^d to R^n, i.e. 
        Mapping x=[x1,x2,...,xd] to y=[f1(x), f2(x), ..., fn(x)]
        based on a regular x-grid. In this implementation, I have split
        the forward mapping into two separate functions (you can combine them
        if you want):
        
            a. fake_autostructure_and_postprocessing works out the quantities
                interest (in this case energies). 
                
            b. fake error function takes the output from the autostructure run
                and computes the deviation from observed values. 
            
    2. Use the python package emcee to generate a MCMC sample of the 
        x vectors in R^d from the posterior distribution, p_post, 
        where 
        
            p_post(x) ~ p_prior(x) p_likelihood(x)
            
    3. Construct a multivariable histogram from the prior
    
    4. Sample from the histogram.
"""


# =============================================================================
# 1. Construct Interpolators
# =============================================================================

def variable_shift_method(ion, n_rates):
    
    E, E_nist, E_shift = structure(ion, method="shift")
    structure_dr(ion, method="shift")
    
    
    T = postprocessing_rates(ion, E, E_nist, method="shift")[0]
    rates = np.zeros((n_rates, T.size))
    
    plt.figure()
    for i in range(rates.shape[0]):
        #shifts = (np.random.rand(*E.shape) * 2 - 1)*2
        shifts = np.full(E.shape, i*0.01)
        rates[i,:] = postprocessing_rates(ion, E, E_nist, method="shift", shift=shifts)[1]
        plt.semilogx(T, rates[i, :])
    plt.show()
    
    std = np.std(rates, axis=0)
    plt.figure()
    plt.loglog(T, std)

    
    
        

