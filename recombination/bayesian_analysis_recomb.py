#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:03:40 2018

@author: kyle
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner 
from recombination_methods import State
from bayesian_methods import log_uniform_prior, log_likelihood, log_posterior, interpolators_from_scratch
import os
import time
from graphing import graph_rates_from_file
from sklearn.svm import SVC

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
    

def make_distribution(atom, seq, shell, nist_cutoff=0.05, n_lambdas=2, n_walkers=10, n_steps=1000, prior_shape="uniform", likelihood_shape="uniform"):
    ion = State(atom, seq, shell)
    
    T = ion.dielectronic_recomb()[0]
    n_points = len(T)
    
    # Interval endpoints for each input component
    x_bnd = np.array([[0.8,1.2],[0.8,1.2]])
    
    # Resolution in each dimension
    x_res = np.array([5,5])
    
    # NIST data
    y_nist = ion.structure(lambdas=[])[1]
    # Construct interpolators
    n_energies = len(y_nist)
    err_interpolators, erg_interpolators, rate_interpolators = interpolators_from_scratch(ion, x_bnd, x_res, n_energies)
    
    
    # =============================================================================
    # Run MCMC Code
    # =============================================================================
    #
    # Likelihood is based on error*100% error in each component 
    # 
    y_bnd = np.zeros((n_energies, n_lambdas))
    for i in range(n_energies):
        y_bnd[i, :] = -1, 1
        
    y_bnd *= nist_cutoff
    
    
    #
    # Specify starting points for each Markov chain (in a tight ball around optimum)
    #
    pos = [np.array([1 for i in range(n_lambdas)]) + 1e-4*np.random.randn(n_lambdas) for i in range(n_walkers)]
    
    #
    # Initialize the sampler
    #
    sampler = emcee.EnsembleSampler(n_walkers, n_lambdas, log_posterior, 
                                    args=(err_interpolators, x_bnd, y_bnd, prior_shape, likelihood_shape))
    #
    # Run the MCMC routine
    #
    sampler.run_mcmc(pos, n_steps);
    
    #
    # The sampler.chain has shape (n_walkers, n_steps, n_dim)
    # Reshape array of samples (but throw away initial sample)
    #
    samples = sampler.chain[:, 50:, :].reshape((-1, n_lambdas))
    
    
    # -----------------------------------------------------------------------------
    # Visualize the posterior density
    # -----------------------------------------------------------------------------
    
    
    corner.corner(samples, labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
    plt.show()
    
    # -----------------------------------------------------------------------------
    # Plug the Monte Carlo samples from the posterior into the function and comput
    # a histogram for the output. Make sure the point [0,0,0] has nonzero density   
    # -----------------------------------------------------------------------------
    
    n_samples = samples.shape[0]
    
    
    rate_samples = np.zeros((n_samples,n_points))
    E_samples = np.zeros((n_samples, n_energies))
    
    for i in range(n_samples):
        for j in range(n_energies):
            E_samples[i,j] = erg_interpolators[j](samples[i]) 
        for j in range(n_points):
            rate_samples[i,j] = rate_interpolators[j](samples[i]) 
        
    
       
    corner.corner(E_samples[:,1:], labels=[f"$E_{i+1}$" for i in range(n_energies-1)], truths=[0 for i in range(n_energies-1)])
    plt.show()
    
    
    # =============================================================================
    # Compute histogram
    # =============================================================================
    """
    H, edges = np.histogramdd(samples, bins=(20,20,20), normed=True)
    CH = np.cumsum(H.ravel())
    plt.plot(CH)
    plt.show()
    """
    
    # =============================================================================
    # Sample from histogram
    # =============================================================================
    """
    y = sample_from_histogram(H, edges, 10000)
    
    corner.corner(y, labels=["$\lambda_1$", "$\lambda_2$"], truths=[1,1] )
    plt.show()
    """
    
    
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"    
    fname = direc + f"rates_{shell}_{int(nist_cutoff*100)}_{prior_shape}P_{likelihood_shape}L.npy"
    #np.savetxt(fname, rate_samples, header=" ".join([str(x) for x in T]), comments="")
    out_data = np.array([T, rate_samples])
    np.save(fname, arr=out_data)
    return samples, E_samples, rate_samples
    

if __name__ == "__main__":
    
    # Uncomment this section to use command line arguments to pass atom, seq, shell
    """
    if (len(sys.argv) < 3):
        sys.exit("Correct usage: python bayesian_analysis.py <atom> <isoelectronic_sequence> <core-ex shell (e.g. 2-2)>")
    
    atom = sys.argv[1]
    seq = sys.argv[2]
    shell = sys.argv[3]
    if len(str(shell)) != 3:
        sys.exit("Correct usage: python bayesian_analysis.py <atom> <isoelectronic_sequence> <core-ex shell (e.g. 2-2)>")
    """

    start = time.time()
    
    atom = "o"
    seq = "be"
    shell = "2-2"
    nist_cutoff=0.05
    prior_shape="gaussian"
    likelihood_shape="gaussian"
    ion = State(atom, seq, shell)
    
    make_distribution(atom=atom, seq=seq, shell=shell, nist_cutoff=nist_cutoff, prior_shape=prior_shape, likelihood_shape=likelihood_shape)
    
    end = time.time()
    print(f"Runtime: {int(end-start)}s")
