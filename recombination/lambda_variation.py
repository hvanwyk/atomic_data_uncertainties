#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 16:20:47 2019

@author: kyle
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
import emcee
import corner 
from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import log_posterior, interpolators, lambdas_grid
import time
from graphing import graph_rates_from_file
from scipy.optimize import minimize
import multiprocessing as mp

def energy_grid(ion, x_ravel, x_res):
    #
    # Generate error data 
    # 
    n = structure(ion=ion)[0].shape[0]
    n_points = x_ravel.shape[0]
    err = np.empty((n_points,n))
    erg = np.empty((n_points, n))
    
    
    for i in range(n_points):
        x = x_ravel[i,:]
        if ion.isoelec_seq == "he":
            x=np.r_[1.0,x]
        potential = np.random.choice([-1, 1])
        data = structure(ion=ion, lambdas=x, potential=potential)
        err[i,:] = data[2]
        erg[i, :] = data[0]
    
    Err = [np.reshape(err[:,j], x_res) for j in range(n)]
    Erg = [np.reshape(erg[:,j], x_res) for j in range(n)]

    return Err, Erg

def rates_grid(ion, x_ravel, x_res, parallel=False):
    
    E, E_nist, shift = structure(ion)
    structure_dr(ion)
    T = postprocessing_rates(ion, E, E_nist)[0]
    n_rates = T.size
    n_points = x_ravel.shape[0]
    
    rates = np.empty((n_points, n_rates))
    if parallel:
        n_procs = mp.cpu_count()
        pool = mp.Pool(processes=n_procs)
        if ion.isoelec_seq == "he":
            rates = [pool.apply(get_rate, args=(ion, np.r_[1.0,x])) for x in x_ravel]
        else:
            rates = [pool.apply(get_rate, args=(ion, x)) for x in x_ravel]
        rates = np.array(rates)
    
    else:
        for i in range(n_points):
            x = x_ravel[i,:]
            if ion.isoelec_seq == "he":
                x=np.r_[1.0,x]
            E, E_nist, delta_E = structure(ion, lambdas=x)
            structure_dr(ion, lambdas=x)
            rates[i, :] = postprocessing_rates(ion, E, E_nist, lambdas=x)[1]
    
    
    Rates = [np.reshape(rates[:,j], x_res) for j in range(n_rates)]
    return Rates

def lambda_distribution(ion, x_bnd, x_res, nist_cutoff=0.05, n_lambdas=2, n_walkers=10, n_steps=1000, 
                      prior_shape="uniform", likelihood_shape="uniform", plot=True, outfile=None):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    # NIST data
    y_nist = structure(ion, lambdas=[])[1]
    
    # Construct interpolators
    n_energies = y_nist.size
    
    Err, Erg = energy_grid(ion, x_ravel, x_res)
    
    err_interpolators = interpolators(X_1D, Err)
    
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
    lambda_samples = sampler.chain[:, 50:, :].reshape((-1, n_lambdas))
    
    
    # -----------------------------------------------------------------------------
    # Visualize the posterior density
    # -----------------------------------------------------------------------------
    
    if plot:
        plt.figure()
        corner.corner(lambda_samples, labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
        #plt.show()
        plt.savefig("lambdas.jpg")
        
    
    # -----------------------------------------------------------------------------
    # Plug the Monte Carlo samples from the posterior into the function and comput
    # a histogram for the output. Make sure the point [0,0,0] has nonzero density   
    # -----------------------------------------------------------------------------
    

    
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
    
    if outfile is not None:
        np.save(outfile, arr=lambda_samples)
    
    return lambda_samples
    
def energy_distribution(ion, lambda_samples, x_bnd, x_res, plot=True, outfile=None):
    
    y_nist = structure(ion, lambdas=[])[1]
    n_energies = y_nist.size
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    Err, Erg = energy_grid(ion, x_ravel, x_res)
    energy_interpolators = interpolators(X_1D, Erg)
    
    n_samples = lambda_samples.shape[0]
    E_samples = np.zeros((n_samples, n_energies))

    for i in range(n_samples):
        for j in range(n_energies):
            E_samples[i,j] = energy_interpolators[j](lambda_samples[i])         
        
    if plot:
        plt.figure()
        corner.corner(E_samples[:,1:], labels=[f"$E_{i+1}$" for i in range(n_energies-1)], truths=[0 for i in range(n_energies-1)])
        #plt.show()
        plt.savefig("Energies.jpg")
    
    return E_samples

"""
def xsec_data(lambda_samples):
    E, E_nist, E_shift = structure(ion)
    np.random.shuffle(lambda_samples)
    lambdas = lambda_samples[:100,:]
    xsec_samples = np.zeros((lambdas.shape[0], E.size))
    for i in range(lambdas.shape[0]):
        structure_xsec(ion)
        xsec_samples[i,:] = postprocessing_rates(ion, lambdas=lambdas[i,:], xsec=True)[1]
"""

def rates_distribution(ion, lambda_samples, x_bnd, x_res, outfile=None):
        
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    E, E_nist, E_shift = structure(ion)
    structure_dr(ion)
    T = postprocessing_rates(ion, E, E_nist)[0]
    n_points = T.size
    n_samples = lambda_samples.shape[0]
    
    Rates = rates_grid(ion, x_ravel, x_res)
    
    rate_interpolators = interpolators(X_1D, Rates)
    
    rate_samples = np.zeros((n_samples,n_points))
    
    for i in range(n_samples):
        for j in range(n_points):
            rate_samples[i,j] = rate_interpolators[j](lambda_samples[i]) 
        
    if outfile:
        np.save(outfile, np.array([T, rate_samples]))
        
    return T, rate_samples
    

def energy_optimization(ion, lambda_samples, x_bnd, x_res):    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    Err, Erg = energy_grid(ion, x_ravel, x_res)
    energy_interpolators = interpolators(X_1D, Erg)
    def func(x): 
        mean = 0
        for inty in energy_interpolators[1:]:
            mean += (inty(x)**2.0)
        return mean
    res = minimize(func, x0=np.array([1.0, 1.0]), bounds=[(0, None), (0, None)])
    lambda_opt = res.x
    E_min = []
    for inty in energy_interpolators[1:]:
        E_min.append(inty(lambda_opt))
        
    return np.array(lambda_opt), np.array(E_min)

def main():
    ion = State(atom, seq, shell)
    
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"    
    file_name_common = f"_{int(nist_cutoff*100)}_{prior_shape}P_{likelihood_shape}L"
    
    # Interval endpoints for each input component
    x_bnd = np.array([[0.8,1.2],[0.8,1.2]])
    
    # Resolution in each dimension
    grid_resolution = 5
    x_res = np.array([grid_resolution, grid_resolution])
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    lambdas_file = direc + "lambdas" + file_name_common+".npy"
    lambda_samples = lambda_distribution(ion, x_bnd=x_bnd, x_res=x_res, nist_cutoff=nist_cutoff, prior_shape=prior_shape, 
                      likelihood_shape=likelihood_shape, outfile=lambdas_file)
    
    """
    energy_distribution(ion, lambda_samples, x_bnd, x_res)
    
    Rates = rates_grid(ion, x_ravel, x_res)    

    rates_file = direc+"rates"+file_name_common+".npy"
    T, rate_samples = rates_distribution(ion, lambda_samples, x_bnd, x_res, outfile=rates_file)
    
    graph_file=direc + "rates"+file_name_common + ".png"
    graph_rates_from_file(ion, infile = rates_file, outfile=graph_file)
    
    end = time.time()
    print(f"Runtime: {int(end-start)}s")
    """
    
    
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
    seq = "li"
    shell = "2-2"
    ion = State(atom, seq, shell)
    
    nist_cutoff=0.01
    prior_shape="uniform"
    likelihood_shape="uniform"
    
    main()
    