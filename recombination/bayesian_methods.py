#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:19:15 2018

@author: kyle
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from recombination_methods import structure, get_rate, structure_dr, postprocessing_rates
import multiprocessing as mp


def log_prior(x, x_bnd, prior_shape="uniform"):
    """
    Natural logarithm of the uniform prior, which encodes our initial belief 
    concerning the distribution of x
    
    Inputs: 
    
        x: double, (d,) array of a single point in x space 
        
        x: double, (d,2) array of lower and upper bounds for the x's
        
    
    Output:
    
        log(pdf) = - sum_{i=1,...,d} log(lmd_max_i-lmd_min_i)
    """
    x_min, x_max = x_bnd[:,0], x_bnd[:,1]
    
    if prior_shape.lower() == "uniform":
        if np.all((x >= x_min)*(x <= x_max)):
            return -np.sum(np.log(x_max-x_min))
        else:
            return -np.infty
        
    elif prior_shape.lower() == "gaussian":
        
        sigma = x_bnd[:,1] - x_bnd[:,0]
        mean = (x_bnd[:,0] + x_bnd[:,1])/2
        
        gauss = np.log((1/(np.sqrt(2*np.pi)*sigma))) - ((x-mean)**2) / (2*(sigma**2))
        
        gauss[np.isnan(gauss)] = -np.infty

        return np.sum(gauss)
    
    

def log_likelihood(x, interpolators, y_bnd, likelihood_shape="uniform"):
    """
    Logarithm of the likelhood function (using uniform bounds on the error).
    
    Inputs:
    
        x: double, (d,) array of a single point in x space
        
        interpolators: RegularGridInterpolator, n-length list of interpolating
            functions.
            
        y_bnd: double, (d,2) lower and upper bounds on the output y
    """
    n = len(interpolators)
    y = np.zeros(n)
    for i, interpolator in zip(range(n),interpolators):
        y[i] = interpolator(x)
        
    if likelihood_shape.lower() == "uniform":
        
        y_min, y_max = y_bnd[:,0], y_bnd[:,1]
        
        if np.all((y>=y_min)*(y<=y_max)):
            return -np.sum(np.log(y_max-y_min))
        else:
            return -np.infty
        
    elif likelihood_shape.lower() == "gaussian":
        
        sigma = y_bnd[:,1] - y_bnd[:,0]
        mean = (y_bnd[:,0] + y_bnd[:,1])/2
        
        gauss = np.log((1/(np.sqrt(2*np.pi)*sigma))) - ((y-mean)**2) / (2*(sigma**2))
        
        gauss[np.isnan(gauss)] = -np.infty
        
        return np.sum(gauss)

    
    
def log_posterior(x, interpolators, x_bnd, y_bnd, prior_shape="uniform", likelihood_shape="uniform"):
    """
    Compute the log of the posterior density, formed from the logs of the
    likelihood and prior density functions. 
    
    Inputs: 
    
        x: double, (d,) array of a single point in x-space
        
        interpolators: RegularGridInterpolator, n-list of regular grid 
            interpolating functions computed using scipy.RegularGridInterpolator
        
        x_bnd: (d, 2) array of lower and upper bounds for the inputs x 
        
        y_bnd: (n, 2) array of lower and upper bounds for the outputs y
        
    
    Output:
    
        log(post) = log(prior) + log(likelihood) 
        
    NOTE: If you want to be a bit more flexible in the prior and likelihood 
        distributions, you can use flags like 'likelihood_type' or 'prior_type'
        
    """
    lp = log_prior(x, x_bnd, prior_shape)
    """
    if not np.isfinite(lp):
        return -np.infty
    else:
    """

    return lp + log_likelihood(x, interpolators, y_bnd, likelihood_shape)

def lambdas_grid(x_bnd, x_res):
    #
    # One dimensional grid in each direction
    # 
    d = x_bnd.shape[0]
    X_1D = []
    for i in range(d):
        X_1D.append(np.linspace(x_bnd[i,0],x_bnd[i,1],x_res[i]))
    
    #
    # Multi-dimensional Grid, using meshgrid (ij indexing)
    #
    X = np.meshgrid(*X_1D, indexing='ij')
        
    #
    # Unravel grid into (n_points, d) array
    # 
    x_ravel = np.array([x.ravel() for x in X]).T
    
    return X_1D, x_ravel

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
        data = structure(ion=ion, lambdas=x)
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
        rates = [pool.apply(get_rate, args=(ion, x)) for x in x_ravel]
        rates = np.array(rates)
    
    else:
        for i in range(n_points):
            x = x_ravel[i,:]
            E, E_nist, delta_E = structure(ion, lambdas=x)
            structure_dr(ion, lambdas=x)
            rates[i, :] = postprocessing_rates(ion, E, E_nist, lambdas=x)[1]
    
    
    Rates = [np.reshape(rates[:,j], x_res) for j in range(n_rates)]
    return Rates
    
def interpolators(X_1D, grid_vals):
    """
    Constructs a list of RegularGridInterpolators, one for each output
    
    """
    #
    # Form the interpolators
    # 
    n = len(grid_vals)
    interpolators = []
    for j in range(n):
        interpolators.append(RegularGridInterpolator(X_1D, grid_vals[j], bounds_error=False))
    
    return interpolators
    
    
def sample_from_histogram(H, edges, n_samples):
    """
    Generate a random sample from a histogram defined over a hypercube.
    We proceed by conditional sampling, i.e.
    
    1. Pick a random bin based on the histogram values.
    2. Choose a random point inside each subcube. 
    
    """
    # 
    # Generate uniformly distributed random variables
    # 
    U = np.random.rand(n_samples)
    
    # 
    # Normalized cumulative sum from 0 to 1  
    # 
    CDF = np.cumsum(H.ravel())
    CDF = CDF/CDF[-1]
    
    #
    # Select bins randomly.
    # 
    bins = np.searchsorted(CDF, U)
    
    #
    # Determine the bin coordinates within the multi-dimensional grid 
    # 
    n_cells = tuple([len(edge)-1 for edge in edges])
    i_bins = np.array(np.unravel_index(bins, n_cells)).T
    
    d = i_bins.shape[1]
    V = np.random.rand(n_samples, d)
    
    X = np.empty((n_samples,d))
    for j in range(d):
        dxj = np.diff(edges[j])
        i_binsj = i_bins[:,j]
        X[:,j] = edges[j][i_binsj] + V[:,j]*dxj[i_binsj]
    return X