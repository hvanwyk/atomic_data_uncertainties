#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 30 13:19:15 2018

@author: kyle
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
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
    
    

def log_likelihood(x, interpolators, y_bnd, likelihood_shape="gaussian"):
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

    
    
def log_posterior(x, interpolators, x_bnd, y_bnd, prior_shape="uniform", likelihood_shape="gaussian"):
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