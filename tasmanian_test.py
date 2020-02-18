#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:59:57 2019

@author: kyle
"""

# Internal modules
from bayesian_methods import lambdas_grid, interpolators

import numpy as np
import Tasmanian
from scipy.interpolate import RegularGridInterpolator
from time import time
import matplotlib.pyplot as plt


def time_grid_interpolators():
    """
    
    """
    dim_lmd = 3
    dim_out = 10
    
    fig, ax = plt.subplots(111)
    
    #
    # Hypercube [-1,1]^dim_lmd
    # 
    x_bnd = np.array([-np.ones(dim_lmd), np.ones(dim_lmd)]).T
    
    #
    # Vary the mesh resolution 
    # 
    for pts_per_dim in [2**i for i in np.arange(1,9)]:
        # Specify resolution
        resolution = [pts_per_dim for dummy in range(dim_lmd)]
        
        # Form the lambda grid
        X_1D, x_ravel = lambdas_grid(x_bnd, resolution)
        
        # Generate output data (random)
        grid_vals = [np.random.rand(*resolution) for dummy in range(dim_out)]
        
        # Form Interpolators
        interpolators = []
        for j in range(dim_out):
            interpolator = RegularGridInterpolator(X_1D, grid_vals[j], 
                                                   bounds_error=False)
            interpolators.append(interpolator)
        
        # Evaluate all interpolators at n points
        n = 10000
        x_eval = np.random.rand(n, dim_lmd)
        tic = time()
        for interpolator in interpolators:
            interpolator(x_eval)
        toc = time() - tic
        
        print(toc)
        
        tic = time()
        for i in range(n):
            for interpolator in interpolators:
                interpolator(x_eval[i])
        toc = time()-tic 
               
        print(' ', toc)
    
    
def time_sparse_grid_interpolators():
    """
    """
    print(dir(Tasmanian.SparseGrid))
    dim_lmd = 2
    dim_out = 10
    order = 1  # Piecewise linear
    depth = 10  # 
    grid = Tasmanian.SparseGrid()
    grid.makeGlobalGrid(dim_lmd,dim_out, depth, 'iptotal', 'leja')
    x = grid.getNeededPoints()
    n = grid.getNumPoints()
    print(x.shape)
    y = np.random.rand(n,dim_out)
    grid.loadNeededPoints(y)
    
    n = 100000
    x_eval = np.random.rand(n,dim_lmd)
    tic = time()
    for i in range(n):
        grid.evaluate(np.array(x_eval[i,:]))
    toc = time()-tic            
    print(toc)
    
    
if __name__=='__main__':
    #time_grid_interpolators()
    time_sparse_grid_interpolators()
