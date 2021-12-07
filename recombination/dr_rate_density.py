#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 12:38:01 2021

@author: hans-werner
"""
import numpy as np
import matplotlib.pyplot as plt
import pickle
import corner 
import sys

from scipy.interpolate import CubicSpline, PPoly
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV


class DRRateDensity(object):
    """
    Class for estimating densities of the dielectronic recombination rates 
    as a function of temperature using kernel density estimators
    """
    def __init__(self, ion, seq, shell):
        pass
    
        self.__bandwidths = []
        
    
    def fit(self, T, r, local=True, verbose=False, lowest_rate=1e-50, 
            bandwidth_range=10**np.linspace(-3,1,11)):
        """
        Fit the distribution 
        
        Inputs:
            
            T: double, m-vector of temperatures
            
            r: double, (m, n) array of rates 
            
            local: 
        """
        m = len(T)
        n = r.shape[1]
        
        if verbose:
            print('Fitting Kernel Density Estimators to Rates Data:')
            print('------------------------------------------------')
            print(f'Temperature Range: [{np.min(T)},{np.max(T)}]')
            print('Number of Temperature Points: ', m)
            print('Number of Samples: ', n)
            
        # 
        # Replace 0 rate with small number
        #
        if verbose:
            print(f'WARNING: {np.sum(r==0)} rates are zero')
            print(f'Replacing zero rates with default min, r={lowest_rate}.')
        r[r==0] = 1e-50  
    
        #
        # Convert the temperatures and rates to logarithmic scale
        #
        logT = np.log(T)
        logr = np.log(r)
        
        #
        # Fit 
        # 
        if local:
            #
            # Fit piecewise cubic spline to data 
            #
            p = CubicSpline(logT,logr)
            
            if verbose:
                print('Fitting cubic splines to data')
            
            #
            # Fit densities to coefficients on each sub-interval
            # 
            kde = []
            for i in range(m-1):
                
                if verbose:
                    print(f'.. {i}th interval [{T[i],T[i+1]}]')
                    
                # Get coefficients for the ith sub-interval
                ci = p.c[:,i,:]    
                
                #
                # Determine the optimal bandwidth using a grid search              
                #
                grid = GridSearchCV(KernelDensity(kernel='gaussian'),
                                    {'bandwidth':bandwidth_range})
                grid.fit(ci.T)
                bw = grid.best_params_['bandwidth']
                
                #
                # Define
                # 
                kdei = KernelDensity(bandwidth=bw)
                kdei.fit(ci.T)
        else:
            #
            # Fit global sum of exponentials to data
            # 
            pass 
        
    
    def sample(self, T, n):
        """
        Generate n samples from the rate distribution at temperature T
        
        Inputs:
            
            T: double, temperature (in Kelvin) within range
            
            n: int, number of samples. 
        """
        pass


    def plot(self):
        """
        Plots the density function over the given range of temperatures
        """
        pass


    def bandwidth(self):
        """
        Returns the bandwidths used to estimate the coefficient distributions
        """
        return self.__bandwidths
    
    