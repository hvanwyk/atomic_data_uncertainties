#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Tasmanian
from scipy.interpolate import RegularGridInterpolator

"""
Interpolation Module
====================
Used to construct vector-valued interpolators over hyper-rectangles. An 
interpolator is a mapping of the form 

    x -> p(x), 
    
where x = [x1,...,xd] and p(x) = [p1(x),....pn(x)]. 

Created on Tue Apr 27 11:00:57 2021

@author: hans-werner
"""
class Interpolator:
    """
    """
    def __init__(self, domain, function, n_outputs, output_labels=None, 
                 resolution='auto', tol=1e-6):
        """
        Constructor
        
        Inputs:
            
            domain: (d, 2) array with the boundaries of the hyper-rectangle, 
                i.e. 
            
                domain = [x1_min, x1_max]
                         [x2_min, x2_max]
                            :       :
                         [xd_min, xd_max]
    
            function: mapping from domain to R^n used to evaluate the true 
                function.
            
            n_outputs: int, number of outputs
            
            output_labels: str/int, list of n_output output labels
            
            interpolation_type: str, type of interpolant to use
                'global-polynomial', 'local-polynomial', 'wavelet'
                
            interpolation_degree: int, 
            
            tensorization:  
                        
            resolution: int/str, mesh resolution for the interpolator
                
        """
        # Dimension of hypercube
        d = domain.shape[0]
        
        
        if type(resolution) is int:
            # Uniform resolution in all directions. 
            p = Tasmanian.makeLocalPolynomialGrid()
        pass
    
    def add_labels(self, inputs=None, outputs=None):
        """
        Labels inputs and outputs
        """
        pass
    
    
    def labels(self):
        """
        Get 
        """
        pass
    
    
    def shape(self,):
        """
        Determine 
        """
        pass
    
    
    def evaluate(self):
        """
        """
        pass
    
    
    def data(self, inputs=True, outputs=True):
        """
        Access the data used to construct the interpolator.
        """
        pass
