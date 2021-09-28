#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import Tasmanian
from scipy.interpolate import RegularGridInterpolator
import numpy as np

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
    def __init__(self, function, domain, n_outputs, output_labels=None, 
                 input_labels=None, interpolation_type='local-polynomial',
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
            
            input_labels: str/int, list of descriptions for input variables
            
            interpolation_type: str, type of interpolant to use
                'global-polynomial', 'local-polynomial', 'wavelet'
                
            interpolation_degree: int, 
            
            tensorization: str, way in which one-dimensional components are
                combined, 'sparse', 'full'.
                        
            resolution: int/str, mesh resolution for the interpolator
                
        """
        #
        # Input/Output Dimensions
        # 
        self.__n_inputs = domain.shape[0]  # dimension of hypercube
        self.__n_outputs = n_outputs  # number of outputs
        
        #
        # Construct the interpolant
        # 
        if interpolation_type=='global-polynomial':
            #
            # Global polynomials
            # 
            pass
        elif interpolation_type=='local-polynomial':
            #
            # Piecwise polynomials
            #
            pass
        elif interpolation_type=='wavelet':
            #
            # Wavelet interpolants
            # 
            pass
        if type(resolution) is int:
            # Uniform resolution in all directions. 
            grid = Tasmanian.makeLocalPolynomialGrid()
            
        elif type(resolution) is list or type(resolution) is np.ndarray:
            #
            # Different resolution in each direction
            # 
            pass
        elif type(resolution) is str:
            #
            # Automatically determine the resolution
            # 
            assert resolution=='auto', \
                'Unknown string for variable "resolution" '
        
        #
        # Store interpolation grid
        # 
        self.__grid = grid
        
        
    def add_labels(self, inputs=None, outputs=None):
        """
        Labels for inputs and outputs
        
        Inputs:
            
            inputs: str, list of length n_inputs containing the 
                labels for the input variables
                
            outputs: list of length n_outputs containing the
                labels for the output variables.
        """
        self.__labels = {}
        n_inputs,n_outputs = self.shape()  # Number of inputs and outputs
        
        #
        # Add labels for the input variables
        # 
        if inputs is not None: 
            # Check that the dimensions match 
            assert len(inputs)==n_inputs, 'The number of input labels'+\
                'does not match the number of inputs.'
                
        # Store input labels            
        self.__labels['inputs'] = inputs
            
        #
        # Add labels for the output variables
        # 
        if outputs is not None:
            # Check that dimensions match
            assert len(outputs)==n_outputs, 'The number of output labels'+\
                'does not match the number of outputs'
        
        # Store output labels
        self.__labels['outputs'] = outputs    
            
            
    def labels(self):
        """
        Returns the labels 
        """
        return self.__labels
    
    
    def shape(self):
        """
        Returns the number of inputs and outputs
        """
        return self.__n_inputs, self.__n_outputs
    
    
    def evaluate(self,x):
        """
        Evaluate the interpolant at an array of points
        
        Inputs:
            
            x: double, (n_points, n_inputs)
        """
        pass
    
    
    def data(self, inputs=True, outputs=True):
        """
        Access the data used to construct the interpolator.
        """
        pass
