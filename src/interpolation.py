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
    def __init__(self, domain, function):
        """
        Constructor
        """
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
    
    
    def eval(self):
        """
        """
        pass
    
    
    def data(self, inputs=True, outputs=True):
        """
        Access the data used to construct the interpolator.
        """
        pass
