#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 28 20:42:50 2019

@author: kyle
"""

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt

def f(x,y):
    return (x**2 + y**2)

x = np.linspace(-10, 10, 101)
y = np.linspace(-10, 10, 101)

data = f(*np.meshgrid(x, y))
print(x.shape)
interp = RegularGridInterpolator((x,y), data)

X,Y = np.meshgrid(x,y)
plt.contourf(X,Y, interp.values)
