#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 11:34:39 2021

@author: hans-werner
"""

import matplotlib.pyplot as plt
import numpy as np



# Generate the data (100 samples at each of 10 temperature points)
Y = np.random.randn(100,10)

# Specify the temperatures
temperatures = np.array([10**i for i in range(10)])

# Plot the violinplot 
plt.violinplot(Y,positions=temperatures, widths=temperatures)

# Set the scale to 'log'
plt.xscale('log')
plt.show()

"""
Y = np.random.randn(100,10)

ax = plt.violinplot(Y,positions=[10**i for i in range(10)],
                    widths=[0.5*10**i for i in range(10)]) 
plt.xscale('log')
plt.show()

"""