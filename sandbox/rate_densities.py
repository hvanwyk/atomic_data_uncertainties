#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 16:22:03 2021

@author: hans-werner
"""
import pickle
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, BarycentricInterpolator
from sklearn.neighbors import KernelDensity

data = pickle.load(open('be-like-oxygen.pkl','rb'))

# Interpolate one line
t = np.log(data['be_like_o_pos']['T'])
y = np.log(data['be_like_o_pos']['rate_samples'])


"""
y_cs = CubicSpline(t,y.T)
y_gp = BarycentricInterpolator(t,y.T)

print(y_cs.c.shape)
print(y.shape)
#kde = KernelDensity(kernel='gaussian',bandwidth=0.5)

#kde.fit(y)
#r_range = np.linspace(np.min(y[:,0]),np.max(y[:,0]),101)

#np.exp(kde.score_samples(r_range))
tt = np.linspace(t[0],t[-1],100)

plt.plot(t,y.T,'k.',tt,y_cs(tt),'r', linewidth=0.1)
plt.show()
"""