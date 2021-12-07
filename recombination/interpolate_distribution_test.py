#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec  5 23:06:36 2021

@author: loch
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 22:21:05 2021

@author: loch
"""


import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
from mpl_toolkits.mplot3d import Axes3D

import scipy
from scipy import stats
#from lmfit.models import SkewedGaussianModel, LognormalModel
from scipy.optimize import curve_fit

if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")

import pickle

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['text.usetex'] = True #Let TeX do the typsetting
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',r'\sansmath']
#Force sans-serif math mode
plt.rcParams['font.family'] = 'sans-serif' # ... for regular text
plt.rcParams['font.sans-serif'] = 'Helvetica' # Choose a nice font here



#%%
#Fitting rate coefficients for each lambda pair

#def drfit(c1,c2,c3,c4,c5,c6,e1,e2,e3,e4,e5,e6,Te):
#    return c1*exp(-e1/te)+c2*exp(-e2/te)+c3*exp(-e3/te)+c4*exp(-e4/te)+c5*exp(-e4/te)+c6*exp(-e6/te)

def drfit(te,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5):
    rate=np.zeros(len(te))
    for i in range(len(te)):
        rate[i]=(c1*np.exp(-e1/te[i])+c2*np.exp(-e2/te[i])+c3*np.exp(-e3/te[i])+c4*np.exp(-e4/te[i])+c5*np.exp(-e5/te[i]))*te[i]**(-1.5)
    return rate


def drcalc(te,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5):
    rate=(c1*np.exp(-e1/te)+c2*np.exp(-e2/te)+c3*np.exp(-e3/te)+c4*np.exp(-e4/te)+c5*np.exp(-e5/te))*te**(-1.5)
    return rate


#%%    
data=pickle.load(open('be_like_o_fe_0.4_1.6_40grid_nist_0.03_uniform_gaussian.pkl','rb'))
   
best_vals_store=np.zeros((len(data['be_like_o_pos']['rate_samples'][:,0]),10))
te=data['be_like_o_pos']['T']/11604.
init_vals=[1.e-10,1.e-10,1.e-10,1.e-10,1.e-10,0.1,0.2,0.5,1.0,10.0]

npts=100
fig1, axs1 = plt.subplots(1,1)
for i in range(npts):
    rate_lambda_pair=data['be_like_o_pos']['rate_samples'][i,:]
    best_vals, covar = curve_fit(drfit, te, rate_lambda_pair, p0=init_vals,method='lm')
    best_vals_store[i,:]=best_vals
    c1=best_vals_store[i,0]
    c2=best_vals_store[i,1]
    c3=best_vals_store[i,2]
    c4=best_vals_store[i,3]
    c5=best_vals_store[i,4]
    e1=best_vals_store[i,5]
    e2=best_vals_store[i,6]
    e3=best_vals_store[i,7]
    e4=best_vals_store[i,8]
    e5=best_vals_store[i,9]

    result=drfit(te,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5)
    axs1.set_xlabel('Electron temperature (eV)')
    axs1.set_ylabel('Rate coefficient')
    axs1.loglog(te,rate_lambda_pair,'o', color='green')
    axs1.loglog(te,result)

plt.show()

#%%    
    
result_all=np.zeros(npts)
te_interp=1000.0
for i in range(npts):
    c1=best_vals_store[i,0]
    c2=best_vals_store[i,1]
    c3=best_vals_store[i,2]
    c4=best_vals_store[i,3]
    c5=best_vals_store[i,4]
    e1=best_vals_store[i,5]
    e2=best_vals_store[i,6]
    e3=best_vals_store[i,7]
    e4=best_vals_store[i,8]
    e5=best_vals_store[i,9]
    result_all[i]=drcalc(te_interp,c1,c2,c3,c4,c5,e1,e2,e3,e4,e5)


fig2, axs2 = plt.subplots(1,1)
axs1.set_ylabel('Number of values')
axs1.set_xlabel('Rate coefficient')
axs2.hist(np.log10(data['be_like_o_pos']['rate_samples'][:,14]),bins=1000,density=True, ls='dotted',alpha=0.5,color='red',label='689 eV')
axs2.hist(np.log10(data['be_like_o_pos']['rate_samples'][:,15]),bins=1000,density=True, ls='dotted',alpha=0.5,color='blue',label='1379 eV')
axs2.hist(np.log10(result_all),bins=10,density=True, ls='dotted',alpha=0.5,color='green',label='1000 eV')
axs2.legend()
plt.show()
   