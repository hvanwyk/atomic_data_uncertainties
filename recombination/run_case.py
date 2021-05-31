#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 11:09:08 2021

@author: loch
"""
import numpy as np
import sys
import matplotlib.pyplot as plt
import emcee
import corner 
import scipy.stats as stats
import math

from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import log_posterior, interpolators, lambdas_grid
import time
from graphing import graph_rates_from_file
from scipy.optimize import minimize
import multiprocessing as mp
from mpl_toolkits import mplot3d
from lambda_variation import energy_grid, rates_grid,  lambda_distribution, energy_distribution, rates_distribution, energy_optimization


def run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps):
  
   lambdas_file = direc + "lambdas" + file_name_common+".npy"   
   print('Evaluating lambda distribution functions using MCMC')
   lambda_samples, Err, Erg, err_interpolators,y_nist = lambda_distribution(ion,cent_pot=cent_pot,x_bnd=x_bnd, x_res=x_res, nist_cutoff=nist_cutoff, 
                      prior_shape=prior_shape, n_walkers=n_walkers, n_steps=n_steps,
                      likelihood_shape=likelihood_shape, outfile=lambdas_file)
    


#   xv=X_1D[0][:]
#   yv=X_1D[1][:]
#   xvals,yvals=np.meshgrid(xv,yv)
#   indx=1

#   zvals=Erg[indx][:][:]

#   nist_vals= np.full((grid_resolution,grid_resolution),y_nist[indx])

#   ax = plt.axes(projection='3d')
##ax.set_zlim3d(0.6,0.8)
#   ax.set_xlabel('lambda_2s')
#   ax.set_ylabel('lambda_2p')
#   ax.set_zlabel('Energy')
#   ax.plot_wireframe(xvals,yvals,zvals,color='r',label='Input lambda grid')
#   ax.plot_wireframe(xvals,yvals,nist_vals,color='b',label='NIST Value')
#   ax.plot_wireframe(xvals,yvals,nist_vals*(1.0+nist_cutoff),color='g',label='NIST Value+3%')
#   ax.plot_wireframe(xvals,yvals,nist_vals*(1.0-nist_cutoff),color='g',label='NIST Value-3%')
#   ax.legend()


   print('Calculating Energy distribution functions')
   E_samples, Err, Erg=energy_distribution(ion, lambda_samples, x_bnd, x_res, cent_pot)

#   fig, axs = plt.subplots(1, 2, sharey=True, tight_layout=True)
#   axs[0].set_xlabel('E_1')
#   axs[0].set_ylabel('Frequency')
#   axs[0].hist(E_samples[:,1], bins=20,density=True)
#   mu_1 = y_nist[1]
#   sigma_1 = mu_1*nist_cutoff
#   x_1 = np.linspace(mu_1 - 3*sigma_1, mu_1 + 3*sigma_1, 100)
#   axs[0].plot(x_1, stats.norm.pdf(x_1, mu_1, sigma_1))

#   axs[1].set_xlabel('E_4')
#   axs[1].hist(E_samples[:,2], bins=20,density=True)
#   mu_2 = y_nist[2]
#   sigma_2 = mu_2*nist_cutoff
#   x_2 = np.linspace(mu_2 - 3*sigma_2, mu_2 + 3*sigma_2, 100)
#   axs[1].plot(x_2, stats.norm.pdf(x_2, mu_2, sigma_2))
 
#Rates = rates_grid(ion, x_ravel, x_res)    
    
   rates_file = direc+"rates"+file_name_common+".npy"
   print('Calculating rate distributions using cent_pot=',cent_pot)
   T, rate_samples = rates_distribution(ion, lambda_samples, x_bnd, x_res, cent_pot,outfile=rates_file)
    
   graph_file=direc + "rates"+file_name_common + ".png"
   graph_rates_from_file(ion, infile = rates_file, outfile=graph_file)

   rate_avg = np.empty(19)
   rate_std = np.empty(19)
   rate_percent = np.empty(19)

   for i in range(len(rate_avg)):
       rate_avg[i]=np.average(rate_samples[:,i])
       rate_std[i]=np.std(rate_samples[:,i])
       if(rate_avg[i] > 0.0):
#        print('i=',i)
           rate_percent[i]=rate_std[i]/rate_avg[i]*100
  
        
#   fig, axs3 = plt.subplots(2,1)
#   axs3[0].set_xscale("log")
#   axs3[0].set_yscale("log")
#   axs3[0].plot(T,rate_avg)
#   axs3[0].errorbar(T,rate_avg,yerr=rate_std*2.0,fmt='ro',ecolor='black')
#   axs3[1].set_xscale("log")
#   axs3[1].plot(T,rate_percent)
    

#   fig, axs4 = plt.subplots(1, 2, sharey=True, tight_layout=True)
#   axs4[0].hist(rate_samples[:,1], bins=100,density=True)
#   axs4[1].hist(rate_samples[:,18], bins=100,density=True)

   data = {"T": T,
    "nist_vals": y_nist,
    "lambda_samples": lambda_samples,
    "E_samples":E_samples,
    "rate_samples":rate_samples,
    "rate_avg":rate_avg,
    "rate_std":rate_std,
    "rate_percent":rate_percent}
  

   return data
