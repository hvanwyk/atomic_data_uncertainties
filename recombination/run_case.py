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


def run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps,emax,up_dir,nist_shift=False):
  
   print('NIST shift=',nist_shift)
   lambdas_file = direc + "lambdas" + file_name_common+".npy"   
   print('Evaluating lambda distribution functions using MCMC')
   lambda_samples, Err, Erg, err_interpolators,y_nist = lambda_distribution(ion,up_dir,cent_pot=cent_pot,x_bnd=x_bnd, x_res=x_res, nist_cutoff=nist_cutoff, 
                      prior_shape=prior_shape, n_walkers=n_walkers, n_steps=n_steps,
                      likelihood_shape=likelihood_shape, outfile=lambdas_file,emax=emax)
    

   print('Calculating Energy distribution functions')
   E_samples, Err, Erg=energy_distribution(ion, up_dir,emax, lambda_samples, x_bnd, x_res, cent_pot)

    
   rates_file = direc+"rates"+file_name_common+".npy"
   print('Calculating rate distributions using cent_pot=',cent_pot,'emax=',emax,'nist_shift=',nist_shift)
   T, rate_samples, rates = rates_distribution(ion, up_dir,emax, lambda_samples, x_bnd, x_res, cent_pot,outfile=rates_file,nist_shift=nist_shift)
    

   rate_avg = np.empty(19)
   rate_std = np.empty(19)
   rate_percent = np.empty(19)

   for i in range(len(rate_avg)):
       rate_avg[i]=np.average(rate_samples[:,i])
       rate_std[i]=np.std(rate_samples[:,i])
       if(rate_avg[i] > 0.0):
#        print('i=',i)
           rate_percent[i]=rate_std[i]/rate_avg[i]*100
  


   data = {"T": T,
    "nist_vals": y_nist,
    "Erg":Erg,
    "X_1D":X_1D,
    "lambda_samples": lambda_samples,
    "E_samples":E_samples,
    "Err":Err,
    "rate_samples":rate_samples,
    "rate_avg":rate_avg,
    "rate_std":rate_std,
    "rate_percent":rate_percent}
  

   return data
