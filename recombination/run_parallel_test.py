#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 20 20:29:34 2021

@author: loch
"""

#%%
import numpy as np
import sys
import matplotlib.pyplot as plt
import corner 
from scipy import stats


from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import log_posterior, interpolators, lambdas_grid
import time
import pickle
import os
from utilities import create_directories, get_nist_energy, read_levels, compare_to_nist
import multiprocessing as mp
from functools import partial
from run_case_par import run_case_par

from lambda_variation import energy_grid, rates_grid,  lambda_distribution, energy_distribution, rates_distribution, energy_optimization
from utilities import create_directories, get_nist_energy, read_levels, compare_to_nist


#Set the iso-electronic sequence and which ions are to be calculated
# Also set which core excitation is to be calculated
seq = "be"
atoms=['o','fe']
shell = "2-2"
    
#Set the % uncertainty to be used for the NIST energy comparison by the Bayesian method.
#Set the distribution function used for the NIST comparison
#Set the prio distribution for the initial lambda distribution functions
nist_cutoff=0.03
prior_shape="uniform"
likelihood_shape="gaussian"
  
#Set this to true for the post-processor to shift series limits to NIST energies.
nist_shift=True
        
#Set the number of walker and number of steps per walker for the Markov-chain Monte-Carlo calculation
n_walkers=100
n_steps=3000

#Lambda range and grid
x_bnd = np.array([[0.4,1.6],[0.4,1.6]])
grid_resolution = 40

up_dir=res=os.getcwd()


#filename="be_like_all_pt4_1pt6_nist_shifts.pkl"
filename="output.pkl"




#%%

data={"seq":seq,"shell":shell,"nist_cutoff":nist_cutoff,"prior_shape":prior_shape,"likelihood_shape":likelihood_shape,
                 "n_walkers":n_walkers,"n_steps":n_steps}

    
start = time.time()
print('Core excitation : ' , shell)
print('nist cutoff = ',nist_cutoff)
print('prior_shape = ',prior_shape)
print('likelihood shape = ',likelihood_shape)
print('number of walkers = ', n_walkers)
print('number of steps per walker = ', n_steps)
#print('EMAX=',emax)
print('up_dir=',up_dir)
print(seq,'-like','atoms:',atoms)
print('Lambda range=',x_bnd)
print('Number of Lambda grid points',grid_resolution)
print('NIST shift selected?=',nist_shift)

print(' ')
print(' ')

sz=len(atoms)
cent_pot=1
partial_func = partial(run_case_par, seq,shell,nist_cutoff, prior_shape,likelihood_shape, cent_pot,n_walkers,n_steps,up_dir,nist_shift,grid_resolution,x_bnd)
with mp.Pool(processes=sz) as pool:
#   data_pos = pool.map_async(partial_func,iterable=args_iter).get()
   data_pos = pool.map(partial_func,iterable=atoms)
   print('Finished positive')

#cent_pot=-1
#partial_func = partial(run_case_par, seq,shell,nist_cutoff, prior_shape,likelihood_shape, cent_pot,n_walkers,n_steps,up_dir,nist_shift)
#with mp.Pool(processes=sz) as pool:
##   data_pos = pool.map_async(partial_func,iterable=args_iter).get()
#   data_neg = pool.map(partial_func,iterable=atoms)
#   print('Finished negative')



for i in range(len(atoms)):
    atom=atoms[i]
    pos_name=seq+'_like_'+atom+'_pos'
    neg_name=seq+'_like_'+atom+'_neg'
    tmp_pos={
       "T":     data_pos[i][0],
      "nist_vals":      data_pos[i][1],
      "Erg":            data_pos[i][2],
      "X_1D":           data_pos[i][3],
      "lambda_samples": data_pos[i][4],
      "E_samples":      data_pos[i][5],
      "Err":            data_pos[i][6],
      "accept":         data_pos[i][7],
      "autocorrel":     data_pos[i][8],
      "rate_samples":   data_pos[i][9],
      "rates":          data_pos[i][10]
    }
#    tmp_neg={
#       "T":     data_neg[i][0],
#      "nist_vals":      data_neg[i][1],
#      "Erg":            data_neg[i][2],
#      "X_1D":           data_neg[i][3],
#      "lambda_samples": data_neg[i][4],
#      "E_samples":      data_neg[i][5],
#      "Err":            data_neg[i][6],
#      "accept":         data_neg[i][7],
#      "autocorrel":     data_neg[i][8],
#      "rate_samples":   data_neg[i][9],
#      "rates":          data_neg[i][10]
#    }
    tmp1={pos_name:tmp_pos}
#    tmp2={neg_name:tmp_neg}
    data.update(tmp1)
#    data.update(tmp2)
 
    
   
#%%
pos_name=seq+'_like_'+atom+'_pos'
tmp1={pos_name:data_pos}
data.update(tmp1)
print('updating dictionary')
end = time.time()
print(f"Runtime: {int(end-start)}s")
print(' ')
print(' ')


#%%
rate_samples_combined=np.concatenate((data['be_like_o_pos']['rate_samples'],data['be_like_o_neg']['rate_samples']))
fig3, axs3 = plt.subplots(1,1)
title_pos='Rate coefficients'
axs3.set_yscale("log")
axs3.set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs3.set_xlabel('Electron Temperature index')
axs3.set_title(title_pos)
axs3.boxplot(rate_samples_combined,showfliers=False)

fig3, axs3 = plt.subplots(1,1)
title_pos='Rate coefficients'
axs3.set_yscale("log")
axs3.set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs3.set_xlabel('Electron Temperature index')
axs3.set_title(title_pos)
axs3.boxplot(data_pos[2][:][:],showfliers=False)

