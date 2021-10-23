#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 23 18:27:50 2021

@author: loch
"""
import numpy as np
import sys

from recombination_methods import State
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import lambdas_grid
from lambda_variation import energy_grid, rates_grid,  lambda_distribution, energy_distribution, rates_distribution, energy_optimization
from utilities import get_nist_energy

def run_case_par(seq,shell,nist_cutoff, prior_shape,likelihood_shape, cent_pot,n_walkers,n_steps,up_dir,nist_shift,grid_resolution,x_bnd,atom):
    
    ion = State(atom, seq, shell)
    nist_vals = get_nist_energy(up_dir + f"/NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist")
    emax=nist_vals[-1]*1.1  
    x_res = np.array([grid_resolution, grid_resolution])
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)

 
    print('Evaluating lambda distribution functions using MCMC for',atom)
    lambda_samples, Err, Erg, err_interpolators,y_nist,accept,autocorrel = lambda_distribution(ion,up_dir,cent_pot=cent_pot,x_bnd=x_bnd, x_res=x_res, nist_cutoff=nist_cutoff, 
                      prior_shape=prior_shape, n_walkers=n_walkers, n_steps=n_steps,
                      likelihood_shape=likelihood_shape,emax=emax)

    print('Calculating Energy distribution functions for',atom)
    E_samples, Err, Erg=energy_distribution(ion, up_dir,emax, lambda_samples, x_bnd, x_res, cent_pot)

    print('Calculating rate distributions for',atom,' using cent_pot=',cent_pot,'emax=',emax,'nist_shift=',nist_shift)
    T, rate_samples,rates = rates_distribution(ion, up_dir,emax, lambda_samples, x_bnd, x_res, cent_pot,nist_shift=nist_shift)
    
    print('Finished',seq,'-like',atom)
    result=[T,nist_vals,Erg,X_1D,lambda_samples,E_samples,Err,accept,autocorrel,rate_samples,rates]
    return result