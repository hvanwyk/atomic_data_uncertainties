#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 14:36:35 2019

@author: kyle
"""

import sys
import numpy as np
import pandas as pd
if ".." not in sys.path:
    sys.path.append("..")
if "../recombination" not in sys.path:
    sys.path.append("../recombination")
from state import State
from bayesian_methods import log_posterior, interpolators, lambdas_grid
from read_adf04 import rates_dataframe, read_adf04
from utilities import get_nist_energy, read_levels, compare_to_nist
from excitation import run_r_matrix, orbitals
import emcee
import corner
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D



def make_energy_grid(ion, x_ravel, x_res, n_lambdas=2, nmax=3, potential_type=1, include_einstein=False):
    # Generate error data 
    # 
    df_nist = get_nist_energy(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")
    y_nist = df_nist["E"].values
    orbs = orbitals(ion, nmax)

    if include_einstein:
        run_r_matrix(ion, lambdas=[1.0]*len(orbs), born_only=True, potential_type=potential_type)
        einstein = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/born/adf04ic")
        einstein = einstein[einstein["final"]==1.0]
        n = y_nist.shape[0] + einstein["A"].values.shape[0]
    else:
        n = y_nist.shape[0]
        
    n_points = x_ravel.shape[0]
    err = np.empty((n_points,n))
    erg = np.empty((n_points, n))
    
    for i in range(n_points):
        x = x_ravel[i,:]
        if n_lambdas < 5:
            x=np.r_[1.0,x,[1.0]*(len(orbs)-n_lambdas-1)]
        else:
            x = np.r_[1.0, x]
        run_r_matrix(ion, lambdas=x, potential_type=potential_type, born_only=True)
        df_comp, ground = read_levels(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/str/LEVELS")
        y_comp = df_comp["E"].values
        if include_einstein:
            einstein = rates_dataframe("isoelectronic/he-like/o6/born/adf04ic")
            einstein = einstein[einstein["final"]==1.0]
        err[i,:] = compare_to_nist(df_comp, df_nist)
        erg[i, :] = y_comp
    
    Err = [np.reshape(err[:,j], x_res) for j in range(n)]
    Erg = [np.reshape(erg[:,j], x_res) for j in range(n)]

    return Err, Erg

def make_rates_grid(ion, x_ravel, x_res, n_lambdas=2, nmax=3, potential_type=1):
    
    orbs = orbitals(ion, nmax)
    run_r_matrix(ion, lambdas=[1.0]*len(orbs), nmax=nmax, potential_type=potential_type)
    df_base = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/{ion.isoelec_seq}like_ks20#{ion.species}{ion.ion_charge}.dat")
    df_base.to_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_dataframe_template_all.csv", index=None)
    df_base = df_base[df_base["final"]==1]
    df_base.to_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_dataframe_template_ground.csv", index=None)
    
    T = df_base.columns[2:]
    n_rates = df_base.values.shape[0]
    n_temperatures = T.size
    n_points = x_ravel.shape[0]
    
    rates = np.empty((n_points, n_rates, n_temperatures))

    for i in range(n_points):
        x = x_ravel[i,:]
        if n_lambdas < 5:
            x=np.r_[1.0, x, [1.0]*(len(orbs)-n_lambdas-1)]
        else:
            x = np.r_[1.0,x]
        run_r_matrix(ion, lambdas=x, nmax=nmax, potential_type=potential_type)
        df = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/{ion.isoelec_seq}like_ks20#{ion.species}{ion.ion_charge}.dat")
        df = df[df["final"]==1]
        rates[i, :, :] = df.values[:, 2:]
    
    
    Rates = [[np.reshape(rates[:,j, k], x_res) for k in range(n_temperatures)] for j in range(n_rates)]
    return np.array(Rates), df_base


def make_lambda_distribution(ion, x_bnd, x_res, n_lambdas=2, nmax=3, n_walkers=10, n_steps=10000, nist_cutoff=0.05, potential_type=1):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    # NIST data
    df_nist = get_nist_energy(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")
    y_nist = df_nist["E"].values
    
    # Construct interpolators
    n_energies = y_nist.size
    
    Err, Erg = make_energy_grid(ion, x_ravel, x_res, potential_type=potential_type)
    
    err_interpolators = interpolators(X_1D, Err)
    
    
    y_bnd = np.zeros((n_energies, 2))
    for i in range(n_energies):
        y_bnd[i, :] = -1, 1
        
    y_bnd *= nist_cutoff
    
    #
    # Specify starting points for each Markov chain (in a tight ball around optimum)
    #
    pos = [np.array([1 for i in range(n_lambdas)]) + 1e-4*np.random.randn(n_lambdas) for i in range(n_walkers)]
    
    #
    # Initialize the sampler
    #
    sampler = emcee.EnsembleSampler(n_walkers, n_lambdas, log_posterior, 
                                    args=(err_interpolators, x_bnd, y_bnd))
    #
    # Run the MCMC routine
    #
    sampler.run_mcmc(pos, n_steps);
    
    #
    # The sampler.chain has shape (n_walkers, n_steps, n_dim)
    # Reshape array of samples (but throw away initial sample)
    #
    lambda_samples = sampler.chain[:, 50:, :].reshape((-1, n_lambdas))
    
    corner.corner(lambda_samples, labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
    plt.show()
    
    
    
    return lambda_samples

def make_rates_distribution(ion, lambda_samples, x_bnd, x_res, n_lambdas=2, potential_type=1):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    n_samples = lambda_samples.shape[0]
    
    Rates, df_base = make_rates_grid(ion=ion, x_ravel=x_ravel, x_res=x_res, n_lambdas=n_lambdas, potential_type=potential_type)
    
    T = df_base.columns[2:]
    n_rates = df_base.values.shape[0]
    n_temperatures = T.size
    
    rate_interpolators = []
    for j in range(n_rates):
        rate_interpolators.append(interpolators(X_1D, Rates[j]))
    
    rate_samples = np.zeros((n_samples, n_rates, n_temperatures))
    
    for i in range(n_samples):
        for j in range(n_rates):
            for k in range(n_temperatures):
                rate_samples[i,j,k] = rate_interpolators[j][k](lambda_samples[i]) 
        

        
    return T, rate_samples, df_base

    
if __name__ == "__main__":
    
    atom = "o"
    seq = "he"
    shell = "1-2"
    
    ion = State(atom, seq, shell)
    
    nmax = 3
    
    # Interval endpoints for each input component
    n_lambdas = 2
    x_bnd = []
    for i in range(n_lambdas):
        x_bnd.append([0.8,1.2])
    x_bnd = np.array(x_bnd)
    
    
    
    
    
    # Interval endpoints for each input component
    grid_size_energies = 5
    x_res_energies = np.array([grid_size_energies]*n_lambdas)
    #lambdas = make_lambda_distribution(ion=ion, x_bnd=x_bnd, x_res=x_res_energies, n_lambdas=n_lambdas, nist_cutoff=0.05, potential_type=-1)
    
    """
    # Resolution in each dimension for rates interpolation grid
    grid_size_rates = 2
    x_res_rates = np.array([grid_size_rates]*n_lambdas)
    make_rates_distribution(ion=ion, lambda_samples=lambdas, x_bnd=x_bnd, x_res=x_res_rates, n_lambdas=n_lambdas)
    
    """
