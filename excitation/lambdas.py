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
from read_adf04 import rates_dataframe, read_adf04_np
from utilities import get_nist_energy, read_levels, compare_to_nist
from excitation import run_r_matrix, orbitals
import emcee
import corner
import matplotlib.pyplot as plt

def make_energy_grid(ion, x_ravel, x_res, nmax=3, include_einstein=False):
    # Generate error data 
    # 
    y_nist = get_nist_energy(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")
    orbs = orbitals(ion, nmax)

    if include_einstein:
        run_r_matrix(ion, lambdas=[1.0]*len(orbs), born_only=True)
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
        x=np.r_[1.0,x,[1.0]*(len(orbs)-3)]
        #potential = np.random.choice([-1, 1])
        potential = 1
        run_r_matrix(ion, lambdas=x, potential_type=potential, born_only=True)
        y_comp,ground = read_levels(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/str/LEVELS")
        if include_einstein:
            einstein = rates_dataframe("isoelectronic/he-like/o6/born/adf04ic")
            einstein = einstein[einstein["final"]==1.0]
        err[i,:] = compare_to_nist(y_comp, y_nist)
        erg[i, :] = y_comp
    
    Err = [np.reshape(err[:,j], x_res) for j in range(n)]
    Erg = [np.reshape(erg[:,j], x_res) for j in range(n)]

    return Err, Erg

def make_rates_grid(ion, x_ravel, x_res, nmax=3):
    
    orbs = orbitals(ion, nmax)
    run_r_matrix(ion, lambdas=[1.0]*len(orbs), nmax=nmax)
    df_base = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")
    df_base = df_base[df_base["final"]==1]
    
    T = df_base.columns[3:]
    n_rates = df_base.values.shape[0]
    n_temperatures = T.size
    n_points = x_ravel.shape[0]
    
    rates = np.empty((n_points, n_rates, n_temperatures))

    for i in range(n_points):
        x = x_ravel[i,:]
        x=np.r_[1.0, x, [1.0]*(len(orbs)-3)]
        run_r_matrix(ion, lambdas=x, nmax=nmax)
        df = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")
        df = df[df["final"]==1]
        rates[i, :, :] = df.values[:, 3:]
    
    
    Rates = [[np.reshape(rates[:,j, k], x_res) for k in range(n_temperatures)] for j in range(n_rates)]
    return Rates, df_base


def make_lambda_distribution(ion, x_bnd, x_res, n_lambdas=2, n_walkers=10, n_steps=1000, nist_cutoff=0.05, outfile=None):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    # NIST data
    y_nist = get_nist_energy(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")
    
    # Construct interpolators
    n_energies = y_nist.size
    
    Err, Erg = make_energy_grid(ion, x_ravel, x_res)
    
    err_interpolators = interpolators(X_1D, Err)
    
    
    y_bnd = np.zeros((n_energies, n_lambdas))
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
    
    if outfile is not None:
        np.save(outfile, arr=lambda_samples)
    
    return lambda_samples

def make_rates_distribution(ion, lambda_samples, x_bnd, x_res, save=True):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    n_samples = lambda_samples.shape[0]
    
    Rates, df_base = make_rates_grid(ion, x_ravel, x_res)
    
    T = df_base.columns[3:]
    n_rates = df_base.values.shape[0]
    n_temperatures = T.size
    
    rate_interpolators = []
    for j in range(n_rates):
        rate_interpolators += interpolators(X_1D, Rates[:,j,:])
    
    rate_samples = np.zeros((n_samples, n_rates, n_temperatures))
    
    for i in range(n_samples):
        for j in range(n_rates):
            for k in range(n_temperatures):
                rate_samples[i,j,k] = rate_interpolators[j][k](lambda_samples[i]) 
        
    if save:
        np.save(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates.npy", np.array([T, rate_samples]))
        df_base.to_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/base_adf04.csv")
        
    return T, rate_samples, df_base

"""
def rate_samples(ion, lambda_samples, n_samples=50, nmax=3, save=True):
    
    orbs = orbitals(ion, nmax)
        
    sample_size = lambda_samples.shape[0]
    for i in range(n_samples):
        lambdas = np.r_[1.0, lambda_samples[np.random.randint(0, high=sample_size), :], [1.0]*(len(orbs)-3)]
        run_r_matrix(ion, lambdas=lambdas, nmax=nmax)
        if i==0:
            df = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")
            df = df[df["final"]==1.0]
            for j in range(lambdas.size):
                df[f"lambda_{orbs[j]}"] = lambdas[j]
        else:
            data = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")
            data = data[data["final"]==1.0]
            for j in range(lambdas.size):
                data[f"lambda_{orbs[j]}"] = lambdas[j]
            df = pd.concat((df, data))
    
    df.sort_values(by=[f"lambda_{i+1}" for i in range(len(orbs))] + ["final","initial"], inplace=True)
    if save:
        df.to_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rate_coefficients.csv")
    
    return df

        

def graph_rate_coefficients(ion, datafile):
    
    df = pd.read_csv(datafile)
    by_row_ind = df.groupby(df.index)
    df_avg = by_row_ind.mean()
    
    df.sort_values(by=["initial"], inplace=True)
    levels = read_adf04_np(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")[0]
    
    
    
    
    for i in range(n_rates):
        fig, ax = plt.subplots(2,1)
        for j in range(n_samples):
            ax[0].plot(T[:-1], rates[j, i, :-1])
            ax[0].set_xlabel("Temperature (K)")
            ax[0].set_ylabel("Excitation Rate Coeff")
            ax[0].set_title(f"Excitation Rate - {levels[levels['config']]}")
"""
    

if __name__ == "__main__":
    
    atom = "o"
    seq = "he"
    shell = "1-2"
    
    ion = State(atom, seq, shell)
    
    nmax = 3
    
    orbs = orbitals(ion, nmax)
    lambdas = [1.0]*len(orbs)
    
    df = rates_dataframe("adf04")
    df = df[df["final"]==1.0]
    print(df.values[:,3:])  
    """
    df.sort_values(by=["final", "initial"], inplace=True)
    df2 = df.copy()
    
    df_comb = pd.concat((df,df2))
    df_comb = pd.concat((df_comb, df))
    
    by_row_ind = df_comb.groupby(df_comb.index)
    df_means = by_row_ind.mean()
    
    print(df_comb)
    
    # Interval endpoints for each input component
    x_bnd = np.array([[0.8,1.2],[0.8,1.2]])
    
    # Resolution in each dimension
    grid_resolution = 5
    x_res = np.array([grid_resolution, grid_resolution])
    
    lambdas = make_lambda_distribution(ion, x_bnd, x_res)
    print(lambdas.shape)
    """
    
    