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



def make_energy_grid(ion, x_ravel, x_res, n_lambdas=2, nmax=3, include_einstein=False):
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
        if n_lambdas < 5:
            x=np.r_[1.0,x,[1.0]*(len(orbs)-n_lambdas-1)]
        else:
            x = np.r_[1.0, x]
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

def make_rates_grid(ion, x_ravel, x_res, n_lambdas=2, nmax=3):
    
    orbs = orbitals(ion, nmax)
    run_r_matrix(ion, lambdas=[1.0]*len(orbs), nmax=nmax)
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
        run_r_matrix(ion, lambdas=x, nmax=nmax)
        df = rates_dataframe(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/{ion.isoelec_seq}like_ks20#{ion.species}{ion.ion_charge}.dat")
        df = df[df["final"]==1]
        rates[i, :, :] = df.values[:, 2:]
    
    
    Rates = [[np.reshape(rates[:,j, k], x_res) for k in range(n_temperatures)] for j in range(n_rates)]
    return np.array(Rates), df_base


def make_lambda_distribution(ion, x_bnd, x_res, n_lambdas=2, n_walkers=10, n_steps=10000, nist_cutoff=0.05, outfile=None):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
    # NIST data
    y_nist = get_nist_energy(f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}_n={nmax}.nist")
    
    # Construct interpolators
    n_energies = y_nist.size
    
    Err, Erg = make_energy_grid(ion, x_ravel, x_res)
    
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
    
    if outfile is not None:
        np.save(outfile, arr=lambda_samples)
    
    return lambda_samples

def make_rates_distribution(ion, lambda_samples, x_bnd, x_res, n_lambdas=2, save=True):
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    n_samples = lambda_samples.shape[0]
    
    Rates, df_base = make_rates_grid(ion=ion, x_ravel=x_ravel, x_res=x_res, n_lambdas=n_lambdas)
    
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
        
    if save:
        np.save(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates.npy", np.array([T, rate_samples]))
        
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
    
    levels, T, _1, _2 = read_adf04("adf04")
    
    w_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(1)1( 1.0)'), "#"].values[0])-1
    x_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 2.0)'), "#"].values[0])-1
    y_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 1.0)'), "#"].values[0])-1
    z_num = int(levels.loc[(levels['config']=='1S1 2S1') & (levels['(2S+1)L( 2J)'] == '(3)0( 1.0)'), "#"].values[0])-1
    
    # Resolution in each dimension for lambdas interpolation grid
    grid_size_energies = 5
    x_res_energies = np.array([grid_size_energies]*n_lambdas)
    
    X_1D, x_ravel = lambdas_grid(x_bnd, x_res_energies)
    
    lam_grid = np.array(x_ravel)

    """
    Err, Erg = make_energy_grid(ion, x_ravel, x_res_energies)
    
    err_grid = np.array(Err)
    
    for i in range(1, len(Err)):
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter(x_ravel[:, 0], x_ravel[:, 1], err_grid[i,:,:])
        ax.set_xlabel("2s Lambda")
        ax.set_ylabel("2p Lambda")
        ax.set_zlabel(f"E$_{i}$ Error")
        fig.tight_layout()
        fig.savefig(f"Err{i}_grid.eps")
    """
    
    # Resolution in each dimension for rates interpolation grid
    grid_size_rates = 2
    x_res_rates = np.array([grid_size_rates]*n_lambdas)
    
<<<<<<< HEAD
    rates_grid, df_base = make_rates_grid(ion, x_ravel, x_res_rates)
        
    fig_X = plt.figure()
    
    ax = fig_X.add_subplot(211, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[x_num,0, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel("X A-Value")
    
    ax = fig_X.add_subplot(212, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[x_num,10, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel(f"X Epsilon at T={T[9]}")
    
    fig_X.tight_layout()
    fig_X.savefig(f"X_grid.eps")
    
    
    
    fig_Z = plt.figure()
    
    ax = fig_Z.add_subplot(211, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[z_num,0, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel("Z A-Value")
    
    ax = fig_Z.add_subplot(212, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[z_num,10, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel(f"Z Epsilon at T={T[9]}")
    
    fig_Z.tight_layout()
    fig_Z.savefig(f"Z_grid.eps")
    
    
    
    fig_W = plt.figure()
    
    ax = fig_W.add_subplot(211, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[w_num,0, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel("W A-Value")
    
    ax = fig_W.add_subplot(212, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[w_num,10, :,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel(f"W Epsilon at T={T[9]}")
    
    fig_W.tight_layout()
    fig_W.savefig(f"W_grid.eps")
    
    """
    lambdas = make_lambda_distribution(ion=ion, x_bnd=x_bnd, x_res=x_res_energies, n_lambdas=n_lambdas)
    
    
    make_rates_distribution(ion=ion, lambda_samples=lambdas, x_bnd=x_bnd, x_res=x_res_rates, n_lambdas=n_lambdas)
    """    
=======
    
>>>>>>> c5f04a32d40bbe0cf6f76c520a256c36d1059214
