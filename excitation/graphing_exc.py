#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:42:20 2020

@author: kyle
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_adf04 import read_adf04
from state import State

def graph_rates(ion, file):
    T, rate_samples = np.load(file)
    rate_samples = np.abs(rate_samples)
    n_samples, n_rates, n_temperatures = rate_samples.shape
    finite_T = T[:-1].values.astype(np.float)
    
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    for j in range(n_rates):
        
        avg = np.mean(rate_samples[:,j,:-1], axis=0)
        std = np.std(rate_samples[:,j,:-1], axis=0)
        
        fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        for sample in rate_samples[::100, j, :-1]:
            ax[0].semilogx(finite_T, sample)
            ax[0].set_xlabel("Temperature (K)")
            ax[0].set_ylabel("Rate Coeff.")
            ax[0].set_title(f"Excitation Rate: {levels['config'].values[j+1]} {levels['(2S+1)L( 2J)'].values[j+1]} to Ground")
            
            ax[1].semilogx(finite_T, (std/avg)*100)
            ax[1].set_xlabel("Temperature (K)")
            ax[1].set_ylabel("% Deviation")
            ax[1].set_title("Relative Uncertainty")
            
            ax[2].hist(rate_samples[:,j,-1])
            ax[2].set_xlabel("Rate Coefficient")
            ax[2].set_ylabel("Count")
            ax[2].set_title("Infinite Energy Rate")
            
        fig.tight_layout()
        fig.savefig(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/excitation_graph_{j+2}_to_1.eps")


def plot_G_ratio(ion, rates_file):
    T, rate_samples = np.load(rates_file)
    finite_T = T[:-1].values.astype(np.float)
    rate_samples = np.abs(rate_samples)
    n_samples, n_rates, n_temperatures = rate_samples.shape

    n_temperatures -= 1
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    w_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(1)1( 1.0)'), "#"].values[0])
    x_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 2.0)'), "#"].values[0])
    y_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 1.0)'), "#"].values[0])
    z_num = int(levels.loc[(levels['config']=='1S1 2S1') & (levels['(2S+1)L( 2J)'] == '(3)0( 1.0)'), "#"].values[0])
    
    w = rate_samples[:, w_num-2, :-1]
    x = rate_samples[:, x_num-2, :-1]
    y = rate_samples[:, y_num-2, :-1]
    z = rate_samples[:, z_num-2, :-1]
    
    G = (x+y+z)/w
    
    G_avg = np.mean(G, axis=0)
    G_err = np.zeros(n_temperatures)
    
    for i in range(n_temperatures):
    
        df = pd.DataFrame(np.c_[x[:,i], y[:,i], z[:,i], x[:,i]+y[:,i]+z[:,i], w[:,i]], columns=["X","Y","Z","X+Y+Z","W"])
        
        cov = df.cov()
        
        G_err[i] = np.mean(G[:,i])*np.sqrt((1/(np.mean(x[:,i]+y[:,i]+z[:,i]))**2)*(cov["X"]["X"] + cov["Y"]["Y"] + cov["Z"]["Z"]) + 
                        cov["W"]["W"]/(np.mean(w[:,i])**2) - 2*cov["W"]["X+Y+Z"]/(np.mean(x[:,i]+y[:,i]+z[:,i])*np.mean(w[:,i])))
    
    plt.figure()
    plt.errorbar(finite_T, G_avg, yerr=G_err)
    plt.xscale("log")
    plt.title(f"G Ratio vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+")
    plt.xlabel("Temperature (K)")
    plt.ylabel("G Ratio")


def plot_R_ratio(ion, rates_file):
    T, rate_samples = np.load(rates_file)
    finite_T = T[:-1].values.astype(np.float)
    rate_samples = np.abs(rate_samples)
    n_samples, n_rates, n_temperatures = rate_samples.shape
    
    n_temperatures-=1
    
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    x_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 2.0)'), "#"].values[0])
    y_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 1.0)'), "#"].values[0])
    z_num = int(levels.loc[(levels['config']=='1S1 2S1') & (levels['(2S+1)L( 2J)'] == '(3)0( 1.0)'), "#"].values[0])
    
    x = rate_samples[:, x_num-2, :-1]
    y = rate_samples[:, y_num-2, :-1]
    z = rate_samples[:, z_num-2, :-1]
    R = (x+y)/z
    
    R_avg = np.mean(R, axis=0)
    R_err = np.zeros(n_temperatures)
    
    for i in range(n_temperatures):
    
        df = pd.DataFrame(np.c_[x[:,i], y[:,i], x[:,i]+y[:,i], z[:,i]], columns=["X","Y","X+Y","Z"])
        
        cov = df.cov()
        
        R_err[i] = np.mean(R[:,i])*np.sqrt((1/(np.mean(x[:,i]+y[:,i]))**2)*(cov["X"]["X"] + cov["Y"]["Y"]) + 
                        cov["Z"]["Z"]/(np.mean(z[:,i])**2) - 2*cov["Z"]["X+Y"]/(np.mean(x[:,i]+y[:,i])*np.mean(z[:,i])))
    
    plt.figure()
    plt.errorbar(finite_T, R_avg, yerr=R_err)
    plt.xscale("log")
    plt.title(f"R Ratio vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+")
    plt.xlabel("Temperature (K)")
    plt.ylabel("R Ratio")
    print(df.corr())