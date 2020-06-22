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
import sys
if ".." not in sys.path:
    sys.path.append("..")
from state import State

def graph_rates(ion, file):
    T, rate_samples_tf, rate_samples_hartree = np.load(file)
    rate_samples_tf = np.abs(rate_samples_tf)
    rate_samples_hartree = np.abs(rate_samples_hartree)
    
    n_samples, n_rates, n_temperatures = rate_samples_tf.shape
    finite_T = T[1:-1].values.astype(np.float)
    
    print(rate_samples_hartree.shape)
    print(rate_samples_tf.shape)
    """
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    for j in range(n_rates):
        
        avg = np.mean(rate_samples_tf[:,j,1:-1], axis=0)
        std = np.std(rate_samples_tf[:,j,1:-1], axis=0)
        
        fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        for sample in rate_samples_tf[::1000, j, 1:-1]:
            ax[0].semilogx(finite_T, sample)
            ax[0].set_xlabel("Temperature (K)")
            ax[0].set_ylabel("Epsilon")
            ax[0].set_title(f"Excitation Rate: {levels['config'].values[j+1]} {levels['(2S+1)L( 2J)'].values[j+1]} to Ground")
            
            ax[1].semilogx(finite_T, (std/avg)*100)
            ax[1].set_xlabel("Temperature (K)")
            ax[1].set_ylabel("% Deviation")
            ax[1].set_title("Relative Uncertainty")
            
        ax[2].hist(rate_samples_tf[:,j,0])
        ax[2].set_xlabel("A-Value")
        ax[2].set_ylabel("Count")
        ax[2].set_title("Einstein A-Value")
            
        fig.tight_layout()
        fig.savefig(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/excitation_graph_tf_{j+2}_to_1.eps")
        plt.close()
    
    n_samples, n_rates, n_temperatures = rate_samples_hartree.shape
    finite_T = T[1:-1].values.astype(np.float)
    
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    for j in range(n_rates):
        
        avg = np.mean(rate_samples_hartree[:,j,1:-1], axis=0)
        std = np.std(rate_samples_hartree[:,j,1:-1], axis=0)
        
        fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        for sample in rate_samples_hartree[::1000, j, 1:-1]:
            ax[0].semilogx(finite_T, sample)
            ax[0].set_xlabel("Temperature (K)")
            ax[0].set_ylabel("Epsilon")
            ax[0].set_title(f"Excitation Rate: {levels['config'].values[j+1]} {levels['(2S+1)L( 2J)'].values[j+1]} to Ground")
            
            ax[1].semilogx(finite_T, (std/avg)*100)
            ax[1].set_xlabel("Temperature (K)")
            ax[1].set_ylabel("% Deviation")
            ax[1].set_title("Relative Uncertainty")
            
        ax[2].hist(rate_samples_hartree[:,j,0])
        ax[2].set_xlabel("A-Value")
        ax[2].set_ylabel("Count")
        ax[2].set_title("Einstein A-Value")
            
        fig.tight_layout()
        fig.savefig(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/excitation_graph_hartree_{j+2}_to_1.eps")
        plt.close()

    """
if __name__ == "__main__":
    
    atom="o"
    seq="be"
    shell="2-2"
    
    ion = State(atom, seq, shell)
    graph_rates(ion, f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rate_samples.npy")
    
    
    
    
    """
    T, rates_tf = np.load(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_thomasfermi.npy")
    T, rates_h = np.load(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_hartree.npy")
    
    T = T[1:-1].astype(np.float)
    
    avg_tf = np.mean(rates_tf, axis=0)[5, 1:-1]
    err_tf = np.std(rates_tf, axis=0)[5, 1:-1]
    
    avg_h = np.mean(rates_h, axis=0)[5, 1:-1]
    err_h = np.std(rates_h, axis=0)[5, 1:-1]
    
    plt.figure()
    
    plt.plot(T, avg_h, "k")
    plt.fill_between(T, avg_h-err_h, avg_h+err_h, color="blue", label="1 sigma")
    plt.fill_between(T, avg_h-2*err_h, avg_h-err_h, color="deepskyblue", label="2 sigma")
    plt.fill_between(T, avg_h+err_h, avg_h+2*err_h, color="deepskyblue")
    plt.fill_between(T, avg_h-3*err_h, avg_h-2*err_h, color="lightskyblue", label="3 sigma")
    plt.fill_between(T, avg_h+2*err_h, avg_h+3*err_h, color="lightskyblue")
    plt.title("7-1 Transition Epsilon")
    
    
    plt.errorbar(T, avg_tf, yerr=err_tf, label="Thomas-Fermi")
    plt.errorbar(T, avg_h, yerr=err_h, label="Hartree")
    plt.title("Effect of Central Field Potential on 2-1 Epsilon")
    
    
    plt.xlabel("Temperature (K)")
    plt.ylabel("Epsilon")
    plt.legend()
    plt.xscale('log')
    
    plt.savefig("7-1 Error Plot - Hartree.eps")
    #plt.savefig("2-1 Thomas-Fermi vs. Hartree.eps")
    """
    