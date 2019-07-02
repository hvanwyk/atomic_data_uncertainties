#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:55:15 2018

@author: kyle
"""

import numpy as np
import pandas as pd
from recombination_methods import State
import matplotlib.pyplot as plt
import os


def experimental_fit(c, E, T):
    k_b = 8.617e-5
    fit = np.zeros(T.shape)
    for i in range(len(c)):
        fit += c[i]*np.exp(-E[i]/(k_b*T))
    fit /= (T**1.5)
    return fit

def graph_rates_from_file(ion, infile, outfile, graph_every=100, x_range=None, y_range=None, err_type="std", experimental=False):

    data = np.load(infile)
    direc = "/".join(outfile.split("/")[:-1])+"/"
    T = data[0]
    rates = data[1]
    
    avg = np.mean(rates, axis=0)
    if err_type=="std":
        err = np.std(rates, axis=0)
    elif err_type=="max":
        err = (np.max(rates, axis=0) - np.min(rates, axis=0))/2
        
    graphed = rates[::graph_every]
    
        
  
    fig, ax = plt.subplots(2, 1, figsize=(6,8))
    for rate in graphed:
        ax[0].plot(T, rate)
    ax[0].set_title(f"Dielectronic Recombination of {ion.species.capitalize()}{ion.ion_charge}")
    ax[0].set_xlabel("Temperature (K)")
    ax[0].set_ylabel("DR Rate (cm^3 s^-1)")
    ax[0].set_xscale("log")
    #ax[0].set_yscale("log")
    
    ax[1].plot(T, 100*err/avg)
    ax[1].set_title("Percent Uncertainty")
    ax[1].set_xlabel("Temperature (K)")
    ax[1].set_ylabel("% Uncertainty Relative to Mean")    
    ax[1].set_xscale("log")
    #ax[1].set_yscale("log")
    
    for a in ax:
        if (x_range != None):
            a.set_xlim(x_range[0], x_range[1])
        if (y_range != None):
            a.set_ylim(y_range[0], y_range[1])
    fig.tight_layout()
    fig.savefig(outfile)
    
    
    #Used when plotting experimental data
    dir_contents = os.listdir(direc)
    if "experimental_coefficients.dat" in dir_contents:
        plt.figure()
        temps = np.linspace(100, 1000000, 100000)
        coeff = pd.read_csv(direc+"experimental_coefficients.dat", delimiter=", ")
        c = coeff["c"].values
        E = coeff["E"].values
        plt.plot(temps, experimental_fit(c, E, temps), label="Experiment", color="orange")
        plt.errorbar(x=T, y=avg, yerr=err, label="Theory")
        #plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Temperature (K)")
        plt.ylabel("DR Rate (cm^3 s^-1)")
        plt.title(f"Dielectronic Recombination of {ion.species.capitalize()}{ion.ion_charge}+")
        plt.legend()
        plt.savefig(direc + "/experiment_vs_theory.png")
    
    plt.figure()
    T_hist_choice = 1
    plt.hist(rates[:, T_hist_choice], 100)
    plt.title(f"T={T[T_hist_choice]}K Rates Histogram")
    plt.xlabel("DR Rate (cm^3 s^-1)")
    plt.ylabel("Count")
    plt.savefig(direc + f"histogram_{T_hist_choice}")

def graph_xsec(ion, energy, xsec, outfile):
    plt.figure()
    plt.title(f"DR Cross Section of {ion.species}{ion.ion_charge}")
    plt.xlabel("Energy (Ryd)")
    plt.ylabel("Cross Section")
    plt.xscale("log")
    #plt.xlim(5e3, 1e5)
    
    plt.plot(energy, xsec)
    
    plt.savefig(outfile)

def graph_experimental(ion, direc):
    temps = np.linspace(100, 1e8, 100000)
    coeff = pd.read_csv(direc+"experimental_coefficients.dat", delimiter=", ")
    c = coeff["c"].values
    E = coeff["E"].values
    plt.plot(temps, experimental_fit(c, E, temps), label="Experiment")
    #plt.yscale("log")
    plt.xscale("log")
    plt.xlabel("Temperature (K)")
    plt.ylabel("DR Rate (cm^3 s^-1)")
    plt.title(f"Dielectronic Recombination of {ion.species.capitalize()}{ion.ion_charge}+")
    plt.legend()
    #plt.savefig(direc + "/experiment_vs_theory.png")

if __name__ == "__main__":
    
    atom = "fe"
    seq = "be"
    shell = "2-2"
    ion = State(atom, seq, shell)
    nist_cutoff=0.05
    prior_shape="gaussian"
    likelihood_shape="gaussian"
    ion = State(atom, seq, shell)
    
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"    
    file_name_common = f"_{int(nist_cutoff*100)}_{prior_shape}P_{likelihood_shape}L"
    
    # Be-like oxygen
    #c = np.array([2.02e-7, 1.51e-5, 5.2e-5, 3.46e-4, 6.77e-4, 4.09e-3, 6.76e-3])
    #E = np.array([0.01, 0.06, 0.54, 2.24, 6.93, 15.67, 21.28])
    
    #Be-like carbon
    #c = np.array([1.12e-5, 2.14e-5, 3.77e-5, 4.4e-3])
    #E = np.array([0.39, 1.33, 3.25, 12.05])
    
    # Be-like nitrogen
    #c = np.array([4.22e-6, 4.8e-5, 1.93e-4, 7.14e-3])
    #E = np.array([0.23, 0.69, 2.92, 14.83])
    
    rates_file = direc+"rates"+file_name_common+".npy"
    graph_file=direc + "rates"+file_name_common + ".png"
    graph_rates_from_file(ion, infile = rates_file, outfile=graph_file)