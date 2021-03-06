#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 18 12:03:40 2018

@author: kyle
"""
import numpy as np
import os
import matplotlib.pyplot as plt
import corner
import emcee
 
from recombination_methods import State, structure, structure_dr, postprocessing_rates
from bayesian_methods import log_posterior, interpolators
from lambda_variation import lambda_distribution, energy_optimization
from graphing import graph_experimental, graph_rates_from_file
import time

"""
Set of codes for implementing the energy shift method for computing uncertainties
on the Dielectronic recombination rates at low temperatures. 
"""


# =============================================================================
# 1. Construct Interpolators
# =============================================================================

def shift_grid(x_bnd, x_res):
    d = x_bnd.shape[0]
    X_1D = []
    for i in range(d):
        X_1D.append(np.linspace(x_bnd[i,0],x_bnd[i,1],x_res[i]))
    
    #
    # Multi-dimensional Grid, using meshgrid (ij indexing)
    #
    X = np.meshgrid(*X_1D, indexing='ij')
        
    #
    # Unravel grid into (n_points, d) array
    # 
    x_ravel = np.array([x.ravel() for x in X]).T
    
    return X_1D, x_ravel

def rates_grid(ion, x_ravel, x_res):
    
    E, E_nist, E_shift = structure(ion, method="shift")
    structure_dr(ion, method="shift")
    T = postprocessing_rates(ion, E, E_nist, method="shift")[0]
    n_rates = T.size
    n_points = x_ravel.shape[0]
    
    rates = np.zeros((n_points, n_rates))
    
    for i in range(n_points):
        x = x_ravel[i,:]
        rates[i, :] = postprocessing_rates(ion, E, E_nist, method="shift", shift=x)[1]
    
    Rates = [np.reshape(rates[:,j], x_res) for j in range(n_rates)]
    return Rates

    
def variable_shift_method(ion, n_points, max_shift, res):
    
    E, E_nist, E_shift = structure(ion, method="shift")
    structure_dr(ion, method="shift")
    x_res = np.full(E.shape, res)
    T = postprocessing_rates(ion, E, E_nist, method="shift", shift=np.zeros(E.shape))[0]
    n_rates = T.size
    rates = np.zeros((n_points, n_rates))
    
    x_bnd = np.zeros((E.size, 2))
    for i in range(x_bnd.shape[0]):
        x_bnd[i,:] = -max_shift, max_shift
    
    X_1D, x_ravel = shift_grid(x_bnd, x_res)
    
    Rates = rates_grid(ion, x_ravel, x_res)
    
    rate_interpolators = interpolators(X_1D, Rates)
    
    plt.figure()
    for i in range(n_points):
        shifts = (np.random.rand(*E.shape) * 2 - 1)*max_shift / 13.6
        for j in range(n_rates):
            rates[i, j] = rate_interpolators[j](shifts)
        #rates[i,:] = postprocessing_rates(ion, E, E_nist, method="shift", shift=shifts)[1]
        plt.semilogx(T, rates[i, :])
    plt.show()
    
    avg = np.mean(rates, axis=0)
    
    err = np.std(rates, axis=0)
    #err = (np.max(rates, axis=0)-np.min(rates, axis=0))/2
    plt.figure()
    plt.errorbar(T, avg, yerr=err)
    
def simple_shift_method(ion, n_samples, max_shift, rates_file, xsec_file="", ECORIC=0, NMIN=3, NMAX=15, LMIN=0, LMAX=7,
                        compute_xsec=False, EWIDTH=0.0001, NBIN=1000, EMIN=0.0, EMAX=2.0):
    """
    Calculate the Dielectronic Radiation cross-sections and rates by means of the simple shift method. 
    
    
    Inputs:
    
        ion: State, specified by atom, isoelectronic sequence, and shell
    
        n_samples: int, sample size (number of lambdas)
    
        max_shift: double, maximum amount to shift by
    
        rates_file: str, file name for storing DR rates
    
        xsec_file: str, file name for storing cross-sections
    
        ECORIC: str, AS namelist input, correction energy added to continuum energy
    
        NMIN, NMAX: int, DRR NAMELIST input, min/max n-value of Rydberg electron
    
        LMIN, LMAX: int, DRR NAMELIST input, min/max l-value of Rydberg electron
        
        compute_xsec: bool, whether cross-sections should be computed 
    
        EWIDTH: double, AS namelist input,
            >0 Convolute the (bin) energy-averaged cross sections with a Gaussian
                of FWHM EWIDTH UNITS.
            =-0.5 (actually any non-integer negative number) convolute with a
                Maxwellian distribution.
            =0.0 Convolute with a cooler distribution,
    
        NBIN: int, AS namelist input, number of bin energies (number of bins=NBIN-1). 
            The bin width, (EMAX-EMIN)/(NBIN-1) should be smaller than the convolution 
            width and larger than the natural width of the resonances.

        EMIN, EMAX: double, AS NAMELIST, specify the energy range (in Rydbergs)
    
    """
    #
    # Run Structure Calculation to get energies and LEVELS file
    # 
    E = structure(ion, method="shift")
    
    #
    # 
    # 
    structure_dr(ion, method="shift", ECORIC=ECORIC, NMIN=NMIN, NMAX=NMAX, LMIN=LMIN, LMAX=LMAX)
    
    #
    # Run ADASDR postprocessor to produce temperatures at which 
    # 
    T = postprocessing_rates(ion, E, method="shift", shift=np.zeros(E.shape))[0]
    n_rates = T.size
    rates = np.zeros((n_samples, n_rates))
    
    
    energy = postprocessing_rates(ion, E, method="shift", shift=np.zeros(E.shape), compute_xsec=True,
                EWIDTH=EWIDTH, NBIN=NBIN, EMIN=EMIN, EMAX=EMAX)[0]
    n_points = energy.size
    xsec = np.zeros((n_samples, n_points))
    
    for i in range(n_samples):
        #
        # Vary shifts randomly around reference energies
        # 
        shifts = (np.random.rand(*E.shape) * 2 - 1)*max_shift / 13.6
        rates[i,:] = postprocessing_rates(ion, E, method="shift", shift=shifts)[1]
        if compute_xsec:
            xsec[i,:] = postprocessing_rates(ion, E, method="shift", shift=shifts, compute_xsec=True, 
                EWIDTH=EWIDTH, NBIN=NBIN, EMIN=EMIN, EMAX=EMAX)[1]
    np.save(rates_file, np.array([T,rates]))
    
    if compute_xsec:
        np.save(xsec_file, np.array([energy, xsec]))

    
def graph_rates_shift(rates_file, xsec_file="", ECORIC=0):
    """
    Plot temperature vs DR rate and energy vs cross section graphs for given
    
    Inputs:
    
        rates_file: str, location of file storing rates
        
        xsec_file: str, location of cross-sections file.
        
    
    """
    T, rates = np.load(rates_file)    
    avg = np.mean(rates, axis=0)
    err = np.std(rates, axis=0)
    fig, ax = plt.subplots(2, 1, figsize=(6,8))
    for rate in rates:
        ax[0].semilogx(T, rate)
    ax[0].set_xlabel("Temperature (K)")
    ax[0].set_ylabel("DR Rate (cm^3 s^-1)")
    ax[0].set_yscale("log")
    ax[0].set_title(f"Dielectronic Recombination of {atom.capitalize()}{ion.ion_charge}+" + 
      ("" if ECORIC==0 else f" - {ECORIC} Ryd Below Threshold"))
    
    ax[1].plot(T, err)
    ax[1].set_xlabel("Temperature (K)")
    ax[1].set_ylabel("Std. Dev. (cm^3 s^-1)")    
    ax[1].set_xscale("log")
    #ax[1].set_yscale("log")
    ax[1].set_title("Uncertainty")
    
    fig.tight_layout()
    fig.savefig(".".join(rates_file.split(".")[:-1]) + ".png")
    plt.show()
    
    if xsec_file != "":
        plt.figure()
        plt.title(f"DR Cross Section of {ion.species}{ion.ion_charge}")
        plt.xlabel("Energy (Ryd)")
        plt.ylabel("Cross Section")
        plt.xscale("log")
        
        energy, xsec = np.load(xsec_file)
        
        for x in xsec:
            plt.plot(energy, x)
        
        plt.savefig(".".join(xsec_file.split(".")[:-1]) + ".png")
    
    
    dir_contents = os.listdir(direc)
    
    if "experimental_coefficients.dat" in dir_contents:
        plt.figure()
        
        plt.errorbar(T, avg, yerr=err, label="Theory")
        #plt.yscale("log")
        graph_experimental(ion, direc)
        plt.title(f"DR of {atom.capitalize()}{ion.ion_charge}+ --- Max Shift={max_shift} eV" + 
                           ("" if ECORIC==0 else f" - {ECORIC*13.6} eV Below Threshold"))
        plt.ylim((0, None))
        plt.xlim(1e2, 1e7)
        plt.savefig(direc + f"experiment_vs_theory_shift_{max_shift}" +  ("" if ECORIC==0 else f"_ECORIC{ECORIC}") + ".png")
        plt.show()


atom = "o"
seq = "b"
shell = "2-2"


ion = State(atom, seq, shell)
print(ion)

max_shift = 1.5
n_points = 50

ECORIC = 0
NMIN = 3
NMAX = 15
LMIN = 0
LMAX = 7
direc = f"results/isoelectronic/{seq}/{atom}{ion.ion_charge}/"
rates_file = direc + f"rates_shift_{max_shift}" + ("" if ECORIC==0 else f"_ECORIC{ECORIC}") + (
        "" if (NMIN==3 and NMAX==15) else f"_N_{NMIN}_{NMAX}") + (
        "" if (LMIN==0 and LMAX==7) else f"_L_{LMIN}_{LMAX}") +".npy"

xsec_file = direc + f"xsec_shift_{max_shift}" + ("" if ECORIC==0 else f"_ECORIC{ECORIC}") + (
        "" if (NMIN==3 and NMAX==15) else f"_N_{NMIN}_{NMAX}") + (
        "" if (LMIN==0 and LMAX==7) else f"_L_{LMIN}_{LMAX}") +".npy"
        
simple_shift_method(ion, n_points, max_shift, rates_file, xsec_file, ECORIC=ECORIC, compute_xsec=True, NMIN=NMIN, NMAX=NMAX, LMIN=LMIN, LMAX=LMAX)
graph_rates_shift(rates_file, xsec_file, ECORIC=ECORIC)

#shift_method_xsec(ion, 100, max_shift, xsec_file, ECORIC=ECORIC, NMIN=NMIN, NMAX=NMAX, LMIN=LMIN, LMAX=LMAX)
#graph_xsec(xsec_file, )


        
