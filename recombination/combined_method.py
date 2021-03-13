#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:16:24 2019

@author: kyle
"""

import numpy as np
import os
import time
import sys
if ".." not in sys.path:
    sys.path.append("..")
from utilities import create_directories
from state import State
from recombination_methods import structure, structure_dr, postprocessing_rates
from lambda_variation import lambda_distribution
from graphing import graph_rates_from_file

def get_them_rates_boi(atom):
    ion = State(atom, seq, shell)
        
    x_bnd = np.array([[0.8, 1.2], [0.8,1.2]])
    x_res = np.array([5,5])
    prior_shape = "gaussian"
    likelihood_shape = "gaussian"
    
    create_directories(ion, method="combined")
    
    direc = f"results/isoelectronic/{seq}/{atom}{ion.ion_charge}/combined_method/"
    
    os.system("cp " + f"results/isoelectronic/{seq}/{atom}{ion.ion_charge}/experimental_coefficients.dat " +
              direc+"experimental_coefficients.dat")
    lambdas_file = direc + "lambdas.npy"
    
    n_samples = 50
    
    if False: #"lambdas.npy" in os.listdir(direc):
        lambda_samples = np.load(lambdas_file)
    else:
        lambda_samples = lambda_distribution(ion, x_bnd=x_bnd, x_res=x_res, nist_cutoff=nist_cutoff, prior_shape=prior_shape,
                                         likelihood_shape=likelihood_shape, outfile=lambdas_file)
        
    lambda_samples = lambda_samples[np.random.randint(0, lambda_samples.shape[0], size=n_samples), :]

    max_shift = 0.2
    
    E, E_nist, delta_E = structure(ion, method="combined")
    structure_dr(ion, method="combined")
    T = postprocessing_rates(ion, E, E_nist, method="combined")[0]
    n_rates = T.size
    rates = np.zeros((n_samples, n_rates))
    
    for i in range(n_samples):
        potential = np.random.choice([-1,1])
        x = lambda_samples[i,:]
        if seq == "he":
            x = np.r_[1.0,x]
        E, E_nist, delta_E = structure(ion, method="combined", lambdas=x, potential=potential)
        structure_dr(ion, method="combined", lambdas=x, potential=potential)
    
        shifts = (np.random.rand(*E.shape) * 2 - 1)*max_shift / 13.6
        """
        shifts = np.zeros(*E.shape)
        for i in range(len(shifts)-5):
            shifts[i] = (np.random.random()*2 - 1) * 0.2 / 13.6
        for i in range(len(shifts)-5, len(shifts)):
            shifts[i] = (np.random.random()*2 - 1) * 1.5 / 13.6
        """
        shifts[0] = 0.0
        
        rates[i,:] = postprocessing_rates(ion, E, E_nist, method="combined", shift=shifts)[1]
    
    
    rates_file = direc + "rates.npy"
    np.save(rates_file, np.array([T,rates]))

    graphs = direc + f"graphs_{n_samples}_samples.png"
    
    graph_rates_from_file(ion, rates_file, graphs, graph_every=1)
    
    
if __name__ == "__main__":
    
    start = time.time()
    
    shell = "1-2" #core excitation shells
    atom = "o"
    seq = "he" #isoelectronic sequence
    nist_cutoff = 0.01

    
    get_them_rates_boi(atom)
    
        
                
        
        
        
