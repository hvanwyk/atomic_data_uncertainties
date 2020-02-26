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
    T, rate_samples = np.load(file)
    rate_samples = np.abs(rate_samples)
    n_samples, n_rates, n_temperatures = rate_samples.shape
    finite_T = T[1:-1].values.astype(np.float)
    
    levels = read_adf04(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")[0]
    
    for j in range(n_rates):
        
        avg = np.mean(rate_samples[:,j,1:-1], axis=0)
        std = np.std(rate_samples[:,j,1:-1], axis=0)
        
        fig, ax = plt.subplots(3, 1, figsize=(10,10), dpi=300)
        ax[0].grid()
        ax[1].grid()
        ax[2].grid()
        for sample in rate_samples[::1000, j, 1:-1]:
            ax[0].semilogx(finite_T, sample)
            ax[0].set_xlabel("Temperature (K)")
            ax[0].set_ylabel("Epsilon")
            ax[0].set_title(f"Excitation Rate: {levels['config'].values[j+1]} {levels['(2S+1)L( 2J)'].values[j+1]} to Ground")
            
            ax[1].semilogx(finite_T, (std/avg)*100)
            ax[1].set_xlabel("Temperature (K)")
            ax[1].set_ylabel("% Deviation")
            ax[1].set_title("Relative Uncertainty")
            
        ax[2].hist(rate_samples[:,j,0])
        ax[2].set_xlabel("A-Value")
        ax[2].set_ylabel("Count")
        ax[2].set_title("Einstein A-Value")
            
        fig.tight_layout()
        fig.savefig(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/excitation_graph_{j+2}_to_1.eps")

if __name__ == "__main__":
    
    atom="o"
    seq="he"
    shell="1-2"
    
    ion = State(atom, seq, shell)
    graph_rates(ion, "rates.npy")