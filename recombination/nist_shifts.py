#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 31 16:33:41 2018

@author: kyle
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from recombination import State, structure

def nist_shift(seq, atom, lambdas):
    
    structure(seq, atom, lambdas)
    
    initial = State(atom, seq)
    direc = "results/isoelectronic/{}/{}{}".format(seq, atom, initial.ion_charge)
    os.chdir(direc)

    with open("LEVELS", "r") as levels:
        lines = levels.read().splitlines()
        shift = np.zeros(len(lines)-2)
        for i in range(1, len(lines)-1):
            line = lines[i].split()
            shift[i-1] = float(line[-1])
            
    with open("../../../../NIST/isoelectronic/{}/{}{}.nist".format(seq, atom, initial.ion_charge)) as nist:
            lines = nist.read().splitlines()
            for i in range(1, len(lines)):
                line = lines[i].split()
                temp = shift[i-1]
                shift[i-1] -= float(line[-1])
                shift[i-1] /= temp
                
    os.chdir("../../../..")
    return (shift*100)
    
    
seq = "li"
atom = "c"

grid_size = 100
cutoff = 5

pts = np.linspace(0.8, 1.2, grid_size)
num_lams = len(nist_shift(seq, atom, [])) - 1

species = State(atom, seq)
ion = "{}{}".format(atom, species.ion_charge)

with open("results/isoelectronic/{}/{}/{}_lambdas".format(seq, ion, ion), "w") as output:
    for i in range(num_lams):
        output.write(" Lambda {} | %Shift {} |".format(i+1, i+1))
    output.write("\n")
    output.write("-"*45 + "\n")
    for i in range(grid_size):
        for j in range(grid_size):
            lam_f = [pts[i], pts[j]]
            lam_str = [str(x) for x in lam_f]
            temp = nist_shift(seq, atom, lam_str)
            if np.ndarray.max(temp[1:]) <= cutoff:
                for i in range(1, num_lams+1):
                    output.write("{:>9.4f}  {:>9.4f}  ".format(lam_f[i-1], temp[i]))
                output.write("\n")
              
                
                    
"""
fig, ax = plt.subplots(2, 1)
ax[0].hist(shifts[:, 1], bins=50)
ax[0].set_title("% Deviation from NIST for 2s of O4")
ax[0].set_xlabel("% Deviation of Calculation from NIST")

ax[1].hist(shifts[:, 2], bins=50)
ax[1].set_title("% Deviation from NIST for 2p of O4")
ax[1].set_xlabel("% Deviation of Calculation from NIST")

fig.tight_layout()
fig.savefig("O4 NIST shifts")
"""

