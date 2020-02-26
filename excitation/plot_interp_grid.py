#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 14:23:19 2020

@author: kyle
"""

import sys
if ".." not in sys.path:
    sys.path.append("..")
from state import State
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_adf04 import read_adf04
from bayesian_methods import lambdas_grid
from lambdas import make_energy_grid, make_rates_grid


atom = "o"
seq = "he"
shell = "1-2"

ion = State(atom, seq, shell)

nmax = 3

rates_grid = np.load(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rate_grid.npy")
# Resolution in each dimension for lambdas interpolation grid
grid_size = rates_grid.shape[2]


# Interval endpoints for each input component
n_lambdas = 2
x_bnd = []
for i in range(n_lambdas):
    x_bnd.append([0.8,1.2])
x_bnd = np.array(x_bnd)

x_res = np.array([grid_size]*n_lambdas)

X_1D, x_ravel = lambdas_grid(x_bnd, x_res)

lam_grid = np.array(x_ravel)



levels, T, _1, _2 = read_adf04("adf04")

w_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(1)1( 1.0)'), "#"].values[0])-2
x_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 2.0)'), "#"].values[0])-2
y_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 1.0)'), "#"].values[0])-2
z_num = int(levels.loc[(levels['config']=='1S1 2S1') & (levels['(2S+1)L( 2J)'] == '(3)0( 1.0)'), "#"].values[0])-2

"""

Err, Erg = make_energy_grid(ion, x_ravel, x_res)

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
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.scatter(x_ravel[:, 0], x_ravel[:, 1], err_grid[i,:,:])
    ax.set_xlabel("2s Lambda")
    ax.set_ylabel("2p Lambda")
    ax.set_zlabel(f"E$_{i}$")
    fig.tight_layout()
    fig.savefig(f"E{i}_grid.eps")


"""

fig_X = plt.figure()

ax = fig_X.add_subplot(211, projection="3d")
ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[x_num,0, :,:])
ax.set_xlabel("2s Lambda")
ax.set_ylabel("2p Lambda")
ax.set_zlabel("X A-Value")

ax = fig_X.add_subplot(212, projection="3d")
ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[x_num,1, :,:])
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
ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[z_num,1, :,:])
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
ax.scatter(x_ravel[:, 0], x_ravel[:, 1], rates_grid[w_num,1, :,:])
ax.set_xlabel("2s Lambda")
ax.set_ylabel("2p Lambda")
ax.set_zlabel(f"W Epsilon at T={T[9]}")

fig_W.tight_layout()
fig_W.savefig(f"W_grid.eps")

