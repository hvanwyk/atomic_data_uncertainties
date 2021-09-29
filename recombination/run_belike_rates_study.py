#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 17:59:04 2021

@author: lochstu
"""

import numpy as np
import sys
import matplotlib.pyplot as plt
import corner 
from scipy import stats


from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import log_posterior, interpolators, lambdas_grid
import time
import pickle
import os
from utilities import create_directories, get_nist_energy, read_levels, compare_to_nist
from run_case import run_case

#atoms=['b','c','n','o','f','ne','na','mg','al','si','p','s','cl','k','ca','sc','ti','v','cr','mn','fe']
atoms=['o']
seq = "be"
shell = "2-2"
    
nist_cutoff=0.03
prior_shape="uniform"
likelihood_shape="gaussian"
  
#Set this to true for the post-processor to shift series limits to NIST energies.
nist_shift=True
    
    
# Interval endpoints for each input component
x_bnd = np.array([[0.4,1.6],[0.4,1.6]])
    
# Resolution in each dimension
grid_resolution = 40
x_res = np.array([grid_resolution, grid_resolution])
    
X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
n_walkers=100
n_steps=3000

up_dir=res=os.getcwd()


#filename="be_like_all_pt4_1pt6_nist_shifts.pkl"
filename="be_like_o_fine_lambda_mesh.pkl"


#%%

data={"seq":seq,"shell":shell,"x_bnd":x_bnd,"nist_cutoff":nist_cutoff,"prior_shape":prior_shape,"likelihood_shape":likelihood_shape,
                 "x_bnd":x_bnd,"x_res":x_res,"X_1D":X_1D,"x_ravel":x_ravel,"grid_resolution":grid_resolution,"n_walkers":n_walkers,"n_steps":n_steps}

for i in range(len(atoms)):
    start = time.time()
    atom = atoms[i]
    ion = State(atom, seq, shell)
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"    
    file_name_common = f"_{int(nist_cutoff*100)}_{prior_shape}P_{likelihood_shape}L"
    nist_vals = get_nist_energy(up_dir + f"/NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist")
    emax=nist_vals[-1]*1.1  
    print('Calculation for ',seq,'-like',atom)
    print('Core excitation : ' , shell)
    print('nist cutoff = ',nist_cutoff)
    print('prior_shape = ',prior_shape)
    print('likelihood shape = ',likelihood_shape)
    print('Lambda range', x_bnd)
    print('number of walkers = ', n_walkers)
    print('number of steps per walker = ', n_steps)
    print('EMAX=',emax)
    print('up_dir=',up_dir)


    cent_pot=1
    data_pos = run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps,emax,up_dir,nist_shift=nist_shift)
    print(' ')
    cent_pot=-1
    data_neg = run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps,emax,up_dir,nist_shift=nist_shift)
    pos_name=seq+'_like_'+atom+'_pos'
    neg_name=seq+'_like_'+atom+'_neg'
    tmp1={pos_name:data_pos}
    tmp2={neg_name:data_neg}
    data.update(tmp1)
    data.update(tmp2)
    print('updating dictionary')
    end = time.time()
    print(f"Runtime: {int(end-start)}s")
    print(' ')
    print(' ')


a_file = open(filename, "wb")
pickle.dump(data, a_file)
a_file.close()

#%%
data=pickle.load(open('be_like_o_fine_lambda_mesh.pkl','rb'))

Lambda_1=data['X_1D'][0][:]
Lambda_2=data['X_1D'][1][:]
Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)

fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
surf = ax1.plot_wireframe(Lambda_1,Lambda_2, data['be_like_o_pos']['rates'][2][:][:],color='black',linewidth=2)

n_lambdas=2
fig2=corner.corner(data['be_like_o_pos']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])


rate_samples_combined=np.concatenate((data['be_like_o_pos']['rate_samples'],data['be_like_fe_neg']['rate_samples']))

fig3, axs3 = plt.subplots(1,1)
title_pos='Rate coefficients'
axs3.set_yscale("log")
axs3.set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs3.set_xlabel('Electrson Temperature index')
axs3.set_title(title_pos)
axs3.boxplot(rate_samples_combined,showfliers=False)
axs3.boxplot(data['be_like_o_pos']['rate_samples'],showfliers=False)
axs3.boxplot(data['be_like_o_neg']['rate_samples'],showfliers=False)
axs3.legend()
plt.show()

indx=3
fig2, ax2 = plt.subplots(1,1)
ax2.set_ylabel('Number of counts')
ax2.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax2.set_xscale("log")
ax2.hist(data['be_like_o_pos']['rate_samples'][:,indx],bins=10000,density=True, ls='dotted',alpha=0.5,label='TF')
#ax2.hist(data['be_like_o_neg']['rate_samples'][:,indx],bins=100,density=True, ls='dashed',edgecolor='black',alpha=0.5,label='STO',color='blue')
ax2.legend()



fig3 = plt.figure()
ax = fig3.gca(projection='3d')
surf = ax.plot_wireframe(Lambda_1,Lambda_2, data['be_like_o_pos']['rates'][indx][:][:],color='black',linewidth=2)
ax.scatter3D(data['be_like_o_pos']['lambda_samples'][:,1], data['be_like_o_pos']['lambda_samples'][:,0], data['be_like_o_pos']['rate_samples'][:,indx])

fig4=corner.corner(data['be_like_o_neg']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])

