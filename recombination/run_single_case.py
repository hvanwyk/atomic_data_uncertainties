#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 25 09:31:02 2021

@author: lochstu
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 17 10:51:07 2021

@author: loch
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 11:09:08 2021

@author: loch
"""
import numpy as np
import sys
import matplotlib.pyplot as plt

from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")
from bayesian_methods import log_posterior, interpolators, lambdas_grid
import time
from run_case import run_case


atom = "o"
seq = "li"
shell = "2-2"
ion = State(atom, seq, shell)
    
nist_cutoff=0.03
prior_shape="uniform"
likelihood_shape="gaussian"
  
ion = State(atom, seq, shell)
    
direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"    
file_name_common = f"_{int(nist_cutoff*100)}_{prior_shape}P_{likelihood_shape}L"
    
# Interval endpoints for each input component
x_bnd = np.array([[0.6,1.4],[0.6,1.4]])
    
# Resolution in each dimension
grid_resolution = 20
x_res = np.array([grid_resolution, grid_resolution])
    
X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    

n_walkers=100
n_steps=3000

print('Calculation for ',seq,'-like',atom)
print('Core excitation : ' , shell)
print('nist cutoff = ',nist_cutoff)
print('prior_shape = ',prior_shape)
print('likelihood shape = ',likelihood_shape)
print('Lambda range', x_bnd)
print('number of walkers = ', n_walkers)
print('number of steps per walker = ', n_steps)

#%%

start = time.time()

cent_pot=1
data_li_like_o_pos = run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps)

cent_pot=-1
data_li_like_o_neg = run_case(atom,seq,shell,ion,nist_cutoff, prior_shape,likelihood_shape,direc,file_name_common, x_bnd, x_res, X_1D, grid_resolution, cent_pot,n_walkers,n_steps)

end = time.time()
print(f"Runtime: {int(end-start)}s")


#%%


#Li-like O
rate_avg_pos=data_li_like_o_pos['rate_avg']
rate_avg_neg=data_li_like_o_neg['rate_avg']
rate_std_pos=data_li_like_o_pos['rate_std']
rate_std_neg=data_li_like_o_neg['rate_std']
rate_percent_pos=data_li_like_o_pos['rate_percent']
rate_percent_neg=data_li_like_o_neg['rate_percent']
rate_samples_pos=data_li_like_o_pos['rate_samples']
rate_samples_neg=data_li_like_o_neg['rate_samples']
rate_samples_combined=np.concatenate((rate_samples_pos,rate_samples_neg))
rate_avg_combined = np.average(rate_samples_combined,axis=0)
rate_std_combined = np.std(rate_samples_combined,axis=0)
rate_percent_combined = rate_std_combined/rate_avg_combined*100.

T=data_li_like_o_pos['T']
T_ev=T/11604.

indx=4
print('Results for T(K)=' ,T[indx], 'K :  T(eV)= ', np.around(T_ev[indx],decimals=2), ' eV')
print('Average: TF       = ', rate_avg_pos[indx])
print('Average: STO      = ', rate_avg_neg[indx])
print('Average: TF + STO = ', rate_avg_combined[indx])
print('Standard Deviation: TF       = ', rate_std_pos[indx])
print('Standard Deviation: STO      = ', rate_std_neg[indx])
print('Standard Deviation: TF + STO = ', rate_std_combined[indx])
print('Percent uncertainty: TF       = ', np.around(rate_percent_pos[indx],decimals=3))
print('Percent uncertainty: STO      = ', np.around(rate_percent_neg[indx],decimals=3))
print('Percent uncertainty: TF + STO = ', np.around(rate_percent_combined[indx],decimals=3))


indx=4
title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig1, ax1 = plt.subplots(1,1)
ax1.set_ylabel('Number of counts')
ax1.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax1.set_title(title)
ax1.hist(rate_samples_pos[:,indx],bins=100,density=True, ls='dotted',edgecolor='blue',alpha=0.5,label='TF',color='red')
ax1.hist(rate_samples_neg[:,indx],bins=100,density=True, ls='dashed',edgecolor='black',alpha=0.5,label='STO',color='blue')
ax1.legend()

indx=10
print('Results for T(K)=' ,T[indx], 'K :  T(eV)= ', np.around(T_ev[indx],decimals=2), ' eV')
print('Average: TF       = ', rate_avg_pos[indx])
print('Average: STO      = ', rate_avg_neg[indx])
print('Average: TF + STO = ', rate_avg_combined[indx])
print('Standard Deviation: TF       = ', rate_std_pos[indx])
print('Standard Deviation: STO      = ', rate_std_neg[indx])
print('Standard Deviation: TF + STO = ', rate_std_combined[indx])
print('Percent uncertainty: TF       = ', np.around(rate_percent_pos[indx],decimals=3))
print('Percent uncertainty: STO      = ', np.around(rate_percent_neg[indx],decimals=3))
print('Percent uncertainty: TF + STO = ', np.around(rate_percent_combined[indx],decimals=3))



indx=10
title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig2, ax2 = plt.subplots(1,1)
ax2.set_ylabel('Number of counts')
ax2.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax2.set_title(title)
ax2.hist(rate_samples_pos[:,indx],bins=100,density=True, ls='dotted',edgecolor='blue',alpha=0.5,label='TF',color='red')
ax2.hist(rate_samples_neg[:,indx],bins=100,density=True, ls='dashed',edgecolor='black',alpha=0.5,label='STO',color='blue')
ax2.legend()


fig3, axs3 = plt.subplots(2,1)
title_pos='Rate coefficient for ' + seq + '-like' + atom +' : TF and STO potentials'
axs3[0].set_xscale("log")
axs3[0].set_yscale("log")
axs3[0].set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs3[1].set_xlabel('Electron Temperature (K)')
axs3[1].set_ylabel('Percentage uncertainty')
axs3[0].set_title(title_pos)
axs3[0].plot(T,rate_avg_pos,'r',label='TF')
axs3[0].errorbar(T,rate_avg_pos,yerr=rate_std_pos*2.0,fmt='ro',ecolor='r')
axs3[0].plot(T,rate_avg_neg,'b',label='STO')
axs3[0].errorbar(T,rate_avg_neg,yerr=rate_std_neg*2.0,fmt='bo',ecolor='b')
axs3[0].plot(T,rate_avg_combined,'b',label='TF+STO')
axs3[0].errorbar(T,rate_avg_combined,yerr=rate_std_combined*2.0,fmt='go',ecolor='b')
axs3[1].set_xscale("log")
axs3[1].set_yscale("log")
axs3[1].plot(T,rate_percent_pos,'r',label='TF')
axs3[1].plot(T,rate_percent_neg,'b',label='STO')
axs3[1].plot(T,rate_percent_combined,'g',label='TF+STO')
axs3[0].legend()
axs3[1].legend()



#%%

np.savez('li_like_o',data_li_like_o_pos=data_li_like_o_pos,data_li_like_o_neg=data_li_like_o_neg)


