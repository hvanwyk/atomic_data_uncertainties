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


#atoms=['b','c','n','o','f','ne','na','mg','al','si','k','ti','cr','fe']
#atoms=['ca','sc','ti','v','cr','mn']
atoms=['o','fe']

seq = "be"
shell = "2-2"
    
nist_cutoff=0.03
prior_shape="uniform"
likelihood_shape="gaussian"
  
#Set this to true for the post-processor to shift series limits to NIST energies.
nist_shift=True
    
    
# Interval endpoints for each input component
x_bnd = np.array([[0.6,1.4],[0.6,1.4]])
    
# Resolution in each dimension
grid_resolution = 20
x_res = np.array([grid_resolution, grid_resolution])
    
X_1D, x_ravel = lambdas_grid(x_bnd, x_res)
    
n_walkers=100
n_steps=3000

up_dir=res=os.getcwd()


filename="be_like_o_fe_0.6_1.4_nist_shifts.pkl"


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
#Make plots for O details of each calculation


#Energies compared with NIST range
data=pickle.load(open('be_like_o_fe_0.4_1.6_nist_shifts.pkl','rb'))
indx=3
Lambda_1=data['X_1D'][0][:]
Lambda_2=data['X_1D'][1][:]
Egy_vals_neg=data['be_like_o_neg']['Erg'][indx][:][:]
#Egy_vals_pos=data['be_like_o_pos']['Erg'][indx][:][:]
Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)

nist_vals=np.empty([grid_resolution,grid_resolution])
nist_vals_max=np.empty([grid_resolution,grid_resolution])
nist_vals_min=np.empty([grid_resolution,grid_resolution])
nist_vals.fill(data['be_like_o_pos']['nist_vals'][indx])
nist_vals_max.fill(data['be_like_o_pos']['nist_vals'][indx]*(1.0+nist_cutoff))
nist_vals_min.fill(data['be_like_o_pos']['nist_vals'][indx]*(1.0-nist_cutoff))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('Lambda 2s')
ax.set_ylabel('Lambda 2p')
ax.set_zlabel('Energy (Ryd)')
surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_neg,color='black',linewidth=2)
#surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=1)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals,color='r',linewidth=0.7)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals_max,color='g',linewidth=0.7)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals_min,color='b',linewidth=0.7)
plt.show()


#Lambda distribution function
n_lambdas=2
fig=corner.corner(data['be_like_o_neg']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
plt.show()

#Energy distribution function
n_energies=10
fig2=corner.corner(data['be_like_o_pos']['E_samples'][:,1:], labels=[f"$E_{i+1}$" for i in range(n_energies-1)], truths=[0 for i in range(n_energies-1)])


#%%
#Make plots for Fe details of each calculation


#Energies compared with NIST range
#data=pickle.load(open('be_like_o_fe_0.4_1.6_nist_shifts.pkl','rb'))
indx=3
Lambda_1=data['X_1D'][0][:]
Lambda_2=data['X_1D'][1][:]
Egy_vals_neg=data['be_like_fe_neg']['Erg'][indx][:][:]
#Egy_vals_pos=data['be_like_o_pos']['Erg'][indx][:][:]
Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)

nist_vals=np.empty([grid_resolution,grid_resolution])
nist_vals_max=np.empty([grid_resolution,grid_resolution])
nist_vals_min=np.empty([grid_resolution,grid_resolution])
nist_vals.fill(data['be_like_fe_pos']['nist_vals'][indx])
nist_vals_max.fill(data['be_like_fe_pos']['nist_vals'][indx]*(1.0+nist_cutoff))
nist_vals_min.fill(data['be_like_fe_pos']['nist_vals'][indx]*(1.0-nist_cutoff))

fig = plt.figure()
ax = fig.gca(projection='3d')
ax.set_xlabel('Lambda 2s')
ax.set_ylabel('Lambda 2p')
ax.set_zlabel('Energy (Ryd)')
surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_neg,color='black',linewidth=2)
#surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=1)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals,color='r',linewidth=0.7)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals_max,color='g',linewidth=0.7)
surf = ax.plot_wireframe(Lambda_1,Lambda_2, nist_vals_min,color='b',linewidth=0.7)
plt.show()


#Lambda distribution function
n_lambdas=2
fig=corner.corner(data['be_like_fe_neg']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
plt.show()

#Energy distribution function
n_energies=10
fig2=corner.corner(data['be_like_fe_pos']['E_samples'][:,1:], labels=[f"$E_{i+1}$" for i in range(n_energies-1)], truths=[0 for i in range(n_energies-1)])

#%%
#Make plots for final rates with uncertainties

#Be-like O
atom='o'
rate_avg_pos=data['be_like_o_pos']['rate_avg']
rate_avg_neg=data['be_like_o_neg']['rate_avg']
rate_std_pos=data['be_like_o_pos']['rate_std']
rate_std_neg=data['be_like_o_neg']['rate_std']
rate_percent_pos=data['be_like_o_pos']['rate_percent']
rate_percent_neg=data['be_like_o_neg']['rate_percent']
rate_samples_pos=data['be_like_o_pos']['rate_samples']
rate_samples_neg=data['be_like_o_neg']['rate_samples']
T=data['be_like_o_pos']['T']


rate_samples_combined=np.concatenate((rate_samples_pos,rate_samples_neg))
#rate_avg_combined = np.average(rate_samples_combined,axis=0)
rate_std_combined = np.std(rate_samples_combined,axis=0)
rate_mode_combined = stats.mode(rate_samples_combined,axis=0)
rate_median_combined = np.median(rate_samples_combined,axis=0)

#The following is to maek an average that does not include zeros.
#This can be a problem for low charge states

rate_samples_2=rate_samples_combined
rate_samples_2[rate_samples_2 == 0] = np.nan
rate_avg_combined =np.nanmean(rate_samples_2, axis=0)
rate_percent_combined = (rate_std_combined)/rate_avg_combined*100.



T_ev=T/11604.

indx=3
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
print(' ')

title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig1, ax1 = plt.subplots(1,1)
ax1.set_ylabel('Number of counts')
ax1.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax1.set_xscale("log")
ax1.set_title(title)
ax1.hist(rate_samples_pos[:,indx],bins=10000,density=True, ls='dotted',edgecolor='blue',alpha=0.5,label='TF',color='red')
ax1.hist(rate_samples_neg[:,indx],bins=10000,density=True, ls='dashed',edgecolor='black',alpha=0.5,label='STO',color='blue')
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
print(' ')


indx=10
title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig2, ax2 = plt.subplots(1,1)
ax2.set_ylabel('Number of counts')
ax2.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax2.set_title(title)
ax2.set_xscale("log")
ax2.hist(rate_samples_pos[:,indx],bins=100,density=True, ls='dotted',edgecolor='blue',alpha=0.5,label='TF',color='red')
ax2.hist(rate_samples_neg[:,indx],bins=100,density=True, ls='dashed',edgecolor='black',alpha=0.5,label='STO',color='blue')
ax2.legend()


rate_samples=rate_samples_combined[:,0:]
str_T=T[0:]
unc_low=(np.quantile(rate_samples,0.5,axis=0)-np.quantile(rate_samples,0.25,axis=0))/np.quantile(rate_samples,0.5,axis=0)
unc_high=(np.quantile(rate_samples,0.75,axis=0)-np.quantile(rate_samples,0.5,axis=0))/np.quantile(rate_samples,0.5,axis=0)
fig3, axs3 = plt.subplots(2,1)
title_pos='Rate coefficient for ' + seq + '-like' + atom +' : TF and STO potentials'
axs3[0].set_xscale("log")
axs3[0].set_yscale("log")
axs3[0].set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs3[1].set_xlabel('Electron Temperature (K)')
axs3[1].set_ylabel('Percentage uncertainty')
axs3[0].set_title(title_pos)
#axs3[0].plot(T,rate_avg_pos,'r',label='TF')
#axs3[0].errorbar(T,rate_avg_pos,yerr=rate_std_pos*2.0,fmt='ro',ecolor='r')
#axs3[0].plot(T,rate_avg_neg,'b',label='STO')
#axs3[0].errorbar(T,rate_avg_neg,yerr=rate_std_neg*2.0,fmt='bo',ecolor='b')
axs3[0].plot(T,rate_median_combined,'r',label='Median')
axs3[0].plot(T,rate_median_combined-rate_median_combined*unc_low,'g',label='1st quartile')
axs3[0].plot(T,rate_median_combined+rate_median_combined*unc_high,'b',label='3rd quartile')
#axs3[0].plot(T,rate_mode_combined[0][0][:],'g',label='TF+STO mode')
#axs3[0].errorbar(T,rate_avg_combined,yerr=rate_std_combined*1.0,fmt='go',ecolor='b')
axs3[1].set_xscale("log")
axs3[1].set_yscale("log")
axs3[1].plot(T,unc_low*100,'g',label='Uncert. 1st quartile')
axs3[1].plot(T,unc_high*100,'b',label='Uncert. 3rd quartile')
axs3[1].plot(T,rate_percent_combined,'black',label='stddev/avg')
axs3[0].legend()
axs3[1].legend()


seq='Be'
#data=pickle.load(open('be_like_o_fe_0.4_1.6_nist_shifts.pkl','rb'))
fig4, axs4 = plt.subplots(1,1)
title_pos='Rate coefficient for ' + seq + '-like '  + atom 
#axs4.set_xscale("log")
axs4.set_yscale("log")
axs4.set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs4.set_xlabel('Electron Temperature index')
axs4.set_title(title_pos)
#axs4.set_xticks(T)
#axs4.set_xticks([T[y] for y in range(len(T))])
#axs4.boxplot(rate_samples[:,1:],positions=T[1:],widths=2.0)
axs4.boxplot(rate_samples,showfliers=False)
#axs4.plot(T,rate_median_combined)
#axs4.boxplot(T[1],rate_samples[3,:])
#ax4.set_xticklabels(string(T))
axs4.legend()
plt.show()

median=np.median(rate_samples,axis=0)
iqr=stats.iqr(rate_samples,axis=0)
uncert=iqr/median*100.

#%% Fe-plots

atom='fe'
rate_avg_pos=data['be_like_fe_pos']['rate_avg']
rate_avg_neg=data['be_like_fe_neg']['rate_avg']
rate_std_pos=data['be_like_fe_pos']['rate_std']
rate_std_neg=data['be_like_fe_neg']['rate_std']
rate_percent_pos=data['be_like_fe_pos']['rate_percent']
rate_percent_neg=data['be_like_fe_neg']['rate_percent']
rate_samples_pos=data['be_like_fe_pos']['rate_samples']
rate_samples_neg=data['be_like_fe_neg']['rate_samples']
T=data['be_like_fe_pos']['T']

rate_samples_combined=np.concatenate((rate_samples_pos,rate_samples_neg))
#rate_avg_combined = np.average(rate_samples_combined,axis=0)
rate_std_combined = np.std(rate_samples_combined,axis=0)
#rate_percent_combined = rate_std_combined/rate_avg_combined*100.

rate_samples_2=rate_samples_combined
rate_samples_2[rate_samples_2 == 0] = np.nan
rate_avg_combined =np.nanmean(rate_samples_2, axis=0)
rate_percent_combined = (rate_std_combined)/rate_avg_combined*100.



T_ev=T/11604.

indx=3
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
print(' ')

title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig1, ax1 = plt.subplots(1,1)
ax1.set_ylabel('Number of counts')
ax1.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax1.set_xscale("log")
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
print(' ')

title='Results for Be-like ' + atom + ' for T='+T.astype(str)[indx]+' K or '+(T_ev.astype(str)[indx])+' eV'
fig2, ax2 = plt.subplots(1,1)
ax2.set_ylabel('Number of counts')
ax2.set_xlabel('Rate coefficient (cm^3 s^-1)')
ax2.set_title(title)
ax2.set_xscale("log")
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
#axs3[0].plot(T,rate_avg_pos,'r',label='TF')
#axs3[0].errorbar(T,rate_avg_pos,yerr=rate_std_pos*2.0,fmt='ro',ecolor='r')
#axs3[0].plot(T,rate_avg_neg,'b',label='STO')
#axs3[0].errorbar(T,rate_avg_neg,yerr=rate_std_neg*2.0,fmt='bo',ecolor='b')
axs3[0].plot(T,rate_avg_combined,'b',label='TF+STO avg')
axs3[0].plot(T,rate_median_combined,'r',label='TF+STO median')
axs3[0].plot(T,rate_mode_combined[0][0][:],'g',label='TF+STO mode')
#axs3[0].errorbar(T,rate_avg_combined,yerr=rate_std_combined*1.0,fmt='go',ecolor='b')
axs3[1].set_xscale("log")
axs3[1].set_yscale("log")
axs3[1].plot(T,rate_percent_pos,'r',label='TF')
axs3[1].plot(T,rate_percent_neg,'b',label='STO')
axs3[1].plot(T,rate_percent_combined,'g',label='TF+STO')
axs3[0].legend()
axs3[1].legend()



seq='Be'
#data=pickle.load(open('be_like_o_fe_0.4_1.6_nist_shifts.pkl','rb'))
rate_samples=rate_samples_combined[:,0:]
str_T=T[0:]
fig4, axs4 = plt.subplots(1,1)
title_pos='Rate coefficient for ' + seq + '-like '  + atom 
#axs4.set_xscale("log")
axs4.set_yscale("log")
axs4.set_ylabel('Rate Coefficeint (cm^3 s^-1)')
axs4.set_xlabel('Electron Temperature index')
axs4.set_title(title_pos)
#axs4.set_xticks(T)
#axs4.set_xticks([T[y] for y in range(len(T))])
#axs4.boxplot(rate_samples[:,1:],positions=T[1:],widths=2.0)
axs4.boxplot(rate_samples,showfliers=False)
#axs4.plot(T,rate_median_combined)
#axs4.boxplot(T[1],rate_samples[3,:])
#ax4.set_xticklabels(string(T))
axs4.legend()
plt.show()