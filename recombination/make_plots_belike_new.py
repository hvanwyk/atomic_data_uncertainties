#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 25 19:54:05 2021

@author: loch
"""

import numpy as np
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator,FormatStrFormatter,MaxNLocator
from mpl_toolkits.mplot3d import Axes3D

import corner 
from scipy import stats


from recombination_methods import State, structure, structure_dr, postprocessing_rates, get_rate
if "../src/" not in sys.path:
    sys.path.append("../src/")
if ".." not in sys.path:
    sys.path.append("..")

import pickle

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.unicode_minus'] = False

plt.rcParams['text.usetex'] = True #Let TeX do the typsetting
plt.rcParams['text.latex.preamble'] = [r'\usepackage{sansmath}',r'\sansmath']
#Force sans-serif math mode
plt.rcParams['font.family'] = 'sans-serif' # ... for regular text
plt.rcParams['font.sans-serif'] = 'Helvetica' # Choose a nice font here


#%%
def adjacent_values(vals, q1, q3):
    upper_adjacent_value = q3 + (q3 - q1) * 1.5
    upper_adjacent_value = np.clip(upper_adjacent_value, q3, vals[-1])

    lower_adjacent_value = q1 - (q3 - q1) * 1.5
    lower_adjacent_value = np.clip(lower_adjacent_value, vals[0], q1)
    return lower_adjacent_value, upper_adjacent_value

#%%
data=pickle.load(open('be_like_o_fe_0.4_1.6_40grid_nist_0.03_uniform_gaussian.pkl','rb'))

#%%
#Plots for O4+
Lambda_1=data['be_like_o_pos']['X_1D'][0][:]
Lambda_2=data['be_like_o_pos']['X_1D'][1][:]
Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)

lambda_labels=["$\lambda_{2s}$","$\lambda_{2p}$"]

n_lambdas=2
#fig0=corner.corner(data['be_like_o_pos']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
fig0=corner.corner(data['be_like_o_pos']['lambda_samples'], labels=[lambda_labels[i] for i in range(n_lambdas)], label_kwargs={"fontsize": 14}, truths=[1 for i in range(n_lambdas)])
plt.savefig('lambda_distr_Be_like_O.eps')

egy_labels=["$E_{2s2p (^3P_0)}$","$E_{2s2p (^3P_1)}$","$E_{2s2p (^3P_2)}$","$E_{2s2p (^1P_1)}$","$E_{2p^2 (^3P_0)}$","$E_{2p^2 (^3P_1)}$","$E_{2p^2 (^3P_2)}$","$E_{2p^2 (^1D_2)}$","$E_{2p^2 (^1S_0)}$"]
n_energies=10
fig1=corner.corner(data['be_like_o_pos']['E_samples'][:,1:], labels=[egy_labels[i] for i in range(n_energies-1)], label_kwargs={"fontsize": 14}, quantiles=(0.25, 0.5, 0.75),show_titles='True',truths=[0 for i in range(n_energies-1)])
plt.savefig('energy_distr_Be_like_O.eps')

indx=2
grid_resolution=40
nist_cutoff=0.03
nist_vals=np.empty([grid_resolution,grid_resolution])
nist_vals_max=np.empty([grid_resolution,grid_resolution])
nist_vals_min=np.empty([grid_resolution,grid_resolution])
nist_vals.fill(data['be_like_o_pos']['nist_vals'][indx])
nist_vals_max.fill(data['be_like_o_pos']['nist_vals'][indx]*(1.0+nist_cutoff))
nist_vals_min.fill(data['be_like_o_pos']['nist_vals'][indx]*(1.0-nist_cutoff))
Egy_vals_pos=data['be_like_o_pos']['Erg'][indx][:][:]
Egy_vals_pos[Egy_vals_pos>0.85]='Nan'
Egy_vals_pos[Egy_vals_pos<0.65]='Nan'

         
fig2 = plt.figure(figsize=(14,10))
ax2 = fig2.gca(projection='3d')
ax2.set_xlabel('$\lambda_{2s}$', fontsize=24)
ax2.set_ylabel('$\lambda_{2p}$', fontsize=24)
ax2.set_zlabel('Energy (Ryd)', fontsize=18)
ax2.set_zlim(0.65, 0.85)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=2.0)
#surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=1)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals,color='r',linewidth=0.7)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals_max,color='g',linewidth=0.7)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals_min,color='b',linewidth=0.7)
plt.show()
plt.savefig('energy_grid_Be_like_O.eps')


rate_samples_combined=np.concatenate((data['be_like_o_pos']['rate_samples'],data['be_like_o_neg']['rate_samples']))
labels=np.round(np.log10(data['be_like_o_pos']['T'][3:]),2)

##Boxplot
#fig3, axs3 = plt.subplots(1,1)
#title_pos='Rate coefficients'
#axs3.set_yscale("log")
#axs3.set_ylabel('Rate Coefficient (cm$^3$ $s^{-1}$)', fontsize=16)
#axs3.set_xlabel('Log$_{10}$ (Electron Temperature (K))', fontsize=16)
#axs3.set_title(title_pos)
#axs3.boxplot(rate_samples_combined[:,3:],showfliers=False,labels=labels)
##axs3.boxplot(data['be_like_o_pos']['rate_samples'][:,3:],showfliers=False)
##axs3.boxplot(data['be_like_o_neg']['rate_samples'][:,3:],showfliers=False)
#axs3.legend()
#plt.show()

quartile1, medians, quartile3 = np.percentile(rate_samples_combined[:,:], [25, 50, 75], axis=0)
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(rate_samples_combined[:,:], quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)

IQR = stats.iqr(rate_samples_combined[:,:], axis=0,interpolation = 'midpoint')
median = np.median(rate_samples_combined, axis=0)
Q1 = np.percentile(rate_samples_combined[:,:],  25, axis=0, interpolation = 'midpoint')
Q3 = np.percentile(rate_samples_combined[:,:],  75, axis=0, interpolation = 'midpoint')
upper_whisker = np.min(rate_samples_combined[:,:],axis=0)
lower_whisker = np.max(rate_samples_combined[:,:],axis=0)
Uncert=IQR[3:]/2./median[3:]*100


fig4, axs4 = plt.subplots(1,1)
axs4.set_ylabel('Log$_{10}$(Rate Coeff. (cm$^3$ s$^{-1}$))',fontsize=12)
axs4.set_xlabel('Log$_{10}$(Electron Temperature (K))',fontsize=12)
#axs4.set_yscale("log")
#axs4.set_xscale("log")
axs4.violinplot(np.log10(rate_samples_combined[:,3:]), showmeans=False, showmedians=True,
        showextrema=True,positions=np.log10(data['be_like_o_pos']['T'][3:]),widths=0.3)
axs4.scatter(np.log10(data['be_like_o_pos']['T'][3:]), np.log10(medians[3:]), marker='o', color='green', s=20, zorder=3)
axs4.vlines(np.log10(data['be_like_o_pos']['T'][3:]), np.log10(quartile1[3:]), np.log10(quartile3[3:]), color='k', linestyle='-', lw=2)
axs4.legend()
plt.show()
plt.savefig('violin_plot_Be_like_O.eps')

#indx=3
#fig5 = plt.figure()
#ax5 = fig5.gca(projection='3d')
#surf = ax5.plot_wireframe(Lambda_1,Lambda_2, np.log10(data['be_like_o_pos']['rates'][indx][:][:]),color='black',linewidth=0.5)
#ax5.scatter3D(data['be_like_o_pos']['lambda_samples'][:,1], data['be_like_o_pos']['lambda_samples'][:,0], np.log10(data['be_like_o_pos']['rate_samples'][:,indx]))
##surf = ax1.plot_wireframe(Lambda_1,Lambda_2, data['be_like_o_neg']['rates'][indx][:][:],color='red',linewidth=0.5)
#ax5.set_xlabel('Lambda 2s')
#ax5.set_ylabel('Lambda 2p')
#ax5.set_zlabel('log$_{10}$(Rate Coefficient cm$^3$ s$^{-1}$)')
#ax5.legend()
#plt.show()
#plt.savefig('DR_rate_grid_plus_points_Be_like_O_low_Te.eps')
#
#

#Note red is TFDA and blue is Hartree potential with STO
indx=3
fig6, ax6 = plt.subplots(1,1)
ax6.set_ylabel('Number of counts',fontsize=16)
ax6.set_xlabel('log$_{10}$(Rate coefficient (cm$^3$ s${^-1}$))',fontsize=16)
#ax6.set_xscale("log")
ax6.hist(np.log10(data['be_like_o_pos']['rate_samples'][:,indx]),bins=1000,density=True, ls='dotted',alpha=0.5,color='red')
ax6.hist(np.log10(data['be_like_o_neg']['rate_samples'][:,indx]),bins=1000,density=True, ls='dashed',alpha=0.5,color='blue')
ax6.legend()
plt.show()
plt.savefig('DR_rate_distr_Be_like_O_low_Te.eps')


indx=10
fig7, ax7 = plt.subplots(1,1)
ax7.set_ylabel('Number of counts',fontsize=16)
ax7.set_xlabel('Rate coefficient (cm$^3$ s${^-1}$)',fontsize=16)
#ax4.set_xscale("log")
ax7.hist(data['be_like_o_pos']['rate_samples'][:,indx],bins=100,density=True, ls='dotted',alpha=0.5,color='red')
ax7.hist(data['be_like_o_neg']['rate_samples'][:,indx],bins=100,density=True, ls='dashed',alpha=0.5,color='blue')
ax7.legend()
plt.show()
plt.savefig('DR_rate_distr_Be_like_O_high_Te.eps')



#%%
#Plots for Fe22+

Lambda_1=data['be_like_fe_pos']['X_1D'][0][:]
Lambda_2=data['be_like_fe_pos']['X_1D'][1][:]
Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)

lambda_labels=["$\lambda_{2s}$","$\lambda_{2p}$"]
n_lambdas=2
#fig0=corner.corner(data['be_like_fe_pos']['lambda_samples'], labels=[f"$lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
fig0=corner.corner(data['be_like_fe_pos']['lambda_samples'], labels=[lambda_labels[i] for i in range(n_lambdas)], label_kwargs={"fontsize": 14}, truths=[1 for i in range(n_lambdas)])
plt.savefig('lambda_distr_Be_like_Fe.eps')

egy_labels=["$E_{2s2p (^3P_0)}$","$E_{2s2p (^3P_1)}$","$E_{2s2p (^3P_2)}$","$E_{2s2p (^1P_1)}$","$E_{2p^2 (^3P_0)}$","$E_{2p^2 (^3P_1)}$","$E_{2p^2 (^3P_2)}$","$E_{2p^2 (^1D_2)}$","$E_{2p^2 (^1S_0)}$"]
n_energies=10
fig1=corner.corner(data['be_like_fe_pos']['E_samples'][:,1:], labels=[egy_labels[i] for i in range(n_energies-1)], label_kwargs={"fontsize": 14}, quantiles=(0.25, 0.5, 0.75),show_titles='True',truths=[0 for i in range(n_energies-1)])
plt.savefig('energy_distr_Be_like_Fe.eps')

indx=2
grid_resolution=40
nist_cutoff=0.03
nist_vals=np.empty([grid_resolution,grid_resolution])
nist_vals_max=np.empty([grid_resolution,grid_resolution])
nist_vals_min=np.empty([grid_resolution,grid_resolution])
nist_vals.fill(data['be_like_fe_pos']['nist_vals'][indx])
nist_vals_max.fill(data['be_like_fe_pos']['nist_vals'][indx]*(1.0+nist_cutoff))
nist_vals_min.fill(data['be_like_fe_pos']['nist_vals'][indx]*(1.0-nist_cutoff))
Egy_vals_pos=data['be_like_fe_pos']['Erg'][indx][:][:]
   
fig2 = plt.figure(figsize=(14,10))
ax2 = fig2.gca(projection='3d')
ax2.set_xlabel('$\lambda_{2s}$', fontsize=24)
ax2.set_ylabel('$\lambda_{2p}$', fontsize=24)
ax2.set_zlabel('Energy (Ryd)', fontsize=18)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=2.0)
#surf = ax.plot_wireframe(Lambda_1,Lambda_2, Egy_vals_pos,color='black',linewidth=1)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals,color='r',linewidth=0.7)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals_max,color='g',linewidth=0.7)
surf = ax2.plot_wireframe(Lambda_1,Lambda_2, nist_vals_min,color='b',linewidth=0.7)
plt.show()
plt.savefig('energy_grid_Be_like_Fe.eps')


rate_samples_combined=np.concatenate((data['be_like_fe_pos']['rate_samples'],data['be_like_fe_neg']['rate_samples']))
labels=np.round(np.log10(data['be_like_fe_pos']['T']),2)

##Boxplot
#fig3, axs3 = plt.subplots(1,1)
#title_pos='Rate coefficients'
#axs3.set_yscale("log")
#axs3.set_ylabel('Rate Coefficient (cm$^3$ $s^{-1}$)', fontsize=16)
#axs3.set_xlabel('Log$_{10}$ (Electron Temperature (K))', fontsize=16)
#axs3.set_title(title_pos)
#axs3.boxplot(rate_samples_combined[:,3:],showfliers=False,labels=labels)
##axs3.boxplot(data['be_like_fe_pos']['rate_samples'][:,3:],showfliers=False)
##axs3.boxplot(data['be_like_fe_neg']['rate_samples'][:,3:],showfliers=False)
#axs3.legend()
#plt.show()

quartile1, medians, quartile3 = np.percentile(rate_samples_combined[:,:], [25, 50, 75], axis=0)
whiskers = np.array([
    adjacent_values(sorted_array, q1, q3)
    for sorted_array, q1, q3 in zip(rate_samples_combined[:,:], quartile1, quartile3)])
whiskers_min, whiskers_max = whiskers[:, 0], whiskers[:, 1]
inds = np.arange(1, len(medians) + 1)

IQR = stats.iqr(rate_samples_combined[:,:], axis=0,interpolation = 'midpoint')
median = np.median(rate_samples_combined, axis=0)
Q1 = np.percentile(rate_samples_combined[:,:],  25, axis=0, interpolation = 'midpoint')
Q3 = np.percentile(rate_samples_combined[:,:],  75, axis=0, interpolation = 'midpoint')
upper_whisker = np.min(rate_samples_combined[:,:],axis=0)
lower_whisker = np.max(rate_samples_combined[:,:],axis=0)
Uncert=IQR/2./median*100


# fig4, axs4 = plt.subplots(2,1)
# axs4[0].set_ylabel('Rate Coeff. (cm$^3$ s$^{-1}$)',fontsize=12)
# #axs4[0].set_xlabel('Log$_{10}($Electron Temperature (K))',fontsize=16)
# #axs4[0].set_title(title_pos)
# axs4[0].xaxis.set_tick_params(direction='out')
# axs4[0].xaxis.set_ticks_position('bottom')
# axs4[0].set_xticks(np.arange(1, len(labels) + 1))
# axs4[0].set_xticklabels(labels)
# axs4[0].violinplot(np.log10(rate_samples_combined[:,3:]), showmeans=False, showmedians=True,
#         showextrema=True)
# axs4[0].scatter(inds[3:]-3, np.log10(medians[3:]), marker='o', color='green', s=30, zorder=3)
# axs4[0].vlines(inds[3:]-3, np.log10(quartile1[3:]), np.log10(quartile3[3:]), color='k', linestyle='-', lw=5)
# axs4[0].legend()
# axs4[1].set_ylabel('(IQR/2)/median*100',fontsize=12)
# axs4[1].set_xlabel('Log$_{10}$(Electron temperature (K))',fontsize=16)
# #axs4[1].set_yscale("log")
# axs4[1].xaxis.set_tick_params(direction='out')
# axs4[1].xaxis.set_ticks_position('bottom')
# axs4[1].set_xticks(np.arange(1, len(labels) + 1))
# axs4[1].set_xticklabels(labels)
# axs4[1].plot(inds[3:]-3,Uncert,'r')
# axs4[1].legend()
# plt.show()
# plt.savefig('violin_plot_Be_like_Fe.eps')

fig4, axs4 = plt.subplots(1,1)
axs4.set_ylabel('Log$_{10}$(Rate Coeff. (cm$^3$ s$^{-1}$))',fontsize=12)
axs4.set_xlabel('Log$_{10}($Electron Temperature (K))',fontsize=12)
#axs4.set_yscale("log")
#axs4.set_xscale("log")
axs4.violinplot(np.log10(rate_samples_combined), showmeans=False, showmedians=True,
        showextrema=True,positions=np.log10(data['be_like_o_pos']['T']),widths=0.3)
axs4.scatter(np.log10(data['be_like_o_pos']['T']), np.log10(medians), marker='o', color='green', s=5, zorder=3)
axs4.vlines(np.log10(data['be_like_o_pos']['T']), np.log10(quartile1), np.log10(quartile3), color='k', linestyle='-', lw=2)
axs4.legend()
plt.show()
plt.savefig('violin_plot_Be_like_Fe.eps')

#
#indx=3
#fig5 = plt.figure()
#ax5 = fig5.gca(projection='3d')
#surf = ax5.plot_wireframe(Lambda_1,Lambda_2, np.log10(data['be_like_fe_pos']['rates'][indx][:][:]),color='black',linewidth=0.5)
#ax5.scatter3D(data['be_like_fe_pos']['lambda_samples'][:,1], data['be_like_fe_pos']['lambda_samples'][:,0], np.log10(data['be_like_fe_pos']['rate_samples'][:,indx]))
##surf = ax1.plot_wireframe(Lambda_1,Lambda_2, data['be_like_fe_neg']['rates'][indx][:][:],color='red',linewidth=0.5)
#ax5.set_xlabel('$\lambda_{2s}$')
#ax5.set_ylabel('$\lambda_{2p$}$')
#ax5.set_zlabel('log$_{10}$(Rate Coefficient cm$^3$ s$^{-1}$)')
#plt.show()
#plt.savefig('DR_rate_grid_plus_points_Be_like_Fe_low_Te.eps')
#
#




#Note red is TFDA and blue is Hartree potential with STO
indx=3
fig6, ax6 = plt.subplots(1,1)
ax6.set_ylabel('Number of counts',fontsize=16)
ax6.set_xlabel('log$_{10}$(Rate coefficient (cm$^3$ s${^-1}$))',fontsize=16)
#ax6.set_xscale("log")
ax6.hist(np.log10(data['be_like_fe_pos']['rate_samples'][:,indx]),bins=100,density=True, ls='dotted',alpha=0.5,color='red')
ax6.hist(np.log10(data['be_like_fe_neg']['rate_samples'][:,indx]),bins=100,density=True, ls='dashed',alpha=0.5,color='blue')
ax6.legend()
plt.show()
plt.savefig('DR_rate_distr_Be_like_Fe_low_Te.eps')


indx=10
fig7, ax7 = plt.subplots(1,1)
ax7.set_ylabel('Number of counts',fontsize=16)
ax7.set_xlabel('log$_{10}$(Rate coefficient (cm$^3$ s${^-1}$))',fontsize=16)
#ax4.set_xscale("log")
ax7.hist(data['be_like_fe_pos']['rate_samples'][:,indx],bins=100,density=True, ls='dotted',alpha=0.5,color='red')
ax7.hist(data['be_like_fe_neg']['rate_samples'][:,indx],bins=100,density=True, ls='dashed',alpha=0.5,color='blue')
ax7.legend()
plt.show()
plt.savefig('DR_rate_distr_Be_like_Fe_high_Te.eps')

#%%Iso-electronic plots
# Be-like
data_iso_belike=pickle.load(open('be_like_all_0.4_1.6_40grid_nist_0.03_uniform_gaussian.pkl','rb'))


rate_samples_combined_b =np.concatenate((data_iso_belike['be_like_b_pos']['rate_samples'],data_iso_belike['be_like_b_neg']['rate_samples']))
rate_samples_combined_c =np.concatenate((data_iso_belike['be_like_c_pos']['rate_samples'],data_iso_belike['be_like_c_neg']['rate_samples']))
rate_samples_combined_n =np.concatenate((data_iso_belike['be_like_n_pos']['rate_samples'],data_iso_belike['be_like_n_neg']['rate_samples']))
rate_samples_combined_o =np.concatenate((data_iso_belike['be_like_o_pos']['rate_samples'],data_iso_belike['be_like_o_neg']['rate_samples']))
rate_samples_combined_f =np.concatenate((data_iso_belike['be_like_f_pos']['rate_samples'],data_iso_belike['be_like_f_neg']['rate_samples']))
rate_samples_combined_ne=np.concatenate((data_iso_belike['be_like_ne_pos']['rate_samples'],data_iso_belike['be_like_ne_neg']['rate_samples']))
rate_samples_combined_mg=np.concatenate((data_iso_belike['be_like_mg_pos']['rate_samples'],data_iso_belike['be_like_mg_neg']['rate_samples']))
rate_samples_combined_na=np.concatenate((data_iso_belike['be_like_na_pos']['rate_samples'],data_iso_belike['be_like_na_neg']['rate_samples']))
rate_samples_combined_al=np.concatenate((data_iso_belike['be_like_al_pos']['rate_samples'],data_iso_belike['be_like_al_neg']['rate_samples']))
rate_samples_combined_si=np.concatenate((data_iso_belike['be_like_si_pos']['rate_samples'],data_iso_belike['be_like_si_neg']['rate_samples']))
rate_samples_combined_p =np.concatenate((data_iso_belike['be_like_p_pos']['rate_samples'],data_iso_belike['be_like_p_neg']['rate_samples']))
rate_samples_combined_s =np.concatenate((data_iso_belike['be_like_s_pos']['rate_samples'],data_iso_belike['be_like_s_neg']['rate_samples']))
rate_samples_combined_cl=np.concatenate((data_iso_belike['be_like_cl_pos']['rate_samples'],data_iso_belike['be_like_cl_neg']['rate_samples']))
rate_samples_combined_k =np.concatenate((data_iso_belike['be_like_k_pos']['rate_samples'],data_iso_belike['be_like_k_neg']['rate_samples']))
rate_samples_combined_ca=np.concatenate((data_iso_belike['be_like_ca_pos']['rate_samples'],data_iso_belike['be_like_ca_neg']['rate_samples']))
rate_samples_combined_sc=np.concatenate((data_iso_belike['be_like_sc_pos']['rate_samples'],data_iso_belike['be_like_sc_neg']['rate_samples']))
rate_samples_combined_ti=np.concatenate((data_iso_belike['be_like_ti_pos']['rate_samples'],data_iso_belike['be_like_ti_neg']['rate_samples']))
rate_samples_combined_v =np.concatenate((data_iso_belike['be_like_v_pos']['rate_samples'],data_iso_belike['be_like_v_neg']['rate_samples']))
rate_samples_combined_cr=np.concatenate((data_iso_belike['be_like_cr_pos']['rate_samples'],data_iso_belike['be_like_cr_neg']['rate_samples']))
rate_samples_combined_mn=np.concatenate((data_iso_belike['be_like_mn_pos']['rate_samples'],data_iso_belike['be_like_mn_neg']['rate_samples']))
rate_samples_combined_fe=np.concatenate((data_iso_belike['be_like_fe_pos']['rate_samples'],data_iso_belike['be_like_fe_neg']['rate_samples']))


quartile1_b, medians_b, quartile3_b = np.percentile(rate_samples_combined_b[:,:], [25, 50, 75], axis=0)
quartile1_c, medians_c, quartile3_c = np.percentile(rate_samples_combined_c[:,:], [25, 50, 75], axis=0)
quartile1_n, medians_n, quartile3_n = np.percentile(rate_samples_combined_n[:,:], [25, 50, 75], axis=0)
quartile1_o, medians_o, quartile3_o = np.percentile(rate_samples_combined_o[:,:], [25, 50, 75], axis=0)
quartile1_f, medians_f, quartile3_f = np.percentile(rate_samples_combined_f[:,:], [25, 50, 75], axis=0)
quartile1_ne, medians_ne, quartile3_ne = np.percentile(rate_samples_combined_ne[:,:], [25, 50, 75], axis=0)
quartile1_mg, medians_mg, quartile3_mg = np.percentile(rate_samples_combined_mg[:,:], [25, 50, 75], axis=0)
quartile1_na, medians_na, quartile3_na = np.percentile(rate_samples_combined_na[:,:], [25, 50, 75], axis=0)
quartile1_al, medians_al, quartile3_al = np.percentile(rate_samples_combined_al[:,:], [25, 50, 75], axis=0)
quartile1_si, medians_si, quartile3_si = np.percentile(rate_samples_combined_si[:,:], [25, 50, 75], axis=0)
quartile1_p, medians_p, quartile3_p = np.percentile(rate_samples_combined_p[:,:], [25, 50, 75], axis=0)
quartile1_s, medians_s, quartile3_s = np.percentile(rate_samples_combined_s[:,:], [25, 50, 75], axis=0)
quartile1_cl, medians_cl, quartile3_cl = np.percentile(rate_samples_combined_cl[:,:], [25, 50, 75], axis=0)
quartile1_k, medians_k, quartile3_k = np.percentile(rate_samples_combined_k[:,:], [25, 50, 75], axis=0)
quartile1_ca, medians_ca, quartile3_ca = np.percentile(rate_samples_combined_ca[:,:], [25, 50, 75], axis=0)
quartile1_sc, medians_sc, quartile3_sc = np.percentile(rate_samples_combined_sc[:,:], [25, 50, 75], axis=0)
quartile1_ti, medians_ti, quartile3_ti = np.percentile(rate_samples_combined_ti[:,:], [25, 50, 75], axis=0)
quartile1_v, medians_v, quartile3_v = np.percentile(rate_samples_combined_v[:,:], [25, 50, 75], axis=0)
quartile1_cr, medians_cr, quartile3_cr = np.percentile(rate_samples_combined_cr[:,:], [25, 50, 75], axis=0)
quartile1_mn, medians_mn, quartile3_mn = np.percentile(rate_samples_combined_mn[:,:], [25, 50, 75], axis=0)
quartile1_fe, medians_fe, quartile3_fe = np.percentile(rate_samples_combined_fe[:,:], [25, 50, 75], axis=0)





fig10, axs10 = plt.subplots(5,1)
axs10[0].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[0].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[0].set_yscale("log")
axs10[0].set_xscale("log")
axs10[0].set_ylim(0.01,100.)
axs10[0].set_xlim(10.0,1.e10)
axs10[1].set_xlim(10.0,1.e10)
axs10[2].set_xlim(10.0,1.e10)
axs10[3].set_xlim(10.0,1.e10)
axs10[4].set_xlim(10.0,1.e10)
axs10[1].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[1].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[1].set_yscale("log")
axs10[1].set_xscale("log")
axs10[2].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[2].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[2].set_yscale("log")
axs10[2].set_xscale("log")
axs10[3].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[3].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[3].set_yscale("log")
axs10[3].set_xscale("log")
axs10[4].set_ylabel('(Q3-M)/M',fontsize=8)
axs10[4].set_xlabel('Electron temperature (K)',fontsize=12)
axs10[4].set_yscale("log")
axs10[4].set_xscale("log")
axs10[4].set_ylim(0.001,1.)
axs10[0].plot(data_iso_belike['be_like_b_pos']['T'],(quartile3_b-medians_b)/medians_b,color='black',linestyle='solid',label='B Q3-median')
#axs10[0].plot(data_iso_belike['be_like_b_pos']['T'],(medians_b-quartile1_b)/medians_b,color='black',linestyle='dashed',label='B Median-Q1')
axs10[0].plot(data_iso_belike['be_like_c_pos']['T'],(quartile3_c-medians_c)/medians_c,color='red',linestyle='solid',label='C Q3-median')
#axs10[0].plot(data_iso_belike['be_like_c_pos']['T'],(medians_c-quartile1_c)/medians_c,color='red',linestyle='dashed',label='C Median-Q1')
axs10[0].plot(data_iso_belike['be_like_n_pos']['T'],(quartile3_n-medians_n)/medians_n,color='blue',linestyle='solid',label='N Q3-median')
#axs10[0].plot(data_iso_belike['be_like_n_pos']['T'],(medians_n-quartile1_n)/medians_n,color='blue',linestyle='dashed',label='N Median-Q1')
axs10[0].plot(data_iso_belike['be_like_o_pos']['T'],(quartile3_o-medians_o)/medians_o,color='green',linestyle='solid',label='O Q3-median')
#axs10[0].plot(data_iso_belike['be_like_o_pos']['T'],(medians_o-quartile1_o)/medians_o,color='green',linestyle='dashed',label='O Median-Q1')
axs10[1].plot(data_iso_belike['be_like_f_pos']['T'],(quartile3_f-medians_f)/medians_f,color='black',linestyle='solid',label='F Q3-median')
#axs10[1].plot(data_iso_belike['be_like_f_pos']['T'],(medians_f-quartile1_f)/medians_f,color='black',linestyle='dashed',label='F Median-Q1')
axs10[1].plot(data_iso_belike['be_like_ne_pos']['T'],(quartile3_ne-medians_ne)/medians_ne,color='red',linestyle='solid',label='Ne Q3-median')
#axs10[1].plot(data_iso_belike['be_like_ne_pos']['T'],(medians_ne-quartile1_ne)/medians_ne,color='red',linestyle='dashed',label='Ne Median-Q1')
axs10[1].plot(data_iso_belike['be_like_mg_pos']['T'],(quartile3_mg-medians_mg)/medians_mg,color='blue',linestyle='solid',label='Mg Q3-median')
#axs10[1].plot(data_iso_belike['be_like_mg_pos']['T'],(medians_mg-quartile1_mg)/medians_mg,color='blue',linestyle='dashed',label='Mg Median-Q1')
axs10[1].plot(data_iso_belike['be_like_na_pos']['T'],(quartile3_na-medians_na)/medians_na,color='green',linestyle='solid',label='Na Q3-median')
#axs10[1].plot(data_iso_belike['be_like_na_pos']['T'],(medians_na-quartile1_na)/medians_na,color='green',linestyle='dashed',label='Na Median-Q1')
axs10[2].plot(data_iso_belike['be_like_al_pos']['T'],(quartile3_al-medians_al)/medians_al,color='black',linestyle='solid',label='Al Q3-median')
#axs10[2].plot(data_iso_belike['be_like_al_pos']['T'],(medians_al-quartile1_al)/medians_al,color='black',linestyle='dashed',label='Al Median-Q1')
axs10[2].plot(data_iso_belike['be_like_si_pos']['T'],(quartile3_si-medians_si)/medians_si,color='red',linestyle='solid',label='Si Q3-median')
#axs10[2].plot(data_iso_belike['be_like_si_pos']['T'],(medians_si-quartile1_si)/medians_si,color='red',linestyle='dashed',label='Si Median-Q1')
axs10[2].plot(data_iso_belike['be_like_p_pos']['T'],(quartile3_p-medians_p)/medians_p,color='blue',linestyle='solid',label='P Q3-median')
#axs10[2].plot(data_iso_belike['be_like_p_pos']['T'],(medians_p-quartile1_p)/medians_p,color='blue',linestyle='dashed',label='P Median-Q1')
axs10[2].plot(data_iso_belike['be_like_s_pos']['T'],(quartile3_s-medians_s)/medians_s,color='green',linestyle='solid',label='S Q3-median')
#axs10[2].plot(data_iso_belike['be_like_s_pos']['T'],(medians_s-quartile1_s)/medians_s,color='green',linestyle='dashed',label='S Median-Q1')
axs10[3].plot(data_iso_belike['be_like_cl_pos']['T'],(quartile3_cl-medians_cl)/medians_cl,color='black',linestyle='solid',label='Cl Q3-median')
#axs10[3].plot(data_iso_belike['be_like_cl_pos']['T'],(medians_cl-quartile1_cl)/medians_cl,color='black',linestyle='dashed',label='Cl Median-Q1')
axs10[3].plot(data_iso_belike['be_like_k_pos']['T'],(quartile3_k-medians_k)/medians_k,color='red',linestyle='solid',label='K Q3-median')
#axs10[3].plot(data_iso_belike['be_like_k_pos']['T'],(medians_k-quartile1_k)/medians_k,color='red',linestyle='dashed',label='K Median-Q1')
axs10[3].plot(data_iso_belike['be_like_ca_pos']['T'],(quartile3_ca-medians_ca)/medians_ca,color='blue',linestyle='solid',label='Ca Q3-median')
#axs10[3].plot(data_iso_belike['be_like_ca_pos']['T'],(medians_ca-quartile1_ca)/medians_ca,color='blue',linestyle='dashed',label='Ca Median-Q1')
axs10[3].plot(data_iso_belike['be_like_sc_pos']['T'],(quartile3_sc-medians_sc)/medians_sc,color='green',linestyle='solid',label='Sc Q3-median')
#axs10[3].plot(data_iso_belike['be_like_sc_pos']['T'],(medians_sc-quartile1_sc)/medians_sc,color='green',linestyle='dashed',label='Sc Median-Q1')
axs10[4].plot(data_iso_belike['be_like_ti_pos']['T'],(quartile3_ti-medians_ti)/medians_ti,color='black',linestyle='solid',label='Ti Q3-median')
#axs10[4].plot(data_iso_belike['be_like_ti_pos']['T'],(medians_ti-quartile1_ti)/medians_ti,color='black',linestyle='dashed',label='Ti Median-Q1')
axs10[4].plot(data_iso_belike['be_like_v_pos']['T'],(quartile3_v-medians_v)/medians_v,color='red',linestyle='solid',label='Ca Q3-median')
#axs10[4].plot(data_iso_belike['be_like_v_pos']['T'],(medians_v-quartile1_v)/medians_v,color='red',linestyle='dashed',label='Ca Median-Q1')
axs10[4].plot(data_iso_belike['be_like_cr_pos']['T'],(quartile3_cr-medians_cr)/medians_cr,color='blue',linestyle='solid',label='Cr Q3-median')
#axs10[4].plot(data_iso_belike['be_like_cr_pos']['T'],(medians_cr-quartile1_cr)/medians_cr,color='blue',linestyle='dashed',label='Cr Median-Q1')
axs10[4].plot(data_iso_belike['be_like_mn_pos']['T'],(quartile3_mn-medians_mn)/medians_mn,color='green',linestyle='solid',label='Mn Q3-median')
#axs10[4].plot(data_iso_belike['be_like_mn_pos']['T'],(medians_mn-quartile1_mn)/medians_mn,color='green',linestyle='dashed',label='Mn Median-Q1')
axs10[4].plot(data_iso_belike['be_like_fe_pos']['T'],(quartile3_fe-medians_fe)/medians_fe,color='purple',linestyle='solid',label='Fe Q3-median')
#axs10[4].plot(data_iso_belike['be_like_fe_pos']['T'],(medians_fe-quartile1_fe)/medians_fe,color='purple',linestyle='dashed',label='Fe Median-Q1')
axs10[0].legend(prop={"size":6})
axs10[1].legend(prop={"size":6})
axs10[2].legend(prop={"size":6})
axs10[3].legend(prop={"size":6})
axs10[4].legend(prop={"size":6})
plt.show()
plt.savefig('isoelectronic_Be_like_Q3.eps')


fig10, axs10 = plt.subplots(5,1)
axs10[0].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[0].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[0].set_yscale("log")
axs10[0].set_xscale("log")
axs10[0].set_xlim(10.0,1.e10)
axs10[1].set_xlim(10.0,1.e10)
axs10[2].set_xlim(10.0,1.e10)
axs10[3].set_xlim(10.0,1.e10)
axs10[4].set_xlim(10.0,1.e10)
axs10[0].set_ylim(0.01,1.1)
axs10[1].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[1].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[1].set_yscale("log")
axs10[1].set_xscale("log")
axs10[2].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[2].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[2].set_yscale("log")
axs10[2].set_xscale("log")
axs10[3].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[3].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[3].set_yscale("log")
axs10[3].set_xscale("log")
axs10[4].set_ylabel('(M-Q1)/M',fontsize=8)
axs10[4].set_xlabel('Electron temperature (K)',fontsize=12)
axs10[4].set_yscale("log")
axs10[4].set_xscale("log")
axs10[4].set_ylim(0.001,1.)
#axs10[0].plot(data_iso_belike['be_like_b_pos']['T'],(quartile3_b-medians_b)/medians_b,color='black',linestyle='solid',label='B Q3-median')
axs10[0].plot(data_iso_belike['be_like_b_pos']['T'],(medians_b-quartile1_b)/medians_b,color='black',linestyle='dashed',label='B Median-Q1')
#axs10[0].plot(data_iso_belike['be_like_c_pos']['T'],(quartile3_c-medians_c)/medians_c,color='red',linestyle='solid',label='C Q3-median')
axs10[0].plot(data_iso_belike['be_like_c_pos']['T'],(medians_c-quartile1_c)/medians_c,color='red',linestyle='dashed',label='C Median-Q1')
#axs10[0].plot(data_iso_belike['be_like_n_pos']['T'],(quartile3_n-medians_n)/medians_n,color='blue',linestyle='solid',label='N Q3-median')
axs10[0].plot(data_iso_belike['be_like_n_pos']['T'],(medians_n-quartile1_n)/medians_n,color='blue',linestyle='dashed',label='N Median-Q1')
#axs10[0].plot(data_iso_belike['be_like_o_pos']['T'],(quartile3_o-medians_o)/medians_o,color='green',linestyle='solid',label='O Q3-median')
axs10[0].plot(data_iso_belike['be_like_o_pos']['T'],(medians_o-quartile1_o)/medians_o,color='green',linestyle='dashed',label='O Median-Q1')
#axs10[1].plot(data_iso_belike['be_like_f_pos']['T'],(quartile3_f-medians_f)/medians_f,color='black',linestyle='solid',label='F Q3-median')
axs10[1].plot(data_iso_belike['be_like_f_pos']['T'],(medians_f-quartile1_f)/medians_f,color='black',linestyle='dashed',label='F Median-Q1')
#axs10[1].plot(data_iso_belike['be_like_ne_pos']['T'],(quartile3_ne-medians_ne)/medians_ne,color='red',linestyle='solid',label='Ne Q3-median')
axs10[1].plot(data_iso_belike['be_like_ne_pos']['T'],(medians_ne-quartile1_ne)/medians_ne,color='red',linestyle='dashed',label='Ne Median-Q1')
#axs10[1].plot(data_iso_belike['be_like_mg_pos']['T'],(quartile3_mg-medians_mg)/medians_mg,color='blue',linestyle='solid',label='Mg Q3-median')
axs10[1].plot(data_iso_belike['be_like_mg_pos']['T'],(medians_mg-quartile1_mg)/medians_mg,color='blue',linestyle='dashed',label='Mg Median-Q1')
#axs10[1].plot(data_iso_belike['be_like_na_pos']['T'],(quartile3_na-medians_na)/medians_na,color='green',linestyle='solid',label='Na Q3-median')
axs10[1].plot(data_iso_belike['be_like_na_pos']['T'],(medians_na-quartile1_na)/medians_na,color='green',linestyle='dashed',label='Na Median-Q1')
#axs10[2].plot(data_iso_belike['be_like_al_pos']['T'],(quartile3_al-medians_al)/medians_al,color='black',linestyle='solid',label='Al Q3-median')
axs10[2].plot(data_iso_belike['be_like_al_pos']['T'],(medians_al-quartile1_al)/medians_al,color='black',linestyle='dashed',label='Al Median-Q1')
#axs10[2].plot(data_iso_belike['be_like_si_pos']['T'],(quartile3_si-medians_si)/medians_si,color='red',linestyle='solid',label='Si Q3-median')
axs10[2].plot(data_iso_belike['be_like_si_pos']['T'],(medians_si-quartile1_si)/medians_si,color='red',linestyle='dashed',label='Si Median-Q1')
#axs10[2].plot(data_iso_belike['be_like_p_pos']['T'],(quartile3_p-medians_p)/medians_p,color='blue',linestyle='solid',label='P Q3-median')
axs10[2].plot(data_iso_belike['be_like_p_pos']['T'],(medians_p-quartile1_p)/medians_p,color='blue',linestyle='dashed',label='P Median-Q1')
#axs10[2].plot(data_iso_belike['be_like_s_pos']['T'],(quartile3_s-medians_s)/medians_s,color='green',linestyle='solid',label='S Q3-median')
axs10[2].plot(data_iso_belike['be_like_s_pos']['T'],(medians_s-quartile1_s)/medians_s,color='green',linestyle='dashed',label='S Median-Q1')
#axs10[3].plot(data_iso_belike['be_like_cl_pos']['T'],(quartile3_cl-medians_cl)/medians_cl,color='black',linestyle='solid',label='Cl Q3-median')
axs10[3].plot(data_iso_belike['be_like_cl_pos']['T'],(medians_cl-quartile1_cl)/medians_cl,color='black',linestyle='dashed',label='Cl Median-Q1')
#axs10[3].plot(data_iso_belike['be_like_k_pos']['T'],(quartile3_k-medians_k)/medians_k,color='red',linestyle='solid',label='K Q3-median')
axs10[3].plot(data_iso_belike['be_like_k_pos']['T'],(medians_k-quartile1_k)/medians_k,color='red',linestyle='dashed',label='K Median-Q1')
#axs10[3].plot(data_iso_belike['be_like_ca_pos']['T'],(quartile3_ca-medians_ca)/medians_ca,color='blue',linestyle='solid',label='Ca Q3-median')
axs10[3].plot(data_iso_belike['be_like_ca_pos']['T'],(medians_ca-quartile1_ca)/medians_ca,color='blue',linestyle='dashed',label='Ca Median-Q1')
#axs10[3].plot(data_iso_belike['be_like_sc_pos']['T'],(quartile3_sc-medians_sc)/medians_sc,color='green',linestyle='solid',label='Sc Q3-median')
axs10[3].plot(data_iso_belike['be_like_sc_pos']['T'],(medians_sc-quartile1_sc)/medians_sc,color='green',linestyle='dashed',label='Sc Median-Q1')
#axs10[4].plot(data_iso_belike['be_like_ti_pos']['T'],(quartile3_ti-medians_ti)/medians_ti,color='black',linestyle='solid',label='Ti Q3-median')
axs10[4].plot(data_iso_belike['be_like_ti_pos']['T'],(medians_ti-quartile1_ti)/medians_ti,color='black',linestyle='dashed',label='Ti Median-Q1')
#axs10[4].plot(data_iso_belike['be_like_v_pos']['T'],(quartile3_v-medians_v)/medians_v,color='red',linestyle='solid',label='Ca Q3-median')
axs10[4].plot(data_iso_belike['be_like_v_pos']['T'],(medians_v-quartile1_v)/medians_v,color='red',linestyle='dashed',label='Ca Median-Q1')
#axs10[4].plot(data_iso_belike['be_like_cr_pos']['T'],(quartile3_cr-medians_cr)/medians_cr,color='blue',linestyle='solid',label='Cr Q3-median')
axs10[4].plot(data_iso_belike['be_like_cr_pos']['T'],(medians_cr-quartile1_cr)/medians_cr,color='blue',linestyle='dashed',label='Cr Median-Q1')
#axs10[4].plot(data_iso_belike['be_like_mn_pos']['T'],(quartile3_mn-medians_mn)/medians_mn,color='green',linestyle='solid',label='Mn Q3-median')
axs10[4].plot(data_iso_belike['be_like_mn_pos']['T'],(medians_mn-quartile1_mn)/medians_mn,color='green',linestyle='dashed',label='Mn Median-Q1')
#axs10[4].plot(data_iso_belike['be_like_fe_pos']['T'],(quartile3_fe-medians_fe)/medians_fe,color='purple',linestyle='solid',label='Fe Q3-median')
axs10[4].plot(data_iso_belike['be_like_fe_pos']['T'],(medians_fe-quartile1_fe)/medians_fe,color='purple',linestyle='dashed',label='Fe Median-Q1')
axs10[0].legend(prop={"size":6})
axs10[1].legend(prop={"size":6})
axs10[2].legend(prop={"size":6})
axs10[3].legend(prop={"size":6})
axs10[4].legend(prop={"size":6})
plt.show()
plt.savefig('isoelectronic_Be_like_Q1.eps')




##Plots to try to figure out stange Q3-Median=0 result for Be-like Cr at the 6th Te index
###Boxplot
#fig3, axs3 = plt.subplots(1,1)
#title_pos='Rate coefficients'
#axs3.set_yscale("log")
#axs3.set_ylabel('Rate Coefficient (cm$^3$ $s^{-1}$)', fontsize=16)
#axs3.set_xlabel('Log$_{10}$ (Electron Temperature (K))', fontsize=16)
#axs3.set_title(title_pos)
##axs3.boxplot(rate_samples_combined_cr[:,:],showfliers=False)
#axs3.violinplot(rate_samples_combined_cr, showmeans=False, showmedians=True,showextrema=True,widths=0.3)
##axs3.boxplot(data['be_like_fe_pos']['rate_samples'][:,3:],showfliers=False)
##axs3.boxplot(data['be_like_fe_neg']['rate_samples'][:,3:],showfliers=False)
#axs3.legend()
#plt.show()
#
#
#
#indx=6
#fig7, ax7 = plt.subplots(1,1)
#ax7.set_ylabel('Number of counts',fontsize=16)
#ax7.set_xlabel('log$_{10}$(Rate coefficient (cm$^3$ s${^-1}$))',fontsize=16)
##ax7.set_yscale("log")
#ax7.hist(data_iso_belike['be_like_cr_pos']['rate_samples'][:,indx],bins=20,density=True, ls='dotted',alpha=0.5,color='red')
#ax7.hist(data_iso_belike['be_like_cr_neg']['rate_samples'][:,indx],bins=20,density=True, ls='dashed',alpha=0.5,color='blue')
#ax7.legend()
#plt.show()
#
#
#Lambda_1=data_iso_belike['be_like_cr_pos']['X_1D'][0][:]
#Lambda_2=data_iso_belike['be_like_cr_pos']['X_1D'][1][:]
#Lambda_1, Lambda_2 = np.meshgrid(Lambda_1,Lambda_2)
#
#indx=6
#fig5 = plt.figure()
#ax5 = fig5.gca(projection='3d')
#surf = ax5.plot_wireframe(Lambda_1,Lambda_2, np.log10(data_iso_belike['be_like_cr_pos']['rates'][indx][:][:]),color='black',linewidth=0.5)
##surf = ax5.plot_wireframe(Lambda_1,Lambda_2, np.log10(data_iso_belike['be_like_cr_pos']['rates'][indx][:][:]),color='black',linewidth=0.5)
#ax5.scatter3D(data_iso_belike['be_like_cr_pos']['lambda_samples'][:,1], data_iso_belike['be_like_cr_pos']['lambda_samples'][:,0], np.log10(data_iso_belike['be_like_cr_pos']['rate_samples'][:,indx]))
#ax5.set_xlabel('$\lambda_{2s}$')
#ax5.set_ylabel('$\lambda_{2p$}$')
#ax5.set_zlabel('log$_{10}$(Rate Coefficient cm$^3$ s$^{-1}$)')
#plt.show()
#
#
#lambda_labels=["$\lambda_{2s}$","$\lambda_{2p}$"]
#n_lambdas=2
##fig0=corner.corner(data['be_like_o_pos']['lambda_samples'], labels=[f"$\lambda_{i+1}$" for i in range(n_lambdas)], truths=[1 for i in range(n_lambdas)])
#fig0=corner.corner(data_iso_belike['be_like_cr_pos']['lambda_samples'], labels=[lambda_labels[i] for i in range(n_lambdas)], label_kwargs={"fontsize": 14}, truths=[1 for i in range(n_lambdas)])
#
#
#fig0=corner.corner(data_iso_belike['be_like_cr_neg']['lambda_samples'], labels=[lambda_labels[i] for i in range(n_lambdas)], label_kwargs={"fontsize": 14}, truths=[1 for i in range(n_lambdas)])
#
#
#np.quantile(rate_samples_combined_cr,0.75,axis=0)
#np.quantile(rate_samples_combined_cr,0.5,axis=0)
#
#
#np.quantile(data_iso_belike['be_like_cr_neg']['rate_samples'],0.5,axis=0)
#np.quantile(data_iso_belike['be_like_cr_neg']['rate_samples'],0.75,axis=0)
#np.quantile(data_iso_belike['be_like_cr_pos']['rate_samples'],0.5,axis=0)
#np.quantile(data_iso_belike['be_like_cr_pos']['rate_samples'],0.75,axis=0)
#

#%% Iso-electronic plots
# Li-like

data_iso_lilike=pickle.load(open('li_like_all_0.4_1.6_40grid_nist_0.02_uniform_gaussian.pkl','rb'))


rate_samples_combined_be =np.concatenate((data_iso_lilike['li_like_be_pos']['rate_samples'],data_iso_lilike['li_like_be_neg']['rate_samples']))
rate_samples_combined_b =np.concatenate((data_iso_lilike['li_like_b_pos']['rate_samples'],data_iso_lilike['li_like_b_neg']['rate_samples']))
rate_samples_combined_c =np.concatenate((data_iso_lilike['li_like_c_pos']['rate_samples'],data_iso_lilike['li_like_c_neg']['rate_samples']))
rate_samples_combined_n =np.concatenate((data_iso_lilike['li_like_n_pos']['rate_samples'],data_iso_lilike['li_like_n_neg']['rate_samples']))
rate_samples_combined_o =np.concatenate((data_iso_lilike['li_like_o_pos']['rate_samples'],data_iso_lilike['li_like_o_neg']['rate_samples']))
rate_samples_combined_f =np.concatenate((data_iso_lilike['li_like_f_pos']['rate_samples'],data_iso_lilike['li_like_f_neg']['rate_samples']))
rate_samples_combined_ne=np.concatenate((data_iso_lilike['li_like_ne_pos']['rate_samples'],data_iso_lilike['li_like_ne_neg']['rate_samples']))
rate_samples_combined_mg=np.concatenate((data_iso_lilike['li_like_mg_pos']['rate_samples'],data_iso_lilike['li_like_mg_neg']['rate_samples']))
rate_samples_combined_na=np.concatenate((data_iso_lilike['li_like_na_pos']['rate_samples'],data_iso_lilike['li_like_na_neg']['rate_samples']))
rate_samples_combined_al=np.concatenate((data_iso_lilike['li_like_al_pos']['rate_samples'],data_iso_lilike['li_like_al_neg']['rate_samples']))
rate_samples_combined_si=np.concatenate((data_iso_lilike['li_like_si_pos']['rate_samples'],data_iso_lilike['li_like_si_neg']['rate_samples']))
rate_samples_combined_p =np.concatenate((data_iso_lilike['li_like_p_pos']['rate_samples'],data_iso_lilike['li_like_p_neg']['rate_samples']))
rate_samples_combined_s =np.concatenate((data_iso_lilike['li_like_s_pos']['rate_samples'],data_iso_lilike['li_like_s_neg']['rate_samples']))
rate_samples_combined_cl=np.concatenate((data_iso_lilike['li_like_cl_pos']['rate_samples'],data_iso_lilike['li_like_cl_neg']['rate_samples']))
rate_samples_combined_k =np.concatenate((data_iso_lilike['li_like_k_pos']['rate_samples'],data_iso_lilike['li_like_k_neg']['rate_samples']))
rate_samples_combined_ca=np.concatenate((data_iso_lilike['li_like_ca_pos']['rate_samples'],data_iso_lilike['li_like_ca_neg']['rate_samples']))
rate_samples_combined_sc=np.concatenate((data_iso_lilike['li_like_sc_pos']['rate_samples'],data_iso_lilike['li_like_sc_neg']['rate_samples']))
rate_samples_combined_ti=np.concatenate((data_iso_lilike['li_like_ti_pos']['rate_samples'],data_iso_lilike['li_like_ti_neg']['rate_samples']))
rate_samples_combined_v =np.concatenate((data_iso_lilike['li_like_v_pos']['rate_samples'],data_iso_lilike['li_like_v_neg']['rate_samples']))
rate_samples_combined_cr=np.concatenate((data_iso_lilike['li_like_cr_pos']['rate_samples'],data_iso_lilike['li_like_cr_neg']['rate_samples']))
rate_samples_combined_mn=np.concatenate((data_iso_lilike['li_like_mn_pos']['rate_samples'],data_iso_lilike['li_like_mn_neg']['rate_samples']))
rate_samples_combined_fe=np.concatenate((data_iso_lilike['li_like_fe_pos']['rate_samples'],data_iso_lilike['li_like_fe_neg']['rate_samples']))



quartile1_be, medians_be, quartile3_be = np.percentile(rate_samples_combined_be[:,:], [25, 50, 75], axis=0)
quartile1_b, medians_b, quartile3_b = np.percentile(rate_samples_combined_b[:,:], [25, 50, 75], axis=0)
quartile1_c, medians_c, quartile3_c = np.percentile(rate_samples_combined_c[:,:], [25, 50, 75], axis=0)
quartile1_n, medians_n, quartile3_n = np.percentile(rate_samples_combined_n[:,:], [25, 50, 75], axis=0)
quartile1_o, medians_o, quartile3_o = np.percentile(rate_samples_combined_o[:,:], [25, 50, 75], axis=0)
quartile1_f, medians_f, quartile3_f = np.percentile(rate_samples_combined_f[:,:], [25, 50, 75], axis=0)
quartile1_ne, medians_ne, quartile3_ne = np.percentile(rate_samples_combined_ne[:,:], [25, 50, 75], axis=0)
quartile1_mg, medians_mg, quartile3_mg = np.percentile(rate_samples_combined_mg[:,:], [25, 50, 75], axis=0)
quartile1_na, medians_na, quartile3_na = np.percentile(rate_samples_combined_na[:,:], [25, 50, 75], axis=0)
quartile1_al, medians_al, quartile3_al = np.percentile(rate_samples_combined_al[:,:], [25, 50, 75], axis=0)
quartile1_si, medians_si, quartile3_si = np.percentile(rate_samples_combined_si[:,:], [25, 50, 75], axis=0)
quartile1_p, medians_p, quartile3_p = np.percentile(rate_samples_combined_p[:,:], [25, 50, 75], axis=0)
quartile1_s, medians_s, quartile3_s = np.percentile(rate_samples_combined_s[:,:], [25, 50, 75], axis=0)
quartile1_cl, medians_cl, quartile3_cl = np.percentile(rate_samples_combined_cl[:,:], [25, 50, 75], axis=0)
quartile1_k, medians_k, quartile3_k = np.percentile(rate_samples_combined_k[:,:], [25, 50, 75], axis=0)
quartile1_ca, medians_ca, quartile3_ca = np.percentile(rate_samples_combined_ca[:,:], [25, 50, 75], axis=0)
quartile1_sc, medians_sc, quartile3_sc = np.percentile(rate_samples_combined_sc[:,:], [25, 50, 75], axis=0)
quartile1_ti, medians_ti, quartile3_ti = np.percentile(rate_samples_combined_ti[:,:], [25, 50, 75], axis=0)
quartile1_v, medians_v, quartile3_v = np.percentile(rate_samples_combined_v[:,:], [25, 50, 75], axis=0)
quartile1_cr, medians_cr, quartile3_cr = np.percentile(rate_samples_combined_cr[:,:], [25, 50, 75], axis=0)
quartile1_mn, medians_mn, quartile3_mn = np.percentile(rate_samples_combined_mn[:,:], [25, 50, 75], axis=0)
quartile1_fe, medians_fe, quartile3_fe = np.percentile(rate_samples_combined_fe[:,:], [25, 50, 75], axis=0)





fig10, axs10 = plt.subplots(5,1)
axs10[0].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[0].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[0].set_yscale("log")
axs10[0].set_xscale("log")
axs10[0].set_ylim(0.01,1.1)
axs10[1].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[1].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[1].set_yscale("log")
axs10[1].set_xscale("log")
axs10[2].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[2].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[2].set_yscale("log")
axs10[2].set_xscale("log")
axs10[2].set_ylim(0.001,1.)
axs10[3].set_ylabel('(M-Q1)/M',fontsize=8)
#axs10[3].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[3].set_yscale("log")
axs10[3].set_xscale("log")
axs10[4].set_ylabel('(M-Q1)/M',fontsize=8)
axs10[4].set_xlabel('Electron temperature (K)',fontsize=12)
axs10[4].set_yscale("log")
axs10[4].set_xscale("log")
axs10[4].set_ylim(0.001,1.)
#axs10[0].plot(data_iso_lilike['li_like_be_pos']['T'],(quartile3_be-medians_be)/medians_be,color='purple',linestyle='solid',label='Be Q3-median')
axs10[0].plot(data_iso_lilike['li_like_be_pos']['T'],(medians_be-quartile1_be)/medians_be,color='purple',linestyle='dashed',label='Be Median-Q1')
#axs10[0].plot(data_iso_lilike['li_like_b_pos']['T'],(quartile3_b-medians_b)/medians_b,color='black',linestyle='solid',label='B Q3-median')
axs10[0].plot(data_iso_lilike['li_like_b_pos']['T'],(medians_b-quartile1_b)/medians_b,color='black',linestyle='dashed',label='B Median-Q1')
#axs10[0].plot(data_iso_lilike['li_like_c_pos']['T'],(quartile3_c-medians_c)/medians_c,color='red',linestyle='solid',label='C Q3-median')
axs10[0].plot(data_iso_lilike['li_like_c_pos']['T'],(medians_c-quartile1_c)/medians_c,color='red',linestyle='dashed',label='C Median-Q1')
#axs10[0].plot(data_iso_lilike['li_like_n_pos']['T'],(quartile3_n-medians_n)/medians_n,color='blue',linestyle='solid',label='N Q3-median')
axs10[0].plot(data_iso_lilike['li_like_n_pos']['T'],(medians_n-quartile1_n)/medians_n,color='blue',linestyle='dashed',label='N Median-Q1')
#axs10[0].plot(data_iso_lilike['li_like_o_pos']['T'],(quartile3_o-medians_o)/medians_o,color='green',linestyle='solid',label='O Q3-median')
axs10[0].plot(data_iso_lilike['li_like_o_pos']['T'],(medians_o-quartile1_o)/medians_o,color='green',linestyle='dashed',label='O Median-Q1')
#axs10[1].plot(data_iso_lilike['li_like_f_pos']['T'],(quartile3_f-medians_f)/medians_f,color='black',linestyle='solid',label='F Q3-median')
axs10[1].plot(data_iso_lilike['li_like_f_pos']['T'],(medians_f-quartile1_f)/medians_f,color='black',linestyle='dashed',label='F Median-Q1')
#axs10[1].plot(data_iso_lilike['li_like_ne_pos']['T'],(quartile3_ne-medians_ne)/medians_ne,color='red',linestyle='solid',label='Ne Q3-median')
axs10[1].plot(data_iso_lilike['li_like_ne_pos']['T'],(medians_ne-quartile1_ne)/medians_ne,color='red',linestyle='dashed',label='Ne Median-Q1')
#axs10[1].plot(data_iso_lilike['li_like_mg_pos']['T'],(quartile3_mg-medians_mg)/medians_mg,color='blue',linestyle='solid',label='Mg Q3-median')
axs10[1].plot(data_iso_lilike['li_like_mg_pos']['T'],(medians_mg-quartile1_mg)/medians_mg,color='blue',linestyle='dashed',label='Mg Median-Q1')
#axs10[1].plot(data_iso_lilike['li_like_na_pos']['T'],(quartile3_na-medians_na)/medians_na,color='green',linestyle='solid',label='Na Q3-median')
axs10[1].plot(data_iso_lilike['li_like_na_pos']['T'],(medians_na-quartile1_na)/medians_na,color='green',linestyle='dashed',label='Na Median-Q1')
#axs10[2].plot(data_iso_lilike['li_like_al_pos']['T'],(quartile3_al-medians_al)/medians_al,color='black',linestyle='solid',label='Al Q3-median')
axs10[2].plot(data_iso_lilike['li_like_al_pos']['T'],(medians_al-quartile1_al)/medians_al,color='black',linestyle='dashed',label='Al Median-Q1')
#axs10[2].plot(data_iso_lilike['li_like_si_pos']['T'],(quartile3_si-medians_si)/medians_si,color='red',linestyle='solid',label='Si Q3-median')
axs10[2].plot(data_iso_lilike['li_like_si_pos']['T'],(medians_si-quartile1_si)/medians_si,color='red',linestyle='dashed',label='Si Median-Q1')
#axs10[2].plot(data_iso_lilike['li_like_p_pos']['T'],(quartile3_p-medians_p)/medians_p,color='blue',linestyle='solid',label='P Q3-median')
axs10[2].plot(data_iso_lilike['li_like_p_pos']['T'],(medians_p-quartile1_p)/medians_p,color='blue',linestyle='dashed',label='P Median-Q1')
#axs10[2].plot(data_iso_lilike['li_like_s_pos']['T'],(quartile3_s-medians_s)/medians_s,color='green',linestyle='solid',label='S Q3-median')
axs10[2].plot(data_iso_lilike['li_like_s_pos']['T'],(medians_s-quartile1_s)/medians_s,color='green',linestyle='dashed',label='S Median-Q1')
#axs10[3].plot(data_iso_lilike['li_like_cl_pos']['T'],(quartile3_cl-medians_cl)/medians_cl,color='black',linestyle='solid',label='Cl Q3-median')
axs10[3].plot(data_iso_lilike['li_like_cl_pos']['T'],(medians_cl-quartile1_cl)/medians_cl,color='black',linestyle='dashed',label='Cl Median-Q1')
#axs10[3].plot(data_iso_lilike['li_like_k_pos']['T'],(quartile3_k-medians_k)/medians_k,color='red',linestyle='solid',label='K Q3-median')
axs10[3].plot(data_iso_lilike['li_like_k_pos']['T'],(medians_k-quartile1_k)/medians_k,color='red',linestyle='dashed',label='K Median-Q1')
#axs10[3].plot(data_iso_lilike['li_like_ca_pos']['T'],(quartile3_ca-medians_ca)/medians_ca,color='blue',linestyle='solid',label='Ca Q3-median')
axs10[3].plot(data_iso_lilike['li_like_ca_pos']['T'],(medians_ca-quartile1_ca)/medians_ca,color='blue',linestyle='dashed',label='Ca Median-Q1')
#axs10[3].plot(data_iso_lilike['li_like_sc_pos']['T'],(quartile3_sc-medians_sc)/medians_ca,color='green',linestyle='solid',label='Sc Q3-median')
axs10[3].plot(data_iso_lilike['li_like_sc_pos']['T'],(medians_sc-quartile1_sc)/medians_sc,color='green',linestyle='dashed',label='Sc Median-Q1')
#axs10[4].plot(data_iso_lilike['li_like_ti_pos']['T'],(quartile3_ti-medians_ti)/medians_ti,color='black',linestyle='solid',label='Ti Q3-median')
axs10[4].plot(data_iso_lilike['li_like_ti_pos']['T'],(medians_ti-quartile1_ti)/medians_ti,color='black',linestyle='dashed',label='Ti Median-Q1')
#axs10[4].plot(data_iso_lilike['li_like_v_pos']['T'],(quartile3_ca-medians_ca)/medians_ca,color='red',linestyle='solid',label='Ca Q3-median')
axs10[4].plot(data_iso_lilike['li_like_v_pos']['T'],(medians_v-quartile1_v)/medians_v,color='red',linestyle='dashed',label='Ca Median-Q1')
#axs10[4].plot(data_iso_lilike['li_like_cr_pos']['T'],(quartile3_cr-medians_cr)/medians_cr,color='blue',linestyle='solid',label='Cr Q3-median')
axs10[4].plot(data_iso_lilike['li_like_cr_pos']['T'],(medians_cr-quartile1_cr)/medians_cr,color='blue',linestyle='dashed',label='Cr Median-Q1')
#axs10[4].plot(data_iso_lilike['li_like_mn_pos']['T'],(quartile3_mn-medians_mn)/medians_mn,color='green',linestyle='solid',label='Mn Q3-median')
axs10[4].plot(data_iso_lilike['li_like_mn_pos']['T'],(medians_mn-quartile1_mn)/medians_mn,color='green',linestyle='dashed',label='Mn Median-Q1')
#axs10[4].plot(data_iso_lilike['li_like_fe_pos']['T'],(quartile3_fe-medians_fe)/medians_fe,color='purple',linestyle='solid',label='Fe Q3-median')
axs10[4].plot(data_iso_lilike['li_like_fe_pos']['T'],(medians_fe-quartile1_fe)/medians_fe,color='purple',linestyle='dashed',label='Fe Median-Q1')
axs10[0].legend(prop={"size":6})
axs10[1].legend(prop={"size":6})
axs10[2].legend(prop={"size":6})
axs10[3].legend(prop={"size":6})
axs10[4].legend(prop={"size":6})
plt.show()
plt.savefig('isoelectronic_Li_like_Q1.eps')


fig10, axs10 = plt.subplots(5,1)
axs10[0].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[0].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[0].set_yscale("log")
axs10[0].set_xscale("log")
axs10[0].set_ylim(0.01,1.1)
axs10[1].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[1].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[1].set_yscale("log")
axs10[1].set_xscale("log")
axs10[2].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[2].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[2].set_yscale("log")
axs10[2].set_xscale("log")
axs10[2].set_ylim(0.001,1.)
axs10[3].set_ylabel('(Q3-M)/M',fontsize=8)
#axs10[3].set_xlabel('Electron temperature (K)',fontsize=16)
axs10[3].set_yscale("log")
axs10[3].set_xscale("log")
axs10[4].set_ylabel('(Q3-M)/M',fontsize=8)
axs10[4].set_xlabel('Electron temperature (K)',fontsize=12)
axs10[4].set_yscale("log")
axs10[4].set_xscale("log")
axs10[4].set_ylim(0.001,1.)
axs10[0].plot(data_iso_lilike['li_like_be_pos']['T'],(quartile3_be-medians_be)/medians_be,color='purple',linestyle='solid',label='Be Q3-median')
#axs10[0].plot(data_iso_lilike['li_like_be_pos']['T'],(medians_be-quartile1_be)/medians_be,color='black',linestyle='dashed',label='B Median-Q1')
axs10[0].plot(data_iso_lilike['li_like_b_pos']['T'],(quartile3_b-medians_b)/medians_b,color='black',linestyle='solid',label='B Q3-median')
#axs10[0].plot(data_iso_lilike['li_like_b_pos']['T'],(medians_b-quartile1_b)/medians_b,color='black',linestyle='dashed',label='B Median-Q1')
axs10[0].plot(data_iso_lilike['li_like_c_pos']['T'],(quartile3_c-medians_c)/medians_c,color='red',linestyle='solid',label='C Q3-median')
#axs10[0].plot(data_iso_lilike['li_like_c_pos']['T'],(medians_c-quartile1_c)/medians_c,color='red',linestyle='dashed',label='C Median-Q1')
axs10[0].plot(data_iso_lilike['li_like_n_pos']['T'],(quartile3_n-medians_n)/medians_n,color='blue',linestyle='solid',label='N Q3-median')
#axs10[0].plot(data_iso_lilike['li_like_n_pos']['T'],(medians_n-quartile1_n)/medians_n,color='blue',linestyle='dashed',label='N Median-Q1')
axs10[0].plot(data_iso_lilike['li_like_o_pos']['T'],(quartile3_o-medians_o)/medians_o,color='green',linestyle='solid',label='O Q3-median')
#axs10[0].plot(data_iso_lilike['li_like_o_pos']['T'],(medians_o-quartile1_o)/medians_o,color='green',linestyle='dashed',label='O Median-Q1')
axs10[1].plot(data_iso_lilike['li_like_f_pos']['T'],(quartile3_f-medians_f)/medians_f,color='black',linestyle='solid',label='F Q3-median')
#axs10[1].plot(data_iso_lilike['li_like_f_pos']['T'],(medians_f-quartile1_f)/medians_f,color='black',linestyle='dashed',label='F Median-Q1')
axs10[1].plot(data_iso_lilike['li_like_ne_pos']['T'],(quartile3_ne-medians_ne)/medians_ne,color='red',linestyle='solid',label='Ne Q3-median')
#axs10[1].plot(data_iso_lilike['li_like_ne_pos']['T'],(medians_ne-quartile1_ne)/medians_ne,color='red',linestyle='dashed',label='Ne Median-Q1')
axs10[1].plot(data_iso_lilike['li_like_mg_pos']['T'],(quartile3_mg-medians_mg)/medians_mg,color='blue',linestyle='solid',label='Mg Q3-median')
#axs10[1].plot(data_iso_lilike['li_like_mg_pos']['T'],(medians_mg-quartile1_mg)/medians_mg,color='blue',linestyle='dashed',label='Mg Median-Q1')
axs10[1].plot(data_iso_lilike['li_like_na_pos']['T'],(quartile3_na-medians_na)/medians_na,color='green',linestyle='solid',label='Na Q3-median')
#axs10[1].plot(data_iso_lilike['li_like_na_pos']['T'],(medians_na-quartile1_na)/medians_na,color='green',linestyle='dashed',label='Na Median-Q1')
axs10[2].plot(data_iso_lilike['li_like_al_pos']['T'],(quartile3_al-medians_al)/medians_al,color='black',linestyle='solid',label='Al Q3-median')
#axs10[2].plot(data_iso_lilike['li_like_al_pos']['T'],(medians_al-quartile1_al)/medians_al,color='black',linestyle='dashed',label='Al Median-Q1')
axs10[2].plot(data_iso_lilike['li_like_si_pos']['T'],(quartile3_si-medians_si)/medians_si,color='red',linestyle='solid',label='Si Q3-median')
#axs10[2].plot(data_iso_lilike['li_like_si_pos']['T'],(medians_si-quartile1_si)/medians_si,color='red',linestyle='dashed',label='Si Median-Q1')
axs10[2].plot(data_iso_lilike['li_like_p_pos']['T'],(quartile3_p-medians_p)/medians_p,color='blue',linestyle='solid',label='P Q3-median')
#axs10[2].plot(data_iso_lilike['li_like_p_pos']['T'],(medians_p-quartile1_p)/medians_p,color='blue',linestyle='dashed',label='P Median-Q1')
axs10[2].plot(data_iso_lilike['li_like_s_pos']['T'],(quartile3_s-medians_s)/medians_s,color='green',linestyle='solid',label='S Q3-median')
#axs10[2].plot(data_iso_lilike['li_like_s_pos']['T'],(medians_s-quartile1_s)/medians_s,color='green',linestyle='dashed',label='S Median-Q1')
axs10[3].plot(data_iso_lilike['li_like_cl_pos']['T'],(quartile3_cl-medians_cl)/medians_cl,color='black',linestyle='solid',label='Cl Q3-median')
#axs10[3].plot(data_iso_lilike['li_like_cl_pos']['T'],(medians_cl-quartile1_cl)/medians_cl,color='black',linestyle='dashed',label='Cl Median-Q1')
axs10[3].plot(data_iso_lilike['li_like_k_pos']['T'],(quartile3_k-medians_k)/medians_k,color='red',linestyle='solid',label='K Q3-median')
#axs10[3].plot(data_iso_lilike['li_like_k_pos']['T'],(medians_k-quartile1_k)/medians_k,color='red',linestyle='dashed',label='K Median-Q1')
axs10[3].plot(data_iso_lilike['li_like_ca_pos']['T'],(quartile3_ca-medians_ca)/medians_ca,color='blue',linestyle='solid',label='Ca Q3-median')
#axs10[3].plot(data_iso_lilike['li_like_ca_pos']['T'],(medians_ca-quartile1_ca)/medians_ca,color='blue',linestyle='dashed',label='Ca Median-Q1')
axs10[3].plot(data_iso_lilike['li_like_sc_pos']['T'],(quartile3_sc-medians_sc)/medians_sc,color='green',linestyle='solid',label='Sc Q3-median')
#axs10[3].plot(data_iso_lilike['li_like_sc_pos']['T'],(medians_sc-quartile1_sc)/medians_sc,color='green',linestyle='dashed',label='Sc Median-Q1')
axs10[4].plot(data_iso_lilike['li_like_ti_pos']['T'],(quartile3_ti-medians_ti)/medians_ti,color='black',linestyle='solid',label='Ti Q3-median')
#axs10[4].plot(data_iso_lilike['li_like_ti_pos']['T'],(medians_ti-quartile1_ti)/medians_ti,color='black',linestyle='dashed',label='Ti Median-Q1')
axs10[4].plot(data_iso_lilike['li_like_v_pos']['T'],(quartile3_v-medians_v)/medians_v,color='red',linestyle='solid',label='Ca Q3-median')
#axs10[4].plot(data_iso_lilike['li_like_v_pos']['T'],(medians_v-quartile1_v)/medians_v,color='red',linestyle='dashed',label='Ca Median-Q1')
axs10[4].plot(data_iso_lilike['li_like_cr_pos']['T'],(quartile3_cr-medians_cr)/medians_cr,color='blue',linestyle='solid',label='Cr Q3-median')
#axs10[4].plot(data_iso_lilike['li_like_cr_pos']['T'],(medians_cr-quartile1_cr)/medians_cr,color='blue',linestyle='dashed',label='Cr Median-Q1')
axs10[4].plot(data_iso_lilike['li_like_mn_pos']['T'],(quartile3_mn-medians_mn)/medians_mn,color='green',linestyle='solid',label='Mn Q3-median')
#axs10[4].plot(data_iso_lilike['li_like_mn_pos']['T'],(medians_mn-quartile1_mn)/medians_mn,color='green',linestyle='dashed',label='Mn Median-Q1')
axs10[4].plot(data_iso_lilike['li_like_fe_pos']['T'],(quartile3_fe-medians_fe)/medians_fe,color='purple',linestyle='solid',label='Fe Q3-median')
#axs10[4].plot(data_iso_lilike['li_like_fe_pos']['T'],(medians_fe-quartile1_fe)/medians_fe,color='purple',linestyle='dashed',label='Fe Median-Q1')
axs10[0].legend(prop={"size":6})
axs10[1].legend(prop={"size":6})
axs10[2].legend(prop={"size":6})
axs10[3].legend(prop={"size":6})
axs10[4].legend(prop={"size":6})
plt.show()
plt.savefig('isoelectronic_Li_like_Q3.eps')

