#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 15:43:22 2021

@author: lochstu
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle

data1=pickle.load(open('be_like_0.4_1.6_b_mg.pkl','rb'))
data2=pickle.load(open('be_like_0.4_1.6_k_fe.pkl','rb'))



b_rate_avg_pos=data1['be_like_b_pos']['rate_avg']
b_rate_avg_neg=data1['be_like_b_neg']['rate_avg']
b_rate_std_pos=data1['be_like_b_pos']['rate_std']
b_rate_std_neg=data1['be_like_b_neg']['rate_std']
b_rate_percent_pos=data1['be_like_b_pos']['rate_percent']
b_rate_percent_neg=data1['be_like_b_neg']['rate_percent']
b_rate_samples_pos=data1['be_like_b_pos']['rate_samples']
b_rate_samples_neg=data1['be_like_b_neg']['rate_samples']
b_T=data1['be_like_b_pos']['T']
b_rate_samples_combined=np.concatenate((b_rate_samples_pos,b_rate_samples_neg))
b_rate_avg_combined = np.average(b_rate_samples_combined,axis=0)
b_rate_std_combined = np.std(b_rate_samples_combined,axis=0)
b_rate_percent_combined = b_rate_std_combined/b_rate_avg_combined*100.
b_atom='B'
seq=data1['seq']
b_T_ev=b_T/11604.

c_rate_avg_pos=data1['be_like_c_pos']['rate_avg']
c_rate_avg_neg=data1['be_like_c_neg']['rate_avg']
c_rate_std_pos=data1['be_like_c_pos']['rate_std']
c_rate_std_neg=data1['be_like_c_neg']['rate_std']
c_rate_percent_pos=data1['be_like_c_pos']['rate_percent']
c_rate_percent_neg=data1['be_like_c_neg']['rate_percent']
c_rate_samples_pos=data1['be_like_c_pos']['rate_samples']
c_rate_samples_neg=data1['be_like_c_neg']['rate_samples']
c_T=data1['be_like_c_pos']['T']
c_rate_samples_combined=np.concatenate((c_rate_samples_pos,c_rate_samples_neg))
c_rate_avg_combined = np.average(c_rate_samples_combined,axis=0)
c_rate_std_combined = np.std(c_rate_samples_combined,axis=0)
c_rate_percent_combined = c_rate_std_combined/c_rate_avg_combined*100.
c_atom='C'
c_T_ev=c_T/11604.

n_rate_avg_pos=data1['be_like_n_pos']['rate_avg']
n_rate_avg_neg=data1['be_like_n_neg']['rate_avg']
n_rate_std_pos=data1['be_like_n_pos']['rate_std']
n_rate_std_neg=data1['be_like_n_neg']['rate_std']
n_rate_percent_pos=data1['be_like_n_pos']['rate_percent']
n_rate_percent_neg=data1['be_like_n_neg']['rate_percent']
n_rate_samples_pos=data1['be_like_n_pos']['rate_samples']
n_rate_samples_neg=data1['be_like_n_neg']['rate_samples']
n_T=data1['be_like_n_pos']['T']
n_rate_samples_combined=np.concatenate((n_rate_samples_pos,n_rate_samples_neg))
n_rate_avg_combined = np.average(n_rate_samples_combined,axis=0)
n_rate_std_combined = np.std(n_rate_samples_combined,axis=0)
n_rate_percent_combined = n_rate_std_combined/n_rate_avg_combined*100.
n_atom='N'
seq=data1['seq']
n_T_ev=n_T/11604.

o_rate_avg_pos=data1['be_like_o_pos']['rate_avg']
o_rate_avg_neg=data1['be_like_o_neg']['rate_avg']
o_rate_std_pos=data1['be_like_o_pos']['rate_std']
o_rate_std_neg=data1['be_like_o_neg']['rate_std']
o_rate_percent_pos=data1['be_like_o_pos']['rate_percent']
o_rate_percent_neg=data1['be_like_o_neg']['rate_percent']
o_rate_samples_pos=data1['be_like_o_pos']['rate_samples']
o_rate_samples_neg=data1['be_like_o_neg']['rate_samples']
o_T=data1['be_like_o_pos']['T']
o_rate_samples_combined=np.concatenate((o_rate_samples_pos,o_rate_samples_neg))
o_rate_avg_combined = np.average(o_rate_samples_combined,axis=0)
o_rate_std_combined = np.std(o_rate_samples_combined,axis=0)
o_rate_percent_combined = o_rate_std_combined/o_rate_avg_combined*100.
o_atom='O'
o_T_ev=o_T/11604.


f_rate_avg_pos=data1['be_like_f_pos']['rate_avg']
f_rate_avg_neg=data1['be_like_f_neg']['rate_avg']
f_rate_std_pos=data1['be_like_f_pos']['rate_std']
f_rate_std_neg=data1['be_like_f_neg']['rate_std']
f_rate_percent_pos=data1['be_like_f_pos']['rate_percent']
f_rate_percent_neg=data1['be_like_f_neg']['rate_percent']
f_rate_samples_pos=data1['be_like_f_pos']['rate_samples']
f_rate_samples_neg=data1['be_like_f_neg']['rate_samples']
f_T=data1['be_like_f_pos']['T']
f_rate_samples_combined=np.concatenate((f_rate_samples_pos,f_rate_samples_neg))
f_rate_avg_combined = np.average(f_rate_samples_combined,axis=0)
f_rate_std_combined = np.std(f_rate_samples_combined,axis=0)
f_rate_percent_combined = f_rate_std_combined/f_rate_avg_combined*100.
f_atom='F'
f_T_ev=f_T/11604.

ne_rate_avg_pos=data1['be_like_ne_pos']['rate_avg']
ne_rate_avg_neg=data1['be_like_ne_neg']['rate_avg']
ne_rate_std_pos=data1['be_like_ne_pos']['rate_std']
ne_rate_std_neg=data1['be_like_ne_neg']['rate_std']
ne_rate_percent_pos=data1['be_like_ne_pos']['rate_percent']
ne_rate_percent_neg=data1['be_like_ne_neg']['rate_percent']
ne_rate_samples_pos=data1['be_like_ne_pos']['rate_samples']
ne_rate_samples_neg=data1['be_like_ne_neg']['rate_samples']
ne_T=data1['be_like_ne_pos']['T']
ne_rate_samples_combined=np.concatenate((ne_rate_samples_pos,ne_rate_samples_neg))
ne_rate_avg_combined = np.average(ne_rate_samples_combined,axis=0)
ne_rate_std_combined = np.std(ne_rate_samples_combined,axis=0)
ne_rate_percent_combined = ne_rate_std_combined/ne_rate_avg_combined*100.
ne_atom='Ne'
ne_T_ev=ne_T/11604.

na_rate_avg_pos=data1['be_like_na_pos']['rate_avg']
na_rate_avg_neg=data1['be_like_na_neg']['rate_avg']
na_rate_std_pos=data1['be_like_na_pos']['rate_std']
na_rate_std_neg=data1['be_like_na_neg']['rate_std']
na_rate_percent_pos=data1['be_like_na_pos']['rate_percent']
na_rate_percent_neg=data1['be_like_na_neg']['rate_percent']
na_rate_samples_pos=data1['be_like_na_pos']['rate_samples']
na_rate_samples_neg=data1['be_like_na_neg']['rate_samples']
na_T=data1['be_like_na_pos']['T']
na_rate_samples_combined=np.concatenate((na_rate_samples_pos,na_rate_samples_neg))
na_rate_avg_combined = np.average(na_rate_samples_combined,axis=0)
na_rate_std_combined = np.std(na_rate_samples_combined,axis=0)
na_rate_percent_combined = na_rate_std_combined/na_rate_avg_combined*100.
na_atom='Na'
na_T_ev=na_T/11604.

mg_rate_avg_pos=data1['be_like_mg_pos']['rate_avg']
mg_rate_avg_neg=data1['be_like_mg_neg']['rate_avg']
mg_rate_std_pos=data1['be_like_mg_pos']['rate_std']
mg_rate_std_neg=data1['be_like_mg_neg']['rate_std']
mg_rate_percent_pos=data1['be_like_mg_pos']['rate_percent']
mg_rate_percent_neg=data1['be_like_mg_neg']['rate_percent']
mg_rate_samples_pos=data1['be_like_mg_pos']['rate_samples']
mg_rate_samples_neg=data1['be_like_mg_neg']['rate_samples']
mg_T=data1['be_like_mg_pos']['T']
mg_rate_samples_combined=np.concatenate((mg_rate_samples_pos,mg_rate_samples_neg))
mg_rate_avg_combined = np.average(mg_rate_samples_combined,axis=0)
mg_rate_std_combined = np.std(mg_rate_samples_combined,axis=0)
mg_rate_percent_combined = mg_rate_std_combined/mg_rate_avg_combined*100.
mg_atom='Mg'
mg_T_ev=mg_T/11604.

k_rate_avg_pos=data2['be_like_k_pos']['rate_avg']
k_rate_avg_neg=data2['be_like_k_neg']['rate_avg']
k_rate_std_pos=data2['be_like_k_pos']['rate_std']
k_rate_std_neg=data2['be_like_k_neg']['rate_std']
k_rate_percent_pos=data2['be_like_k_pos']['rate_percent']
k_rate_percent_neg=data2['be_like_k_neg']['rate_percent']
k_rate_samples_pos=data2['be_like_k_pos']['rate_samples']
k_rate_samples_neg=data2['be_like_k_neg']['rate_samples']
k_T=data2['be_like_k_pos']['T']
k_rate_samples_combined=np.concatenate((k_rate_samples_pos,k_rate_samples_neg))
k_rate_avg_combined = np.average(k_rate_samples_combined,axis=0)
k_rate_std_combined = np.std(k_rate_samples_combined,axis=0)
k_rate_percent_combined = k_rate_std_combined/k_rate_avg_combined*100.
k_atom='K'
k_T_ev=k_T/11604.

ca_rate_avg_pos=data2['be_like_ca_pos']['rate_avg']
ca_rate_avg_neg=data2['be_like_ca_neg']['rate_avg']
ca_rate_std_pos=data2['be_like_ca_pos']['rate_std']
ca_rate_std_neg=data2['be_like_ca_neg']['rate_std']
ca_rate_percent_pos=data2['be_like_ca_pos']['rate_percent']
ca_rate_percent_neg=data2['be_like_ca_neg']['rate_percent']
ca_rate_samples_pos=data2['be_like_ca_pos']['rate_samples']
ca_rate_samples_neg=data2['be_like_ca_neg']['rate_samples']
ca_T=data2['be_like_ca_pos']['T']
ca_rate_samples_combined=np.concatenate((ca_rate_samples_pos,ca_rate_samples_neg))
ca_rate_avg_combined = np.average(ca_rate_samples_combined,axis=0)
ca_rate_std_combined = np.std(ca_rate_samples_combined,axis=0)
ca_rate_percent_combined = ca_rate_std_combined/ca_rate_avg_combined*100.
ca_atom='Ca'
ca_T_ev=ca_T/11604.

sc_rate_avg_pos=data2['be_like_sc_pos']['rate_avg']
sc_rate_avg_neg=data2['be_like_sc_neg']['rate_avg']
sc_rate_std_pos=data2['be_like_sc_pos']['rate_std']
sc_rate_std_neg=data2['be_like_sc_neg']['rate_std']
sc_rate_percent_pos=data2['be_like_sc_pos']['rate_percent']
sc_rate_percent_neg=data2['be_like_sc_neg']['rate_percent']
sc_rate_samples_pos=data2['be_like_sc_pos']['rate_samples']
sc_rate_samples_neg=data2['be_like_sc_neg']['rate_samples']
sc_T=data2['be_like_sc_pos']['T']
sc_rate_samples_combined=np.concatenate((sc_rate_samples_pos,sc_rate_samples_neg))
sc_rate_avg_combined = np.average(sc_rate_samples_combined,axis=0)
sc_rate_std_combined = np.std(sc_rate_samples_combined,axis=0)
sc_rate_percent_combined = sc_rate_std_combined/sc_rate_avg_combined*100.
sc_atom='Sc'
sc_T_ev=sc_T/11604.

ti_rate_avg_pos=data2['be_like_ti_pos']['rate_avg']
ti_rate_avg_neg=data2['be_like_ti_neg']['rate_avg']
ti_rate_std_pos=data2['be_like_ti_pos']['rate_std']
ti_rate_std_neg=data2['be_like_ti_neg']['rate_std']
ti_rate_percent_pos=data2['be_like_ti_pos']['rate_percent']
ti_rate_percent_neg=data2['be_like_ti_neg']['rate_percent']
ti_rate_samples_pos=data2['be_like_ti_pos']['rate_samples']
ti_rate_samples_neg=data2['be_like_ti_neg']['rate_samples']
ti_T=data2['be_like_ti_pos']['T']
ti_rate_samples_combined=np.concatenate((ti_rate_samples_pos,ti_rate_samples_neg))
ti_rate_avg_combined = np.average(ti_rate_samples_combined,axis=0)
ti_rate_std_combined = np.std(ti_rate_samples_combined,axis=0)
ti_rate_percent_combined = ti_rate_std_combined/ti_rate_avg_combined*100.
ti_atom='Ti'
ti_T_ev=ti_T/11604.

v_rate_avg_pos=data2['be_like_v_pos']['rate_avg']
v_rate_avg_neg=data2['be_like_v_neg']['rate_avg']
v_rate_std_pos=data2['be_like_v_pos']['rate_std']
v_rate_std_neg=data2['be_like_v_neg']['rate_std']
v_rate_percent_pos=data2['be_like_v_pos']['rate_percent']
v_rate_percent_neg=data2['be_like_v_neg']['rate_percent']
v_rate_samples_pos=data2['be_like_v_pos']['rate_samples']
v_rate_samples_neg=data2['be_like_v_neg']['rate_samples']
v_T=data2['be_like_v_pos']['T']
v_rate_samples_combined=np.concatenate((v_rate_samples_pos,v_rate_samples_neg))
v_rate_avg_combined = np.average(v_rate_samples_combined,axis=0)
v_rate_std_combined = np.std(v_rate_samples_combined,axis=0)
v_rate_percent_combined = v_rate_std_combined/v_rate_avg_combined*100.
v_atom='V'
v_T_ev=v_T/11604.

cr_rate_avg_pos=data2['be_like_cr_pos']['rate_avg']
cr_rate_avg_neg=data2['be_like_cr_neg']['rate_avg']
cr_rate_std_pos=data2['be_like_cr_pos']['rate_std']
cr_rate_std_neg=data2['be_like_cr_neg']['rate_std']
cr_rate_percent_pos=data2['be_like_cr_pos']['rate_percent']
cr_rate_percent_neg=data2['be_like_cr_neg']['rate_percent']
cr_rate_samples_pos=data2['be_like_cr_pos']['rate_samples']
cr_rate_samples_neg=data2['be_like_cr_neg']['rate_samples']
cr_T=data2['be_like_cr_pos']['T']
cr_rate_samples_combined=np.concatenate((cr_rate_samples_pos,cr_rate_samples_neg))
cr_rate_avg_combined = np.average(cr_rate_samples_combined,axis=0)
cr_rate_std_combined = np.std(cr_rate_samples_combined,axis=0)
cr_rate_percent_combined = cr_rate_std_combined/cr_rate_avg_combined*100.
cr_atom='Cr'
cr_T_ev=cr_T/11604.

mn_rate_avg_pos=data2['be_like_mn_pos']['rate_avg']
mn_rate_avg_neg=data2['be_like_mn_neg']['rate_avg']
mn_rate_std_pos=data2['be_like_mn_pos']['rate_std']
mn_rate_std_neg=data2['be_like_mn_neg']['rate_std']
mn_rate_percent_pos=data2['be_like_mn_pos']['rate_percent']
mn_rate_percent_neg=data2['be_like_mn_neg']['rate_percent']
mn_rate_samples_pos=data2['be_like_mn_pos']['rate_samples']
mn_rate_samples_neg=data2['be_like_mn_neg']['rate_samples']
mn_T=data2['be_like_mn_pos']['T']
mn_rate_samples_combined=np.concatenate((mn_rate_samples_pos,mn_rate_samples_neg))
mn_rate_avg_combined = np.average(mn_rate_samples_combined,axis=0)
mn_rate_std_combined = np.std(mn_rate_samples_combined,axis=0)
mn_rate_percent_combined = mn_rate_std_combined/mn_rate_avg_combined*100.
mn_atom='Mn'
mn_T_ev=mn_T/11604.

fe_rate_avg_pos=data2['be_like_fe_pos']['rate_avg']
fe_rate_avg_neg=data2['be_like_fe_neg']['rate_avg']
fe_rate_std_pos=data2['be_like_fe_pos']['rate_std']
fe_rate_std_neg=data2['be_like_fe_neg']['rate_std']
fe_rate_percent_pos=data2['be_like_fe_pos']['rate_percent']
fe_rate_percent_neg=data2['be_like_fe_neg']['rate_percent']
fe_rate_samples_pos=data2['be_like_fe_pos']['rate_samples']
fe_rate_samples_neg=data2['be_like_fe_neg']['rate_samples']
fe_T=data2['be_like_fe_pos']['T']
fe_rate_samples_combined=np.concatenate((fe_rate_samples_pos,fe_rate_samples_neg))
fe_rate_avg_combined = np.average(fe_rate_samples_combined,axis=0)
fe_rate_std_combined = np.std(fe_rate_samples_combined,axis=0)
fe_rate_percent_combined = fe_rate_std_combined/fe_rate_avg_combined*100.
fe_atom='Fe'
fe_T_ev=fe_T/11604.


fig3, axs3 = plt.subplots(1,1)
title_pos='Uncertainty coefficient for ' + seq + '-like' 
axs3.set_xscale("log")
axs3.set_yscale("log")
axs3.set_xlabel('Electron Temperature (K)')
axs3.set_ylabel('Percentage uncertainty')
axs3.set_title(title_pos)
axs3.plot(b_T,b_rate_percent_combined,label=b_atom)
axs3.plot(c_T,c_rate_percent_combined,label=c_atom)
axs3.plot(n_T,n_rate_percent_combined,label=n_atom)
axs3.plot(o_T,o_rate_percent_combined,label=o_atom)
axs3.plot(f_T,f_rate_percent_combined,label=f_atom)
axs3.plot(ne_T,ne_rate_percent_combined,label=ne_atom)
axs3.plot(na_T,na_rate_percent_combined,label=na_atom)
axs3.plot(mg_T,mg_rate_percent_combined,label=mg_atom)
axs3.plot(k_T,k_rate_percent_combined,label=k_atom)
axs3.plot(ca_T,ca_rate_percent_combined,label=ca_atom)
axs3.plot(sc_T,sc_rate_percent_combined,label=sc_atom)
axs3.plot(ti_T,ti_rate_percent_combined,label=ti_atom)
axs3.plot(v_T,v_rate_percent_combined,label=v_atom)
axs3.plot(cr_T,cr_rate_percent_combined,label=cr_atom)
axs3.plot(mn_T,mn_rate_percent_combined,label=mn_atom)
axs3.plot(fe_T,fe_rate_percent_combined,label=fe_atom)
axs3.legend()


