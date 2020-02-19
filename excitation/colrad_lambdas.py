#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 10:49:29 2020

@author: kyle
"""

import numpy as np
import pandas as pd
import sys
import re
if ".." not in sys.path:
    sys.path.append("..")
from state import State
from read_adf04 import read_adf04, rates_dataframe
from excitation import run_r_matrix, orbitals
if "../../ColRadPy/" not in sys.path:
    sys.path.append("../../ColRadPy/")
from colradpy import colradpy
import matplotlib.pyplot as plt


def write_adf04(filename, ion, rates=[], adf04_template=f"isoelectronic/he-like/o6/adas/adf04", nmax=3):
    if adf04_template is None:
        lmb = [1.0]*len(orbitals(ion, nmax)) 
        run_r_matrix(ion, lambdas=lmb, nmax=nmax)
        adf04_template = f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04"
    hdr = read_adf04(adf04_template)[3]
    hdr = re.sub("\(  \)", "(2S)", hdr)
    df_all = pd.read_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_dataframe_template_all.csv")
    df_ground = pd.read_csv(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/rates_dataframe_template_ground.csv")
    df_ground[list(df_ground.columns.values[2:])] = rates
    
    df_all[df_all["final"]==1] = df_ground
    with open(filename, "w") as file:
        file.write(hdr)
        for i in range(df_all.values.shape[0]):
            line = list(df_all.values[i, :])
            line[0] = str(int(line[0]))
            line[1] = str(int(line[1]))
            
            file.write(" "*(4-len(line[0])))
            file.write(str(int(line[0])))
            file.write(" "*(4-len(line[1])))
            file.write(str(int(line[1])))
            
            for j in range(2, len(line)):
                line[j] = "{:.2e}".format(line[j])
                if line[j][0]=="-" and j!=2:
                    file.write(f"{line[j]}".replace("e", ""))
                else:
                    file.write(f" {line[j]}".replace("e", ""))
            file.write("\n")
        file.write("  -1\n")
        file.write("  -1  -1\n")
        
        
if __name__ == "__main__":
    
    atom = "o"
    seq = "he"
    shell = "1-2"
    
    ion = State(atom, seq, shell)
    T, data = np.load("rates_gridres_5.npy")
    
   
    adf04 = "test_adf04"
    metastable_levels = np.array([0])
    temperature_arr = np.geomspace(100, 1000, 50)
    density_arr  = np.geomspace(1e1, 1e18, 50)
    
    o = colradpy(adf04, metastable_levels, temperature_arr, density_arr, use_recombination=False, 
              use_recombination_three_body=False, use_ionization=False, suppliment_with_ecip=False)
    o.solve_cr()
    

    
    fig1, ax1 = plt.subplots(2,1)
    fig2, ax2 = plt.subplots(2,1)
    
    g_samples = []
    r_samples = []

    for rates in data[::500, :, :]:
        write_adf04("test_adf04", ion, rates=rates, adf04_template=f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/adas/adf04")
        adf04 = "test_adf04"
        metastable_levels = np.array([0])
        temperature_arr = np.geomspace(100, 1000, 50)
        density_arr  = np.geomspace(1e1, 1e18, 50)
        
        o = colradpy(adf04, metastable_levels, temperature_arr, density_arr, use_recombination=False, 
                  use_recombination_three_body=False, use_ionization=False, suppliment_with_ecip=False)
        o.solve_cr()
        
        levels = read_adf04(adf04)[0]
        w_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(1)1( 1.0)'), "#"].values[0])-1
        x_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 2.0)'), "#"].values[0])-1
        y_num = int(levels.loc[(levels['config']=='1S1 2P1') & (levels['(2S+1)L( 2J)'] == '(3)1( 1.0)'), "#"].values[0])-1
        z_num = int(levels.loc[(levels['config']=='1S1 2S1') & (levels['(2S+1)L( 2J)'] == '(3)0( 1.0)'), "#"].values[0])-1
                               
        
        x_pec_ind = np.where( (o.data['processed']['pec_levels'][:,0] == x_num) & (o.data['processed']['pec_levels'][:,1] == 0))[0]
        y_pec_ind = np.where( (o.data['processed']['pec_levels'][:,0] == y_num) & (o.data['processed']['pec_levels'][:,1] == 0))[0]
        z_pec_ind = np.where( (o.data['processed']['pec_levels'][:,0] == z_num) & (o.data['processed']['pec_levels'][:,1] == 0))[0]
        w_pec_ind = np.where( (o.data['processed']['pec_levels'][:,0] == w_num) & (o.data['processed']['pec_levels'][:,1] == 0))[0]
        
        met_ind = 0
        dens_ind = 0
        
        x1 = o.data['processed']['pecs'][x_pec_ind, met_ind, :, dens_ind].ravel()
        y1 = o.data['processed']['pecs'][y_pec_ind, met_ind, :, dens_ind].ravel()
        z1 = o.data['processed']['pecs'][z_pec_ind, met_ind, :, dens_ind].ravel()
        w1 = o.data['processed']['pecs'][w_pec_ind, met_ind, :, dens_ind].ravel()
        
        g = (x1+y1+z1)/w1
        g_samples.append(g)
        temperatures = o.data['user']['temp_grid']
        
        ax1[0].plot(temperatures, g)
        
        
        temp_ind = 30
        x2 = o.data['processed']['pecs'][x_pec_ind, met_ind, temp_ind, :].ravel()
        y2 = o.data['processed']['pecs'][y_pec_ind, met_ind, temp_ind, :].ravel()
        z2 = o.data['processed']['pecs'][z_pec_ind, met_ind, temp_ind, :].ravel()
        
        r = z2/(x2+y2)
        r_samples.append(r)
        densities = o.data['user']['dens_grid']
        
        ax2[0].plot(densities, r)
        
        
    ax1[0].set_xscale('log')
    #ax1[0].set_yscale('log')
    ax1[0].set_title(f"G Ratio vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+ at a Density of {densities[dens_ind]} cm$^{-3}$")
    ax1[0].set_xlabel("Electron Temperature (eV)")
    ax1[0].set_ylabel("G Ratio")
        
    g_samples = np.array(g_samples)
    g_avg = np.mean(g_samples, axis=0)
    g_err = np.std(g_samples, axis=0)
    
    ax1[1].plot(temperatures, (g_err/g_avg)*100)
    ax1[1].set_xscale('log')
    ax1[1].set_title("% Error in G Ratio")
    ax1[1].set_xlabel("Electron Temperature (eV)")
    ax1[1].set_ylabel("% Error")

    fig1.tight_layout()
    
    fig1.savefig("G Ratio vs. T.eps")
    
    
    ax2[0].set_xscale('log')
    #ax2[0].set_yscale('log')
    ax2[0].set_title(f"R Ratio vs. Density for {ion.species.capitalize()}{ion.ion_charge}+ at Electron Temperature of {temperatures[temp_ind]:.2f} eV")
    ax2[0].set_xlabel("Density (cm$^{-3}$)")
    ax2[0].set_ylabel("R Ratio")
    
    r_samples = np.array(r_samples)
    r_avg = np.mean(r_samples, axis=0)
    r_err = np.std(r_samples, axis=0)
    
    ax2[1].plot(densities, (r_err/r_avg)*100)
    ax2[1].set_xscale('log')
    ax2[1].set_title("% Error in R Ratio")
    ax2[1].set_xlabel("Density (cm$^{-3}$)")
    ax2[1].set_ylabel("% Error")
    
    fig2.tight_layout()
    
    fig2.savefig("R Ratio vs. Density.eps")
    
    
    
    