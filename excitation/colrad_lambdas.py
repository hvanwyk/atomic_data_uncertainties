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
    

    
    fig_G, ax_G = plt.subplots(2,1)
    fig_R, ax_R = plt.subplots(2,1)
    fig_X_temp, ax_X_temp = plt.subplots(2, 1)
    fig_Z_temp, ax_Z_temp = plt.subplots(2, 1)
    fig_X_dens, ax_X_dens = plt.subplots(2, 1)
    fig_Z_dens, ax_Z_dens = plt.subplots(2, 1)
    
    g_samples = []
    r_samples = []
    x_samples_temp = []
    x_samples_dens = []
    z_samples_temp = []
    z_samples_dens = []

    for rates in data[:, :, :]:
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
        x_samples_temp.append(x1)
        z_samples_temp.append(z1)
        temperatures = o.data['user']['temp_grid']
        
        ax_G[0].plot(temperatures, g)
        ax_X_temp[0].plot(temperatures, x1)
        ax_Z_temp[0].plot(temperatures, z1)
        
        
        temp_ind = 30
        x2 = o.data['processed']['pecs'][x_pec_ind, met_ind, temp_ind, :].ravel()
        y2 = o.data['processed']['pecs'][y_pec_ind, met_ind, temp_ind, :].ravel()
        z2 = o.data['processed']['pecs'][z_pec_ind, met_ind, temp_ind, :].ravel()
        
        r = z2/(x2+y2)
        r_samples.append(r)
        x_samples_dens.append(x2)
        z_samples_dens.append(z2)
        densities = o.data['user']['dens_grid']
        
        ax_R[0].plot(densities, r)
        ax_X_dens[0].plot(densities, x2)
        ax_Z_dens[0].plot(densities, z2)
        
        
    ax_G[0].set_xscale('log')
    #ax_G[0].set_yscale('log')
    ax_G[0].set_title(f"G Ratio vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+ at N$_e$={densities[dens_ind]} cm$^{-3}$")
    ax_G[0].set_xlabel("Electron Temperature (eV)")
    ax_G[0].set_ylabel("G Ratio")
        
    g_samples = np.array(g_samples)
    g_avg = np.mean(g_samples, axis=0)
    g_err = np.std(g_samples, axis=0)
    
    ax_G[1].plot(temperatures, (g_err/g_avg)*100)
    ax_G[1].set_xscale('log')
    ax_G[1].set_title("% Error in G Ratio")
    ax_G[1].set_xlabel("Electron Temperature (eV)")
    ax_G[1].set_ylabel("% Error")

    fig_G.tight_layout()
    
    fig_G.savefig("G Ratio vs. T.eps")
    
    
    ax_R[0].set_xscale('log')
    #ax_R[0].set_yscale('log')
    ax_R[0].set_title(f"R Ratio vs. Density for {ion.species.capitalize()}{ion.ion_charge}+ at T$_e$={temperatures[temp_ind]:.2f}eV")
    ax_R[0].set_xlabel("Density (cm$^{-3}$)")
    ax_R[0].set_ylabel("R Ratio")
    
    r_samples = np.array(r_samples)
    r_avg = np.mean(r_samples, axis=0)
    r_err = np.std(r_samples, axis=0)
    
    ax_R[1].plot(densities, (r_err/r_avg)*100)
    ax_R[1].set_xscale('log')
    ax_R[1].set_title("% Error in R Ratio")
    ax_R[1].set_xlabel("Density (cm$^{-3}$)")
    ax_R[1].set_ylabel("% Error")
    
    fig_R.tight_layout()
    fig_R.savefig("R Ratio vs. Density.eps")
    
    
    ax_X_temp[0].set_xscale('log')
    ax_X_temp[0].set_title(f"X PEC vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+ at N$_e$={densities[dens_ind]}cm$^{-3}$")
    ax_X_temp[0].set_xlabel("Electron Temperature (eV)")
    ax_X_temp[0].set_ylabel("X PEC")
    
    x_samples_temp = np.array(x_samples_temp)
    x_avg_temp = np.mean(x_samples_temp, axis=0)
    x_err_temp = np.std(x_samples_temp, axis=0)
    ax_X_temp[1].plot(temperatures, (x_err_temp/x_avg_temp)*100)
    ax_X_temp[1].set_xscale('log')
    ax_X_temp[1].set_title("% Error in X Line PEC")
    ax_X_temp[1].set_xlabel("Electron Temperature (eV)")
    ax_X_temp[1].set_ylabel("% Error")
    
<<<<<<< HEAD
<<<<<<< HEAD
    fig_X_temp.tight_layout()
    fig_X_temp.savefig("X PEC vs. Temperature.eps")
    
    ax_Z_temp[0].set_xscale('log')
    ax_Z_temp[0].set_title(f"Z PEC vs. Temperature for {ion.species.capitalize()}{ion.ion_charge}+ at N$_e$={densities[dens_ind]}cm$^{-3}$")
    ax_Z_temp[0].set_xlabel("Electron Temperature (eV)")
    ax_Z_temp[0].set_ylabel("Z PEC")
    
    z_samples_temp = np.array(z_samples_temp)
    z_avg_temp = np.mean(z_samples_temp, axis=0)
    z_err_temp = np.std(z_samples_temp, axis=0)
    ax_Z_temp[1].plot(temperatures, (z_err_temp/z_avg_temp)*100)
    ax_Z_temp[1].set_xscale('log')
    ax_Z_temp[1].set_title("% Error in Z Line PEC")
    ax_Z_temp[1].set_xlabel("Electron Temperature (eV)")
    ax_Z_temp[1].set_ylabel("% Error")
    
    fig_Z_temp.tight_layout()
    fig_Z_temp.savefig("Z PEC vs. Temperature.eps")

    
    
    
    ax_X_dens[0].set_xscale('log')
    ax_X_dens[0].set_title(f"X PEC vs. Density for {ion.species.capitalize()}{ion.ion_charge}+ at T$_e$={temperatures[temp_ind]:.2f}eV")
    ax_X_dens[0].set_xlabel("Density (cm$^{-3}$)")
    ax_X_dens[0].set_ylabel("X PEC")
    
    x_samples_dens = np.array(x_samples_dens)
    x_avg_dens = np.mean(x_samples_dens, axis=0)
    x_err_dens = np.std(x_samples_dens, axis=0)
    ax_X_dens[1].plot(densities, (x_err_dens/x_avg_dens)*100)
    ax_X_dens[1].set_xscale('log')
    ax_X_dens[1].set_title("% Error in X Line PEC")
    ax_X_dens[1].set_xlabel("Density (cm$^{-3}$)")
    ax_X_dens[1].set_ylabel("% Error")
    
    fig_X_dens.tight_layout()
    fig_X_dens.savefig("X PEC vs. Density.eps")

    
    ax_Z_dens[0].set_xscale('log')
    ax_Z_dens[0].set_title(f"Z PEC vs. Density for {ion.species.capitalize()}{ion.ion_charge}+ at T$_e$={temperatures[temp_ind]:.2f}eV")
    ax_Z_dens[0].set_xlabel("Density (cm$^{-3}$)")
    ax_Z_dens[0].set_ylabel("Z PEC")
    
    z_samples_dens = np.array(z_samples_dens)
    z_avg_dens = np.mean(z_samples_dens, axis=0)
    z_err_dens = np.std(z_samples_dens, axis=0)
    ax_Z_dens[1].plot(densities, (z_err_dens/z_avg_dens)*100)
    ax_Z_dens[1].set_xscale('log')
    ax_Z_dens[1].set_title("% Error in Z Line PEC")
    ax_Z_dens[1].set_xlabel("Density (cm$^{-3}$)")
    ax_Z_dens[1].set_ylabel("% Error")
    
    fig_Z_dens.tight_layout()
    fig_Z_dens.savefig("Z PEC vs. Density.eps")
=======
    
>>>>>>> c5f04a32d40bbe0cf6f76c520a256c36d1059214
=======
    
>>>>>>> c5f04a32d40bbe0cf6f76c520a256c36d1059214
