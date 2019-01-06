#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:55:15 2018

@author: kyle
"""

import numpy as np
from recombination import State, dr_fit_func
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

seq = "be" #isoelectronic sequence
atom = "o"
shell="2-2" #a-b, recombining from shell a to shell b
nist_cutoff=0.05 #fractional deviation of energies from NIST values allowed 

ion = State(atom, seq, shell)
rates_file = f"rate_samples_{shell}_{int(100*nist_cutoff)}.dat"

def graph_rates_from_file(file, graph_every=100, x_range=None, y_range=None, animate=False):
    direc = "results/isoelectronic/{}/{}{}/".format(ion.isoelec_seq, ion.species, ion.ion_charge)
    
    data = np.genfromtxt(direc + file)
    T = data[0, :]
    rates = data[1:,:]
    
    err = np.std(rates, axis=0)
        
    graphed = rates[::100]
    if animate:
        fig, ax = plt.subplots()
        lines = []
        line, = ax.plot(T, graphed[0], "r", animated=True)
        
        def init():
            ax.set_title("Dielectronic Recombination of {}{}".format(ion.species, ion.ion_charge))
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("DR Rate (cm^3 s^-1)")
            ax.set_xscale("log")
            ax.set_yscale("log")
            if (x_range != None):
                ax.set_xlim(x_range[0], x_range[1])
            if (y_range != None):
                ax.set_ylim(y_range[0], y_range[1])
            line.set_xdata(T)
            return line,
        
        def update(i):
            line, = ax.plot(T, graphed[i,:], animated=True)
            lines.append(line)
            return lines
        
        ani = FuncAnimation(fig, update, frames=graphed.shape[0], init_func=init, blit=True)
            
        
        
     
        plt.show()
        ani.save(direc + "/rates_animation.mp4")
    
    else:
        fig, ax = plt.subplots(2, 1)
        for rate in graphed:
            ax[0].plot(T, rate)
        ax[0].set_title("Dielectronic Recombination of {}{}".format(ion.species, ion.ion_charge))
        ax[0].set_xlabel("Temperature (K)")
        ax[0].set_ylabel("DR Rate (cm^3 s^-1)")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        
        ax[1].plot(T, err)
        ax[1].set_title("Standard Devitation")
        ax[1].set_xlabel("Temperature (K)")
        ax[1].set_ylabel("Std. Dev. (cm^3 s^-1)")    
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")
        
        for a in ax:
            if (x_range != None):
                a.set_xlim(x_range[0], x_range[1])
            if (y_range != None):
                a.set_ylim(y_range[0], y_range[1])
        fig.tight_layout()
        fig.savefig(direc + "/rates_graph.png")
        
        
        #Used when plotting experimental data
        """
        plt.figure()
        temps = np.linspace(100, 100000000, 10000)
        plt.plot(temps, vogle_o4_f(temps), label="Experiment")
        plt.errorbar(x=T, y=avg, yerr=err, label="Theory")
        plt.yscale("log")
        plt.xscale("log")
        plt.xlabel("Temperature (K)")
        plt.ylabel("DR Rate (cm^3 s^-1)")
        plt.title("Dielectronic Recombination of Oxygen 4+")
        plt.legend()
        plt.savefig(direc + "/experiment_vs_theory.png")
        
        out = np.transpose(np.array([T, avg, err]))
        np.savetxt(direc + "/avg_rate.dat", out, header = "T, rate, absolute error")
        """
        

def graph_fits_from_file(ion, file):
    direc = "results/isoelectronic/{}/{}{}/".format(ion.isoelec_seq, ion.species, ion.ion_charge)
    coeff = np.genfromtxt(direc+file)
    n_points, n_coeff = np.shape(coeff)
    T = np.linspace(100, 1e8, 100000)
    plt.figure()
    for i in range(0, n_points, 1000):
        plt.semilogx(T, dr_fit_func(T, *coeff[i,:]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Rate Coefficient")
    plt.title("Fitted DR Rates")
    plt.savefig(direc + "{0}{1} Fit Plot.png".format(ion.species, ion.ion_charge))
    

if __name__ == "__main__":
    graph_rates_from_file(rates_file, x_range=(1e2, 1e6), y_range=(1e-15, 1e-10), animate=False)


