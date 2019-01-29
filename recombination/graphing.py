#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 20 11:55:15 2018

@author: kyle
"""

import numpy as np
from recombination_methods import State, dr_fit_func
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

def experimental_fit(c, E, T):
    k = 8.617e-5
    fit = np.zeros(T.shape)
    for i in range(len(c)):
        fit += c[i]*np.exp(-E[i]/(k*T))
    fit /= (T**1.5)
    return fit

def graph_rates_from_file(ion, nist_cutoff=0.05, prior_shape="uniform", likelihood_shape="uniform", 
                          graph_every=10, x_range=None, y_range=None, err_type="std", animate=False, experimental=False):
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"

    infile = direc + f"rates_{ion.shell}_{int(100*nist_cutoff)}_{prior_shape}P_{likelihood_shape}L.npy"
    outfile =  direc + f"rates_{ion.shell}_{int(100*nist_cutoff)}_{prior_shape}P_{likelihood_shape}L.png"
    #data = np.genfromtxt(infile)
    data = np.load(infile)
    T = data[0]
    rates = data[1]
    
    avg = np.mean(rates, axis=0)
    if err_type=="std":
        err = np.std(rates, axis=0)
    elif err_type=="max":
        err = np.max(rates, axis=0) - np.min(rates, axis=0)
        
    graphed = rates[::100]
    if animate:
        fig, ax = plt.subplots()
        lines = []
        line, = ax.plot(T, graphed[0], "r", animated=True)
        
        def init():
            ax.set_title(f"Dielectronic Recombination of {ion.species}{ion.ion_charge}")
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
        
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)
        ani = FuncAnimation(fig, update, frames=graphed.shape[0], init_func=init, blit=True)
        
        ani.save(filename="f.mp4", writer=writer)
    
    else:
        fig, ax = plt.subplots(2, 1)
        for rate in graphed:
            ax[0].plot(T, rate)
        ax[0].set_title(f"Dielectronic Recombination of {ion.species}{ion.ion_charge}")
        ax[0].set_xlabel("Temperature (K)")
        ax[0].set_ylabel("DR Rate (cm^3 s^-1)")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        
        ax[1].plot(T, err)
        if err_type=="std":
            ax[1].set_title("Standard Devitation")
        elif err_type=="max":
            ax[1].set_title("Max-Min")
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
        fig.savefig(outfile)
        
        
        #Used when plotting experimental data
        if experimental==True:
            plt.figure()
            temps = np.linspace(100, 100000000, 1000000)
            plt.plot(temps, experimental_fit(c, E, temps), label="Experiment")
            plt.errorbar(x=T, y=avg, yerr=err, label="Theory")
            plt.yscale("log")
            plt.xscale("log")
            plt.xlabel("Temperature (K)")
            plt.ylabel("DR Rate (cm^3 s^-1)")
            plt.title("Dielectronic Recombination of Oxygen 4+")
            plt.legend()
            plt.savefig(direc + "/experiment_vs_theory.png")
        
        
        plt.figure()
        T_hist_choice = 1
        plt.hist(rates[:, T_hist_choice], 100)
        plt.title(f"T={T[T_hist_choice]}K Rates Histogram")
        plt.xlabel("DR Rate (cm^3 s^-1)")
        plt.ylabel("Count")
        plt.savefig(direc + f"histogram_{T_hist_choice}")
        

def graph_fits_from_file(ion, file):
    direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/"
    coeff = np.genfromtxt(direc+file)
    n_points, n_coeff = np.shape(coeff)
    T = np.linspace(100, 1e8, 100000)
    plt.figure()
    for i in range(0, n_points, 1000):
        plt.semilogx(T, dr_fit_func(T, *coeff[i,:]))
    plt.xlabel("Temperature (K)")
    plt.ylabel("Rate Coefficient")
    plt.title("Fitted DR Rates")
    plt.savefig(direc + f"{ion.species}{ion.ion_charge} Fit Plot.png")
    

def r_ratio(fname, animate=False):
    data = np.genfromtxt(fname, skip_header=2)
    
    if animate:
        fig, ax = plt.subplots()
        x_dat, y_dat = [], []
        ln, = ax.plot([], [], "ro", animated=True)
        
        def init():
            ax.set_xlim((15e-19, 1e-17))
            ax.set_ylim((2e-15, 1e-14))
            
            return ln,
        
        def update(i):
            x_dat.append(float(data[i,0]))
            y_dat.append(float(data[i,3]) + float(data[i,6]))
            ln.set_data(x_dat, y_dat)
            return ln,
        
        anim = FuncAnimation(fig, update, frames=data.shape[0], init_func=init, blit=True)
        
        plt.show()
        anim.save("test3.mp4", writer="ffmpeg")
        
    else:
        x_dat = data[:,0]
        y_dat = data[:,3] + data[:,6]
        plt.scatter(x_dat, y_dat)
        plt.xlim((1e-19, 2e-17))
        plt.ylim((0,1e-14))
        plt.show()

if __name__ == "__main__":
    
    atom = "o"
    seq = "be" #isoelectronic sequence
    shell="2-2" #a-b, recombining from shell a to shell b
    ion = State(atom, seq, shell)

    nist_cutoff=0.05 #fractional deviation of energies from NIST values allowed 
    prior_shape="gaussian"
    likelihood_shape="gaussian"
    
    c = np.array([2.02e-7, 1.51e-5, 5.2e-5, 3.46e-4, 6.77e-4, 4.09e-3, 6.76e-3])
    E = np.array([0.01, 0.06, 0.54, 2.24, 6.93, 15.67, 21.28])
    
    graph_rates_from_file(ion, nist_cutoff=nist_cutoff, prior_shape=prior_shape, likelihood_shape=likelihood_shape, x_range=(1e3, 1e7), 
                          y_range=(1e-13, 1e-9), err_type="max", animate=False, experimental=True)