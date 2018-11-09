#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 12:01:25 2018

@author: kyle
"""

import os
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
       

# Recombining from 2 shell to 2 shell
shell = "2-2"

# Class used to hold information about recombining ion
class State:
   
    nuclear_charge = 2
    configuration = "" 
    species = "he"
    ionized = True
    isoelec_seq = "he"

    def __init__(self, species, isoelec_seq):
        self.species = species
        self.isoelec_seq = isoelec_seq
        with open("../adf00/{}.dat".format(self.isoelec_seq)) as seq_file:
            lines = seq_file.read().splitlines()
            self.seq_num_electrons = abs(int(lines[0].split()[1]))
            line = lines[1].split()
            for i in range(len(line)):
                if line[i] == "(":
                    self.seq_config = line[i-1][:2]
                    break
        with open("../adf00/{}.dat".format(self.species.lower()), "r") as adf00:
            lines = adf00.read().splitlines()
            self.nuclear_charge = abs(int(lines[0].split()[1]))
            self.ion_charge = self.nuclear_charge - self.seq_num_electrons
            state = lines[self.ion_charge + 1].split()
            for i in range(len(state)):
                if (state[i] == "("):
                    self.configuration = '  '.join(state[2:i]).split()
                    break
                elif (i == len(state) - 1):
                    self.configuration =  '  '.join(state[2:]).split()
            self.nist_ground = 0
            for i in range(self.ion_charge+1, self.nuclear_charge+1):
                self.nist_ground += float(lines[i].split()[1].replace('d', 'e'))
            
            self.nist_ground /= -13.605698 

"""
Run the autostructure code and generate LEVELS file with energies.
seq - base atom of the isoelectronic sequence
atom - nucleus of interest
lambdas - array of lambda parameters
"""

def structure(seq, atom, lambdas, absolute_ground=False):
    initial = State(atom, seq)
    
    # checks if directory exists, creates it if not
    if not os.access("results/isoelectronic/{}".format(seq), os.F_OK):
        os.mkdir("results/isoelectronic/{}".format(seq))
    direc = "results/isoelectronic/{}/{}{}".format(seq, atom, initial.ion_charge)
    if not os.access(direc, os.F_OK):
        os.mkdir(direc)
       
    asdeck_file = "{}{}_das_{}_str".format(atom,initial.ion_charge, shell)
    os.system("cp asdeck/structure/{0}-like_str {1}/{2}".format(
            seq, direc, asdeck_file))
    os.chdir(direc)
    
    with open(asdeck_file, "a+") as asdeckin:
        NZION = -initial.nuclear_charge
        PRINT = 'FORM'
        asdeckin.write(" &SMINIM  NZION={} NLAM={} PRINT='{}' &END\n".format(NZION, len(lambdas), PRINT))
        lam = [str(lambd) for lambd in lambdas]
        asdeckin.write("  " + ' '.join(lam) + "\n")
        MENG = -15
        EMIN = 0
        EMAX = 2
        asdeckin.write(" &SRADCON  MENG={} EMIN={} EMAX={} &END\n\n".format(MENG, EMIN, EMAX))

    os.system("./../../../../asdeck.x < " + asdeck_file)
    
    
    levels = np.transpose(np.genfromtxt("LEVELS", skip_header=1))
    y = levels[-1][:len(levels[-1])-1]
    ground = levels[-1][-1]
            
    nist = np.transpose(np.genfromtxt("../../../../NIST/isoelectronic/{}/{}{}.nist".format(seq, atom, initial.ion_charge),
                                      skip_header=1))
    y_nist = nist[-1]
    
    if absolute_ground == True:
        y += ground
        y_nist += initial.nist_ground
    
    shift = (y - y_nist) / (1 + y_nist)  

    os.chdir("../../../..")
    
    return y, y_nist, shift
    
# Run autostructure and postprocessing, shifting to NIST energies
def dielectronic_recomb(seq, atom, lambdas=[], xsec=False):
            
    y, y_nist, shift = structure(seq, atom, lambdas)
    initial = State(atom, seq)
    ion = "{}{}".format(atom, initial.ion_charge)

    if not os.access("results/isoelectronic/{}".format(seq), os.F_OK):
        os.mkdir("results/{}".format(atom))
    direc = "results/isoelectronic/{}/{}".format(seq, ion)
    if not os.access(direc, os.F_OK):
        os.mkdir(direc)
    
    asdeck_file = "{}_das_{}_n".format(ion, shell)
    os.system("cp asdeck/dr/{0}-like.dr {1}/{2}".format(seq, direc, asdeck_file))
    os.chdir(direc)
    
    with open("adasin", "w") as adasin:
        with open("LEVELS", "r") as levels:
            lines = levels.read().splitlines()
            NTAR1=1
            #NTAR2=len(lines)-2
            NTAR2 = len(y)
            NECOR=NTAR2
            adasin.write("/IC/\n")
            adasin.write(" &ONE NTAR1={} NTAR2={} COREX=\'{}\' &END\n".format(NTAR1, NTAR2, shell))
            adasin.write(" &TWO NECOR={} ".format(NECOR))
            
            if xsec == True:
                EWIDTH = 0.01
                NBIN = 1000
                EMIN = 0.0
                EMAX = 1.6
                adasin.write("EWIDTH={} NBIN={} EMIN={} EMAX={} ".format(EWIDTH, NBIN, EMIN, EMAX))
                
            adasin.write("&END\n")
            
            for i in range(1, len(lines)-1):
                line = lines[i].split()
                adasin.write(" " + " ".join(line[0:4]) + "\n")
           
            y_str = [str(x) for x in y]
            nist_str = [str(x) for x in y_nist]
            adasin.write(" ".join(y_str))
            adasin.write("\n")
            adasin.write(" ".join(nist_str))
    
    with open(asdeck_file, "a") as asdeckin:       
        """
        Write the DRR namelist.
        NMIN, NMAX - n shells used
        LMIN, LMAX - l orbitals used per n shell
        """
        NMIN = 3
        NMAX = 3
        LMIN = 0
        LMAX = 1
        LCON = 5
        asdeckin.write(" &DRR    NMIN={} NMAX={} LMIN={} LMAX={} LCON={} &END\n".format(NMIN, NMAX, LMIN, LMAX, LCON))
       
        """
        Write the SMINIM namelist, including lambda parameters
        NZION - Z of recombining atom
        PRINT - output file formatting
        NLAM - # of lambda parameters used
        """
        NZION = -initial.nuclear_charge
        PRINT = 'FORM'
        lam = [str(lambd) for lambd in lambdas]
        asdeckin.write(" &SMINIM  NZION={} NLAM={} PRINT='{}' &END\n".format(NZION, len(lambdas), PRINT))
        asdeckin.write("  " + ' '.join(lam) + "\n")
        
        MENG = -15
        EMIN = 0
        EMAX = 2
        asdeckin.write(" &SRADCON  MENG={} EMIN={} EMAX={} &END\n\n".format(MENG, EMIN, EMAX))
        
        JND = 14
        LMAX = 15
        NMAX = 15
        asdeckin.write("JND={} LMAX={} NMAX={}\n".format(JND, LMAX, NMAX))
        asdeckin.write("    16   20   25   35   45   55   70  100  140  200  300  450  700  999")
        
    os.system("./../../../../asdeck.x < " + asdeck_file)
    
    os.system("cp oic o1")
    os.system("./../../../../adasdr.x < adasin")

    if xsec == True:
        with open("XDRTOT", "r") as xdrtot:
            n_points = sum(1 for line in xdrtot) - 1
            
            E = np.zeros(n_points)
            cross_sec = np.zeros(n_points)
            
            xdrtot.seek(0)
            xdrtot.readline()
            for j, line in enumerate(xdrtot):
                dat = line.split()
                E[j] = float(dat[0])
                cross_sec[j] = float(dat[1])
       
        os.chdir("../../../..")

        return(E, cross_sec)
        
    else:
        with open("adasout", "r") as f:
            string = f.read()
            f.seek(string.find("T(K)"))
            f.readline()
            f.readline()
            
            start = f.tell()
            count = 0
            line = f.readline()
            while(line[0] != "C"):
                count += 1
                line = f.readline()
                
            
            T = np.zeros(count)
            rate = np.zeros(count)
            
            f.seek(start)
            for i in range(count):
                line = f.readline().split()
                try:
                    T[i] = float(line[0])
                    rate[i] = float(line[1])
                except:
                    pass
        
        os.chdir("../../../..")
        return(T, rate)
"""  
Reads in lambda parameters from a {}_lambdas.dat file and runs structure and 
dielectronic_recomb for given isoelectronic sequence and nucleus for each
set of lambdas. Returns temperature array and 2D array of rates for each 
T value and lambda set. Writes to output file if write_data is True.
"""


def get_recombination_rates(seq, atom, xsec=False):
    
    recombining = State(atom, seq)
    ion = "{}{}".format(atom, recombining.ion_charge)
    
    lambdas = np.genfromtxt("results/isoelectronic/{}/{}/{}_lambdas.dat".format(seq, ion, ion), skip_header=1)
    
    n_reps, n_lam = np.shape(lambdas)
    
    if xsec == False:
        T0, rate0= dielectronic_recomb(seq, atom, [])
        n_points = len(T0)
        rates = np.zeros((n_reps, n_points))
        col_width = "15"
        precision = "6"
        
        with open("results/isoelectronic/{}/{}/recomb_rates.dat".format(seq, ion), "w") as output:
    
            output.write("{}  {}  {}  {}\n".format(atom, recombining.nuclear_charge, recombining.ion_charge, n_reps))
            output.write("{string:>{width}}".format(string="T (K)", width=col_width))
            
            for i in range(n_reps):
                
                lam = [str(x) for x in lambdas[i]]
                output.write("{string:>{width}}".format(string=str(lam), width=col_width))
                rates[i,:] = dielectronic_recomb(seq, atom, lam)[1]
            
            output.write("\n")
            for i in range(n_points):
                output.write("{num:{width}.{dec}}".format(num=T0[i], width=col_width, dec=precision))
                for j in range(n_reps):
                    output.write("{num:{width}.{dec}}".format(num=rates[j, i], width=col_width, dec=precision))
                output.write("\n")
            
        return (T0, rates)
    
    else:
        structure(seq, atom, [])
        E0, cross_sec0 = dielectronic_recomb(seq, atom, [], xsec=True)
        n_points = len(E0)
        cross_sec = np.zeros((n_reps, n_points))
        
        col_width = "15"
        precision = "6"
        
        with open("results/isoelectronic/{}/{}/cross_sections.dat".format(seq, ion), "w") as output:
    
            output.write("{}  {}  {}  {}\n".format(atom, recombining.nuclear_charge, recombining.ion_charge, n_reps))
            output.write("{string:>{width}}".format(string="E (Ryd)", width=col_width))
            
            for i in range(n_reps):
                
                lam = [str(x) for x in lambdas[i]]
                output.write("{string:>{width}}".format(string=str(lam), width=col_width))
                structure(seq, atom, lam)
                cross_sec[i,:] = dielectronic_recomb(seq, atom, lam, xsec=True)[1]
            
            output.write("\n")
            for i in range(n_points):
                output.write("{num:{width}.{dec}}".format(num=E0[i], width=col_width, dec=precision))
                for j in range(n_reps):
                    output.write("{num:{width}.{dec}}".format(num=cross_sec[j, i], width=col_width, dec=precision))
                output.write("\n")
            
        return (E0, cross_sec)

def graph_from_file(seq, atom, T_min=100, T_max=1000000, xsec=False):
    species = State(atom, seq)
    ion = "{}{}".format(atom, species.ion_charge)
    
    if xsec == True:
        with open("results/isoelectronic/{0}/{1}/cross_sections.dat".format(seq, ion, ion), "r") as data:
            data.readline()
            data.readline()
            n_reps = len(data.readline().split())-1
            
            data.seek(0)
            n_points = sum(1 for line in data) - 2
            E = np.zeros(n_points)
            cross_sections = np.zeros((n_reps, n_points))
            
            data.seek(0)
            data.readline()
            data.readline()
            for j, line in enumerate(data):
                nums = line.split()
                E[j] = float(nums[0])
                for i in range(n_reps):
                    cross_sections[i, j] = float(nums[i+1])
                    
            avg = np.mean(cross_sections, axis=0)  
            err = np.std(cross_sections, axis=0)
            
            err /= avg
            
        fig, ax = plt.subplots(2, 1)
    
        for i in range(n_reps):
            ax[0].plot(E, cross_sections[i, :])
        ax[0].set_title("DR Cross Section of " + ion)
        ax[0].set_xlabel("Energy (Ryd)")
        ax[0].set_ylabel("Total DR X-section")
        #ax[0].set_xscale("log")
        #ax[0].set_yscale("log")
        
        ax[1].plot(E, err)
        ax[1].set_title("Std. Dev.")
        ax[1].set_xlabel("Energy (Ryd)")
        ax[1].set_ylabel("Std. Dev.")    
        #ax[1].set_xscale("log")
        #ax[1].set_yscale("log")

        """
        ax[2].errorbar(x=E, y=avg, yerr=err)
        ax[2].set_title("Average DR Cross Section")
        ax[2].set_xlabel("Energy (Ryd)")
        ax[2].set_ylabel("Total DR X-section")
        ax[2].set_xscale("log")
        ax[2].set_yscale("log")
        """
        
        fig.tight_layout()
        fig.savefig("results/isoelectronic/{0}/{1}/{2} Cross Section Graph".format(seq, ion, ion))
        plt.close(fig)
    else:
        
        with open("results/isoelectronic/{0}/{1}/recomb_rates.dat".format(seq, ion, ion), "r") as data:
            lines = data.read().splitlines()
            min_line = 2
            max_line = len(lines)
            for i in range(2, len(lines)):
                nums = lines[i].split()
                if (float(nums[0]) <= T_min):
                    min_line = i
                if (float(nums[0]) >= T_max):
                    max_line = i-1
                    break
                
            n_points = max_line - min_line + 1
            n_reps = len(lines[2].split()) - 1
            T = np.zeros(n_points)
            rates = np.zeros((n_reps, n_points))
            
            for j in range(n_points):
                nums = lines[min_line + j].split()
                T[j] = float(nums[0])
                for i in range(n_reps):
                    rates[i, j] = float(nums[i+1])
                    
            avg = np.mean(rates, axis=0)  
            err = np.std(rates, axis=0)
            
            frac_err = err / avg
            
        
        fig, ax = plt.subplots(2, 1)
    
        for i in range(n_reps):
            ax[0].plot(T, rates[i, :])
        ax[0].set_title("Dielectronic Recombination of " + ion)
        ax[0].set_xlabel("Temperature (K)")
        ax[0].set_ylabel("DR Rate (cm^3 s^-1)")
        ax[0].set_xscale("log")
        ax[0].set_yscale("log")
        
        ax[1].plot(T, frac_err)
        ax[1].set_title("Fractional Standard Devitation")
        ax[1].set_xlabel("Temperature (K)")
        ax[1].set_ylabel("Std. Dev. / Rate")    
        ax[1].set_xscale("log")
        ax[1].set_yscale("log")

        
        fig.tight_layout()
        fig.savefig("results/isoelectronic/{0}/{1}/{2} DR".format(seq, ion, ion))
        plt.close(fig)
        
        fig, ax = plt.subplots(2, 1)
    
        for i in range(n_reps):
            ax[0].plot(T, rates[i, :])
        ax[0].set_title("Dielectronic Recombination of " + ion)
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

        
        fig.tight_layout()
        fig.savefig("results/isoelectronic/{0}/{1}/{2} DR- Absolute Error".format(seq, ion, ion))
        plt.close(fig)
        
        plt.figure()
        plt.errorbar(x=T, y=avg, yerr=err)
        plt.title("Average DR Rate")
        plt.xlabel("Temperature (K)")
        plt.ylabel("DR Rate (cm^3 s^-1)")
        plt.xscale("log")
        #plt.yscale("log")
        plt.savefig("results/isoelectronic/{0}/{1}/{2} Errorbars".format(seq, ion, ion))
        
        
        

def xsec_graph(seq, atom, lam):
    
    E, cross_sec = dielectronic_recomb(seq, atom, lam, xsec=True)
    state = State(atom, seq)
    ion = "{}{}".format(atom, state.ion_charge)
    
    plt.figure()
    plt.plot(E, cross_sec)  
    plt.xlabel("Energy (Ryd)")
    plt.ylabel("Total DR X-Section (MB)")
    plt.title("Total DR Cross Section of {}{}".format(atom, state.ion_charge))
    
    plt.savefig("results/isoelectronic/{0}/{1}/{2}_xsec_graph".format(seq, ion, ion))
    
def fit_rate_fortran(seq, atom, T, rate):
    state = State(atom, seq)
    os.chdir("results/isoelectronic/{}/{}{}".format(seq, atom, state.ion_charge))
    
    n_metas = 1
    
    header = " {} {} {}\n{}".format(state.nuclear_charge, state.seq_num_electrons, n_metas, len(T))
    np.savetxt("cfin", np.transpose([T, rate]), header=header, comments='')
    
    os.system("./../../../../fits/cfit.x < cfin")
    os.chdir("../../../..")

def dr_fit_func(T, A1, A2, A3, A4, A5, A6, k1, k2, k3, k4, k5, k6):
    y = (A1*np.exp(-k1/T) + A2*np.exp(-k2/T) + A3*np.exp(-k3/T) + A4*np.exp(-k4/T) + A5*np.exp(-k5/T) + A6*np.exp(-k6/T))/(T**1.5)
    return y
  

def fit_rate(T, rate, graph=False):
    popt, popv = curve_fit(dr_fit_func, T, rate, p0=[1e7, 1e7, 1e7, 1e7, 1e7, 1e7, 1e3, 1e4, 1e5, 1e5, 1e5, 1e5],
                           maxfev=100000)
    start = 0
    if graph == True:
        plt.semilogx(T[start:], dr_fit_func(T[start:], *popt))
        plt.scatter(T[start:], rate[start:], c="r")
    return popt

def fits_from_file(seq, atom, graph=False):
    state = State(atom, seq)
    os.chdir("results/isoelectronic/{}/{}{}".format(seq, atom, state.ion_charge))
    
    data = np.transpose(np.genfromtxt("recomb_rates.dat", skip_header=2))
    T = data[0]
    rates = data[1:]
    n_reps = len(rates) // 100
    coeff = np.zeros((n_reps, 12))
    
    plt.figure()
    x = [i for i in range(1,13)]
    for i in range(n_reps):
        try:
            coeff[i,:] = fit_rate(T, rates[i*100,:])
            plt.plot(x, coeff[i,:])
        except:
            pass
    np.savetxt("fit_coeff.dat", np.transpose(coeff))
    os.chdir("../../../..")


