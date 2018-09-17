#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 12:01:25 2018

@author: kyle
"""

import os
import numpy as np
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

"""
Run the autostructure code and generate LEVELS file with energies.
seq - base atom of the isoelectronic sequence
atom - nucleus of interest
lambdas - array of lambda parameters
"""

def structure(seq, atom, lambdas):
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
        asdeckin.write("  " + ' '.join(lambdas) + "\n")
    
        MENG = -15
        EMIN = 0
        EMAX = 2
        asdeckin.write(" &SRADCON  MENG={} EMIN={} EMAX={} &END\n\n".format(MENG, EMIN, EMAX))

    os.system("./../../../../asdeck.x < " + asdeck_file)

    os.chdir("../../../..")
    
# Run autostructure and postprocessing, shifting to NIST energies
def dielectronic_recomb(seq, atom, lambdas):
            
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
            NTAR2=len(lines)-2
            NECOR=NTAR2
            adasin.write("/IC/\n")
            adasin.write(" &ONE NTAR1={} NTAR2={} COREX=\'{}\' &END\n".format(NTAR1, NTAR2, shell))
            adasin.write(" &TWO NECOR={} &END\n".format(NECOR))
            for i in range(1, len(lines)-1):
                line = lines[i].split()
                adasin.write(" " + line[0] + " " + line[1] + "  " + line[2] + " " + line[3] + "\n")
            for i in range(1, len(lines)-1):
                line = lines[i].split()
                adasin.write(line[len(line)-1] + " ")
            adasin.write("\n")
        with open("../../../../NIST/isoelectronic/{}/{}.nist".format(seq, ion)) as nist:
            lines = nist.read().splitlines()
            for i in range(1, len(lines)):
                line = lines[i].split()
                adasin.write(line[len(line)-1] + " ")

    
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
        asdeckin.write(" &SMINIM  NZION={} NLAM={} PRINT='{}' &END\n".format(NZION, len(lambdas), PRINT))
        asdeckin.write("  " + ' '.join(lambdas) + "\n")
        
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
    
    os.system("mv oic o1")
    os.system("./../../../../adasdr.x < adasin")

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


def get_recombination_rates(seq, atom, write_data):
    
    recombining = State(atom, seq)
    ion = "{}{}".format(atom, recombining.ion_charge)
    
    structure(seq, atom, [])
    T0, rate0= dielectronic_recomb(seq, atom, [])
        
    with open("results/isoelectronic/{}/{}/{}_lambdas.dat".format(seq, ion, ion), "r") as lambda_file:
        lines = lambda_file.read().splitlines()
        num_lam = len(lines[2].split()) // 2
        lambdas = np.zeros((len(lines)-2, num_lam))
        for i in range(2, len(lines)):
            for j in range(num_lam):
                lambdas[i-2, j] = lines[i].split()[2*j] 
    
    n_reps = len(lambdas[:,0])
    n_points = len(T0)
    rates = np.zeros((n_reps, n_points))
    col_width = "15"
    precision = "6"
    
    if (write_data == True):
        with open("results/isoelectronic/{}/{}/recomb_rates.dat".format(seq, ion), "w") as output:
    
            output.write("{}  {}  {}  {}\n".format(atom, recombining.nuclear_charge, recombining.ion_charge, n_reps))
            output.write("{string:>{width}}".format(string="T (K)", width=col_width))
            
            for i in range(n_reps):
                
                lam = [str(x) for x in lambdas[i]]
                output.write("{string:>{width}}".format(string=str(lam), width=col_width))
                structure(seq, atom, lam)
                T, rates[i,:] = dielectronic_recomb(seq, atom, lam)
            
            output.write("\n")
            for i in range(n_points):
                output.write("{num:{width}.{dec}}".format(num=T[i], width=col_width, dec=precision))
                for j in range(n_reps):
                    output.write("{num:{width}.{dec}}".format(num=rates[j, i], width=col_width, dec=precision))
                output.write("\n")
        
    return (T, rates)

def graphing(seq, atom, T_min, T_max):
    species = State(atom, seq)
    ion = "{}{}".format(atom, species.ion_charge)
    
    with open("results/isoelectronic/{0}/{1}/recomb_rates.dat".format(seq, ion, ion), "r") as data:
        lines = data.read().splitlines()
        
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
        
    
    fig, ax = plt.subplots(3, 1)

    for i in range(n_reps):
        ax[0].plot(T, rates[i, :])
    ax[0].set_title("Dielectric Recombination of " + ion)
    ax[0].set_xlabel("Temperature (K)")
    ax[0].set_ylabel("Rate")
    ax[0].set_xscale("log")
    ax[0].set_yscale("log")
    
    ax[1].errorbar(x=T, y=avg, yerr=err)
    ax[1].set_title("Average Recombination Rate")
    ax[1].set_xlabel("Temperature (K)")
    ax[1].set_ylabel("Recombination Rate")
    ax[1].set_xscale("log")
    ax[1].set_yscale("log")
    
    ax[2].plot(T, err)
    ax[2].set_title("Std. Dev.")
    ax[2].set_xlabel("Temperature (K)")
    ax[2].set_ylabel("Std. Dev.")    
    ax[2].set_xscale("log")
    ax[2].set_yscale("log")

    fig.tight_layout()
    fig.savefig("results/isoelectronic/{0}/{1}/{2} Recombination Graph".format(seq, ion, ion))
    plt.close(fig)
        
seq = "li"
atom = "c"
#get_recombination_rates(seq, atom, True)
#graphing(seq, atom, 1000, 1000000)
