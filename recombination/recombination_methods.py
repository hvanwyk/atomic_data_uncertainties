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
       


class State:
    """
    Stores attributes of the initial ion state, before 
    recombination/ionization.
    """
    
    nuclear_charge = 2
    configuration = "" 
    ionized = True

    def __init__(self, species="c", isoelec_seq="li", shell="2-2"):
        self.species = species
        self.isoelec_seq = isoelec_seq
        self.shell = shell
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
        
    
    def structure(self, lambdas=[],PRINT='FORM', MENG=-15, EMIN=0, EMAX=2, E_absolute=False):
        
        """
        Run the autostructure code and generate LEVELS file with energies.
        lambdas - array of lambda parameters
        E_absolute - if True, use absolute energies instead of relative to ground level.
        """
        
        # checks if directory exists, creates it if not
        if not os.access("results/isoelectronic/{}".format(self.isoelec_seq), os.F_OK):
            os.mkdir("results/isoelectronic/{}".format(self.isoelec_seq))
        direc = "results/isoelectronic/{}/{}{}".format(self.isoelec_seq, self.species, self.ion_charge)
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
           
        asdeck_file = "{}{}_das_{}_str".format(self.species, self.ion_charge, self.shell)
        os.system("cp asdeck/structure/{0}-like_str {1}/{2}".format(
                self.isoelec_seq, direc, asdeck_file))
        os.chdir(direc)
        
        with open(asdeck_file, "a+") as asdeckin:
            asdeckin.write(" &SMINIM  NZION={} NLAM={} PRINT='{}' &END\n".format(-self.nuclear_charge, len(lambdas), PRINT))
            lam = [str(lambd) for lambd in lambdas]
            asdeckin.write("  " + ' '.join(lam) + "\n")
            asdeckin.write(" &SRADCON  MENG={} EMIN={} EMAX={} &END\n\n".format(MENG, EMIN, EMAX))
    
        os.system("./../../../../asdeck.x < " + asdeck_file)
        
        
        levels = np.transpose(np.genfromtxt("LEVELS", skip_header=1))
        y = levels[-1][:len(levels[-1])-1]
        ground = levels[-1][-1]
                
        nist = np.transpose(np.genfromtxt("../../../../NIST/isoelectronic/{}/{}{}.nist".format(self.isoelec_seq, self.species, self.ion_charge),
                                          skip_header=1))
        y_nist = nist[-1]
        
        if E_absolute == True:
            y += ground
            y_nist += self.nist_ground
        
        shift = (y - y_nist) / (1 + y_nist)  
    
        os.chdir("../../../..")
        
        return y, y_nist, shift
    
    def dielectronic_recomb(self, lambdas=[], PRINT='FORM', MENG=-15, EMIN=0, EMAX=2, NTAR1=1, NMIN=3, NMAX=3, LMIN=0, LMAX=1, LCON=5,
                            xsec=False, EWIDTH_XSEC=0.01, NBIN_XSEC=1000, EMIN_XSEC=0.0, EMAX_XSEC=2.0, E_absolute=False):
            
        """
        Run DR autostructure and postprocessing, shifting to NIST energies
        """

        y, y_nist, shift = self.structure(lambdas=lambdas, PRINT=PRINT, MENG=MENG, EMIN=EMIN, EMAX=EMAX, E_absolute=E_absolute)
        ion = "{}{}".format(self.species, self.ion_charge)
        
        asdeck_file = "{}_das_{}_n".format(ion, self.shell)
        direc = "results/isoelectronic/{}/{}{}".format(self.isoelec_seq, self.species, self.ion_charge)
        os.system("cp asdeck/dr/{0}-like.dr {1}/{2}".format(self.isoelec_seq, direc, asdeck_file))
        os.chdir(direc)
        
        with open("adasin", "w") as adasin:
            with open("LEVELS", "r") as levels:
                lines = levels.read().splitlines()
                NTAR2 = len(y)
                NECOR=NTAR2
                adasin.write("/IC/\n")
                adasin.write(" &ONE NTAR1={} NTAR2={} COREX=\'{}\' &END\n".format(NTAR1, NTAR2, self.shell))
                adasin.write(" &TWO NECOR={} ".format(NECOR))
                
                if xsec == True:
                    adasin.write("EWIDTH={} NBIN={} EMIN={} EMAX={} ".format(EWIDTH_XSEC, NBIN_XSEC, EMIN_XSEC, EMAX_XSEC))
                    
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
            asdeckin.write(" &DRR    NMIN={} NMAX={} LMIN={} LMAX={} LCON={} &END\n".format(NMIN, NMAX, LMIN, LMAX, LCON))
           
            """
            Write the SMINIM namelist, including lambda parameters
            NZION - Z of recombining atom
            PRINT - output file formatting
            NLAM - # of lambda parameters used
            """
            lam = [str(lambd) for lambd in lambdas]
            asdeckin.write(" &SMINIM  NZION={} NLAM={} PRINT='{}' &END\n".format(-self.nuclear_charge, len(lambdas), PRINT))
            asdeckin.write("  " + ' '.join(lam) + "\n")
            
            asdeckin.write(" &SRADCON  MENG={} EMIN={} EMAX={} &END\n\n".format(MENG, EMIN, EMAX))
            
        os.system("./../../../../asdeck.x < " + asdeck_file)
        
        os.system("cp oic o1")
        os.system("./../../../../adasdr.x < adasin")
    
        if xsec == True:
            data = np.transpose(np.genfromtxt("XDRTOT", skip_header=1))
            E = data[0,:]
            cross_sec = data[1,:]
           
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


def dr_fit_func(T, A1, A2, A3, A4, A5, k1, k2, k3, k4, k5):
    y = (A1*np.exp(-k1/T) + A2*np.exp(-k2/T) + A3*np.exp(-k3/T) + A4*np.exp(-k4/T) + A5*np.exp(-k4/T))/(T**1.5)
    return y
  

def fit_rate(T, rate, graph=False):
    popt, popv = curve_fit(dr_fit_func, T, rate, p0=[1e7, 1e7, 1e7, 1e7, 1e7, 1e3, 1e4, 1e5, 1e5, 1e5],
                           maxfev=100000)
    start = 0
    if graph == True:
        plt.semilogx(T[start:], dr_fit_func(T[start:], *popt))
        plt.scatter(T[start:], rate[start:], c="r")
    return popt



