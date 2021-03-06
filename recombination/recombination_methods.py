#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 12:01:25 2018
@author: kyle
"""

import os
import sys
import numpy as np
if ".." not in sys.path:
    sys.path.append("..")
from utilities import create_directories, get_nist_energy, read_levels, compare_to_nist
from state import State
    
def structure(ion, method="lambdas", lambdas=[], potential=1, MENG=-15, EMIN=0, EMAX=2, E_absolute=False):
    """
    Run the autostructure code and generate LEVELS file with energies.
    
    Inputs:
    
        method: str, specify uncertainty propagation method ('lambdas', 'combined') 
        
        lambdas: list, array of lambda parameters
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are Schmidt 
                orthogonalized. 
                
            NZION <0: nl-depdendent Hartree potential evaluated with Slater type 
                orbitals (NOT Schmidt orthogonalized).
                
        MENG: 
        
        EMIN, EMAX: 
        
        E_absolute: bool, if True, use absolute energies instead of relative to 
            ground level.
            
        
    """
    
    # checks if directory exists, creates it if not
    direc = create_directories(ion, method)
    up_dir = "../../../../../"
    
    if method == "lambdas" and lambdas != []:
        direc += "_".join([str(x) for x in lambdas])
        up_dir = "../../../../../../"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        
    asdeck_file = f"{ion.species}{ion.ion_charge}_das_{ion.shell}_str"
    os.system(f"cp asdeck/structure/{ion.isoelec_seq}-like_str {direc}/{asdeck_file}")
    os.chdir(direc)
    
    with open(asdeck_file, "a+") as asdeckin:
        asdeckin.write(f" &SMINIM  NZION={np.sign(potential)*ion.nuclear_charge} NLAM={len(lambdas)} PRINT='FORM' &END\n")
        lam = [str(lambd) for lambd in lambdas]
        asdeckin.write("  " + ' '.join(lam) + "\n")
        asdeckin.write(f" &SRADCON  MENG={MENG} EMIN={EMIN} EMAX={EMAX} &END\n\n")
        
    os.system("./" + up_dir + "asdeck.x < " + asdeck_file)

    df_comp, ground = read_levels("LEVELS") 
    
    if method=="lambdas" or method=="combined":
        df_nist = get_nist_energy(up_dir + f"/NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist")
        y_nist = df_nist["E"].values
    if E_absolute == True:
        df_comp["E"] += ground
        if method == "lambdas" or method=="combined":
            df_nist["E"] += ion.nist_ground
    
    y = df_comp["E"].values
    y_nist = df_nist["E"].values
    
    os.remove("oic")
    os.remove("olg")
    os.remove("ols")
    os.remove("TERMS")

    os.chdir(up_dir)
    
    if method=="lambdas" or method=="combined":
        y_shift = compare_to_nist(df_comp, df_nist)
        return np.array([y]).flatten(), np.array([y_nist]).flatten(), np.array([y_shift]).flatten()
    else:
        return np.array([y]).flatten()

def structure_dr(ion, method="lambdas", lambdas=[], potential=1, NMIN=3, NMAX=15, JND=14, LMIN=0, LMAX=7, 
                 MENG=-15, EMIN=0, EMAX=2, ECORIC=0):
    """
    Structure run for dielectronic radiation
    
    Inputs:
    
        ion: State, 
        
        method: str, 'lambdas', 'combined'
        
        lambdas: list, of lambda tuples 
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are Schmidt 
                orthogonalized. 
                
            NZION <0: nl-depdendent Hartree potential evaluated with Slater type 
                orbitals (NOT Schmidt orthogonalized).
                
        
        NMIN, NMAX: int, DRR NAMELIST, min/max n-value of Rydberg electron.
        
        JND: int, AS NAMELIST,
        
        LMIN, LMAX: int, DRR NAMELIST, min/max l-value of Rydberg electron.
        
        MENG: signed int, MENG is the number of interpolation energies (in Rydbergs). 
            MENG energies follow and the continuum orbitals are calculated at those 
            energies. 

            = 0 (default) uses 0, DE/3, DE, 3*DE (& 8*DE if MAXLT/JT.gt.35/70) 
                where DE=max(TEAPOT, DELTAX) and TEAPOT is the estimated ionization 
                potential and DELTAX is the highest (spectroscopic) term energy.
            > 0 reads MENG final scattered energies following the NAMELIST.
            < 0 internally sets -MENG energies between a user specified range:
        
        EMIN, EMAX: 
        
        ECORIC:  
    """
    direc = create_directories(ion, method)
    up_dir = "../../../../../"
    ion_name = f"{ion.species}{ion.ion_charge}"

    if method == "lambdas" and lambdas != []:
        direc += "_".join([str(x) for x in lambdas])
        up_dir = "../../../../../../"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
    
    asdeck_file = f"{ion_name}_das_{ion.shell}_n"

    os.system(f"cp asdeck/dr/{ion.isoelec_seq}-like.dr {direc}/{asdeck_file}")
    os.chdir(direc)
    
    with open(asdeck_file, "a") as asdeckin:       
        """
        Write the DRR namelist.
        NMIN, NMAX - n shells used
        LMIN, LMAX - l orbitals used per n shell
        """
        asdeckin.write(f" &DRR    NMIN={NMIN} NMAX={NMAX} JND={JND} LMIN={LMIN} LMAX={LMAX} &END\n")
        asdeckin.write("16   20   25   35   45   55   70  100  140  200  300  450  700  999\n")
       
        """
        Write the SMINIM namelist, including lambda parameters
        NZION - Z of recombining atom
        PRINT - output file formatting
        NLAM - # of lambda parameters used
        """
        lam = [str(lambd) for lambd in lambdas]
        asdeckin.write(f" &SMINIM  NZION={np.sign(potential)*ion.nuclear_charge} NLAM={len(lambdas)} PRINT='UNFORM' &END\n")
        asdeckin.write("  " + ' '.join(lam) + "\n")
        
        asdeckin.write(f" &SRADCON  MENG={MENG} EMIN={EMIN} EMAX={EMAX} ")
        if ECORIC != 0:
            asdeckin.write(f"ECORIC={ECORIC} ")
        asdeckin.write("&END\n\n")
    
    os.system("./" + up_dir + "asdeck.x < " + asdeck_file)
    os.system("mv oicu o1u")

    os.chdir(up_dir)
    

def postprocessing_rates(ion, E, E_nist=[], method="lambdas", lambdas=[], shift=[], NTAR1=1, 
                         compute_xsec=False, EWIDTH=0.001, NBIN=1000, EMIN=0.0, EMAX=2.0):
    """
    Write input deck for ADASDR postprocessor and run ADASDR 
    
    Inputs:
    
        ion: State, specify ion
        
        E: double, energy
        
        E_nist: list, 
        
        method: str,
        
        lambdas: list, 
        
        shift: list, 
        
        NTAR1: int, AS NAMELIST input,
        
        compute_xsec: bool, 
        
        EWIDTH: double, AS NAMELIST input,
        
        NBIN: int, AS NAMELIST input,
        
        EMIN, EMAX: double, AS NAMELIST input, 
        
        
    Outputs:
    
        If compute_xsec = True, return 
        
            energy, 
            
            cross-sections
            
        If compute_xsec = False, return
        
            T: temperature, 
            
            rate: DR rate for each point on the temperature grid 
    """
    direc = create_directories(ion, method)
    up_dir = "../../../../../"

    levels_file = "LEVELS"
    
    if method == "lambdas" and lambdas != []:
        direc += "_".join([str(x) for x in lambdas])
        up_dir = "../../../../../../"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
    
    os.chdir(direc)
    
    with open("adasin", "w") as adasin:
        with open(levels_file, "r") as levels:
            lines = levels.read().splitlines()
            NTAR2 = len(E)
            NECOR=NTAR2
            adasin.write("/IC/\n")
            adasin.write(f" &ONE NTAR1={NTAR1} NTAR2={NTAR2} COREX=\'{ion.shell}\' &END\n")
            adasin.write(f" &TWO NECOR={NECOR} ")
            
            if compute_xsec:
                adasin.write(f"EWIDTH={EWIDTH} NBIN={NBIN} EMIN={EMIN} EMAX={EMAX} ")
                
            adasin.write("&END\n")
            
            for i in range(1, len(lines)-1):
                line = lines[i].split()
                adasin.write(" " + " ".join(line[0:4]) + "\n")
            E_str = [str(x) for x in E]
            if method == "lambdas":
                nist_str = [str(x) for x in E_nist]
            else:
                if shift == []:
                    nist_str = [str(x) for x in E]
                else:
                    nist_str = [str(x) for x in E + shift]
            adasin.write(" ".join(E_str))
            adasin.write("\n")
            adasin.write(" ".join(nist_str))
    os.system("./" + up_dir + "adasdr.x < adasin")
    if not compute_xsec:
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
        os.chdir(up_dir)
        return T, rate
    else:
        data = np.transpose(np.genfromtxt("XDRTOT", skip_header=1))
        energy = data[0,:]
        xsec= data[1,:]
        os.chdir(up_dir)
        return energy, xsec
        

    
def get_rate(ion, lambdas):
    E, E_nist, E_shift = structure(ion, lambdas=lambdas)
    structure_dr(ion, lambdas=lambdas)
    return postprocessing_rates(ion, E, E_nist, lambdas=lambdas)[1]
