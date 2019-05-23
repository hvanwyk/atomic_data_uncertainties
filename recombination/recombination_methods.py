#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 10 12:01:25 2018

@author: kyle
"""

import os
import numpy as np
       
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
        with open(f"../adf00/{self.isoelec_seq}.dat", "r") as seq_file:
            lines = seq_file.read().splitlines()
            self.seq_num_electrons = abs(int(lines[0].split()[1]))
            line = lines[1].split()
            for i in range(len(line)):
                if line[i] == "(":
                    self.seq_config = line[i-1][:2]
                    break
        with open(f"../adf00/{self.species.lower()}.dat", "r") as adf00:
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
        
    def __repr__(self):
        formatted = f"{self.isoelec_seq.capitalize()}-like {self.species.capitalize()}:\n"
        formatted += "\t Z: {self.nuclear_charge}\n"
        formatted += "\t Charge: {self.ion_charge}+\n"
        formatted += "\t Ground-state energy (NIST): {self.nist_ground}\n"
        formatted += "\t Core Excitation: {self.shell}\n"
        return formatted
        
    

def create_directories(ion, method="lambdas"):

    if not os.access("results", os.F_OK):
        os.mkdir("results")
        
    if not os.access(f"results/isoelectronic", os.F_OK):
        os.mkdir("results/isoelectronic")
        
    if not os.access(f"results/isoelectronic/{ion.isoelec_seq}", os.F_OK):
        os.mkdir(f"results/isoelectronic/{ion.isoelec_seq}")
    
    if not os.access(f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}", os.F_OK):
        os.mkdir(f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}")
        
    if method == "lambdas":
        direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/lambda_method/"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        return direc
    else:
        direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/shift_method/"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        return direc

def get_nist_energy(ion, main_direc=""):
    fname = main_direc + f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist"
    nist = np.transpose(np.genfromtxt(fname, skip_header=1))
    E_nist = nist[-1]
    return E_nist
    
def structure(ion, method="lambdas", lambdas=[], potential=1, MENG=-15, EMIN=0, EMAX=2, E_absolute=False):
        
    """
    Run the autostructure code and generate LEVELS file with energies.
    lambdas - array of lambda parameters
    E_absolute - if True, use absolute energies instead of relative to ground level.
    """
    
    # checks if directory exists, creates it if not
    direc = create_directories(ion, method)
    up_dir = "../../../../../"
    
    if method == "lambdas" and lambdas != []:
        direc += "_".join([str(x) for x in lambdas])
        up_dir = "../../../../../../"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        
    asdeck_file = "{}{}_das_{}_str".format(ion.species, ion.ion_charge, ion.shell)
    os.system(f"cp asdeck/structure/{ion.isoelec_seq}-like_str {direc}/{asdeck_file}")
    os.chdir(direc)
    
    with open(asdeck_file, "a+") as asdeckin:
        asdeckin.write(f" &SMINIM  NZION={np.sign(potential)*ion.nuclear_charge} NLAM={len(lambdas)} PRINT='FORM' &END\n")
        lam = [str(lambd) for lambd in lambdas]
        asdeckin.write("  " + ' '.join(lam) + "\n")
        asdeckin.write(f" &SRADCON  MENG={MENG} EMIN={EMIN} EMAX={EMAX} &END\n\n")
        
    os.system("./" + up_dir + "asdeck.x < " + asdeck_file)

    levels = np.transpose(np.genfromtxt("LEVELS", skip_header=1))
    y = levels[-1][:len(levels[-1])-1]
    ground = levels[-1][-1]
            
    #nist = np.transpose(np.genfromtxt(up_dir +f"NIST/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist", skip_header=1))
    #y_nist = nist[-1]
    y_nist = get_nist_energy(ion, up_dir)
    if E_absolute == True:
        y += ground
        y_nist += ion.nist_ground
    
    y_shift = (y - y_nist) / (1 + y_nist)  

    os.chdir(up_dir)
    
    return np.array([y]).flatten(), np.array([y_nist]).flatten(), np.array([y_shift]).flatten()

def structure_dr(ion, method="lambdas", lambdas=[], potential=1, NMIN=3, NMAX=15, JND=14, LMIN=0, LMAX=7, MENG=-15, EMIN=0, EMAX=2, ECORIC=0):
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
        asdeckin.write(f" &SMINIM  NZION={np.sign(potential)*ion.nuclear_charge} NLAM={len(lambdas)} PRINT='FORM' &END\n")
        asdeckin.write("  " + ' '.join(lam) + "\n")
        
        asdeckin.write(f" &SRADCON  MENG={MENG} EMIN={EMIN} EMAX={EMAX} ")
        if ECORIC != 0:
            asdeckin.write(f"ECORIC={ECORIC} ")
        asdeckin.write("&END\n\n")
    
    os.system("./" + up_dir + "asdeck.x < " + asdeck_file)
    os.system("cp oic o1")

    os.chdir(up_dir)

def postprocessing_rates(ion, E, E_nist, method="lambdas", lambdas=[], shift=[], NTAR1=1):
    
    direc = create_directories(ion, method)
    up_dir = "../../../../../"

    levels_file = "LEVELS"
    
    if method == "lambdas" and lambdas != []:
        direc += "_".join([str(x) for x in lambdas])
        up_dir = "../../../../../../"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
    """
    elif method == "shift" and shift != []:
        direc += "_".join([str(x) for x in shift])
        up_dir = "../../../../../"
        levels_file = "../LEVELS"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
    """
    
    os.chdir(direc)
    """
    if method == "shift" and shift != []:
        os.system("cp ../o1 o1")
        os.system("cp ../olg olg")
        os.system("cp ../ols ols")
    """
    with open("adasin", "w") as adasin:
        with open(levels_file, "r") as levels:
            lines = levels.read().splitlines()
            NTAR2 = len(E)
            NECOR=NTAR2
            adasin.write("/IC/\n")
            adasin.write(f" &ONE NTAR1={NTAR1} NTAR2={NTAR2} COREX=\'{ion.shell}\' &END\n")
            adasin.write(f" &TWO NECOR={NECOR} ")
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

    
def get_rate(ion, lambdas):
    E, E_nist, E_shift = structure(ion, lambdas=lambdas)
    structure_dr(ion, lambdas=lambdas)
    return postprocessing_rates(ion, E, E_nist, lambdas=lambdas)[1]
