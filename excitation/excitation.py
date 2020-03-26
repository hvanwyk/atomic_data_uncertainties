#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:22:30 2019

@author: kyle
"""

import os
import sys
if ".." not in sys.path:
    sys.path.append("..")
from state import State
            
def create_directories(ion):
    
    if not os.access(f"isoelectronic", os.F_OK):
        os.mkdir("isoelectronic")
        
    if not os.access(f"isoelectronic/{ion.isoelec_seq}-like", os.F_OK):
        os.mkdir(f"isoelectronic/{ion.isoelec_seq}-like")
    
    if not os.access(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}", os.F_OK):
        os.mkdir(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}")
    
    direc = f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/"
    
    return direc    

def orbitals(ion, nmax): 
    #Generate list of orbitals through shell "nmax"
    orb = []
    
    with open(f"configs/{ion.isoelec_seq}_n_{nmax}.config", "r") as config:
        lines = config.readlines()
        for line in lines:
            line = line.split(" ")
            for term in line:
                x = term.rstrip()[:2]
                if x not in orb:
                    orb.append(x)
    return orb

def gen_input(ion, lambdas, nmax=3, max_ex=30, max_nx=70, maxc=50):
    #Create input file for R-matrix run
    direc = create_directories(ion)
    mesh_fine = 0.000325
    mesh_coarse = 0.01
    maxe_ionpot = 4
    rdamp = 0
    adamp = 0
    accel = 0
    
    with open(direc + "input.dat", "w") as file:
        file.write("GENERAL\n")
        file.write(f"2Jmax_ex = {max_ex}\n")
        file.write(f"2Jmax_nx = {max_nx}\n")
        file.write(f"maxc = {maxc}\n")
        file.write(f"mesh_fine = {mesh_fine}\n")
        file.write(f"mesh_coarse = {mesh_coarse}\n")
        file.write(f"rdamp = {rdamp}\n")
        file.write(f"adamp = {adamp}\n")
        file.write(f"accel = {accel}\n")
        
        file.write("\n")
        file.write("CONFIGURATION LIST\n")
        with open(f"configs/{ion.isoelec_seq}_n_{nmax}.config", "r") as config:
            lines = config.readlines()
            for line in lines:
                file.write(line)
        file.write("\n")
        file.write("SCALING PARAMETERS\n")
        
        orbs = orbitals(ion, nmax)
        for i, orb in enumerate(orbs):
            file.write(f"{orb} = {lambdas[i]}\n")

def run_r_matrix(ion, lambdas, nmax=3, max_ex=30, max_nx=70, maxc=50, potential_type=1, born_only=False):
    direc = create_directories(ion)
    gen_input(ion, lambdas, nmax=nmax, max_ex=max_ex, max_nx=max_nx, maxc=maxc)
    if "pp" not in os.listdir(direc):
        os.system("cp ../../r_matrix/bin/parallel_procfile " + direc+"pp")
    if "adas803.pl" not in os.listdir(direc):
        os.system("cp ../../r_matrix/adas803.pl " + direc+"adas803.pl")
    os.chdir(direc)
    
    #delete all the subdirectories before doing new run
    os.system(f"./adas803.pl --delete") 

    #only compute structure and A-values
    if born_only:
        inp_string = f"./adas803.pl --proc=pp --inp input.dat {ion.nuclear_charge}"
        born_string = f"./adas803.pl --proc=pp --born input.dat {ion.nuclear_charge}"
        if potential_type==-1:
            inp_string += " --tf_potential"
            born_string += " --tf_potential"
        os.system(inp_string) 
        os.system(born_string)
        
    #full R-matrix calculation
    else:
        command_line = f"./adas803.pl --proc=pp input.dat {ion.nuclear_charge}"
        if potential_type==-1:
            command_line += " --tf_potential"
        os.system(command_line) 
    os.chdir("../../../")

if __name__ == "__main__":
    
    seq = "be"
    atom = "o"
    
    ion = State(atom, seq)
    
    max_ex = 30
    max_nx = 70
    maxc = 50
    
    
    nmax=3

    orbs = orbitals(ion, nmax)
    lambdas = [1.0]*len(orbs)

    direc = create_directories(ion)
    run_r_matrix(ion, lambdas=lambdas, nmax=nmax, max_ex=max_ex, max_nx=max_nx, maxc=maxc, potential_type=-1)
    