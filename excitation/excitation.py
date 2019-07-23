#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 15:22:30 2019

@author: kyle
"""

import os

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
            
def create_directories(ion):
    
    if not os.access(f"isoelectronic", os.F_OK):
        os.mkdir("isoelectronic")
        
    if not os.access(f"isoelectronic/{ion.isoelec_seq}-like", os.F_OK):
        os.mkdir(f"isoelectronic/{ion.isoelec_seq}-like")
    
    if not os.access(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}", os.F_OK):
        os.mkdir(f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}")
    
    direc = f"isoelectronic/{ion.isoelec_seq}-like/{ion.species}{ion.ion_charge}/"
    
    return direc

def orbitals(ion):
    direc = create_directories(ion)
    
    orb = []
    
    with open("/".join(direc.split("/")[:-2]) + "/config.dat") as config:
        lines = config.readlines()
        for line in lines:
            line = line.split(" ")
            orb += line.rstrip()
    return orb

def gen_input(ion):
    direc = create_directories(ion)
    
    with open(direc + "input.dat") as file:
        file.write("GENERAL\n")
        file.write(f"2Jmax_ex = {max_ex}\n")
        file.write(f"2Jmax_nx = {max_nx}\n")
        file.write(f"maxc = {maxc}\n")
        file.write(f"mesh_fine = {mesh_fine}\n")
        file.write(f"mesh_coarse = {mesh_coarse}\n")
        file.write(f"rdamp = {rdamp}\n")
        file.write(f"adamp = {adamp}\n")
        file.write(f"accel = {accel}\n")
        
        file.write("")
        file.write("CONFIGURATION LIST\n")
        with open(direc+"../config.dat") as config:
            lines = config.readlines()
            for line in lines:
                file.write(line)

if __name__ == "__main__":
    
    seq = "he"
    atom = "o"
    
    ion = State(atom, seq)
    
    max_ex = 31
    max_nx = 71
    maxc = 30
    mesh_fine = 0.000325
    mesh_coarse = 0.01
    maxe_ionpot = 4
    rdamp = 0
    adamp = 0
    accel = 0
    
    create_directories(ion)
    print(orbitals(ion))