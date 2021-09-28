#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 17:15:32 2019

@author: kyle
"""
import utilities

class State:
    """
    Stores attributes of the initial ion state, before 
    recombination/ionization.
    """
    
    nuclear_charge = 2
    configuration = "" 
    ionized = True

    def __init__(self, species="c", isoelec_seq="li", shell="2-2"):
        """
        Constructor
        
        Inputs:
            
            species: str, the atom
            
            isolec_seq: str, ion stage
            
            shell: str, shell considered
        
        """
        
        self.species = species
        self.isoelec_seq = isoelec_seq
        self.shell = shell
        
        # Get root directory
        rootdir= utilities.root_directory()
        
        #
        # Get Sequence-level information from adf00 file
        # 
        
        # Set sequence path
        seq_path = f'{rootdir}/adf00/{self.isoelec_seq}.dat'
        print(seq_path)
        with open(seq_path, "r") as seq_file:
            lines = seq_file.read().splitlines()
            self.seq_num_electrons = abs(int(lines[0].split()[1]))
            line = lines[1].split()
            for i in range(len(line)):
                if line[i] == "(":
                    self.seq_config = line[i-1][:2]
                    break
                
        #
        # Get species-level information from adf00 file
        #
        species_path = f'{rootdir}/adf00/{self.species.lower()}.dat'
        with open(species_path,'r') as adf00: 
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
        formatted += f"\t Z: {self.nuclear_charge}\n"
        formatted += f"\t Charge: {self.ion_charge}+\n"
        formatted += f"\t Ground-state energy (NIST): {self.nist_ground}\n"
        formatted += f"\t Core Excitation: {self.shell}\n"
        return formatted