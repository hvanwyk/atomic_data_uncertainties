#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  1 14:04:01 2018

@author: kyle
"""

import os

shell = "2-2"

seq = "be"
atom = "c"

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

def structure(seq, atom, lambdas):
    initial = State(atom, seq)
    if not os.access("results/isoelectronic/{}".format(seq), os.F_OK):
        os.mkdir("results/isoelectronic/{}".format(seq))
    direc = "results/isoelectronic/{}/{}{}".format(seq, atom, initial.ion_charge)
    if not os.access(direc, os.F_OK):
        os.mkdir(direc)
    
    os.system("cp asdeck.x {}/asdeck.x".format(direc))
    os.system("cp adasdr.x {}/adasdr.x".format(direc))
    
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
        
    with open("n_conv.dat", "a") as n_conv:
        with open("LEVELS", "r") as levels:
            n_conv.write(levels.read())
            n_conv.write("\n\n")

    os.system("./asdeck.x < " + asdeck_file)

    os.chdir("../../..")
    
structure(seq, atom, [])