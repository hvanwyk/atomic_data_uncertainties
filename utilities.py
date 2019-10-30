#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:31:24 2019

@author: kyle
"""

import os
import numpy as np


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
    elif  method == "shift":
        direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/shift_method/"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        return direc
    elif method == "combined":
        direc = f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}/combined_method/"
        if not os.access(direc, os.F_OK):
            os.mkdir(direc)
        return direc
    
def get_nist_energy(fname):
    nist = np.transpose(np.genfromtxt(fname, skip_header=1))
    E_nist = nist[-1]
    return E_nist

def read_levels(levels_file):
    data = np.transpose(np.genfromtxt(levels_file, skip_header=1))
    y = data[-1][:len(data[-1])-1]
    ground = data[-1][-1]
    
    return y, ground

def compare_to_nist(y_comp, y_nist):
    return ((y_comp - y_nist) / (1 + y_nist))