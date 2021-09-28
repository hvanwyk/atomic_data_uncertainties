#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:31:24 2019

@author: kyle
"""

import os
import numpy as np
import pandas as pd
from pathlib import Path

def root_directory():
    """
    Returns the name of the project's root directory as a string
    """
    return Path(__file__).parent.parent
    
def create_directories(ion, method="lambdas"):
    """
    Create Directory based on the ion, and method in the 'results' folder 
    and returns directory location.General format:
    
        results/isoelectronic/isoelectronic-sequence/ion-species + charge/method
        
    Inputs:
    
        ion: State, used to specify iso-electronic sequence, species and charge
        
        method: str, 
        
            'lambdas' - vary orbital scaling parameters/lambdas  
            'shift' - use energy shifts to get the right location  
            'combined' - vary orbitals and use energy shifts
            
        
    Outputs: 
    
        direc: str, name of directory
    """
    
    # Base directory 
    if not os.access("results", os.F_OK):
        os.mkdir("results")
    
    # Base subdirectory
    if not os.access(f"results/isoelectronic", os.F_OK):
        os.mkdir("results/isoelectronic")
        
    # Subdirectory for isoelectronic sequence
    if not os.access(f"results/isoelectronic/{ion.isoelec_seq}", os.F_OK):
        os.mkdir(f"results/isoelectronic/{ion.isoelec_seq}")
    
    # Subdirectory for species and charge
    if not os.access(f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}", os.F_OK):
        os.mkdir(f"results/isoelectronic/{ion.isoelec_seq}/{ion.species}{ion.ion_charge}")
    
    # Subdirectory for methods
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
    """
    Look up NIST energies from file
    
    Inputs:
    
        fname: str, name of file with NIST energy
        
    
    Outputs:
    
        df: Pandas.DataFrame, NIST energies (Ryd) of each level, ordered by 
        configuration, L, S, J
    """
    nist = np.transpose(np.genfromtxt(fname, skip_header=1))
    df = pd.DataFrame(np.abs(nist.T), columns=['2J', 'P', 'S', 'L', 'CF', 'NI', 'E'])
    df = df.sort_values(["CF", "L", "S", "2J"])
    return df


def read_levels(levels_file):
    """
    Read a LEVELS file generated by Autostructure (AS)
    
    Inputs: 
    
        levels_file: str, name of levels file
        
    Output:
    
        df: Pandas.DataFrame, energies (Ryd) of each level, ordered by 
        configuration, L, S, J
        
        ground: energy (Ryd) of ground state computed by Autostructure
    """
    data = np.transpose(np.genfromtxt(levels_file, skip_header=1))
    df = pd.DataFrame(np.abs(data.T[:-1,:]), columns=['2J', 'P', 'S', 'L', 'CF', 'NI', 'E'])
    
    print(df)
    
    df = df.sort_values(["CF", "L", "S", "2J"])
    ground = data[-1][-1]
    
    return df, ground


def compare_to_nist(df_comp, df_nist):
    """
    Compute signed relative error based on NIST value.
    
    Inputs:
        
        df_comp: Pandas.DataFrame, computed energies of each level, ordered by 
        configuration, L, S, J
        
        df_nist: Pandas.DataFrame, NIST energies of each level, ordered by 
        configuration, L, S, J
        
    
    Outputs:
    
        err: double, relative error
    """
    y_comp = df_comp["E"].values
    y_nist = df_nist["E"].values
    
    
    err = ((y_comp - y_nist) / (1+y_nist))
    err[err==0]=1e-30
    return err

if __name__ == "__main__":
    df_nist = get_nist_energy("excitation/NIST/isoelectronic/he/he0_n=3.nist")
    df_comp, ground = read_levels("excitation/LEVELS_backup")
    
    print(df_comp)