#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 14:46:21 2019

@author: kyle
"""

import numpy as np
import os
import matplotlib.pyplot as plt
from recombination_methods import State, structure, structure_dr, postprocessing_rates
from graphing import graph_xsec

seq = "be"
atom = "fe"
ion = State(atom, seq, "2-2")
ECORIC = 0
NMIN = 3
NMAX = 15
LMIN = 0
LMAX = 7
direc = f"results/isoelectronic/{seq}/{atom}{ion.ion_charge}/"

xsec_file = direc + f"xsec" + ("" if ECORIC==0 else f"_ECORIC{ECORIC}") + (
        "" if (NMIN==3 and NMAX==15) else f"_N_{NMIN}_{NMAX}") + (
        "" if (LMIN==0 and LMAX==7) else f"_L_{LMIN}_{LMAX}") +".png"
        
E, E_nist, E_shift = structure(ion)
structure_dr(ion, NMIN=NMIN, NMAX=NMAX, LMIN=LMIN, LMAX=LMAX, ECORIC=ECORIC) 

energy, xsec = postprocessing_rates(ion, E, E_nist, compute_xsec=True, EWIDTH=0.0001, NBIN=10000)
graph_xsec(ion, energy, xsec, xsec_file)

