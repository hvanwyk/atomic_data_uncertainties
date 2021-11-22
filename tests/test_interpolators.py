#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 15 10:59:57 2019

@author: kyle
"""
# Internal modules
from bayesian_methods import lambdas_grid, interpolators
from utilities import root_folder, get_nist_energy

# External modules
import os
import numpy as np
import Tasmanian
from scipy.interpolate import RegularGridInterpolator
from time import time
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import axes3d


from recombination_methods import State, read_levels, get_rate
from recombination_methods import structure, structure_dr, postprocessing_rates

def make_structure_input_deck(folder, lmd, ion, potential=1, MENG=-15, EMIN=0,
                              EMAX=2, debug=False):
    """
    Writes AUTOSTRUCTURE input deck for structure run. 
    
    Inputs:
        
        folder: str, folder in which to store AUTOSTRUCTURE input deck
        
        lmd: tuple/list, (n,) size list of orbital scaling parameters
        
        ion: State, containing information on ion state
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are 
                Schmidt orthogonalized. 
            
            NZION <0: nl-depdendent Hartree potential evaluated with Slater 
                type orbitals (NOT Schmidt orthogonalized).
                
        MENG: signed int, MENG is the number of interpolation energies (in
            Rydbergs). MENG energies follow and the continuum orbitals are 
            calculated at those energies. 
    
            = 0 (default) uses 0, DE/3, DE, 3*DE (& 8*DE if MAXLT/JT.gt.35/70) 
                where DE=max(TEAPOT, DELTAX) and TEAPOT is the estimated ionization 
                potential and DELTAX is the highest (spectroscopic) term energy.
            > 0 reads MENG final scattered energies following the NAMELIST.
            < 0 internally sets -MENG energies between a user specified range 
            
        EMIN, EMAX: double, minimum and maximum energies
                
                
    Outputs:
        
        asdeck_file: str, path to AUTOSTRUCTURE input deck file 
        
        
    Modified: 10/05/2021 @hvanwyk
    """    
    #
    # Copy Autostructure structure input deck template to current folder
    # 
    
    # Name of input file
    asdeck_file = f"{ion.species}{ion.ion_charge}_das_{ion.shell}_str"
    
    # Location of template
    asdeck_path = root_folder()/'recombination'/'asdeck'/'structure'
    
    # Copy file template to local folder
    os.system(f"cp {asdeck_path}/{ion.isoelec_seq}-like_str " +
              f"{folder}/{asdeck_file}")
    
    
    os.chdir(folder)
    
    #
    # Lambda: Add add closed shell orbital 1s at the beginning
    # 
    lmd = lmd.ravel()  # ensure lambda is a vector
    lmd = np.insert(lmd, 0, 1.0)  # add 1.0 in first entry
    
    #
    # Modify input deck to add information on lambda and the potential
    # 
    with open(asdeck_file, "a+") as asdeckin:
        if (potential == 1):
            #
            # Thomas-Fermi DA potential
            # 
            asdeckin.write(" &SMINIM  "+
                           f"NZION={np.sign(potential)*ion.nuclear_charge} "+
                           f"NLAM={len(lmd)} PRINT='FORM' &END\n")
        else:
            #
            # Hartree Potential
            # 
            asdeckin.write(" &SMINIM "+
                           f"NZION={np.sign(potential)*ion.nuclear_charge} " +
                           f"NLAM={len(lmd)} ORTHOG='YES' PRINT='FORM' &END\n")
        
        #
        # Lambda Values
        # 
        asdeckin.write("  " + ' '.join([str(l) for l in lmd]) + "\n")
        
        #
        # Add Energy Ranges
        #
        asdeckin.write(f" &SRADCON  MENG={MENG} " + 
                       f"EMIN={EMIN} EMAX={EMAX} &END\n\n")
        
    return asdeck_file
    

def run_structure(lmd, ion, potential=1, MENG=-15, EMIN=0, EMAX=2, E_abs=False,
                  folder=None, keep_folder=False, clean_folder=True, 
                  return_info=False):
    """
    Mapping from orbital scaling parameters to energies
    
    Inputs:
        
        lmd: double, (1,d) numpy array of orbital scaling parameters, where
            d is the dimension of the parameter space. 
        
        ion: State, containing information on ion state
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are 
                Schmidt orthogonalized. 
            
            NZION <0: nl-depdendent Hartree potential evaluated with Slater 
                type orbitals (NOT Schmidt orthogonalized).
                
        MENG: signed int, MENG is the number of interpolation energies (in
            Rydbergs). MENG energies follow and the continuum orbitals are 
            calculated at those energies. 
    
            = 0 (default) uses 0, DE/3, DE, 3*DE (& 8*DE if MAXLT/JT.gt.35/70) 
                where DE=max(TEAPOT, DELTAX) and TEAPOT is the estimated 
                ionization potential and DELTAX is the highest (spectroscopic) 
                term energy.
            > 0 reads MENG final scattered energies following the NAMELIST.
            < 0 internally sets -MENG energies between a user specified range.
            
        EMIN: double, minimum energy
            
        EMAX: double, maximum energy
            
        E_abs: bool, whether to return the absolute energy or relative to ground
            
        folder: str, full path of folder in which to store results
        
        keep_folder: bool, whether to keep the folder or remove it afterwards.
        
        clean_folder: bool, whether to delete oic, olg, ols, and TERMS files
        
        return_info: bool, whether to return the ground energy and the 
            configurations in addition to the energies.
            
    
    Outputs:
        
        E: double, (1,n) array of energies, where n is the number of energies.
            If E_abs = False, return energies relative to ground, 
            If E_abs = True, return absolute energies.
        
        E_ground: double, ground energy relative to which E is measured
        
        E_names: str, list of tuples of the form (2J,P,S,L,CF,NI) specifying the
            configuration at which the given energy level occurs
    """ 
    current_folder = os.getcwd()
    
    #
    # Create results folder
    # 
    if folder is None:
        folder = current_folder + '/as_temp'
        
    if not os.access(folder, os.F_OK):
        os.mkdir(folder)
            
    #
    # Make Autostructure Input Deck
    #
    asdeck_file = \
        make_structure_input_deck(folder, lmd, ion, potential=potential,
                                  MENG=MENG, EMIN=EMIN, EMAX=EMAX)
    #
    # Run AUTOSTRUCTURE
    #
    
    # Path to AUTOSTRUCTURE
    asdeck_path = str(root_folder()/'recombination/')
    os.system(asdeck_path + "/asdeck.x < " + asdeck_file)
    
    #
    # Extract Energies from LEVELS file
    # 
    data = np.genfromtxt('LEVELS', skip_header=1)
    E = data[:-1,-1]
    E_ground = data[-1,-1]
    E_names = [tuple(l) for l in data[:-1,:-1].astype(int)]
    
    #E, Eg = read_levels(f"{folder}/LEVELS") 
    if E_abs:
        E += E_ground
    
    
    #
    # Clean Up
    # 
    if not keep_folder:
        #
        # Remove folder
        #
        os.chdir(current_folder)
        os.system(f"rm -rf {folder}")
    else:
        
        if clean_folder:
            #
            # Clean up results folder
            #
            os.chdir(f"{folder}")
            
            os.remove("oic")
            os.remove("olg")
            os.remove("ols")
            os.remove("TERMS")
        
        # Change back to current folder
        os.chdir(current_folder)
     
    if return_info:
        return E, E_ground, E_names
    else:
        return E[:,np.newaxis].T


def make_asdr_input_deck(folder, lmd, ion, potential=1, NMIN=3, NMAX=15,
                         JND=14, LMIN=0, LMAX=7, MENG=-15, EMIN=0, EMAX=2, 
                         ECORIC=0, debug=False):
    """
    Make input deck for an AS-DR run (based on a template) in given folder. 
    
    
    Inputs:

        folder: str, folder in which input deck is to be stored
            
        lmd: double, (1,d) numpy array of orbital scaling parameters, where
            d is the dimension of the parameter space.
            
        ion: State, ion 
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are 
                Schmidt orthogonalized. 
                
            NZION <0: nl-depdendent Hartree potential evaluated with Slater type 
                orbitals (NOT Schmidt orthogonalized).
                
        
        NMIN, NMAX: int, DRR NAMELIST, min/max n-value of Rydberg electron.
        
        JND: int, AS NAMELIST,
        
        LMIN, LMAX: int, DRR NAMELIST, min/max l-value of Rydberg electron.
        
        MENG: signed int, MENG is the number of interpolation energies (in
            Rydbergs). MENG energies follow and the continuum orbitals are 
            calculated at those energies. 
    
            = 0 (default) uses 0, DE/3, DE, 3*DE (& 8*DE if MAXLT/JT.gt.35/70) 
                where DE=max(TEAPOT, DELTAX) and TEAPOT is the estimated ionization 
                potential and DELTAX is the highest (spectroscopic) term energy.
            > 0 reads MENG final scattered energies following the NAMELIST.
            < 0 internally sets -MENG energies between a user specified range.
        
        EMIN, EMAX: double, minimum and maximum energies
        
        ECORIC:
            
    
    Outputs:
        
        asdeck_file: str, name of input file for the AS DR run.
    """
    current_folder = os.getcwd()
    
    #
    # Copy Autostructure DR input deck template to current folder
    # 
    
    # Name of input file
    asdeck_file = f"{ion.species}{ion.ion_charge}_das_{ion.shell}_n"
    
    # Location of template
    asdeck_path = root_folder()/'recombination'/'asdeck'/'dr'
    
    # Copy file template to local folder
    os.system(f"cp {asdeck_path}/{ion.isoelec_seq}-like.dr " +
              f"{folder}/{asdeck_file}")
    
    # Change to destination folder   
    os.chdir(folder)
 
    #
    # Lambda: Add add closed shell orbital 1s at the beginning
    # 
    lmd = lmd.ravel()  # ensure lambda is a vector
    lmd = np.insert(lmd, 0, 1.0)  # add 1.0 in first entry
    
    print(lmd)
    
    if debug:
        print('lambda = ', lmd)
        
    #
    # Modify DR input file
    # 
    with open(asdeck_file, "a") as asdeckin:
        """
        Write the DRR namelist.
        NMIN, NMAX - n shells used
        LMIN, LMAX - l orbitals used per n shell
        """
        asdeckin.write(f" &DRR    NMIN={NMIN} NMAX={NMAX} JND={JND} " + 
                       f"LMIN={LMIN} LMAX={LMAX} &END\n")
        
        asdeckin.write("16   20   25   35   45   55   70  100  140  " +
                       "200  300  450  700  999\n")
       
        """
        Write the SMINIM namelist, including lambda parameters
        NZION - Z of recombining atom
        PRINT - output file formatting
        NLAM - # of lambda parameters used
        """
        if (potential == 1):
           #
           # Thomas Fermi DA Potential
           # 
           asdeckin.write(" &SMINIM  " + 
                          f"NZION={np.sign(potential)*ion.nuclear_charge} " +
                          f"NLAM={len(lmd)} PRINT='UNFORM' &END\n")
        else:
            #
            # Hartree Potential
            # 
            asdeckin.write(" &SMINIM  " +
                           f"NZION={np.sign(potential)*ion.nuclear_charge} " +
                           f"NLAM={len(lmd)} ORTHOG='YES' PRINT='UNFORM' &END\n")
         
        #
        # Lambdas 
        #         
        asdeckin.write("  " + ' '.join([str(l) for l in lmd]) + "\n")
        
        #
        # Energy ranges
        # 
        asdeckin.write(f" &SRADCON  MENG={MENG} EMIN={EMIN} EMAX={EMAX} ")
        if ECORIC != 0:
            asdeckin.write(f"ECORIC={ECORIC} ")
        asdeckin.write("&END\n\n")
    
    # Change back to current folder
    os.chdir(current_folder)
    
    return asdeck_file
    

def run_asdr(lmd, ion, potential=1, NMIN=3, NMAX=15, JND=14, LMIN=0, LMAX=7, 
             MENG=-15, EMIN=0, EMAX=2, ECORIC=0, folder=None, 
             keep_folder=False, clean_folder=True, return_info=False,
             debug=True):
    """
    AUTOSTRUCTURE run for dielectronic recombination. Generates files for 
    ADAS-DR postprocessor.
    
    Inputs:
        
        lmd: 
    
        ion: State, 
        
        potential: int, potential type used to specify NZION in AS
            NZION >0: nl-dependent Thomas-Fermi DA potential. The orbitals are 
                Schmidt orthogonalized. 
                
            NZION <0: nl-depdendent Hartree potential evaluated with Slater 
                type orbitals (NOT Schmidt orthogonalized).
                
        
        NMIN, NMAX: int, DRR NAMELIST, min/max n-value of Rydberg electron.
        
        JND: int, AS NAMELIST,
        
        LMIN, LMAX: int, DRR NAMELIST, min/max l-value of Rydberg electron.
        
        MENG: signed int, MENG is the number of interpolation energies 
            (in Rydbergs). MENG energies follow and the continuum orbitals are 
            calculated at those energies. 

            = 0 (default) uses 0, DE/3, DE, 3*DE (& 8*DE if MAXLT/JT.gt.35/70) 
                where DE=max(TEAPOT, DELTAX) and TEAPOT is the estimated ionization 
                potential and DELTAX is the highest (spectroscopic) term energy.
            > 0 reads MENG final scattered energies following the NAMELIST.
            < 0 internally sets -MENG energies between a user specified range:
        
        EMIN, EMAX: double, minimum and maximum energy levels
        
        ECORIC:  
            
    
    Outputs: 
        
        None
    """
    current_folder = os.getcwd()
    
    #
    # Prepare Results Folder (if necessary)
    # 
    
    # Assign default name
    if folder is None:
        folder = current_folder + '/as_temp'
    
    # Make folder if it's not there already
    if not os.access(folder, os.F_OK):
        os.mkdir(folder)
    
    if debug:
        print('Writing Input Deck')
    
    #
    # Make Autostructure DR Input Deck
    #     
    asdeck_file = \
        make_asdr_input_deck(folder, lmd, ion, potential=potential, NMIN=NMIN, 
                             NMAX=NMAX, JND=JND, LMIN=LMIN, LMAX=LMAX, 
                             MENG=MENG, EMIN=EMIN, EMAX=EMAX, ECORIC=ECORIC,
                             debug=debug)
    
    os.chdir(folder)
    
    #
    # Run AUTOSTRUCTURE DR
    #
    
    # Path to AUTOSTRUCTURE
    asdeck_path = str(root_folder().resolve())+'/recombination'
    
    # Run compiled AS code
    os.system(asdeck_path + "/asdeck.x < " + asdeck_file)

    # Change name of the output for the postprocessor
    os.system("mv oicu o1u")
    
    # Go back to current folder
    os.chdir(current_folder)
    

def make_adasdr_input_deck(folder, lmd, ion, energies, ref_energies=None, 
                           NTAR1=1, compute_xsec=False, EWIDTH=0.001, 
                           NBIN=1000, EMIN=0.0, EMAX=2.0):
    """
    Write input deck for ADASDR Postprocessor (requires LEVELS file)
    
    Inputs:
        
        folder: str, path to folder in which to store results
        
        ion: State, specifying ion
        
        energies: double, energies (obtained from AUTOSTRUCTURE run).
        
        ref_energies: double, reference energies (usually NIST) for adjusting 
                
        NTAR1: int, AS NAMELIST input,
        
        compute_xsec: bool, whether to compute the cross-sections 
            (compute_xsec = TRUE) or the total rates (compute_xsec = FALSE)
        
        EWIDTH: double, AS NAMELIST input, width of convolution kernel.
        
        NBIN: int, AS NAMELIST input, number of bins used to compute 
            cross-section.
        
        EMIN, EMAX: double, AS NAMELIST input, minimum and maximum energies
        
        
    Outputs:
    
        If compute_xsec = True, return 
        
            energy, 
            
            cross-sections
            
        If compute_xsec = False, return
        
            T: temperature, 
            
            rate: DR rate for each point on the temperature grid 
        
    
    """
    # Record current folder 
    current_folder = os.getcwd()
    
    # Go to folder in which to store input deck
    os.chdir(folder)
    
    #
    # Determine whether to shift to reference energies
    # 
    nist_shift = ref_energies is not None
    
    #
    # Create adasin file
    # 
    adasdr_pp_file = 'adasin'
    with open(adasdr_pp_file, "w") as adasin:
        #
        # Set NTAR2 and NECOR variables 
        # 
        NTAR2 = len(energies)
        if nist_shift:
            NECOR=NTAR2
        else:
            NECOR=0
           
#            print('NIST shifts are set to',nist_shift,NTAR2,NECOR)
         
        adasin.write("/IC/\n")
        adasin.write(f" &ONE NTAR1={NTAR1} NTAR2={NTAR2} " +
                     f"COREX=\'{ion.shell}\' &END\n")
        adasin.write(f" &TWO NECOR={NECOR} ")
        
        #
        # Set parameters for computing cross-sections
        #  
        if compute_xsec:
            adasin.write(f"EWIDTH={EWIDTH} NBIN={NBIN} "+
                         f"EMIN={EMIN} EMAX={EMAX} ")
            
        adasin.write("&END\n")
        
        #
        # Add Energies and NIST reference energies (if necessary)
        # 
        with open("LEVELS", "r") as levels:
            lines = levels.read().splitlines()
            
            #
            # Add Energies 
            # 
            for i in range(1, len(lines)-1):
                line = lines[i].split()
                adasin.write(" " + " ".join(line[0:4]) + "\n")
            E_str = [str(x) for x in energies]
            
            adasin.write(" ".join(E_str))
            adasin.write("\n")
            
            #
            # Add reference energies if shifting to them
            #                 
            if nist_shift:
                #
                # Reference energies included 
                # 
                adasin.write(" ".join([str(E) for E in ref_energies]))
            else:
                #
                # No shift, use given energies
                # 
                adasin.write(" ".join([str(E) for E in energies]))
                
    
    os.chdir(current_folder)
    return adasdr_pp_file



def parse_adasdr_output(compute_xsec=False):
    """
    Read output file generated by the AUTOSTRUCTURE ADAS-DR Postprocessor and 
    return either
    """
    if compute_xsec: 
        #
        # Compute Cross-Sections (from "XRTOT" output file)
        #  
        data = np.transpose(np.genfromtxt("XDRTOT", skip_header=1))
        
        # First row has energy values
        energy = data[0,:]
        
        # Second row has cross-sections
        xsec= data[1,:]
    
        return energy, xsec

    else:
        #
        # Compute Total DR rates
        # 
        with open("adasout", "r") as f:
            #
            # Figure out how many temperature points there are 
            # 
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
                
            # Initialize temperature and rates vectors
            temperature = np.zeros(count)
            rate = np.zeros(count)
            
            f.seek(start)
            for i in range(count):
                line = f.readline().split()
                try:
                    temperature[i] = float(line[0])
                    rate[i] = float(line[1])
                except:
                    pass
        return temperature, rate
    



def run_adasdr(lmd, ion, energies, ref_energies=None, NTAR1=1, 
               compute_xsec=False, EWIDTH=0.001, NBIN=1000, EMIN=0.0, 
               EMAX=2.0, folder=None, keep_folder=False, clean_folder=True, 
               debug=False):
    """
    Run the ADAS-DR Postprocessor to compute the total dielectronic 
    recombination rates or the cross-sections
    
    Inputs:
    
        folder: str, folder in which to store input files and results
        
        ion: State, specify ion
        
        energies: double, energies computed by AUTOSTRUCTURE 
        
        ref_energies: double, reference (NIST) energies to shift towards during
            postprocessing
        
        NTAR1: int, AS NAMELIST input,
        
        compute_xsec: bool, whether to return cross sections vs energies
            (compute_xsec=True), or total dr rate coefficients vs temperature
            (compute_xsec=False).
        
        EWIDTH: double, AS NAMELIST input, width of convolution kernel.
        
        NBIN: int, AS NAMELIST input, number of bins used to compute 
            cross-section.
            
        EMIN, EMAX: double, AS NAMELIST input, minimum/maximum energies
        
        
    Outputs:
    
        If compute_xsec = True, return 
        
            energy, 
            
            cross-sections
            
        If compute_xsec = False, return
        
            T: temperature, 
            
            rate: DR rate for each point on the temperature grid 
    """
    current_folder = os.getcwd()
    
    #
    # Prepare Results Folder (if necessary)
    # 
    
    # Assign default name
    if folder is None:
        folder = current_folder + '/as_temp'
    
    # Make folder if it's not there already
    if not os.access(folder, os.F_OK):
        os.mkdir(folder)
        
    #
    # Make Post-Processor Input File
    # 
    adasdr_input_file = \
        make_adasdr_input_deck(folder, lmd, ion, energies, 
                               ref_energies=ref_energies, NTAR1=NTAR1, 
                               EWIDTH=EWIDTH, NBIN=NBIN, EMIN=EMIN, EMAX=EMAX)
    
    if debug:
        print('Input File', adasdr_input_file)
    
    os.chdir(folder)
    
    #
    # Run ADAS-DR Postprocessor
    #
    adasdr_path = str(root_folder().resolve())+'/recombination'
    os.system(adasdr_path + "/adasdr.x < " + adasdr_input_file)
    
    
    #
    # Return Information
    # 
    return parse_adasdr_output(compute_xsec)


    
def test01_energy_interpolation_accuracy():
    """
    Show 
    """
    #
    # Information on the Ion
    # 
    atom = "o"
    seq = "be"
    shell = "2-2"
    ion = State(atom, seq, shell)
    domain = np.array([[0.6,1.4],[0.6,1.4]])
    
    #
    # Reference Grid
    # 
    lmd_range = np.linspace(0.6,1.4,11)
    X,Y = np.meshgrid(lmd_range, lmd_range)
    ref_points = np.column_stack([X.ravel(),Y.ravel()])
    
    #
    # Reference solution
    # 
    E_ref = np.concatenate([run_structure(lmd,ion,1) for lmd in ref_points])
    n_outputs = E_ref.shape[1]
    n_inputs = 2  # dimension of parameter space 

    #
    # Regular Grid Interpolators
    #
    n_rgi = 5
    rgi_res = np.array([2**i for i in np.arange(1,n_rgi+1)])
    rgi_err = np.empty((n_rgi,n_outputs))
    for i in range(n_rgi):
        #
        # Generate mesh at given resolution
        #
        lmd_1d, lmd_pts = lambdas_grid(domain, [rgi_res[i],rgi_res[i]]) 
        
        #
        # Compute energies at given mesh points 
        # 
        E_vals = np.concatenate([run_structure(lmd,ion,1) for lmd in lmd_pts])
    
        #
        # Construct list of interpolators - one for each energy
        # 
        for j in range(n_outputs):
            # Reshape data to match inputs in lmd_1d
            data = E_vals[:,j].reshape(rgi_res[i],rgi_res[i])
            
            # Construct interpolators
            rgi = RegularGridInterpolator(lmd_1d,data,bounds_error=False)
            
            # Evaluate interpolators on fine mesh
            E_rgi = rgi(ref_points) 
        
            # Compute the maximum error
            rgi_err[i,j] = np.amax(np.abs(E_rgi-E_ref[:,j]))
    #    
    # Global sparse grids
    #
    n_pres = 5
    sg_res = np.empty(n_pres)
    sg_err = np.empty((n_pres,n_outputs))
    for i in range(n_pres):
        #
        # Form global sparse grid
        #  
        grid = Tasmanian.makeGlobalGrid(n_inputs,n_outputs,i,'level',
                                        'clenshaw-curtis')
        
        # Specify domain 
        grid.setDomainTransform(np.array([[0.6,1.4],[0.6,1.4]]))
        
        # Compute true energies at interpolation points
        GridPoints = grid.getNeededPoints()
        Energies = np.concatenate([run_structure(lmd,ion,1) \
                                   for lmd in GridPoints])
        
        # Construct interpolant
        grid.loadNeededPoints(Energies)
    
        # Number of interpolation points
        sg_res[i] = grid.getNumPoints()
    
        # Evaluate interpolators on fine mesh
        E_sg = grid.evaluateBatch(ref_points)
        
        for j in range(n_outputs):
            # Compare with reference
            sg_err[i,j] = np.amax(np.abs(E_sg[:,j]-E_ref[:,j]))
    #   
    # Plot results
    #  
    fig, axes = plt.subplots(3,3,sharex=True, sharey=True, figsize=(9,6))
    fig.suptitle('Maximum Interpolation Error for Energies \n Comparison: ' +
                 'Sparse Grid vs Piecewise Linear')
    count = 1
    for i in range(3):
        for j in range(3):
            axes[i,j].loglog(rgi_res**2,rgi_err[:,count],'.-',label=r'pwl')
            axes[i,j].loglog(sg_res,sg_err[:,count],'.-',label=r'sg')
            axes[i,j].set_title(f'$E_{count}$')
            if count == 1:
                axes[i,j].legend()        
            
            # y-label
            if j==0:
                axes[i,j].set_ylabel(r'Error')
            
            # x-label
            if i==2:
                axes[i,j].set_xlabel('Number of Points')
            count += 1
            
            axes[i,j].grid()
    fig.savefig('Interpolation_error_energies.png')
      
    
def test02_dr_interpolation_accuracy():
    """
    Test the accuracy of interpolators in reproducing the DR rates
    """        
    #
    # Information on the Ion
    # 
    atom = "o"
    seq = "be"
    shell = "2-2"
    ion = State(atom, seq, shell)
    domain = np.array([[0.6,1.4],[0.6,1.4]])
    
    #
    # Reference Grid
    # 
    lmd_range = np.linspace(0.6,1.4,11)
    X,Y = np.meshgrid(lmd_range, lmd_range)
    ref_points = np.column_stack([X.ravel(),Y.ravel()])
    
    #
    # Reference solution
    # 
    lmd = ref_points[0,:]
    r = get_rate(ion, lmd)
    energies = run_structure(lmd, ion, keep_folder=True)
    nist_path =root_folder()/'recombination'/'NIST'/'isoelectronic'/\
               f'{ion.isoelec_seq}/{ion.species}{ion.ion_charge}.nist'
               
    nist_energies = get_nist_energy(nist_path)
    run_asdr(lmd, ion, debug=True)
    temperature, rates = \
        run_adasdr(lmd, ion, energies, ref_energies=nist_energies, debug=True)
    #E_ref = np.concatenate([run_structure(lmd,ion,1) for lmd in ref_points])
    #n_outputs = E_ref.shape[1]
    #n_inputs = 2  # dimension of parameter space 


def comparison__structure():
    """
    """
    pass


def comparison_b_asdr():
    """
    Manual Test: Compare input files and results for asdr run  
    """ 
    # Record Current Folder
    current_folder = os.getcwd()
    
    # Ion 
    atom = "o"
    seq = "be"
    shell = "2-2"
    ion = State(atom, seq, shell)
    
    CFP = 1  # Central field potential
    
    # Lambda
    lmd = np.array([1.0, 0.6])
    
    # Updir
    up_dir = str(root_folder().resolve())+'/recombination/'
    
    # Location of NIST file
    nist_path = up_dir + \
                f"/NIST/isoelectronic/{ion.isoelec_seq}/"+ \
                f"{ion.species}{ion.ion_charge}.nist"
    
    # Get NIST energies          
    nist_vals = get_nist_energy(nist_path)
    
    # Set Maximum Energy
    EMAX=nist_vals[-1]*1.1  

    os.chdir(up_dir)
    
    print('Running AS')
    # Run Structure
    energies, e_nist, e_shift = structure(up_dir, ion, lambdas=lmd, 
                                          potential=CFP, emax=EMAX)
    
    print('Energies \n',energies)
    print('NIST Energies \n', e_nist)
    print('NIST shift \n', e_shift)
    
    print('Running Legacy ASDR')
    # Run ASDR
    structure_dr(ion, up_dir, potential=CFP, emax=EMAX, lambdas=lmd)
      
    print('Running ADAS-DR Post-Processor')
    temperature, rates = \
        postprocessing_rates(up_dir, ion, energies, e_nist, lambdas=lmd, 
                             compute_xsec=True)
    
    print('temperature \n', temperature)
    print('rates \n', rates)
    
    
    run_structure(lmd, ion, potential=CFP, EMAX=EMAX)
    """
    os.chdir(current_folder)
    print('Running New ASDR')
    run_asdr(lmd, ion)
    """
    
    # Run post-processor
    #temperature, rates = postprocessing_rates(up_dir, ion, energies)
    
    
def comparison_c_adasdr():
    """
    Compare the adasdr inputs with the 
    """
    pass


def time_grid_interpolators():
    """
    
    """
    dim_lmd = 3
    dim_out = 10
    
    fig, ax = plt.subplots(111)
    
    #
    # Hypercube [-1,1]^dim_lmd
    # 
    x_bnd = np.array([-np.ones(dim_lmd), np.ones(dim_lmd)]).T
    
    #
    # Vary the mesh resolution 
    # 
    for pts_per_dim in [2**i for i in np.arange(1,9)]:
        # Specify resolution
        resolution = [pts_per_dim for dummy in range(dim_lmd)]
        
        # Form the lambda grid
        X_1D, x_ravel = lambdas_grid(x_bnd, resolution)
        
        # Generate output data (random)
        grid_vals = [np.random.rand(*resolution) for dummy in range(dim_out)]
        
        # Form Interpolators
        interpolators = []
        for j in range(dim_out):
            interpolator = RegularGridInterpolator(X_1D, grid_vals[j], 
                                                   bounds_error=False)
            interpolators.append(interpolator)
        
        # Evaluate all interpolators at n points
        n = 10000
        x_eval = np.random.rand(n, dim_lmd)
        tic = time()
        for interpolator in interpolators:
            interpolator(x_eval)
        toc = time() - tic
        
        print(toc)
        
        tic = time()
        for i in range(n):
            for interpolator in interpolators:
                interpolator(x_eval[i])
        toc = time()-tic 
               
        print(' ', toc)
    
    
def time_sparse_grid_interpolators():
    """
    """
    print(dir(Tasmanian.SparseGrid))
    dim_lmd = 2
    dim_out = 10
    order = 1  # Piecewise linear
    depth = 10  # 
    grid = Tasmanian.SparseGrid()
    grid.makeGlobalGrid(dim_lmd,dim_out, depth, 'iptotal', 'leja')
    x = grid.getNeededPoints()
    n = grid.getNumPoints()
    print(x.shape)
    y = np.random.rand(n,dim_out)
    grid.loadNeededPoints(y)
    
    n = 100000
    x_eval = np.random.rand(n,dim_lmd)
    tic = time()
    for i in range(n):
        grid.evaluate(np.array(x_eval[i,:]))
    toc = time()-tic            
    print(toc)
    
    
if __name__=='__main__':
    #time_grid_interpolators()
    #time_sparse_grid_interpolators()
    #test01_energy_interpolation_accuracy()
    #test02_dr_interpolation_accuracy()
    comparison_b_asdr()
    """
    atom = "o"
    seq = "be"
    shell = "2-2"
    ion = State(atom, seq, shell)
    lmd = np.array([1,1])
    E = energy_map(lmd,ion,1)
    print(E.shape)
    lmd_range = np.linspace(0.6,1.4,101)
    X,Y = np.meshgrid(lmd_range, lmd_range)
    lmds = np.array([X.ravel(),Y.ravel()]).transpose()
    
    print(len(lmds))
    Egy = []
    count = 0
    for lmd in lmds:
        count += 1
        E = energy_map(lmd,ion,1)
        Egy.append(E)
        
    Egy =np.array(Egy).squeeze()
    print(Egy.shape)
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_wireframe(X, Y, Egy[:,1].reshape(X.shape))
    plt.show()
    
    print('Global Interpolant')
    n_inputs = 2
    n_outputs = 10
    precision = 10
    
    grid = Tasmanian.makeGlobalGrid(n_inputs,n_outputs,precision,'iptotal',
                                    'clenshaw-curtis')
    grid.setDomainTransform(np.array([[0.6,1.4],[0.6,1.4]]))
    GridPoints = grid.getNeededPoints()
    Energies = np.concatenate([energy_map(lmd,ion,1) for lmd in GridPoints])
    grid.loadNeededPoints(Energies)
    
    grid.plotPoints2D()
    plt.show()
    
    grid.plotResponse2D(iOutput=1)
    plt.show()
    
    print('Constructing Surrogate')
    num_inputs = 2
    num_outputs = E.shape[1]
    print(num_outputs)
    depth = 7
    order = -1
    rule = 'localp'
    grid = Tasmanian.makeLocalPolynomialGrid(num_inputs,num_outputs,depth,order,rule)
    grid.setDomainTransform(np.array([[0.6,1.4],[0.6,1.4]]))
    ftol = 1e-2
    f = lambda x, tid: energy_map(x,ion,2)
    budget = 1000
    Tasmanian.constructSurplusSurrogate(f,budget,1,1,grid,ftol,'classic')
    
    XY = np.column_stack([X.ravel(),Y.ravel()])
    Zi = grid.evaluateBatch(XY)
    
    grid.plotPoints2D()
    plt.show()
    
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    ax.plot_wireframe(X, Y, Zi[:,1].reshape(X.shape))
    plt.show()
    
    
    grid.plotResponse2D(iOutput=1)
    plt.colorbar()
    """  
    """
    while grid.getNumNeeded()>0:
        #  
        Tasmanian.loadNeededPoints(f,grid,1)
        
        # Set adaptive refinement strategy
        grid.setSurplusRefinement(ftol,0,'fds')
    """
    
    
    