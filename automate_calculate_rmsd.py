#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 17:42:58 2019

@author: Arriën Symon Rauh
"""
import os
import re
#import subproceess
import argparse
#import xvg
#import numpy as np
import pandas as pd
# Import the impoprtant package
import calculate_rmsd as rmsd


def main(directory, references_dir, output_file, ndx_file, version='n'):
    """
    Automates the RMSD calculation for all simulations
    (all COM distances and starting structures).

    Parameters
    ----------
    directory : string
        Directory containing all the trajectory data.
    references_dir : string
        DESCRIPTION.
    output_file : string
        DESCRIPTION.
    ndx_file : string
        Filename of the GROMACS index file (.ndx extension).
    version : str, optional
        Set the version of GROMACS you want to use for RMSD calculation. 
        The default is 'n'.

    Returns
    -------
    None.

    """
    check_ndx = os.path.isfile(ndx_file)
    rmsd_data = pd.DataFrame()
    # Simulation directory: sim...
    for simulation in os.listdir(directory):
        if os.path.isdir(simulation) and re.match("sim\d{0,4}", simulation):
            # Distance directory; d...
            for dist in os.listdir(simulation):
                loc_dir = os.path.join(simulation, dist)
                if os.path.isdir(loc_dir) and re.match("d\d.\d{1,2}", dist):
                    # Make a string with full path of trajectory file
                    traj_file = "traj_" + simulation + "_" + dist[1:] + ".xtc"
                    traj_file = os.path.join(simulation, dist, traj_file)
                    #
                    ref_dir = os.path.join(references_dir, dist)
#                    print(traj_file, ref_dir)

                    # Check if all files are present
                    check_traj = os.path.isfile(traj_file)
                    check_ref = os.path.isdir(ref_dir)
                    check = (check_traj and check_ref and check_ndx)
#                    print(check_traj, check_ref, check_ndx, check)
                    if check:
                        # Calculate RMSD for all five references
                        new_data = rmsd.calculate_rmsd(traj_file, ref_dir, 
                                                   ndx_file, version='n',
                                                   remove=False)
                        # Add new data to a dataframe
                        rmsd_data = pd.concat([rmsd_data, new_data], ignore_index=True,
                                              axis=0)
                        
                        print("Calculated RMSD for: " + simulation + "-" + dist)
                    else:
                        print("Error for : " + simulation + '-' + dist)
        else:
            print("Error for : " + simulation)
    
    # Write DataFrame to output_file
    rmsd.write_output(rmsd_data, output_file)


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=''' DESCRIPTION ''')
    parser.add_argument("-d", "--main_dir", help=""""Path to a directory containing 
                        all trajectory files (xtc/pdb)""",required=True)
    parser.add_argument("-s", "--references_dir", help="Directory containing ", required=True)
    parser.add_argument("-o", "--output_file", help="...", required=True)
    parser.add_argument("-n", "--index_file", help="...", required=True)
    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help='''To enebale usage with both old and new versions of GroMaCS''')

    args = parser.parse_args()
    directory = args.main_dir
    references_dir = args.references_dir
    output_file = args.output_file
    ndx_file = args.index_file
    version = args.version
    
    main(directory, references_dir, output_file, ndx_file, version)