#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 12:34:31 2019

@author: Arriën Symon Rauh


This script takes a pdb containing a trajectory file () and calculates 
the root mean square deveiation for every frame towards a set of specified refence structures.

Input:

Output:


Example usage:
python3 -t -ref -o -v -r


!!DEPENDENCIES!!
- Needs GROMACS installation.
- Needs xvg.py
- Needs OsTools3.py
"""
import os
import OsTools3 as ot
import subprocess
import argparse
#import numpy as np
import pandas as pd
import xvg



# Dictionary translating the reference number to the reference
# coordinate file and the xvg output file.
references = {1:["sel2.pdb","ref1.xvg"], 2:["sel3.pdb","ref2.xvg"], 
              3:["sel4.pdb","ref3.xvg"], 4:["sel6.pdb","ref4.xvg"], 
              5:["sel7.pdb","ref5.xvg"]}


def gmx_rms(traj_file, reference, xvg_file, ndx_file, version='n'):
    """
    Call the 

    Parameters
    ----------
    traj_file : string
        Filename of the GROMACS trajectory file (.xtc).
    reference : string
        coordinate file of the reference structure.
    xvg_file : string
        Name of the output xvg file (.xvg extension).
    ndx_file : string
        Filename of the GROMACS index file (.ndx extension).
    version : str, optional
        Set the version of GROMACS you want to use for RMSD calculation. 
        The default is 'n'.

    Returns
    -------
    None.

    """ 
    # Compatible with all versions of Gromacs
    if version == 'n':
        cmd = F"echo 0 0 | gmx rms -s {reference} -f {traj_file} -o {xvg_file} -n {ndx_file}"
    else:
        cmd = F"echo 0 0 | g_rms -s {reference} -f {traj_file} -o {xvg_file}  -n {ndx_file}"

    print(cmd)
    # Running shell command in python script. 
    #See https://docs.python.org/2/library/subprocess.html#popen-constructor
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, 
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         close_fds=True)
    p.wait()


def extract_xvg_data(xvg_file, label):
    """
    Extracts the y column from an xvg file (RMSD) and adds it to a dataframe.

    Parameters
    ----------
    xvg_file : string
        An xvg filename from which the data should be extracted.
    label : string
        Label for the column in the DataFrame.

    Returns
    -------
    data : pandas.DataFrame
        The data extracted from the xvg file.
    """
    data = pd.DataFrame(xvg.XvgSingle.from_fname(xvg_file)._array[:,1], 
                        columns=[label])
    
    return data


def calculate_average(index, data, label):
    """
    Calculates the average around a timepoint t.
    Average of timepooint t is calculuted over the measurement at t
    and five before and after (t±5).
    """
    if index == 0:
        start = index
        end = index + 6
        n = 6.0
    elif index == 20000:
        start = index - 5
        end = index + 1
        n = 6.0
    else:
        start = index - 5
        end = index + 6
        n = 11.0
    # Start cumilative sum
    cumsum = 0.0
    for i in range(start, end):
        cumsum+=data.loc[i,label]
    # Calculate the average
    average = cumsum / n
    return average

def extract_pull_force(force_xvg, label):
    """
    Extracts the measured constraint force (force needed to keep the
    system in its contstaints) from the GROMACS output file.

    Parameters
    ----------
    force_xvg : string
        Filename of the GROMACS output file containing the constraint
        force values.
    label : string
        DESCRIPTION.

    Returns
    -------
    force : pandas.DataFrame
        The values of the constraint force for all the frames.

    """
    # Extract data from pullf.xvg file
    data = extract_xvg_data(force_xvg, label)
    # Caclulate averages
    force_averages = []
    for frame in range(0,20010,10): # Number of frames times 10
        force_averages.append(calculate_average(frame, data, label))
    
    force = pd.DataFrame({label:force_averages})
    return force


def write_output(data, output_file):
    data.to_csv(output_file, sep='\t', index=False, float_format='%.7f')


def extract_metadata(traj_file):
    """
    Parses the filename to extract metadata
    """
    file = os.path.split(traj_file)[-1][:-4].split('_')
    simulation = float(file[1][3:-1])
    replicate = float(file[1][-1])
    distance = float(file[2])
    return (simulation, replicate, distance)


def calculate_rmsd(traj_file, references_dir, ndx_file, 
                   output_file="rmsd_data.txt",
                   version='n', remove=False):
    """
    Calculates the RMSD of one trajectory to multiple reference structures adds
    the constraint forces and writes to a file in tsv format.

    Parameters
    ----------
    traj_file : string
         Filename of the GROMACS trajectory file (.xtc).
    references_dir : string
        DESCRIPTION.
    ndx_file : string
        Filename of the GROMACS index file (.ndx extension).
    output_file : string, optional
        The filename for the file to which the RMSD data is written.
        The default is "rmsd_data.txt".
    version : str, optional
        Set the version of GROMACS you want to use for RMSD calculation. 
        The default is 'n'.
    remove : boolean, optional
        Option to remove the data files containing the RMSD to single
        references to prevent cluttering. The default is False.

    Returns
    -------
    None.

    """
    # GroMaCS version handling
    version = 'n'
    # Initialise a DataFrame
    rmsd_data = pd.DataFrame()
    file_path = ot.file_path(traj_file)[0]

    for ref, xvg_file in references.values():
#        print(ref, xvg_file)
        ref_file = os.path.join(references_dir, ref)
        xvg_file = os.path.join(file_path, xvg_file)
#        print(os.path.isfile(ref_file), ot.file_extension(ref_file), ref_file, xvg_file)
        if (ot.file_extension(ref_file) == "pdb") and os.path.isfile(ref_file):
            
            # Calculate RMSD using GroMACS
            # outputs a xvg-file
            print(traj_file, ref_file, xvg_file, ndx_file, version)
            gmx_rms(traj_file, ref_file, xvg_file, ndx_file, version)
            
            # Mine the results from the xvg output file file
            label = ot.filename(xvg_file)
            new_column = extract_xvg_data(xvg_file, label)
            rmsd_data  = pd.concat([rmsd_data, new_column], axis=1)
            if remove:
                #os.remove(xvg_file)  # Remove xvg file
                pass

    # Add Pull-force
    force_xvg = os.path.join(file_path, "pullf.xvg")
    force = extract_pull_force(force_xvg, "force")
    rmsd_data  = pd.concat([rmsd_data, force], axis=1)

    # Add metadata...
    # Will add com_dist, simulation #, replicate number, frame #
    simulation, replicate, distance = extract_metadata(traj_file)
    rmsd_data["com_dist"] = distance
    rmsd_data["sim"] = simulation
    rmsd_data["rep"] = replicate
    rmsd_data["frame"] = list(range(0,2001))

    # Write output to file
    write_output(rmsd_data, os.path.join(file_path, output_file))
    return(rmsd_data)


def main(traj_file, references_dir,  ndx_file, version, remove=False):
    calculate_rmsd(traj_file, references_dir, ndx_file, 
                   version=version, remove=False)


if __name__ == "__main__":

    # Parse command line arguments
    parser = argparse.ArgumentParser(description=''' DESCRIPTION ''')
    parser.add_argument("-t", "--trajectory_file", help="...", required=True)
    parser.add_argument("-ref", "--references_dir", help="...", required=True)
#    parser.add_argument("-o", "--output_file", help="...", required=True)
    parser.add_argument("-n", "--index_file", help="...", required=True)
    parser.add_argument("-v", "--version", dest="version", action="store_true",
                        help='''To enebale usage with both old and new versions of GroMaCS''')
    parser.add_argument("-r", "--remove", dest="remove", action="store_true",
                        help='''Remove xvg files after mining RMSD data''')

    args = parser.parse_args()
    traj_file = args.trajectory_file
    references_dir = args.references_dir
#    output_file = args.output_file
    ndx_file = args.index_file
    version = args.version
    remove = args.remove
    
    main(traj_file, references_dir, ndx_file, version, remove)