#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 22:23:14 2020

@author: Arriën Symon Rauh
"""
import os
import OsTools3 as ot
import sys
import argparse
import pandas as pd
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
import math


distances = [2.5 , 2.55, 2.6 , 2.65, 2.7 , 2.75, 2.8 , 2.85, 2.9 , 2.95, 3.0 ,
             3.1 , 3.2 , 3.3 , 3.4 , 3.5 , 3.6 , 3.7 , 3.8 , 3.9 , 4.0 , 4.2 ,
             4.4 , 4.6 , 4.8 , 5.0]

def read_data(data_file):
    """
    """
    data = pd.read_csv(data_file, sep='\t')
    return data


def blocked(forces, nblocks):
    """
    Performs block avaeraging.

    Parameters
    ----------
    forces : dataFrame
        The data to calcluate an average over.
    nblocks : integer
        Number of blocks to separate the data into.
    """
    size = len(forces)/nblocks
    blocks = []
    for block in range(nblocks):
        a = []
        for ts in forces[int(block*size):int((block+1)*size)]:
            a.append(ts)
        blocks.append(np.average(a))
    blockaverage = np.average(blocks)
    blockstd = np.std(blocks)
    blockvar = np.var(blocks)
    return nblocks, size, blockaverage, blockstd, blockvar


def deltaG_off(pmf):
    
    dGoff = pmf.min() - pmf[pmf.argmin():].max()
    
    return dGoff



def mean_force_dist(forces, r):
    """
    Extract the mean force
    """
    # Calculate mean per distance
    F_mean = blocked(forces, 100)[2]#-np.average(forces)
    # Estimate the error htrough block averaging *NEEDS CHECKING*
    F_mean_error = blocked(forces, 100)[3]# Block averaging
    return [r, F_mean, F_mean_error]


def mean_force(data):
    """
    """
    # Calculate mean per distance
    mean_force_per_distance = []
    for r in np.sort(data.com_dist.unique()):
        forces_r = data.loc[data['com_dist'] == r, 'force']
        mean_force_per_distance.append(mean_force_dist(forces_r, r))
#        print(r, mean_force_per_distance[-1][0])
    return mean_force_per_distance


def mf_w_entropy_corr(mf,kB,T):
    """
    Performs entropy correction.
    """
    mf_corr = []
    for d in mf:
        mf_corr.append([d[0],d[1]+(2.0*T*kB/d[0]),d[2]])
    return mf_corr

def calculate_pmf(mf,mf_corr):
    """
    Calculates the PMF by numerically integrating the mean force.
    """
    pmf_of_d = []
    # Weights??
    w        = []

    for d in mf:
        pmf_of_d.append([d[0],0.0,0.0])
        w.append(0.0)

    pmf_of_d[-1][1] = 0.0
    for i in reversed(range(len(mf_corr)-1)):
        # Weights to acconut for distance between COM-constraints
        hh = 0.5*(mf[i+1][0]-mf[i][0])
        #
        pmf_of_d[i][1] = pmf_of_d[i+1][1] - hh*(mf_corr[i+1][1]+mf_corr[i][1])
#        print(pmf_of_d[i][1], " = ", pmf_of_d[i+1][1]," - " , hh*(mf_corr[i+1][1]+mf_corr[i][1]))

        w[i]   += hh
        w[i+1] += hh

    pmf_of_d[-1][2] = 0.0
#    print(w, mf_corr)
    var_int = math.pow(w[-1]*mf_corr[-1][2],2)

    for i in reversed(range(len(mf_corr)-1)):
        hh = 0.5*(mf[i+1][0]-mf[i][0])
        pmf_of_d[i][2] = var_int + math.pow(hh*mf_corr[i][2],2)
        var_int       += math.pow(w[i]*mf_corr[i][2],2)

    for i in range(len(pmf_of_d)):
        pmf_of_d[i][2] = math.sqrt(pmf_of_d[i][2])

    return pmf_of_d


def output_pmf(pmf, output_file):
    """
    """
    pmf.to_csv(output_file, sep='\t', index=False, float_format='%.7f')
    # with open(output_file,'w') as of:
    #     for i in range(0,len(pmf)):
    #         of.write(F"{pmf[i,0]}'\t'{pmf[i,1]}'\t'{pmf[i,2]}\n")


def plot_pmf(pmf, png_file):
    """
    """
    plt.plot(pmf[:,0],pmf[:,1])
    plt.xlabel("COM distance (nm)")
    plt.ylabel(r"PMF (F $(kJ mol^{-1})$")
    plt.title(F"{png_file.split('.')[-2].split('/')[-1].replace('_', ' ')}")
    plt.savefig(png_file, dpi=300)
    plt.close()


def constraint_force_integration(data_file, output_file, png_file, ret=False):
    """
    Perfoms the constraint force integration.
    """
    if type(data_file) == type(pd.DataFrame()):
        data = data_file
    elif os.path.isfile(data_file):
        data = read_data(data_file)
    try:
         # Calculate mean force for every distance
         mean_forces = mean_force(data)
    except NameError:
        print("Data Error: data input needs to be checked.")
        sys.exit(1)
    # Definitions of constants
    kB = 8.31451e-3
    T = 303
    # Calculate correction of mean forces
    mf_corr = mf_w_entropy_corr(mean_forces,kB,T)
    # Calculate the PMF
    pmf = np.array(calculate_pmf(mean_forces,mf_corr))
    df_pmf = pd.DataFrame({'com_dist':pmf[:,0], 'pmf':pmf[:,1], 'error':pmf[:,2]})

    # Check output directory
    if not os.path.isdir(ot.file_path(output_file)):
        os.makedirs(ot.file_path(output_file))
    # Write PMF data to file
    output_pmf(df_pmf, output_file)
    # Write plot to file
    if png_file:
        plot_pmf(pmf, png_file)
    if ret:
        return df_pmf


def main(data_file, output_file, png_file):
    """
    """
    constraint_force_integration(data_file, output_file, png_file)


#############MAIN##############
if __name__ == "__main__":
    "Do the work"
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=
                                     '''This script prodcues a PMF from pull
                                        forces''')
    parser.add_argument("-i", "--data_file", help="", required=True)
    parser.add_argument("-o", "--output_file", help="", required=False)
    parser.add_argument("-p", "--plot_pmf", dest="plot_pmf", action="store_true",
                        help='''For including a plot of the pmf''')

    args = parser.parse_args()
    data_file = args.data_file
    output_file = args.output_file
    png_file = args.plot_pmf

    main(data_file, output_file, png_file)

