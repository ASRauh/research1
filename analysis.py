#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 29 18:57:02 2020

@author: Arrienrauh

Script to tun all the analysis steps for creating a 2D free energy landscape.

"""
# general modules
import pandas as pd 
import os
import OsTools3 as ot
# import specific analysis modules
import pca_rmsd as pca
import bin_orientations as bo
import calculate_free_energy as cfe
import constraint_force_integration as cfi
import energy_landscape as el


def read_data(data_file):
    """
    """
    data = pd.read_csv(data_file, sep='\t')
    return data


####################################
##########   Load data    ##########
####################################
data = read_data("../data/rmsd/rmsd.txt")

# For the RMSD to native
data_native = read_data("../data/rmsd_native/rmsd_native.txt")


####################################
##########      PCA       ##########
####################################
data = pca.pca_rmsd(data, "../data/rmsd/pca", 10, 5)


####################################
##########   Binning PC   ##########
####################################
counts = bo.bin_orientations(data, "../data/rmsd/energies/conf_counts.txt", 'pc2', 0.01)

# For the RMSD to native
counts_native = bo.bin_orientations(data_native, "../data/rmsd_native/energies/conf_counts.txt", 'native', 0.01)


####################################
##########   Orient.En.   ##########
####################################
probs, orientational_energies = cfe.calculate_per_distance(counts, "../data/rmsd/energies", plot=True)

# For the RMSD to native
probs_native, fes_native = cfe.calculate_per_distance(counts_native, "../data/rmsd_native/energies", plot=True)



####################################
##########      PMF       ##########
####################################
# def calculate_pmf(data, output_file, png_file):
#     """
#     """
#     if not os.path.isdir(ot.file_path(output_file)):
#         os.makedirs(ot.file_path(output_file))
#     #output_file = F"{directory}/pmf_all_data.txt"
#     #png_file = F"{directory}/pmf_all_data.png"
#     # Calculate PMF
#     pmf = cfi.constraint_force_integration(data, output_file, png_file, ret=True)
#     return pmf

# Calculate pmf for all data
pmf = cfi.constraint_force_integration(data, "../data/rmsd/pmf_data/pmf_all_data.txt",
                                       "../data/rmsd/pmf_data/pmf_all_data.png")

# Calculate pmf per class of structures
pmf_capri = pd.DataFrame()
for capri in data.capri.unique():
    pmf_temp = pmf = cfi.constraint_force_integration(data.loc[data['capri']==capri,], 
                                                      F"../data/rmsd/pmf_data/capri/pmf_{capri}.txt",
                                                      F"../data/rmsd/pmf_data/capri/pmf_{capri}.png")
    pmf_capri['com_dist'] = pmf_temp['com_dist']
    pmf_capri[F'{capri}'] = pmf_temp['pmf']
    pmf_capri[F'error_{capri}'] = pmf_temp['error']


####################################
##########   Landscape    ##########
####################################
energy_landscape = el.determine_energy_landscape(pmf, orientational_energies,
                                                 output_file="../data/rmsd/energy_landscape/energy_landscape.png",
                                                 png_file="../data/rmsd/energy_landscape/energy_landscape.txt")



