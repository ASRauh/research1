#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 13 00:19:25 2019

@author: Arriën Symon Rauh
"""
# import scipy.constants as const
import argparse
import os
import OsTools3 as ot
from math import log, exp
import pandas as pd
import matplotlib.pyplot as plt
import scipy.constants as cst

def free_energy(probability, R=8.31451e-3, T=303):
    """
    Calculates a free energy value from a probability.

    Parameters
    ----------
    probability : float
        Probability of a state/orientation.
    R : float, optional
        Gas constant in kJ/molK. The default is 8.31451e-3 kJ/molK.
    T : float, optional
        Temperature in Kelvin. The default is 303 K.

    Returns
    -------
    free_energy : float
        Free energy value in kJ/mol.
    """
    try:
        free_energy = -R * T * log(probability)
    # Incase of log(0)
    except ValueError:
        free_energy = 0
    return free_energy


def probability_from_energy(energy, R=8.31451e-3, T=303):
    """
    Calculates a free energy value from a probability.

    Parameters
    ----------
    free_energy : float
        Free energy value in kJ/mol.
    R : float, optional
        Gas constant in kJ/molK. The default is 8.31451e-3 kJ/molK.
    T : float, optional
        Temperature in Kelvin. The default is 303 K.

    Returns
    -------
    probability : float
        Probability of a state/orientation.
    """
    energy_ratio = (energy / (R*T))
    p = exp(-energy_ratio)
    return p


def calculate_free_energy(probabilities):
    """
    Calculates free energies for a set of probabilities.
    """
    free_energy_data = []
    for probability in probabilities:
        free_energy_data.append(free_energy(probability))
    return free_energy_data


def calculate_probabilities(counts_dist):
    """
    Calculates free energies for a set of counts.
    """
    probabilities_dist = []
    total_count = sum(counts_dist)
    for i in counts_dist:
        probabilities_dist.append(i/total_count)

    return probabilities_dist


def eyring_equation(experimental_value, experimental_sd , T=303, rate=False):
    """
    Two way calculations with the Eyring equations.

    Parameters
    ----------
    experimental_value : float
        Can be either a rate (/s) or a free energy value (kJ/mol).
    experimental_sd : float
        The standard deviation of the experimental value.
    T : float, optional
        The temperature in Kelvin. The default is 303 K.
    rate : Boolean, optional
        Check for what the experimental value is: True means a rate as input
        and False means an energy as input. The default is False.

    Returns
    -------
    dGs OR ks : float
        The output is either an energy value (dGs in kJ/mol) or a rate (ks in /s).
    """
    A =  (cst.Boltzmann/1000 * T) / (cst.Planck/1000)
    experimental_values = [experimental_value - experimental_sd , 
                           experimental_value, 
                           experimental_value + experimental_sd]

    if rate: # Input is a rate 
        dGs = []
        for k in experimental_values:
            dGs_t = -(cst.gas_constant / 1000) * T * log((k / A))
            dGs.append(dGs_t)
        return dGs
    elif not rate: # Thus experimental value is an energy
        ks = []
        for dG in experimental_values:
            ks_t = A * exp(-dG/((cst.gas_constant / 1000) * T))
            ks.append(ks_t)
        return ks


def plot_energies(xcoords, energies, com_dist, output_dir):
    plt.plot(xcoords, energies)
    plt.xlabel("Conformation - RMSD to native (nm)")
    plt.ylabel(r"F $(kJ mol^{-1})$")
    plt.title(F"Free energies of conformation per bin for COM distance {com_dist} nm")
    plt.savefig(os.path.join(output_dir, F"energies_d{com_dist}.png"), dpi=300)
    plt.close()

# "Conformation (pc2) a.u."

def plot_probabilities(xcoords, counts, com_dist, output_dir):
    plt.plot(xcoords, counts)
    plt.xlabel("Conformation - RMSD to native (nm)")
    plt.ylabel("Probability")
    plt.title(F"Conformation count per bin for COM distance {com_dist} nm")
    plt.savefig(os.path.join(output_dir, F"probabilities_d{com_dist}.png"), dpi=300)
    plt.close()


def calculate_per_distance(counts_data, output_dir, plot=False):
    """
    Calcluates the orientational free energies at every COM distance based on 
    the bin counts of the orienational axis.

    Parameters
    ----------
    counts_data : pandas.DataFrame
        All the bin counts resulting from binning of the orientational axis.
    output_dir : string
        A directory to store the data in.
    plot : boolean, optional
        An option to plot the orientational probabilities/free energies.
        The default is False.

    Returns
    -------
    probabilities : pandas.DataFrame
        The caclulated orientational probabilities (normalised counts).
    free_energies : pandas.DataFrame
        The caclulated orientational free energies in kJ/mol.
    """
    probabilities = pd.concat([pd.DataFrame(), counts_data[counts_data.columns[0]]])
    free_energies = pd.concat([pd.DataFrame(), counts_data[counts_data.columns[0]]])

    # Calculate per distance:
    for dist in counts_data.columns[1:]:
        # The probabilities for every conformation
        probabilities[dist] = calculate_probabilities(counts_data[dist])
        # And the corresponding free energy in kJ/mol
        free_energies[dist] = calculate_free_energy(probabilities[dist])
        if plot:
            # Plot the probabilities and free energies
            plot_probabilities(counts_data[counts_data.columns[0]],
                               probabilities[dist], dist, output_dir)
            plot_energies(counts_data[counts_data.columns[0]],
                          free_energies[dist], dist, output_dir)

    # Write output
    probabilities.to_csv(os.path.join(output_dir, F"probabilities.txt"), sep='\t', index=False, float_format='%.7f')
    free_energies.to_csv(os.path.join(output_dir, F"energies.txt"), sep='\t', index=False, float_format='%.7f')
    return probabilities, free_energies


def main(counts_data, output_dir, plot=False):
    """
    """
    calculate_per_distance(counts_data, output_dir, plot=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
                                     '''
                                     This script takes 
                                     ''')
    parser.add_argument("-i", "--counts_data", help="", required=True)
    parser.add_argument("-o", "--output_dir", help="", required=True)

    args = parser.parse_args()
    counts = args.counts_data
    output_dir = args.output_dir

    main(ounts_data, output_dir, plot=False)
