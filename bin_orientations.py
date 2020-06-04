#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 21:50:04 2019

@author: Arriën Symon Rauh


"""
import os
import argparse
import pandas as pd
import numpy as np
import calculate_free_energy as cfe
import matplotlib.pyplot as plt


def read_data(data_file):
    """
    """
    data = pd.read_csv(data_file, sep='\t')
    return data


def define_bin_ranges(data, column, bin_width):
    """
    Calculate the bin ranges based on the bin width.
    """
    # Extract column of interest
    outter_ranges = (min(data[column]), max(data[column]))
    bin_ranges = []
    x = outter_ranges[0]
    y = outter_ranges[1]

    n_bins = int((y - x) / bin_width)+1
    for i in range(n_bins):
        bin_range = (round(x,6), round(x + bin_width,6))
        bin_ranges.append(bin_range)
        x += bin_width
    
    ranges = pd.IntervalIndex.from_tuples(bin_ranges)
    return ranges #, bin_ranges]


def get_counts_per_bin(data, ranges):
    """
    Calculate bin counts.
    """
    categorical_object = pd.cut(data, ranges)
    counts = pd.DataFrame(categorical_object.value_counts().sort_index())
    return counts


def plot_energies(xcoords, energies, com_dist, output_dir):
    plt.plot(xcoords, energies)
    plt.xlabel("Orientation bin #")
    plt.ylabel("F (kT)")
    plt.title(F"Free energies per bin for COM distance {com_dist} nm")
    plt.savefig(os.path.join(output_dir, "energies_d{com_dist}.png"), dpi=300)


def plot_probabilities(xcoords, counts, com_dist, output_dir):
    plt.plot(xcoords, counts)
    plt.xlabel("Orientation (pc2)")
    plt.ylabel("Probability")
    plt.title(F"Probability per bin for COM distance {com_dist} nm")
    plt.savefig(os.path.join(output_dir, "probabilities_d{com_dist}.png"), dpi=300)


def bin_orientations(data, output_file, column,  bin_width=-0.01):
    """
    Bins the orientational axis and calculates the orientational
    counts for every COM distance.

    Parameters
    ----------
    data : pandas.DatFrame
        All data.
    output_file : string
        Filename of the of the output file.
    column : string
        Column label for the column of interest.
    bin_width : float, optional
        The bin width for the binning. The default is -0.01.

    Returns
    -------
    all_counts : pnadas.DatFrame
        Collection of the orientational counts for every COM distance.

    """
    # Load the data from 
    iteration_list = pd.unique(data["com_dist"])
    all_counts = pd.DataFrame()

    # Define the bin ranges
    ranges = define_bin_ranges(data, column, bin_width)
    # Define coordinates
    coords = [ ranges[i].mid for i in range(len(ranges))]
    all_counts[F"{column}_bin"] = coords

    for dist in iteration_list:
        # Extract colomn of data of interest
        indices_of_interest = data["com_dist"] == dist
        data_of_interest= data.loc[indices_of_interest, column]
        # Get counts per bin in a DataFrame
        counts = get_counts_per_bin(data_of_interest, ranges)
        all_counts[dist] = counts[column].values
    all_counts.to_csv(output_file, sep='\t', index=False, float_format='%.7f')
    return all_counts

#    fig = plt.figure()  # an empty figure with no axes


def main(data_file, output_file, column,  bin_width=-0.01):
    data = read_data(data_file)
    bin_orientations(data, output_file, column,  bin_width=-0.01)


if __name__=="__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='''The script retrieves orientational RMSD data
                                                    from different file and performs a principal 
                                                    components analysis (PCA) on it''')
    parser.add_argument("-i", "--data_file", help="", required=True)
    parser.add_argument("-o", "--output_file", help="", required=False)
    parser.add_argument("-col", "--column", help="The column of the data matrix to bin.", required=True)
    parser.add_argument("-w", "--bin_width", help="The width of the bins", required=False)

    args = parser.parse_args()

    data_file = args.data_file
    output_file = args.output_file
    column = args.column
    bin_width = float(args.bin_width)
#    distance = args.distance
#    std = args.std

    # Check if the directory specified under data_dir and results_dir exist.
    # If not make them:
    results_dir = os.path.split(output_file)[0]
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)

    # Call main()
    main(data_file, output_file, column, bin_width)