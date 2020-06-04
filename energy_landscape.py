#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 10:40:16 2020

@author: Arrienrauh
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import calculate_free_energy as cfe


def read_data(data_file):
    """
    """
    data = pd.read_csv(data_file, sep='\t')
    return data


############################################
######     Through probabilities      ######
############################################
###
### Part that shows that the derivations
### works.
###
def calculate_weights(pmf, R=8.31451e-3, T=303):
    """
    """
    weights = []
    ref_index = pmf.shape[0]-1
    for r in range(0, pmf.shape[0]-1):
        # Get the difference in Free energy between distance r and the reference
        df = round(pmf.loc[r,'pmf'] - pmf.loc[ref_index,'pmf'], 7)
        pj_ref = cfe.probability_from_energy(df)
        print(df, pj_ref)
        weights.append(pj_ref)
    weights.append(1.0)
    total = sum(weights)
    print(total)
    # Normalise the weights
    for i in range(len(weights)):
        weights[i] /= total
    return weights


def reweight_probabilities(probabilities, weights):
    """
    Reweights the orientational probability by multiplying 
    with the corresponding weight based on the PMF.
    """
    reweighted_probabilities = np.zeros(probabilities.shape)

    # For all distances
    for r in range(0, probabilities.shape[1]): #-1):
        # And all orientations
        for i in range(probabilities.shape[0]):
            if probabilities.iloc[i,r] != 0.0:
                # Reweight the orientational probability by multiplying
                # with the corresponding weight based on the PMF.
                reweighted_probabilities[i,r] = probabilities.iloc[i,r] * weights[r]
                print(probabilities.iloc[i,r] * weights[r])
            else:
                pass
    return reweighted_probabilities


def reweighted_energies(reweighted_probabilities):
    """
    Calculates the energy value for the reweighted probabilities.
    """
    energies = np.zeros(reweighted_probabilities.shape)

    for r in range(0, reweighted_probabilities.shape[1]): #-1):
        for i in range(reweighted_probabilities.shape[0]):
            if reweighted_probabilities[i,r] != 0.0:
                energies[i,r] = cfe.free_energy(reweighted_probabilities[i,r])
            else:
                energies[i,r] = 100.0
    return energies


def determine_energy_landscape_probabilities(pmf, conf_energies, probabilities, 
                               output_file, png_file, ret=True):
    weights = calculate_weights(pmf)
    reweighted_probabilities = reweight_probabilities(probabilities.iloc[:,1:],
                                                      weights)
    energies = reweighted_energies(reweighted_probabilities)

    # X: com_dist
    x = pmf["com_dist"]
    # Y: conformation axis
    y = conf_energies["0"]
    # Create grid
    X, Y = np.meshgrid(x, y)
    # Z: Free energy
    z = energies
    return()



############################################
######        Through energies        ######
############################################
def reweight_orientational_energies(orientational_energies, weights):
    """
    Reweights the orientational free energies by adding the corresponding value
    from the PMF.

    Parameters
    ----------
    orientational_energies : pandas.DataFrame
        Free energy values for the different orientations.
    weights : pandas.DataFrame
        The free energy values from the PMF.

    Returns
    -------
    reweighted_energies : TYPE
        The reweighted orientational free energies.

    """
    reweighted_energies = np.zeros(orientational_energies.shape)

    # For all distances
    for r in range(0, orientational_energies.shape[1]):
        # And all orientations
        for i in range(orientational_energies.shape[0]):
            # If the energy is not 0
            if orientational_energies.iloc[i,r]!= 0.0:
                # Add the corresponding value of the PMF
                reweighted_energies[i,r] = orientational_energies.iloc[i,r] + weights[r] # check
            else:
                reweighted_energies[i,r] = 100.0 - 45.0
            print(reweighted_energies[i,r])
    return reweighted_energies


def determine_energy_landscape(pmf, orientational_energies, 
                               output_file, png_file, ret=True):
    """
    Determines the free energy landscape for a protein-protein interaction
    based on the PMF and the orientational free energies.

    Parameters
    ----------
    pmf : pandas.DataFrame
        The potential of mean force data.
    conf_energies : pandas.DataFrame
        DESCRIPTION.
    output_file : TYPE
        DESCRIPTION.
    png_file : TYPE
        Filename for the plot of the energy landscape.
    ret : TYPE, optional
        Option to return data. The default is True.

    Returns
    -------
    energy_landscape : pandas.DataFrame
        The data for the free energy landscape.

    """
    energies = reweight_orientational_energies(orientational_energies.iloc[:,1:],
                                               pmf['pmf'])

    # X: com_dist
    x = pmf["com_dist"]
    # Y: conformation axis
    y = conf_energies["0"]
    # Create grid
    x, y = np.meshgrid(x, y)
    # Z: Free energy
    z = energies

    plot_energy_landscape(x, y, z, png_file=png_file, contour=20, filled=False)

    # Convert x, y, z into a dataframe
    tempor = pd.DataFrame()
    tempor['pc2'] = orientational_energies['0']
    energy_landscape = pd.DataFrame(index=tempor['pc2'])
    energy_landscape["pc2"] = conf_energies['0']
    for i in range(len(pmf.com_dist.unique())):
        energy_landscape[pmf.com_dist.unique()[i]] = z[:,i]
    return(energy_landscape)


def plot_energy_landscape(x, y, z, png_file=None, contour=20, filled=False):
    """
    Plots a contour plot of the free energy landscape.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    z : TYPE
        DESCRIPTION.
    png_file : TYPE, optional
        DESCRIPTION. The default is None.
    contour : TYPE, optional
        DESCRIPTION. The default is 20.
    filled : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    if filled:
        plt.contourf(x, y, z, contour, cmap='coolwarm');
    elif not filled:
        plt.contour(x, y, z, contour, cmap='coolwarm');
    plt.colorbar();
    plt.title("2D free energy landscape of the E9 DNase-Im2 complex")
    plt.xlabel("Conformation (PC2) a.u.")
    plt.ylabel("Scaled COM distance a.u.")
    if png_file:
        plt.savefig(png_file, dpi=600)


def main(pmf, conf_energies, probabilities, output_file, png_file, ret=False):
    """
    """
    # Get the spacing between the com distances from PMF (dF)
    spaced_com_dist = get_y_spacings(pmf)
    # 
    dist_to_coord, coord_to_dist = make_dist_to_coord(spaced_com_dist, 
                                                      pmf.com_dist.unique())
    # 
    x = np.array(conf_energies['pc2_bin'])
    
    X, Y = np.meshgrid(x, spaced_com_dist)
    energies = np.array(conf_energies[conf_energies.columns[1:]]).transpose
    
    plot_com_dist_spaced(pmf, spaced_com_dist)
    plot_energy_landscape(X, Y, energies)
    
    energies = pd.DataFrame("")
    #energy landscape




if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description=
                                     '''The script retrieves orientational RMSD data
                                        from different file and performs a principal 
                                        components analysis (PCA) on it''')
    parser.add_argument("-pmf", "--pmf", help="", required=True)
    parser.add_argument("-conf", "--conf_energies", help="", required=True)
    parser.add_argument("-prob", "--probabilities", help="", required=True)
    parser.add_argument("-o", "--output_file", help="", required=True)
    parser.add_argument("-png", "--png_file", help="", required=True)

    args = parser.parse_args()
    pmf = read_data(args.pmf)
    conf_energies = read_data(args.conf_energies)
    probabilities = read_data(args.probabilities)
    output_file = args.output_file
    png_file = args.png_file

    main(pmf, conf_energies, probabilities, output_file, png_file)
