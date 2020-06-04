#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 10:30:00 2019

@author: Arriën Symon Rauh

USAGE:
python pca_final-all.py -df ... -od ... -nc [int] -lab [id, com_dist, origin] -dist [int] 
python pca_final-all.py -df ... -od ... -nc [int] -lab [id, com_dist, origin] -dist [int] -std
"""
import os
import argparse
import pickle
import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pca_plotting as pp


def pca_rmsd(data, results_dir, label, distance, n_components, std=False):
    """
    Performs principal component analysis 

    Parameters
    ----------
    data : pandas.DataFrame
        Full dataset with orientational RMSDs.
    results_dir : string
        Directory to store the results in.
    label : string
        DESCRIPTION.
    distance : integer/boolean
        Option to add distance in the PCA.
    n_components : integer
        Number of principal compontents to calculate.
    std : boolean, optional
        Option to standardise the data. The default is False.

    Returns
    -------
    new_data : pandas.DataFrame
        Full dataset with principal components added.

    """
    columns = data.columns

    #
    filename = os.path.split(data_file)[-1].replace(".txt","")
    
    if distance != 0:
        # Extract the RMSD values
        rmsd_data = data.loc[:,columns[:6]]
        filename += "_distance"
        rmsd_data["COM_dist"] = rmsd_data.loc[:,"COM_dist"] * float(distance) #.apply(lambda x: x * distance)
    else:
        # Extract the RMSD values
        rmsd_data = data.loc[:,columns[:5]]


    if std:
        # Standardise the rmsd data
        rmsd_data_std =  StandardScaler().fit_transform(rmsd_data)
    else:
        rmsd_data_std = rmsd_data

    # PCA column names
    pca_labels = [F"pc{i+1}" for i in range(n_components)]

    # Perform PCA
    pca_rmsd = PCA(n_components=n_components)
    principalComponents_rmsd = pca_rmsd.fit_transform(rmsd_data_std)
    principal_rmsd = pd.DataFrame(data = principalComponents_rmsd,
                                  columns = pca_labels)
    
    #################
    ##### XRAY ######
    #################
    xray_df = pd.read_csv("/Users/Arrien/Documents/Universiteit/master/semester_3-4/research_1/data/xray_rmsd_10.txt", sep='\t')
    xray_pca = pd.DataFrame(data = pca_rmsd.transform(xray_df), columns = pca_labels)
    ######
    

    # Pickle the sklearn PCA object
    with open(os.path.join(results_dir, F"pca_{filename}.pkl"), "wb") as pf:
        pickle.dump(pca_rmsd, pf)    

    
    ### Extract the loadings ###
    loadings = pca_rmsd.components_

    ### Output the PCA data (scores and loaings) ###
    # Scores
    principal_rmsd.to_csv(os.path.join(results_dir, F"pca_scores_{filename}.txt"), sep='\t',
                      index=False, float_format='%.7f')
    # Loadings
    load_labels = columns[:5]
    if distance != 0: #label == "COM_dist":
        load_labels.append(label)
    pd.DataFrame(loadings, columns=load_labels).to_csv(os.path.join(results_dir,
                                                     F"pca_loadings_{filename}.txt"),
                                                     sep='\t', index=False, float_format='%.7f')

    # Output new dataset with PC scores added
    new_data = pd.concat([data.iloc[:,:5], principal_rmsd, data.iloc[:,5:]], axis=1)
    new_data.to_csv(os.path.join(results_dir, F"new_data_{filename}.txt"),
                    sep='\t', index=False, float_format='%.7f')
    
    #######################
    #####  Plotting  ######
    #######################
    
    
    
    
    #######################
    return new_data



def main(data_file, results_dir, label, distance, n_components, std=False):
    """
    

    Parameters
    ----------
    data_file : TYPE
        DESCRIPTION.
    results_dir : TYPE
        DESCRIPTION.
    label : TYPE
        DESCRIPTION.
    distance : TYPE
        DESCRIPTION.
    n_components : TYPE
        DESCRIPTION.
    std : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """
    # Read data from file into pandas dataframe
    data = pd.read_csv(data_file, sep='\t')
    #
    pca_rmsd(data, results_dir, label, distance, n_components, std=False)




if __name__ == "__main__":
    # Parse command line arguments
    parser = argparse.ArgumentParser(description='''The script retrieves orientational RMSD data
                                                    from different file and performs a principal 
                                                    components analysis (PCA) on it''')
    parser.add_argument("-df", "--data_file", help="", required=True)
    parser.add_argument("-od", "--results_dir", help="", required=False)
    parser.add_argument("-nc", "--n_components", help="The number of PCs", required=False)
    parser.add_argument("-lab", "--label", help="Pick what label to use in the plot {id or com_dist or origin}", required=True)
    parser.add_argument("-dist", "--distance", dest="distance", help='''For including the COM distant restraint 
                        as a variable.''', required=False)
    parser.add_argument("-std", "--std", dest="std", action="store_true",
                        help='''Option to standardise the data before performing the PCA''', required=False)

    args = parser.parse_args()

    data_file = args.data_file
    results_dir = args.results_dir
    n_components = int(args.n_components)
    label = args.label
    distance = float(args.distance)
    std = args.std
    
    
    # Check if the directories specified under data_dir and results_dir exist.
    # If not make them:
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    if not n_components:
        n_components = 2
    if not distance:
        distance = 1

    # Call main()
    main(data_file, results_dir, label, distance, n_components, std)