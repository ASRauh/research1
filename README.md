# research1
Author: Arriën Symon Rauh
## Introduction
This directory contains all the scripts used to perform the analysis CG MARTINI pulling simulations.

## List of contents
<ul>
  <li>calculate_rmsd.py: </li>
  <li>automate_calculate_rmsd.py: </li>
  <li>get_traj_backbone.py</li>
  <li>get_traj_backbone.sh</li>
  <li>analysis.py: </li>
  <li>pca_rmsd.py: </li>
  <li>bin_orientations.py: </li>
  <li>calculate_free_energy.py: </li>
  <li>constraint_force_integration.py: </li>
  <li>energy_landscape.py: </li>
</ul>

## Description of contents
### Calculation of the orientational RMSDs
The calculation of the orientational RMSDs was automated in two python scripts. The first script, calculate_rmsd.py, automates the calculations of the orientational RMSDs for one starting structure, using <strong>gmx rms</strong>, iterating over all COM distances. The second script, <strong style="font-size:20px">automate_calculate_rmsd.py</strong>, automates the first script enabling efficient calculation for all starting structures in one go. 
For these scripts to work the following files are needed: first, trajectory files (<i>traj.xtc</i>) containing only the backbone coordinates, extracted from the simulation data using <strong>gmx trjconv</strong>. Second, the topology/coordinate files (<i>pdb/gro</i>) of the backbones of the references. And finally, a tailor-made index file (<i>index.ndx</i>) that works for all starting structures and references. Overall, these scripts produce a file that contains the RMSD data, the mean force, the COM distance, the ID of the starting structure, replicate number and frame number for all the frames of all the simulations in a tab-separated format. (For extraction of the backbone coordinates: <strong>get_traj_backbone.py</strong> & <strong>get_traj_backbone.sh</strong>

### Further analysis
To perform the further analysis steps (PMF calculation, PCA, Calculating orientational energies and Calculation of the free energy landscape) a set of python scripts has been written. The workflow is summarised in <strong>analysis.py</strong>.
<strong>pca_rmsd.py</strong> is used for the principal component analysis. <strong>bin_orientations.py</strong> and <strong>calculate_free_energy.py</strong> are used for binning the orientational axis and calculating the orientational free energies. <strong>constraint_force_integration.py</strong> is used to obtain the potential of mean force through constraint force integration. 
Finally, <strong>energy_landscape.py</strong> is used to reweight the orientational free energies using the potential of mean force and to ultimately obtain the 2D free energy landscape. All scripts have been designed to in principle function as either script or module
