#!/bin/bash

if [ -f /etc/profile.d/modules.sh ]
then
        . /etc/profile.d/modules.sh
        module load null default-ethernet
        module load gromacs/4.0.5
        module load python
        module load modeller
fi

export PYTHONPATH="/usr/local/Cluster-Apps/modeller/modlib/:/usr/local/Cluster-Apps/modeller/lib/x86_64-intel8/python2.5/"

export LD_LIBRARY_PATH="/usr/local/Cluster-Apps/modeller/lib/x86_64-intel8/"

cd $SGE_O_WORKDIR

python ./get_traj_backbone.py
