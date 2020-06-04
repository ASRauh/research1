#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 15:56:09 2020

@author: Arriën Symon Rauh
"""
import os
import shutil
import subprocess

print(os.path.abspath("."))


#root = "/data01/arrien/traj_data/"
#root_suzan = "/data01/ssa630/Target41
index_dir = "./ndx/"

files = ["cg.top","mdout.mdp","topol.tpr","grompp.err","state_prev.cpt","state.cpt","confout.gro","ener.edr","pullf.xvg"]#,"traj.trr","pullf.xvg"]#,"traj.xtc","md.log","ana.log","ana.err"]

def copy_files(src, dst):
    for file in files:
        shutil.copyfile(os.path.join(src, file), os.path.join(dst, file))


def convert_trajectory(source, dest, simulation, struct_num):
    print(dest+"/traj_" + simulation + ".xtc -n")
    cmd = "echo 14 | trjconv -f " + source + "/traj.xtc -o ./" + dest + "/traj_" + simulation + ".xtc -n " + index_dir + "T41_" + struct_num + "_cg.ndx"
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                     close_fds=True)
    p.wait()

for dirr in os.listdir("."):
    if os.path.isdir(dirr):
        for dist in os.listdir(dirr):
            loc_dir = os.path.join(dirr, dist)
            if os.path.isdir(loc_dir):
                # directory in suzan's 
                structure_number = dirr[3:-1]
                new_trajectory = dirr + "_" + dist[1:]
                dest = loc_dir
                # To catch both cases of subdirectories (Target41_... or T41_...)
                source1 = "/home/ssa630/xd/Target41/" + dirr + "/Target41_" + structure_number + "/" + dist + "/md_sol_prod"
                source2 = "/home/ssa630/xd/Target41/" + dirr + "/T41_" + structure_number + "/" + dist + "/md_sol_prod"
                if os.path.isdir(source1):
                    print(source1, dest)
                    #print(loc_dir)
                    #print("/home/ssa630/xd/Target41/" + dirr + "/Target41_" + dirr[:-1] + "/" + dist + "/md_sol_prod")
                    copy_files(source1, dest)
                    print("Copied files from: " + source1 + " to " + dest)
                    convert_trajectory(source1, dest, new_trajectory, structure_number)
                    print("Extracted trajectory of " + dirr + "-" + dist)
                elif os.path.isdir(source2):
                    print(source2, dest)
                    #print(loc_dir)
                    #print("/home/ssa630/xd/Target41/" + dirr + "/Target41_" + dirr[:-1] + "/" + dist + "/md_sol_prod")
                    copy_files(source2, dest)
                    print("Copied files from: " + source1 + " to " + dest)
                    convert_trajectory(source2, dest, new_trajectory, structure_number)
                    print("Extracted trajectory of " + dirr + "-" + dist)
                else:
                    print("Error for : " + dirr)
