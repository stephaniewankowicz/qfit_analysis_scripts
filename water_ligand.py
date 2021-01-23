#Edited by Stephanie Wankowicz
#began: 2020-04-27

import pkg_resources  # part of setuptools
import os.path
import os
import sys
import time
import copy
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from . import Structure
from .structure.base_structure import _BaseStructure


def parse_args():
    p = ArgumentParser(description=__doc__)
    p.add_argument("--dist", type=float, default='3.0',  help="Distance between water residue and close residues")
    p.add_argument("-l", "--ligand", help="Name of Ligand of interest")
    p.add_argument("holo_structure", type=str,
                   help="PDB-file containing structure.")
    p.add_argument("apo_structure", type=str,
                   help="PDB-file containing structure.") #this should be a superposed PDB files
    p.add_argument("holo_pdb_name", type=str, help="name of Holo PDB.")
    p.add_argument("apo_pdb_name", type=str, help="name of Holo PDB.")
    p.add_argument("-write", type=str, help="output PDB file of subset residues")
    p.add_argument("-dir", type=str, help="directory of output.")
    args = p.parse_args()
    return args


def select_resi(ligand_pdb, water, dist):
    neighbors = {}
    for coor in ligand_pdb.coor:
        value = 0
        dist_pdb = np.linalg.norm(water.coor - coor, axis=1)
        #print(dist_pdb)
        for near_residue, near_chain in zip(water.resi[dist_pdb < dist], water.chain[dist_pdb < dist]):
               key = str(near_residue)+" "+near_chain #+" "+str(dist_pdb[near_residue])
               if key not in neighbors:
                  neighbors[key] = value
    return neighbors
    

def ligand_close(holo_structure, apo_structure, dist, ligand, holo_pdb_name, apo_pdb_name):
        pdb = holo_structure.extract('resn', ligand, '==')
        water = holo_structure.extract('resn', 'HOH', '==')
        holo_neighbors = select_resi(pdb, water, dist)
        water = apo_structure.extract('resn', 'HOH', '==')
        apo_neighbors = select_resi(pdb, water, dist)
        print(holo_neighbors)
        print(apo_neighbors)
        with open(apo_pdb_name + '_' + ligand + '_' + str(dist) + '_watercloseresidue.txt', 'w') as file:
          for key,value in apo_neighbors.items():
            residue_id, chain = key.split()
            file.write(chain + ',' + str(residue_id) + "\n")

        with open(holo_pdb_name + '_' + ligand + '_' + str(dist) + '_watercloseresidue.txt', 'w') as file:
          for key,value in holo_neighbors.items():
            residue_id, chain = key.split()
            file.write(chain + ',' + str(residue_id)  + "\n")

def main():
    args = parse_args()
    if not args.dir == None:
       dir = args.dir
    else:
       dir = ''
    try:
        holo_structure = Structure.fromfile(args.holo_structure)
        apo_structure = Structure.fromfile(args.apo_structure)
    except:
        print('I need the apo and holo PBDs!')
        return
    # Load structure and prepare it
    ligand_close(holo_structure, apo_structure, args.dist, args.ligand, args.apo_pdb_name, args.holo_pdb_name)
