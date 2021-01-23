Edited by Stephanie Wankowicz
#began: 2020-04-09

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
    p.add_argument("structure", type=str,
                   help="PDB-file containing structure.")
    p.add_argument("--dist", type=float, default='3.0',  help="Distance between water residue and close residues")
    p.add_argument("--pdb", help="Name of the input PDB.")
    args = p.parse_args()
    return args


def partial_occ_water(pdb, pdb_name):
        water = pdb.extract('resn', 'HOH', '==')
        water_summary = pd.DataFrame()
        water_df = pd.DataFrame()
        n=1
        water_summary.loc[n, 'PDB'] = pdb_name
        water_summary.loc[n, 'average_bfactor'] = np.average(water.b)
        water_summary.loc[n, 'average_occ'] = np.average(water.q)
        m=1
        for chain in np.unique(water.chain):
            #print(chain)
            tmp_water = water.extract('chain', chain, '==')
            for i in tmp_water.resi:
             #print(i)
             #print(tmp_water.extract('resi',i, '=='))
             water_df.loc[m, 'PDB'] = pdb_name
             water_df.loc[m, 'Chain'] = chain
             water_df.loc[m, 'Resi'] = i
             water_df.loc[m, 'bfactor'] = np.average(water.extract('resi',i, '==').b)
             water_df.loc[m, 'occ'] = np.average(water.extract('resi',i, '==').q)
             m+=1
        if len(np.unique(water.altloc)) > 1:
        #   water_summary.loc[n, 'any_altloc'] = 'yes'
        #   for i in np.unique(water.altloc):
        #       print(i)
               print(water.extract('altloc', i, '=='))
        else:
           water_summary.loc[n, 'any_altloc'] = 'no'
           altloc = 0
           print('no alt loc')
        water_summary.to_csv(pdb_name + '_water_summary.csv')
        water_df.to_csv(pdb_name + '_water_df.csv')
        return altloc
        
        
def residue_closest(structure, dist, pdb_name):
        neighbors = {}
        pdb = structure.extract('resn', 'HOH', '!=')
        #print(pdb)
        water = structure.extract('resn', 'HOH', '==')
        for chain in np.unique(water.chain):
            tmp_water = water.extract('chain', chain, '==')
            for i in tmp_water.resi:
                value = str(i) + ',' + chain
                if (len(tmp_water.extract('resi',i, '==').altloc) == 1):
                  dist_pdb = np.linalg.norm(pdb.coor - tmp_water.extract('resi',i, '==').coor, axis=1)
                  #print('np.amin(dist_pdb)')
                  #print(np.amin(dist_pdb))
                else:
                  tmp_water2 = tmp_water.extract('resi',i, '==')
                  for alt in (tmp_water2.altloc):
                      #print(alt)
                      value = value + alt
                      dist_pdb = np.linalg.norm(pdb.coor - tmp_water2.extract('altloc',alt, '==').coor, axis=1)
                      #print(np.amin(dist_pdb))
                if str(np.amin(dist_pdb)) not in neighbors:
                   neighbors[str(np.amin(dist_pdb))] = value
        with open(pdb_name + '_' + str(dist) + '_waterclosestresidue.txt', 'w') as file:
          for key,value in neighbors.items():
            #print(key)
            #print(value)
            #residue_id,chain,dist = key.split()
            file.write(value + ',' + key + "\n")
            
def residues_close(structure, dist, pdb_name):
        neighbors = {}
        pdb = structure.extract('resn','HOH' , '!=') #should this be atom?
        #print(pdb)
        water = structure.extract('resn', 'HOH', '==')
        for chain in np.unique(water.chain):
            tmp_water = water.extract('chain', chain, '==')
            for i in tmp_water.resi:
                #print(tmp_water.extract('resi',i, '==').coor)
                value = str(i) + ',' +chain
                if (len(tmp_water.extract('resi',i, '==').altloc) == 1):
                   dist_pdb = np.linalg.norm(pdb.coor - tmp_water.extract('resi', i , '==').coor, axis=1)
                   for near_residue, near_chain in zip(pdb.resi[dist_pdb < dist], pdb.chain[dist_pdb < dist]):
                    #print(near_residue)
                    #print(dist_pdb)
                    #print(dist_pdb[near_residue])
                    #print(dist_pdb[near_chain])
                    key = str(near_residue)+" "+near_chain +" "+str(dist_pdb[near_residue])
                    #print(key)
                    if key not in neighbors:
                        neighbors[key] = value
        with open(pdb_name + '_' + str(dist) + '_watercloseresidue.txt', 'w') as file:
          for key,value in neighbors.items():
            #print(key)
            #print(value)
            residue_id,chain,dist = key.split()
            file.write(value + ',' + chain + ',' + residue_id + ',' + dist + "\n")

def residues_close_partial(structure, dist, pdb_name):
        neighbors = {}
        pdb = structure.extract('resn', 'HOH', '!=')
        water = structure.extract('resn', 'HOH', '==')
        for chain in np.unique(water.chain):
              tmp_water = water.extract('chain', chain, '==')
              for i in tmp_water.resi:
                #print(tmp_water.extract('resi',i, '==').coor)
                if (tmp_water.extract('resi',i, '==').q != 1):
                  value = str(i) + ',' +chain
                  dist_pdb = np.linalg.norm(pdb.coor - tmp_water.extract('resi',i, '==').coor, axis=1)
                  for near_residue, near_chain in zip(pdb.resi[dist_pdb < dist], pdb.chain[dist_pdb < dist]):
                    key = str(near_residue)+" "+near_chain
                    if key not in neighbors:
                        neighbors[key]=value
        with open(pdb_name + '_' + str(dist) + '_partialwatercloseresidue.txt', 'w') as file:
          for key,value in neighbors.items():
            #print(key)
            #print(value)
            residue_id,chain = key.split()
            file.write(value + ',' + chain + ',' + residue_id + "\n")

def main():
    args = parse_args()
    #print(os.listdir())
    # Load structure and prepare it
    structure = Structure.fromfile(args.structure).reorder() #put H20 on the bottom
    #structure = structure.extract('record', 'HETATM')
    if not args.pdb == None:
       pdb_name = args.pdb
    else:
       pdb_name = ''
    residue_closest(structure, args.dist, pdb_name)
    residues_close(structure, args.dist, args.pdb)
    residues_close_partial(structure, args.dist, args.pdb)
