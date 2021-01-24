import os
import sys
import time
import copy
import numpy as np
import pandas as pd
from argparse import ArgumentParser
from qfit.structure import Structure


def parse_args():
    p = ArgumentParser(description=__doc__)
    p.add_argument("structure", type=str,
                   help="PDB-file containing structure.")
    p.add_argument("ligand", type=str, help="ligand of interest")

   # Output options
    p.add_argument("--pdb", help="Name of the input PDB.")

    args = p.parse_args()
    return args

def get_bfactors(ligand, pdb_name, ligand_name):
        res_b = []
        altloc = []
        for atom in ligand.name:
            print(atom)
            res_b.append(ligand.extract('name', atom, '==').b)
            res_q.append(ligand.extract('name', atom, '==').q)
        for chain in np.unique(ligand.chain):
          if len(np.unique(ligand.extract('chain', chain, '==').altloc)) > 1:
             altloc.append('yes')
          else:
             altloc.append('no')
        B_factor['altloc'] = None
        B_factor['altloc'] = B_factor['altloc'].astype(object)
        B_factor.loc[n, 'altloc'] = [altloc]
        B_factor.loc[n, 'num_atoms'] = len(ligand.q)
        B_factor.loc[n, 'Occ'] = np.array2string(ligand.q, separator=',')
        B_factor.loc[n, 'Average_Occ'] = np.average(ligand.q)
        B_factor.loc[n, 'PDB'] = pdb_name
        B_factor.loc[n, 'Ligand_Name'] = ligand_name
        B_factor.loc[n, 'Max_Bfactor'] = np.amax(res_b)
        B_factor.loc[n, 'Average_Bfactor'] = np.average(res_b)
        n += 1
        B_factor = pd.DataFrame()
        B_factor.to_csv(pdb_name + '_ligand_B_factors.csv', index=False)

def main():
    args = parse_args()
    # Load structure and prepare it
    structure = Structure.fromfile(args.structure).reorder()
    structure = structure.extract('e', 'H', '!=')
    structure = structure.extract('resn', args.ligand, '==')
    if not args.pdb == None:
       pdb_name = args.pdb
    else:
       pdb_name = ''
    get_bfactors(structure, pdb_name, args.ligand)
