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
    p.add_argument("ligand", type=str, help="ligand of interest")

   # Output options
    p.add_argument("--pdb", help="Name of the input PDB.")

    args = p.parse_args()
    return args

def get_occ(structure, ligand, pdb):
        Occup = pd.DataFrame()
        select = structure.extract('resn', ligand, '==')
        occ = select._q
        print(occ)
        Occup.loc[1,'PDB'] = pdb
        Occup.loc[1,'Ligand_Name'] = ligand
        Occup.loc[1,'Max_Occ'] = np.amin(occ)
        Occup.loc[1, 'Average_Occ'] = np.average(occ)
        Occup.to_csv(self.pdb + 'ligand_occupancy.csv', index=False)

def main():
    args = parse_args()
    # Load structure and prepare it
    structure = Structure.fromfile(args.structure).reorder()
    structure = structure.extract('e', 'H', '!=')
    get_occ(structure, args.ligand, ars.pdb)
