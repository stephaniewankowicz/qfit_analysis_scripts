#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from scipy import stats
from figure_functions import *
from rdkit import Chem	
from rdkit.Chem import Descriptors

os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/')
path = os.getcwd()

all_files = glob.glob(path + "/*_ligand_overlap.txt")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=None)
    df['PDB'] = filename[54:58]
    li.append(df)

lig_overlap = pd.concat(li, axis=0, ignore_index=True)
lig_overlap = lig_overlap.rename(columns={0: "Lig1", 1: "Lig2", 2: "Lig3", 3:"Lig4"})
print(lig_overlap.head())

os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ')
#pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()
print(len(AH_pairs.index))
AH_pairs_true_apo = AH_pairs[~AH_pairs['Apo'].isin(lig_overlap['PDB'])]
print(len(AH_pairs_true_apo.index))
AH_pairs_true_apo.to_csv('ligand_supplementary_table1_QCed_updated_200422_trueapo.txt', sep=' ', index=False)

print(AH_pairs.head())

#lig_overlap.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/PDB_names_lig_overlap.csv')
# lig1 = lig_overlap['Lig1'].values.tolist()
# all_ligs = lig1 + lig_overlap['Lig2'].values.tolist()
# all_ligs = all_ligs + lig_overlap['Lig3'].values.tolist()
# print(set(all_ligs))

all_files = glob.glob("/Users/stephaniewankowicz/Downloads/overlapping_lig_sdf/*sdf")
overlap_lig = pd.DataFrame()

n = 1
for filename in all_files:
	overlap_lig.loc[n, 'Ligand_Name'] = filename[56:59]
	suppl = Chem.SDMolSupplier(filename)
	for mol in suppl:
		overlap_lig.loc[n, 'Num_Atoms'] = mol.GetNumAtoms()
		overlap_lig.loc[n, 'Molecular_Weight'] = Descriptors.ExactMolWt(mol)
	n += 1
print(overlap_lig.head())


overlap_lig.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/overlap_lig.csv')


#AH_pairs.to_csv('ligand_supplementary_table1_QCed_updated_200421_manual3.txt', sep=' ')