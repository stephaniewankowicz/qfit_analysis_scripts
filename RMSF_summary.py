from __future__ import division
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from figure_functions import *

#Figure Colors
colors = ["#1b9e77","#d95f02", "#7570b3","#e7298a","#66a61e", "#e6ab02", "#666666"]
sns.set_palette(sns.color_palette(colors))


#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

#___________________________READ IN FILES_______________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()
all_files = glob.glob(path + "/*qfit_RMSF.csv")
li = []
for filename in all_files:
    df = pd.read_csv(filename, sep=',')
    li.append(df)
print(len(all_files))
RMSF = pd.concat(li, axis=0, ignore_index=True)
RMSF['PDB'] = RMSF['PDB_name'].astype(str).str[0:4]
RMSF = RMSF.rename(columns={"Chain": "chain", "resseq":"resi"})
#print(RMSF.head())

RMSF.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_summary.csv', index=False)

RMSF_merged = merge_apo_holo_df(RMSF)
print(RMSF_merged.head())
#RMSF_merged['rmsf_y'] = RMSF_merged['rmsf_y'].clip(upper=1.5)
#RMSF_merged['rmsf_x'] = RMSF_merged['rmsf_x'].clip(upper=1.5)
RMSF_merged = RMSF_merged.dropna(subset=['RMSF_x', 'RMSF_y'])

RMSF_merged['Difference'] = RMSF_merged['RMSF_x'] - RMSF_merged['RMSF_y']

fig = plt.figure()
sns.distplot(RMSF_merged['Difference'], label='Holo-Apo', bins=50, kde=False)
plt.text(-2, 40000, 'Increase RMSF in Apo') 
plt.text(1, 40000, 'Increased RMSF in Holo')
plt.xlabel('Difference in RMSF by Residues (Holo-Apo)')
#plt.legend()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_diff_dist.png')

print('Whole Protien Difference:')
print(RMSF_merged['Difference'].describe())

make_dist_plot_AH(RMSF_merged['RMSF_x'], RMSF_merged['RMSF_y'],'RMSF Entire Protein', 'Number of Residues', 'RMSF across entrie protein', '/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF.png')
print('Difference of RMSF on Side Chains between Bound/Unbound [Entire Protein]')
paired_wilcoxon(RMSF_merged['RMSF_x'], RMSF_merged['RMSF_y'])


#___________________________________RMSF SUBSET POST QFIT_________________________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()

all_files = glob.glob(path + "/*5.0_rmsf_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, header=0)
    li.append(df)

close_RMSF = pd.concat(li, axis=0, ignore_index=True)
close_RMSF['PDB'] = close_RMSF['PDB_name'].astype(str).str[0:4]
close_RMSF = close_RMSF.rename(columns={"Chain": "chain", "resseq":"resi"})
print(close_RMSF.head())
merged_close_RMSF = merge_apo_holo_df(close_RMSF)


make_dist_plot_AH(merged_close_RMSF['RMSF_x'], merged_close_RMSF['RMSF_y'], 'RMSF', 'Number of Residues', 'Bound v. Unbound RMSF (Close Residues)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_RMSF_5A')

merged_close_RMSF['Difference'] = merged_close_RMSF['RMSF_x'] - merged_close_RMSF['RMSF_y']

fig = plt.figure()
sns.distplot(merged_close_RMSF['Difference'], label='Holo-Apo', bins=50, kde=False)
plt.text(-2, 2500, 'Increase RMSF in Apo') 
plt.text(1, 2500, 'Increased RMSF in Holo')
plt.xlabel('Difference in RMSF by Residue in Binding Site Residues (Holo-Apo)')
#plt.legend()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_close_RMSF_dist.png')

print('Binding Site Residue Difference:')
print(merged_close_RMSF['Difference'].describe())

#STATS
print('Difference of RMSF between Holo/Apo qFit [5A]')
paired_wilcoxon(merged_close_RMSF['RMSF_x'], merged_close_RMSF['RMSF_y'])

