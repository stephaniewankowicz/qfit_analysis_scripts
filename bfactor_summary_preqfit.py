#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from scipy import stats
from figure_functions import *	

colors = ["#1b9e77","#d95f02", "#7570b3","#e7298a","#66a61e", "#e6ab02", "#666666"]
sns.set_palette(sns.color_palette(colors))

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
path=os.getcwd()

all_files = glob.glob(path + "/*B_factors.csv")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[62:66]
    li.append(df)

bfactor = pd.concat(li, axis=0, ignore_index=True)
print('number of PDBs:')
print(len(bfactor.index))

#fix input
bfactor['AA'] = bfactor.AA.str.replace('[','')
bfactor['AA'] = bfactor.AA.str.replace(']','')
bfactor['Chain'] = bfactor.Chain.str.replace(']','')
bfactor['Chain'] = bfactor.Chain.str.replace('[','')
bfactor['resseq'] = bfactor.resseq.str.replace('[','')
bfactor['resseq'] = bfactor.resseq.str.replace(']','')
bfactor['chain'] = bfactor.Chain.str.replace("\'", '')
bfactor['resi'] = bfactor['resseq'].astype(int)

#summary
bfactor_summary = pd.DataFrame()
n = 1
for i in bfactor['PDB'].unique():
    tmp = bfactor[bfactor['PDB'] == i]
    bfactor_summary.loc[n, 'PDB'] = i
    bfactor_summary.loc[n, 'Average_Bfactor'] = tmp['Average_Bfactor'].mean()
    n += 1

bfactor_summary = bfactor_summary.merge(AH_key, on='PDB')


#MERGE
merged_bfactor = merge_apo_holo_df(bfactor)
print(merged_bfactor.head())

bfactor_summary.to_csv('bfactor_summary_preqfit.csv', index=False)
bfactor.to_csv('all_bfactor_preqfit.csv', index=False)


#Distribution Plot
make_dist_plot_AH(merged_bfactor['Average_Bfactor_x'], merged_bfactor['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Bound v. Unbound B-factors (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_bfactors_preqfit')

#STATS
print('Difference of Bfactor Entire Protin')
paired_wilcoxon(merged_bfactor['Average_Bfactor_x'], merged_bfactor['Average_Bfactor_y'])


#_________________________SUBSET________________________________
print('5A Subset')
all_files = glob.glob(path + "/*5.0_bfactor_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[62:66]
    li.append(df)

bfactor_subset = pd.concat(li, axis=0, ignore_index=True)
bfactor_subset = bfactor_subset.rename(columns={"Chain": "chain", "resseq":"resi"})
print(bfactor_subset.head())


bfactor_subset_m = merge_apo_holo_df(bfactor_subset)
print('number of pairs:')
print(len(bfactor_subset_m['Holo'].unique()))
#Distribution Plot
make_dist_plot_AH(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Apo v. Holo B-factors (5A)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_bfactors_5A_preqfit')
make_boxenplot_AH(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'], 'B-Factors', 'B-Factors', '', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_bfactors_5A_preqfit_boxen')

#STATS
print('Difference of Bfactor 5A')
paired_wilcoxon(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'])

#Summary
bfactor_sub_summary = pd.DataFrame()
n = 1
for i in bfactor_subset['PDB'].unique():
    #print(i)
    tmp = bfactor_subset[bfactor_subset['PDB'] == i]
    bfactor_sub_summary.loc[n, 'PDB'] = i
    bfactor_sub_summary.loc[n, 'Average_Bfactor'] = tmp['Average_Bfactor'].mean()
    n += 1

fig = plt.figure()
merged_bfactor['Difference'] = merged_bfactor['Average_Bfactor_x'] - merged_bfactor['Average_Bfactor_y']
bfactor_subset_m['Difference'] = bfactor_subset_m['Average_Bfactor_x'] - bfactor_subset_m['Average_Bfactor_y']
sns.kdeplot(merged_bfactor['Difference'], bw=0.02, label='Entire Protein')
sns.kdeplot(bfactor_subset_m['Difference'], bw=0.02, label='Binding Site')
plt.text(-75, 0.05, 'B-Factors higher in Apo') #transform=ax.transAxes, 
plt.text(25, 0.05, 'B-Factors higher in Holo')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_bfactor_difference.png')
print((merged_bfactor['Difference'].describe()))
print((bfactor_subset_m['Difference'].describe()))
bfactor_subset_m.to_csv('bfactor_subset_m.csv')
