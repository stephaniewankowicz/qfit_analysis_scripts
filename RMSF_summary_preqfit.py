from __future__ import division
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from figure_functions import *


#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

#READ IN FILES
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
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


RMSF_merged = merge_apo_holo_df(RMSF)

make_dist_plot_AH(RMSF_merged['RMSF_x'], RMSF_merged['RMSF_y'],'RMSF Entire Protein', 'Number of Residues', 'RMSF across entrie protein', '/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_preqfit', )
print('Difference of RMSF on Side Chains between Bound/Unbound [Entire Protein]')
paired_wilcoxon(RMSF_merged['RMSF_x'], RMSF_merged['RMSF_y'])


# RMSF_summary = pd.DataFrame()
# n = 1
# for i in RMSF['PDB'].unique():
#     tmp = RMSF[RMSF['PDB'] == i]
#     RMSF_summary.loc[n, 'PDB'] = i
#     RMSF_summary.loc[n, 'Num_Residues'] = len(tmp.index)
#     RMSF_summary.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['RMSF']>0].index)
#     if tmp.RMSF.ge(0).any() == True:
#         RMSF_summary.loc[n, 'Alt_Loc'] = 1
#     else:
#         RMSF_summary.loc[n, 'Alt_Loc'] = 0
#     #RMSF_summary.loc[n, 'Apo_Holo'] = tmp['Apo_Holo'].unique()
#     RMSF_summary.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
#     n += 1

# RMSF_summary['per_altloc'] = RMSF_summary['Num_Alt_Loc'] / RMSF_summary['Num_Residues']

#RMSF SUBSET
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
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

print(merged_close_RMSF.head())

print(merged_close_RMSF.columns)
print(merged_close_RMSF['RMSF_x'])
make_dist_plot_AH(merged_close_RMSF['RMSF_x'], merged_close_RMSF['RMSF_y'],'RMSF Entire Protein', 'Number of Residues', 'RMSF across entrie protein', '/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_preqfit_5A')

#STATS
print('Difference of RMSF between Holo/Apo Pre qFit [5A]')
paired_wilcoxon(merged_close_RMSF['RMSF_x'], merged_close_RMSF['RMSF_y'])


# close_RMSF_summary = pd.DataFrame()
# n = 1
# for i in close_RMSF['PDB'].unique():
#     tmp = close_RMSF[close_RMSF['PDB'] == i]
#     close_RMSF_summary.loc[n, 'PDB'] = i
#     close_RMSF_summary.loc[n, 'Num_Residues'] = len(tmp.index)
#     close_RMSF_summary.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['RMSF']>0].index)
#     if tmp.RMSF.ge(0).any() == True:
#         close_RMSF_summary.loc[n, 'Alt_Loc'] = 1
#     else:
#         close_RMSF_summary.loc[n, 'Alt_Loc'] = 0
#     #close_RMSF_summary.loc[n, 'Apo_Holo'] = tmp['Apo_Holo'].unique()
#     close_RMSF_summary.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
#     n += 1

# close_RMSF_summary['per_altloc'] = close_RMSF_summary['Num_Alt_Loc'] / close_RMSF_summary['Num_Residues']

#close_RMSF_summary_merge = merge_apo_holo_df(close_RMSF_summary)
#close_RMSF_summary_merge.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/close_RMSF_summary.csv')


# close_RMSF_summary = close_RMSF_summary.merge(AH_key)
# close_RMSF_sum_holo = close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Holo']
# close_RMSF_sum_apo = close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']
# test = close_RMSF_sum_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
# merged_close_sum_RMSF = test.merge(close_RMSF_sum_apo, left_on='Apo', right_on='PDB')

# merged_close_sum_RMSF['Percent_Holo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_x']/ merged_close_sum_RMSF['Num_Residues_x']
# merged_close_sum_RMSF['Percent_Apo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_y']/ merged_close_sum_RMSF['Num_Residues_y']
# merged_close_sum_RMSF['Apo_Holo_Multi_Diff'] = merged_close_sum_RMSF['Percent_Holo_Close'] - merged_close_sum_RMSF['Percent_Apo_Close']
# merged_close_sum_RMSF.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/close_RMSF_summary.csv')
