
#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from scipy import stats
import matplotlib.patches as mpatches
from figure_functions import *	

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
path=os.getcwd()

#read in summary files
AH_rotamer_summary = pd.read_csv('AH_rotamer_summary.csv')
merged_order_summary_5 = pd.read_csv('merged_order_summary_5.csv')
merged_order_all = pd.read_csv('merged_order_all.csv')
close_RMSF_summary = pd.read_csv('close_RMSF_summary.csv')
#subset down to CA2
tank_apo = ['3u9h']
tank_holo = ['4buv', '4bu8', '4bu5', '3p0q', '4bu6']
tank_all = ['4buv', '4bu8', '4bu5', '3p0q', '4bu6', '3u9h']

#merged_order_all = merged_order_all[merged_order_all['PDB_y']=='3u9h']
#merged_order_all["res"] = merged_order_all["resi"].astype(str) + merged_order_all["chain"]

#read in individual files only starting with specific pdb names
#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/tank/')
path=os.getcwd()

all_files = glob.glob(path + "/*qFit_methyl.out")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[52:56]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))

order_all[order_all.s2ang < 0] = 0
order_all[order_all.s2ortho < 0] = 0
order_all[order_all.s2calc < 0] = 0

print(order_all.head())

merged_order_all = merge_apo_holo_df(order_all)
merged_order_all["res"] = merged_order_all["resi"].astype(str) + merged_order_all["chain"]
merged_order_all["Difference"] = merged_order_all["s2calc_x"] - merged_order_all["s2calc_y"]
#print(merged_order_all.head())
#merged_order_all.to_csv('merged_order_all_subset.csv')

all_files = glob.glob(path + "/*_*_residue_dist.txt")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=None)
    df['PDB'] = filename[52:56]
    df['Ligand'] = filename[57:60]
    li.append(df)

dist_all = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))

dist_all = dist_all.rename(columns={0:'Resn', 1:'Residue', 2:'AltLoc', 3:'Distance'})

dist_all_collapse = pd.DataFrame()
n=1
for i in dist_all['PDB'].unique():
    print(i)
    for lig in dist_all[dist_all['PDB'] == i]['Ligand'].unique():
        print(lig)
        tmp = dist_all[(dist_all['PDB'] == i) & (dist_all['Ligand'] == lig)]
        for res in tmp['Residue'].unique():
            #print(res)
            dist_all_collapse.loc[n, 'PDB'] = i
            dist_all_collapse.loc[n, 'res'] = res
            dist_all_collapse.loc[n,'Ligand'] = lig
            dist_all_collapse.loc[n, 'Distance'] = np.min(dist_all[(dist_all['PDB'] == i) & (dist_all['Residue'] == res) & (dist_all['Ligand'] == lig)]['Distance'])
            n += 1

dist_all_collapse = dist_all_collapse.merge(AH_key, on = ['PDB'])
#print(dist_all_collapse.head())
dist_all_collapse_holo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Holo']
dist_all_collapse_apo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Apo']
test = dist_all_collapse_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
dist_all_collapse_m = test.merge(dist_all_collapse_apo, left_on='Apo', right_on='PDB')
dist_all_collapse_m = dist_all_collapse_m.drop_duplicates()


#print('dist_all_collapse_m')
#print(dist_all_collapse_m.head())
#print('merged_order_all')
#print(merged_order_all.head())
print('merging')
dist_op = dist_all.merge(merged_order_all, left_on=['PDB','Residue'], right_on=['PDB_x','res'])
print(dist_op.head())
dist_op = dist_op.drop_duplicates()
dist_op.to_csv('dist_op.csv')

#figure
fig = plt.figure()
sns.lineplot(x=dist_op['resi'], y=dist_op['s2calc_y'], hue=dist_op['Ligand_x'])
plt.xlabel('')
plt.legend()
plt.ylabel('')
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/tank.png')


fig = plt.figure()
sns.lineplot(x=dist_op['resi'], y=dist_op['Difference'], hue=dist_op['Ligand_x'])
plt.xlabel('')
plt.legend()
plt.ylabel('')
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/tank_difference.png')


fig = plt.figure()
sns.lineplot(x=dist_op['Distance'], y=dist_op['Difference'], hue=dist_op['Ligand_x'])
plt.xlabel('')
plt.legend()
plt.ylabel('')
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/tank_difference_distance.png')




















#residue by residue

