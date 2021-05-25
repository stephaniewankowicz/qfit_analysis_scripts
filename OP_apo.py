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
#reference files
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200211.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

Apo_Holo_Key = create_AH_key(AH_pairs)

#reference files
pairs = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/apo_pairs_final.csv', sep=',')
pairs = pairs.rename(columns={"PDB1": "Apo", "PDB2": "Holo", "Res1": "Apo_Res", "Res2": "Holo_Res"})
AH_pairs = pairs.drop_duplicates()

print(AH_pairs.head())
AH_key = create_AH_key(AH_pairs)

AH_key = AH_key[AH_key['PDB'].isin(Apo_Holo_Key['PDB'])]

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/')
path=os.getcwd()
all_files = glob.glob(path + "/*_methyl.out")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[54:58]
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)
order_all = order_all[order_all['PDB'].isin(AH_key['PDB'])]

print(order_all.head())
order_all['s2calc']= order_all['s2calc'].apply(lambda x: np.where(x < 0, 0, x))
order_all['s2ang']= order_all['s2ang'].apply(lambda x: np.where(x < 0, 0, x))
order_all['s2ortho']= order_all['s2ortho'].apply(lambda x: np.where(x < 0, 0, x))


#MERGE
df = order_all.merge(AH_key, on=['PDB'])
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
test=test.drop_duplicates()
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
order_all_merged = df_merged.drop_duplicates()
order_all_merged.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_merged_apo.csv')

order_all_merged['Difference'] = order_all_merged['s2calc_x'] - order_all_merged['s2calc_y']

fig = plt.figure()
sns.kdeplot((order_all_merged['s2calc_x'] - order_all_merged['s2calc_y']), label='Difference (Ap)', bw=0.02)
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/apo_difference.png')

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
sns.boxenplot(order_all_merged['s2calc_x'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
sns.boxenplot(order_all_merged['s2calc_y'], orient='v', ax=axes[1]).set(xlabel='All', ylabel='% Residues with Alt Loc')
plt.show()
