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
pairs = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/apo_pairs_final.csv', sep=',')
pairs = pairs.rename(columns={"PDB1": "Apo", "PDB2": "Holo", "Res1": "Apo_Res", "Res2": "Holo_Res"})
AH_pairs = pairs.drop_duplicates()

print(AH_pairs.head())
AH_key = create_AH_key(AH_pairs)

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/')
path=os.getcwd()

all_files = glob.glob(path + "/*qFit__B_factors.csv")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[54:58]
    li.append(df)

bfactor = pd.concat(li, axis=0, ignore_index=True)
print(bfactor)
print(len(bfactor.index))
bfactor[bfactor['PDB'].isin(AH_key['PDB'])]
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

#MERGE
df = bfactor.merge(AH_key, on=['PDB'])
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
test=test.drop_duplicates()
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
merged_bfactor = df_merged.drop_duplicates()
merged_bfactor.to_csv('merged_bfactor.csv')
#Distribution Plot
make_dist_plot_AH(merged_bfactor['Average_Bfactor_x'], merged_bfactor['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Bound v. Unbound B-factors (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/Apo_control_bfactors')

#STATS
print('Difference of Bfactor Entire Protin')
paired_wilcoxon(merged_bfactor['Average_Bfactor_x'], merged_bfactor['Average_Bfactor_y'])


#SUBSET
all_files = glob.glob(path + "/*_5.0_bfactor_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[54:58]
    li.append(df)

bfactor_subset = pd.concat(li, axis=0, ignore_index=True)
bfactor_subset = bfactor_subset.rename(columns={"Chain": "chain", "resseq":"resi"})
bfactor_subset[bfactor_subset['PDB'].isin(AH_key['PDB'])]
print('subset:')
print(bfactor_subset.head())
print(len(bfactor_subset.index))
df = bfactor_subset.merge(AH_key, on=['PDB'])
print(df.head())
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
print(test.head())
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
print(df_merged.head())
bfactor_subset_m = df_merged.drop_duplicates()
print(bfactor_subset_m.head())

#Distribution Plot
make_dist_plot_AH(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Apo v. Holo B-factors (5A)', '/Users/stephaniewankowicz/Downloads/qfit_paper/Apo_control_bfactors_5A')
fig = plt.figure()
bfactor_subset_m['Difference'] = bfactor_subset_m['Average_Bfactor_x'] - bfactor_subset_m['Average_Bfactor_y']
merged_bfactor['Difference'] = merged_bfactor['Average_Bfactor_x'] - merged_bfactor['Average_Bfactor_y']
sns.distplot(bfactor_subset_m['Difference'], label='All')
sns.distplot(merged_bfactor['Difference'], label='Binding Site')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/difference_b_apo.png')




#STATS
print('Difference of Bfactor 5A')
paired_wilcoxon(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'])

#Summary
bfactor_sub_summary = pd.DataFrame()
n = 1
for i in bfactor_subset['PDB'].unique():
    tmp = bfactor_subset[bfactor_subset['PDB'] == i]
    bfactor_sub_summary.loc[n, 'PDB'] = i
    bfactor_sub_summary.loc[n, 'Average_Bfactor'] = tmp['Average_Bfactor'].mean()
    n += 1


#permutation test
diff_obs = np.mean(bfactor_subset_m['Average_Bfactor_x']) - np.mean(bfactor_subset_m['Average_Bfactor_y']) #get observed difference
print(diff_obs)
diff_swap = [] #holding all of the differences

for i in range(10000): #repeat 100 times
    holo = bfactor_subset.sample(frac=0.5, replace=False) #create random 'holo dataset'
    apo = bfactor_subset.sample(frac=0.5, replace=False) #create random 'apo dataset'
    diff_swap.append(np.mean(holo['Average_Bfactor']) - np.mean(apo['Average_Bfactor'])) #calculate swap values

print('subset:')
print(sum(i > diff_obs for i in diff_swap))




#permutation test
diff_obs = np.mean(merged_bfactor['Average_Bfactor_x']) - np.mean(merged_bfactor['Average_Bfactor_y']) #get observed difference
print(diff_obs)
diff_swap = [] #holding all of the differences

for i in range(10000): #repeat 10000 times
    holo = bfactor_summary.sample(frac=0.5, replace=False) #create random 'holo dataset'
    apo = bfactor_summary.sample(frac=0.5, replace=False) #create random 'apo dataset'
    diff_swap.append(np.mean(holo['Average_Bfactor']) - np.mean(apo['Average_Bfactor'])) #calculate swap values

print('Full Protein:')
print(sum(i > diff_obs for i in diff_swap))

