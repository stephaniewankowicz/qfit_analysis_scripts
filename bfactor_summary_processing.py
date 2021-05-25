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
import matplotlib.patches as mpatches
from figure_functions import *	

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()

all_files = glob.glob(path + "/*_qFit__B_factors.csv")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[54:58]
    li.append(df)

bfactor = pd.concat(li, axis=0, ignore_index=True)
print(bfactor.head())

#fix input
bfactor['AA'] = bfactor.AA.str.replace('[','')
bfactor['AA'] = bfactor.AA.str.replace(']','')
bfactor['Chain'] = bfactor.Chain.str.replace(']','')
bfactor['Chain'] = bfactor.Chain.str.replace('[','')
bfactor['resseq'] = bfactor.resseq.str.replace('[','')
bfactor['resseq'] = bfactor.resseq.str.replace(']','')
bfactor['chain'] = bfactor.Chain.str.replace("\'", '')
bfactor['resi'] = bfactor['resseq'].astype(int)
#bfactor['resi'].str.replace("'", '')


#summary
bfactor_summary = pd.DataFrame()
n = 1
for i in bfactor['PDB'].unique():
    tmp = bfactor[bfactor['PDB'] == i]
    bfactor_summary.loc[n, 'PDB'] = i
    bfactor_summary.loc[n, 'Average_Bfactor'] = tmp['Average_Bfactor'].mean()
    n += 1

bfactor_summary = bfactor_summary.merge(AH_key, on='PDB')
bfactor.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/bfactor.csv', index=False)
print(bfactor_summary['Average_Bfactor'].describe())
bfactor_summary[bfactor_summary['Average_Bfactor']>30].to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/highbfactor.csv')

fig = plt.figure()
sns.distplot(bfactor_summary['Average_Bfactor'], kde=False, label='Average_Bfactor')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/BFactor_Average.png')

#MERGE
merged_bfactor = merge_apo_holo_df(bfactor)
print(merged_bfactor.head())

#bfactor_summary.to_csv('bfactor_summary.csv', index=False)
#bfactor.to_csv('all_bfactor.csv', index=False)


#Distribution Plot
make_dist_plot_AH(merged_bfactor['Average_Bfactor_x'], merged_bfactor['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Bound v. Unbound B-factors (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_bfactors')

merged_bfactor['Difference'] = merged_bfactor['Average_Bfactor_x'] - merged_bfactor['Average_Bfactor_y']

fig = plt.figure()
sns.distplot(merged_bfactor['Difference'], label='Holo-Apo', bins=50, kde=False)
plt.title('Differences in B-Factors [Entire Protein]')
plt.xlabel('Residue Differences in B-Factors (Holo-Apo)')
plt.text(-600, 80000, 'Increase RMSF in Apo') 
plt.text(150, 80000, 'Increased RMSF in Holo')
plt.legend()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/Bfactor_compare_dist.png')



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
print(bfactor_subset.head())


bfactor_subset_m = merge_apo_holo_df(bfactor_subset)

bfactor_subset_m['Difference'] = bfactor_subset_m['Average_Bfactor_x'] - bfactor_subset_m['Average_Bfactor_y']

fig = plt.figure()
sns.distplot(bfactor_subset_m['Difference'], label='Holo-Apo', bins=50, kde=False)
plt.xlabel('Residue Differences in B-Factors (Holo-Apo)')
plt.title('Differences in B-Factors [Binding Site Residues]')
plt.text(-85, 1500, 'Increase RMSF in Apo') 
plt.text(30, 1500, 'Increased RMSF in Holo')
plt.legend()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/Bfactor_close_diff_dist.png')



#Distribution Plot
make_dist_plot_AH(bfactor_subset_m['Average_Bfactor_x'], bfactor_subset_m['Average_Bfactor_y'], 'B-Factors', 'Number of Residues', 'Apo v. Holo B-factors (Binding Site)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_bfactors_5A')

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

print(bfactor_sub_summary['Average_Bfactor'].describe())
bfactor_sub_summary[bfactor_sub_summary['Average_Bfactor']>20].to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/highbfactor_close.csv')
fig = plt.figure()
sns.distplot(bfactor_sub_summary['Average_Bfactor'], kde=False, label='Average_Bfactor')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/BFactor_Average_close.png')


# #LIGAND B-FACTOR
# os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
# path=os.getcwd()

# all_files = glob.glob(path + "/*_ligand_B_factors.csv")

# li = []
# for filename in all_files:
#     df = pd.read_csv(filename, index_col=None, header=0)
#     df['PDB'] = filename[47:51]
#     li.append(df)

# bfactor_ligand = pd.concat(li, axis=0, ignore_index=True)
# #print(len(all_files))

# fig = plt.figure()
# sns.distplot(bfactor_ligand['Average_Bfactor'], kde=False)
# plt.xlabel('Histogram of Average Ligand B-factors')
# plt.legend()
# plt.ylabel('Number of Structures')
# plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/Ligand_BFactor_Dist.png')

# bfactor_ligand_merged = bfactor_sub_summary.merge(bfactor_ligand, on='PDB')
# bfactor_ligand_merged['ligand_close_bfactors'] = bfactor_ligand_merged['Average_Bfactor_y']/bfactor_ligand_merged['Average_Bfactor_x']

# bfactor_ligand_merged_upper = bfactor_ligand_merged[bfactor_ligand_merged['ligand_close_bfactors']>1.4]
# bfactor_ligand_merged_lower = bfactor_ligand_merged[bfactor_ligand_merged['ligand_close_bfactors']<0.85]

# print(bfactor_ligand_merged.head())
# bfactor_ligand_merged.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/Bfactor_Ligand_Summary.csv', index=False)



#permutation test
diff_obs = np.mean(bfactor_subset_m['Average_Bfactor_x']) - np.mean(bfactor_subset_m['Average_Bfactor_y']) #get observed difference
print(diff_obs)
diff_swap = [] #holding all of the differences

for i in range(10000): #repeat 100 times
    holo = bfactor_subset.sample(frac=0.5, replace=False) #create random 'holo dataset'
    apo = bfactor_subset.sample(frac=0.5, replace=False) #create random 'apo dataset'
    diff_swap.append(np.mean(holo['Average_Bfactor']) - np.mean(apo['Average_Bfactor'])) #calculate swap values

print('Subset:')
print(sum(i > diff_obs for i in diff_swap))


#permutation test
diff_obs = np.mean(merged_bfactor['Average_Bfactor_x']) - np.mean(merged_bfactor['Average_Bfactor_y']) #get observed difference
print(diff_obs)
diff_swap = [] #holding all of the differences

for i in range(10000): #repeat 100 times
    holo = bfactor_summary.sample(frac=0.5, replace=False) #create random 'holo dataset'
    apo = bfactor_summary.sample(frac=0.5, replace=False) #create random 'apo dataset'
    diff_swap.append(np.mean(holo['Average_Bfactor']) - np.mean(apo['Average_Bfactor'])) #calculate swap values

print('Full Protein:')
print(sum(i > diff_obs for i in diff_swap))
