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
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()
AH_key = create_AH_key(AH_pairs)

pd.set_option('display.max_columns', None)

os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()

all_files = glob.glob(path + "/*_methyl.out") #read in full protein files

li = []
pdb_remove =[]

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[54:58]
    if len(df.index) < 30:
    	print(filename[54:58])
    	print(len(df.index))
    	pdb_remove.append(filename[54:58])
    li.append(df)

order_all = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))

order_all = order_all[order_all['PDB'].isin(AH_key['PDB'])]
order_all = order_all[~order_all['PDB'].isin(pdb_remove)]

order_all['s2calc']= order_all['s2calc'].apply(lambda x: np.where(x < 0, 0, x))
order_all['s2ang']= order_all['s2ang'].apply(lambda x: np.where(x < 0, 0, x))
order_all['s2ortho']= order_all['s2ortho'].apply(lambda x: np.where(x < 0, 0, x))


print(order_all.head())
#MERGE
merged_order_all = merge_apo_holo_df(order_all)

print('merged_order_all')
print(merged_order_all.head())

order_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all.csv')
merged_order_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv', index=False)

#SUBSET
merged_order_all_polar = merged_order_all[merged_order_all['resn_x'].isin(['R','N','D','C','Q','E', 'H', 'K', 'S', 'T', 'Y'])]
merged_order_all_nonpolar = merged_order_all[merged_order_all['resn_x'].isin(['V','W','P','F','M','L','I','G','A'])]

#_________All Order Parameter Distribution Plots_________
#make_dist_plot_AH(merged_order_all['s2calc_x'], merged_order_all['s2calc_y'], 's2calc', 'Number of Residues', 'Bound/Unbound OP (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/OP_all_s2calc')
#make_dist_plot_AH(merged_order_all['s2ortho_x'], merged_order_all['s2ortho_y'], 's2ortho', 'Number of Residues', 'Bound/Unbound s2ortho (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/OP_all_s2ortho')
#make_dist_plot_AH(merged_order_all['s2ang_x'], merged_order_all['s2ang_y'], 's2ang', 'Number of Residues', 'Bound/Unbound s2ang (Entire Protein)', '/Users/stephaniewankowicz/Downloads/qfit_paper/OP_all_s2ang')

#_________All Order Parameter A/H Comparison Plot_________
make_boxenplot_AH(merged_order_all['s2calc_x'], merged_order_all['s2calc_y'], 's2calc', 'Number of Residues', 's2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/OP_calc_boxen')


print('Difference of s2calc between Bound/Unbound [Entire Protein]')
paired_wilcoxon(merged_order_all['s2calc_x'], merged_order_all['s2calc_y'])

#__________STATS on Nonpolar_________
print('Difference of s2calc on only Polar Side Chains between Bound/Unbound [Entire Protein]')
paired_wilcoxon(merged_order_all_polar['s2calc_x'], merged_order_all_polar['s2calc_y'])

print('Difference of s2calc on only nonpolar Side Chains between Bound/Unbound [Entire Protein]')
paired_wilcoxon(merged_order_all_nonpolar['s2calc_x'], merged_order_all_nonpolar['s2calc_y'])

print('number of pairs:')
print(len(merged_order_all['Holo'].unique()))

#FIGURE 
fig = plt.figure()
f, axes = plt.subplots(1, 4, sharey=True, sharex=True)

p1 = sns.boxenplot(merged_order_all_polar['s2calc_x'], orient='v', 
ax=axes[0]).set(xlabel='Polar Bound', ylabel='S2calc Order Parameter')
p2 = sns.boxenplot(merged_order_all_polar['s2calc_y'], orient='v', ax=axes[1]).set(xlabel='Polar Unbound', ylabel='')
p3 = sns.boxenplot(merged_order_all_nonpolar['s2calc_x'], orient='v', ax=axes[2]).set(xlabel='NonPolar Bound', ylabel='')
p4 = sns.boxenplot(merged_order_all_nonpolar['s2calc_y'], orient='v', ax=axes[3]).set(xlabel='NonPolar Unbound', ylabel='')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/FullProtein_OP_byResidueType.png')

fig = plt.figure()
merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']
sns.kdeplot(merged_order_all['Difference'], label='Difference (Bound-Unbound)', bw=0.02)
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/difference.png')


#SUBSET TO SOLVENT EXPOSED RESIDUES
#dictionary of max ACC (Tien et al 2013, PLOS ONE)
# RES_MAX_ACC = {'A': 129.0, 'R': 274.0, 'N': 195.0, 'D': 193.0, \
#                'C': 167.0, 'Q': 225.0, 'E': 223.0, 'G': 104.0, \
#                'H': 224.0, 'I': 197.0, 'L': 201.0, 'K': 236.0, \
#                'M': 224.0, 'F': 240.0, 'P': 159.0, 'S': 155.0, \
#                'T': 172.0, 'W': 285.0, 'Y': 263.0, 'V': 174.0} 
# #read in SASA DF
# all_files = glob.glob(path + "/*_qFit_sasa.csv")

# li = []

# for filename in all_files:
#     df = pd.read_csv(filename, sep=',', header=0)
#     df['PDB'] = filename[54:58]
#     if len(df.index) < 30:
#     	print(filename[54:58])
#     	print(len(df.index))
#     li.append(df)
# sasa = pd.concat(li, axis=0, ignore_index=True)

# print(sasa.head())

# # def get_rsa(amino_acid):
# # 	if amino_acid in RES_MAX_ACC:
# #         max_acc = RES_MAX_ACC[amino_acid]  # Find max SA for residue
# #         rsa = float(solvent_acc) / max_acc # Normalize SA
# #         return rsa

# for index, row in sasa.iterrows():
#  	sasa.loc[index,'RSA'] = float(row['ss'])/RES_MAX_ACC[row['aa']]

# sasa['Solvent_Exposed'] = sasa['RSA'].apply(lambda x: np.where(x < 0.20, 'Not Solved Exposed', 'Solvently Exposed'))

# sasa.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/sasa_all.csv', index=False)
# sasa = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/sasa_all.csv')
# print(sasa.head())
# #divide acc column by total surface area of residue
# #selection 20% or 0.2 as the cut off value for being solvently exposed.
# #SUBSET TO ONLY SEQUENCE UNIQUE
# print('SUBSET TO ONLY UNIQUE SEQUENCE:')
# PDB_to_keep = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/PDB_to_keep.csv')
# subset = merged_order_all[merged_order_all['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# make_boxenplot_AH(subset['s2calc_x'], subset['s2calc_y'], 's2calc', 'Number of Residues', 's2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/OP_calc_boxen_subset')

# print('Difference of s2calc between Bound/Unbound [Entire Protein]')
# paired_wilcoxon(subset['s2calc_x'], subset['s2calc_y'])


# print(order_all[order_all['PDB']=='5glm'])

order_all_sum = pd.DataFrame()
n = 1
for i in order_all['PDB'].unique():
  order_all_sum.loc[n,'PDB'] = i
  order_all_sum.loc[n, 'Min_OP'] = order_all[order_all['PDB']==i]['s2calc'].min()
  order_all_sum.loc[n, 'Max_OP'] = order_all[order_all['PDB']==i]['s2calc'].max()
  order_all_sum.loc[n, 'Quartile1'] = order_all[order_all['PDB']==i]['s2calc'].quantile(0.25)
  order_all_sum.loc[n, 'Quartile3'] = order_all[order_all['PDB']==i]['s2calc'].quantile(0.75)
  order_all_sum.loc[n, 'Median'] = order_all[order_all['PDB']==i]['s2calc'].median()
  n += 1 

order_all_sum.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_sum.csv') 
order_all_PDB = order_all['PDB'].sample(n=100)
order_all_ra = order_all[order_all['PDB'].isin(order_all_PDB)]
order_all_melt = pd.melt(order_all_ra, id_vars=['PDB'], value_vars=['s2calc'])
print(order_all_melt.head())
order_all_melt = order_all_melt.dropna()


# #figure
fig = plt.figure()
ax = sns.boxplot(x=order_all_melt['PDB'], y=order_all_melt['value']) #hue=CaM_melt['variable']
# #ax = sns.lineplot(x=CaM['residue'], y=CaM['s2calc_NMR'], linewidth=2)
plt.xlabel('PDB')
# plt.legend()
plt.ylabel('Order Parameter (Full Protein)')
# labels = ax.axes.get_xticklabels()
#ax.axes.set_xticklabels(labels, rotation=45)
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/all_PDBs_OP.png')

fig = plt.figure()
ax = sns.boxplot(x=order_all_melt['PDB'], y=order_all_melt['value']) #hue=CaM_melt['variable']
# #ax = sns.lineplot(x=CaM['residue'], y=CaM['s2calc_NMR'], linewidth=2)
plt.xlabel('PDB')
# plt.legend()
plt.ylabel('Order Parameter (Full Protein)')
# labels = ax.axes.get_xticklabels()
#ax.axes.set_xticklabels(labels, rotation=45)
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/all_PDBs_OP.png')
