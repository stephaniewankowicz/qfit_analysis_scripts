#packages
from __future__ import division
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
from scipy import stats
import matplotlib.patches as mpatches
from figure_functions import *	

pd.set_option('display.max_columns', None)
#REFERENCE FILE
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)


os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/')
#os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
path=os.getcwd()

# all_files = glob.glob(path + "/*rotamer_output.txt")

# li = []
# for filename in all_files:
# 	df = pd.read_csv(filename, index_col=None, header=0, sep=':')
# 	df['PDB'] = filename[54:58]
# 	li.append(df)
# rotamer = pd.concat(li, axis=0, ignore_index=True)
# print(rotamer.head())
# rotamer = rotamer[rotamer['residue']!= 'SUMMARY'].reset_index()
# split = rotamer['residue'].str.split(" ")
# #rotamer.to_csv('rotamer_all.csv')

# for i in range(0,len(rotamer.index)-1):
# 	rotamer.loc[i,'chain'] = split[i][1]
# 	STUPID = str(rotamer.loc[i,'residue'])[3:10]
# 	#print(STUPID)
# 	tmp = []
# 	try:
# 		tmp = (int(''.join(i for i in STUPID if i.isdigit())))
# 	except:
# 		newstr = ''.join((ch if ch in '0123456789.-e' else ' ') for ch in STUPID)
# 		tmp = [float(i) for i in newstr.split()]
# 	rotamer.loc[i,'resi'] = tmp #int("".join(tmp))


# rotamer.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/rotamer_all.csv')

rotamer = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/rotamer_all.csv')

#print(rotamer.head())
#ROTAMER SUBSET
all_files = glob.glob(path + "/*5.0_rotamer_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0)
    df['PDB'] = filename[54:58]
    li.append(df)

subset_rotamer = pd.concat(li, axis=0, ignore_index=True)
print(subset_rotamer.head())
subset_rotamer = subset_rotamer[subset_rotamer['PDB'].isin(AH_key['PDB'])]

print(subset_rotamer.head())

subset_rotamer['chain_resi'] = subset_rotamer['resi'].map(str) + subset_rotamer['chain']
subset_rotamer['altloc'] = subset_rotamer['residue'].astype(str).str[7]
subset_rotamer['AA'] = subset_rotamer['residue'].astype(str).str[8:11]

subset_rotamer = subset_rotamer[subset_rotamer['rotamer'] != 'OUTLIER']

AH_subset_rotamer = subset_rotamer.merge(AH_key, on='PDB')

subset_rotamer_apo = AH_subset_rotamer[AH_subset_rotamer['Apo_Holo'] == 'Apo']
subset_rotamer_holo = AH_subset_rotamer[AH_subset_rotamer['Apo_Holo'] == 'Holo']

test = subset_rotamer_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
AH_subset_rotamer_merge = test.merge(subset_rotamer_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
print(AH_subset_rotamer_merge.head())

def f(x):    
   return 'yes' if x['rotamer_x'] == x['rotamer_y'] else 'no'

AH_subset_rotamer_merge['Matching_Rotamer'] = AH_subset_rotamer_merge.apply(f, axis=1)
print(AH_subset_rotamer_merge.head())
AH_sum = pd.DataFrame()
n=0
for i in AH_subset_rotamer_merge['Holo'].unique():
	tmp = AH_subset_rotamer_merge[AH_subset_rotamer_merge['Holo'] == i]
	AH_sum.loc[n, 'Holo'] = i
	AH_sum.loc[n, 'Apo'] = str(tmp['Apo'].unique())
	AH_sum.loc[n, 'Total_Residues'] = len(tmp.index)
	AH_sum.loc[n, 'Num_Same'] = len(tmp[tmp['Matching_Rotamer'] == 'yes'].index)
	n += 1
AH_sum['Percentage_Same'] = AH_sum['Num_Same'] / AH_sum['Total_Residues']
print(AH_sum.head())

print(AH_sum['Percentage_Same'].describe())
#FIGURE
fig=plt.figure()
sns.distplot(AH_sum['Percentage_Same'], label='% Same', bins=20, kde=False)
#plt.set_title('Rotamer Status within Structures', fontsize=18)
plt.ylabel('Percentage of Residues', fontsize=18)
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/RotamerStatus_by_AH.png')


#AH_rotamer = rotamer_compare(subset_rotamer)
#print(AH_rotamer.head())
# summary_holo = rotamer_summary(subset_rotamer_holo)
# summary_apo = rotamer_summary(subset_rotamer_apo)

# print('summary_apo')
# print(summary_apo.head())

# summary_apo.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/rotamer_summary_apo.csv')
# summary_holo.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/rotamer_summary_holo.csv')
# different_holo = len(summary_holo[summary_holo['Rotamer_Status'] == 'different'].index)
# same_holo = len(summary_holo[summary_holo['Rotamer_Status'] == 'same'].index)

# apo_len = len(summary_apo.index)
# holo_len = len(summary_holo.index)
# different_apo = len(summary_apo[summary_apo['Rotamer_Status'] == 'different'].index)
# same_apo = len(summary_apo[summary_apo['Rotamer_Status'] == 'same'].index)


# #FIGURE
# fig, ax = plt.subplots()
# barwidth = 0.4

# holo = [(different_holo/holo_len), (same_holo/holo_len)]
# apo = [(different_apo/apo_len), (same_apo/apo_len)]
# print(holo, apo)

# r1 = np.arange(len(holo))

# ax.bar(r1 - barwidth/2, holo, width=barwidth, edgecolor='white', label='Holo')
# ax.bar(r1 + barwidth/2, apo, width=barwidth, edgecolor='white', label='Apo')

# #p1 = plt.bar(x[0], different, width=0.7)
# #p2 = plt.bar(x[1], same, width=0.7)
# #p3 = plt.bar(x[2], both, width=0.7)

# ax.set_title('Rotamer Status within Structures', fontsize=18)
# ax.set_ylabel('Percentage of Residues', fontsize=18)
# ax.set_xticks(r1)
# ax.legend()
# ax.set_xticklabels(('Different','Same'), fontsize=12)
# fig.tight_layout()
# plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/RotamerStatus_by_AH.png')


# AH_rotamer = rotamer_compare(subset_rotamer)
# AH_rotamer = AH_rotamer.drop_duplicates()
# print(AH_rotamer.head())

# HA_single = AH_rotamer[(AH_rotamer['Holo_altloc'] == 1) & (AH_rotamer['Apo_altloc'] == 1)]
# HA_multi = AH_rotamer[(AH_rotamer['Holo_altloc'] > 1) | (AH_rotamer['Apo_altloc'] > 1)]

# multi_summary, single_summary= rotamer_AH_summary(HA_multi, HA_single)

# AH_rotamer_summary = pd.concat([multi_summary, single_summary], axis=0)
# print(AH_rotamer_summary.head())

# #CREATE FIGURE
# different = len(AH_rotamer_summary[AH_rotamer_summary['Rotamer'] == 'Different'].index)
# same = len(AH_rotamer_summary[AH_rotamer_summary['Rotamer'] == 'Same'].index)
# both = len(AH_rotamer_summary[AH_rotamer_summary['Rotamer'] == 'Same and Different'].index)
# print(different, same, both)
# print('% Different:')
# print(different/len(AH_rotamer_summary.index))
# print('% Same:')
# print(same/len(AH_rotamer_summary.index))

# print('% Both:')
# print(both/len(AH_rotamer_summary.index))

# fig = plt.figure
# x = range(3)

# p1 = plt.bar(x[0], different, width=0.7)
# p2 = plt.bar(x[1], same, width=0.7)
# p3 = plt.bar(x[2], both, width=0.7)

# plt.title('Rotamer Status')
# plt.ylabel('Number of Residues')
# plt.xticks(x, ('All Different','All Same', 'Same & Different'))
# plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/RotamerStatus_5A'+ '.png')


# AH_rotamer_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/AH_rotamer_summary.csv', index=False)



