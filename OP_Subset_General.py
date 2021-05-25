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

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
#pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()

all_files = glob.glob(path + "/*qFit_methyl.out")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[54:58]
    li.append(df)
order_all = pd.concat(li, axis=0, ignore_index=True)


order_all = order_all[order_all['PDB'].isin(AH_key['PDB'])]
order_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/order_all.csv')

all_files = glob.glob(path + "/*_5.0_order_param_subset.csv")
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[54:58]
    li.append(df)
order_5 = pd.concat(li, axis=0, ignore_index=True)

order_5 = order_5[order_5['PDB'].isin(AH_key['PDB'])]

order_5['s2calc']= order_5['s2calc'].apply(lambda x: np.where(x < 0, 0, x))
order_5['s2ang']= order_5['s2ang'].apply(lambda x: np.where(x < 0, 0, x))
order_5['s2ortho']= order_5['s2ortho'].apply(lambda x: np.where(x < 0, 0, x))

all_files = glob.glob(path + "/*_10.0_order_param_subset.csv")
li = []


for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=0)
    df['PDB'] = filename[54:58]
    li.append(df)
order_10 = pd.concat(li, axis=0, ignore_index=True)
order_10 = order_10[order_10['PDB'].isin(AH_key['PDB'])]

order_10['s2calc']= order_10['s2calc'].apply(lambda x: np.where(x < 0, 0, x))
order_10['s2ang']= order_10['s2ang'].apply(lambda x: np.where(x < 0, 0, x))
order_10['s2ortho']= order_10['s2ortho'].apply(lambda x: np.where(x < 0, 0, x))


order_all_tmp = order_all.merge(AH_key, on = ['PDB'])
order_all_apo = order_all_tmp[order_all_tmp['Apo_Holo'] == 'Apo']
order_all_holo = order_all_tmp[order_all_tmp['Apo_Holo'] == 'Holo']


order_10.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_10.csv')
order_5.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_5.csv')

#MERGE
df = order_5.merge(AH_key, on=['PDB'])
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
merged_order_5 = df_merged.drop_duplicates()


merged_order_all = merge_apo_holo_df(order_all)


df = order_10.merge(AH_key, on=['PDB'])
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
merged_order_10 = df_merged.drop_duplicates()

print('5, 10, all:')
print(len(merged_order_5['Holo'].unique()))
print(len(merged_order_10['Holo'].unique()))
print(len(merged_order_all['Holo'].unique()))

merged_order_5.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_5.csv')
merged_order_10.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_10.csv')
merged_order_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv')


print('WE EXPORTED EVERYTHING!')
#All Order Parameter Distribution Plots
make_dist_plot_AH(merged_order_5['s2calc_x'], merged_order_5['s2calc_y'], '', 'Number of Residues', 'Bound/Unbound within 5A', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2calc_5A')
make_boxenplot_AH(merged_order_10['s2calc_x'], merged_order_10['s2calc_y'], '', 'Number of Residues', 'Bound/Unbound within 10A', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2calc_10A_boxen')



#__________________________STATS____________________________
print('Difference of s2calc on Side Chains with 5A between Bound/Unbound')
paired_wilcoxon(merged_order_5['s2calc_x'], merged_order_5['s2calc_y'])

print('Difference of s2calc on Side Chains with 10A between Bound/Unbound')
paired_wilcoxon(merged_order_10['s2calc_x'], merged_order_10['s2calc_y'])

#Create OP 5A Summary
order5_summary = pd.DataFrame()
n = 1
for i in order_5['PDB'].unique():
    tmp = order_5[order_5['PDB'] == i]
    order5_summary.loc[n, 'PDB'] = i
    order5_summary.loc[n, 'Average_Order5_Calc'] = tmp['s2calc'].mean()
    n += 1


order_5_summary_AH = order5_summary.merge(AH_key, left_on=['PDB'], right_on=['PDB'])
order_5_summary_holo = order_5_summary_AH[order_5_summary_AH['Apo_Holo']=='Holo']
order_5_summary_apo = order_5_summary_AH[order_5_summary_AH['Apo_Holo']=='Apo']
order_5_summary_AH = order_5_summary_AH.drop_duplicates()
print('number of pairs:')
print(len(order_5_summary_AH.index))

order_5_summary_holo.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/holo_order_summary_5.csv')
order_5_summary_apo.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/apo_order_summary_5.csv')


test = order_5_summary_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
merged_order_summary_5 = test.merge(order_5_summary_apo, left_on=['Apo'], right_on=['PDB']) 
merged_order_summary_5 = merged_order_summary_5.drop_duplicates()
merged_order_summary_5.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_summary_5.csv')
print('merged order summary 5:')

#________________________FIGURE________________________
fig = plt.figure()
sns.boxplot(x=merged_order_summary_5["Ligand"], y=merged_order_summary_5["Average_Order5_Calc_x"])

plt.title('Order Parameter by Ligand')
plt.ylabel('Order Parameter')
plt.xlabel('Ligand')
fig.savefig('Ligand_OrderParameter.png')

#____________________________Create OP Not Binding______________________
order_notbinding = pd.DataFrame()
for i in order_all['PDB'].unique():
   tmp_all = order_all[order_all['PDB'] == i]
   tmp_5 = order_5[order_5['PDB']==i] 
   merged = tmp_all.merge(tmp_5.drop_duplicates(), on=['resn','resi', 'chain', 'PDB']) 
   order_notbinding = order_notbinding.append(merged, ignore_index=True)

print('_____________order_notbinding_________________:')
df = order_notbinding.merge(AH_key, on=['PDB'])
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
merged_order_notbinding = df_merged.drop_duplicates()
print(merged_order_notbinding.head())

#____________________________Create OP Far______________________

order_far = pd.DataFrame()
for i in order_all['PDB'].unique():
   tmp_all = order_all[order_all['PDB'] == i]
   tmp_10 = order_10[order_10['PDB']==i] 
   merged = tmp_all.merge(tmp_10.drop_duplicates(), on=['resn','resi', 'chain', 'PDB'])#, 
   order_far = order_far.append(merged, ignore_index=True)


print('_____________order_far_________________:')
print(len(order_far['PDB'].unique()))
df = order_far.merge(AH_key, on=['PDB'])
#print(df.head())
df_holo = df[df['Apo_Holo'] == 'Holo']
df_apo = df[df['Apo_Holo'] == 'Apo']
#print(df_apo[df_apo['PDB']=='5e9m'])
test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
merged_order_far = df_merged.drop_duplicates()
print(merged_order_far.head())
#merged_order_10 = merge_apo_holo_df(order_10)

#order_far = order_far.merge(AH_key, left_on = ['PDB'], right_on=['PDB'])
#order_far_holo = order_far[order_far['Apo_Holo'] == 'Holo']
#order_far_apo = order_far[order_far['Apo_Holo'] == 'Apo']

#test = order_far_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
#merged_order_far = test.merge(order_far_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
#merged_order_far = merged_order_far.drop_duplicates()
merged_order_far.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_far.csv', index=False)
#merged_order_far = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_far.csv')
make_dist_plot_AH(merged_order_far['s2calc_x_x'], merged_order_far['s2calc_x_y'], 's2calc', 'Number of Residues', 'Bound v. Unbound s2calc (Further than 10A)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2calc_>10A')
#make_dist_plot_AH(merged_order_far['s2ortho_x_x'], merged_order_far['s2ortho_x_y'], 's2ortho', 'Number of Residues', 'Bound v. Unbound s2ortho (Further than 10A)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2ortho_>10A')

print('Difference of s2calc on Side Chains >10 A between Bound/Unbound [Entire Protein]')
paired_wilcoxon(merged_order_far['s2calc_x_x'], merged_order_far['s2calc_x_y'])

order_all = order_all[order_all['resn'] != 0]
plot_order = order_all.groupby(by=["resn"])["s2calc"].mean().sort_values().index


#OP by AA type
fig = plt.figure()
sns.boxplot(x=order_all["resn"], y=order_all["s2calc"], order=plot_order)

plt.title('Order Parameter by Amino Acid Type')
plt.ylabel('Order Parameter')
plt.xlabel('Amino Acids')
fig.savefig('AA_OrderParameter.png')

#merged_order_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv')
print(merged_order_all.head())
merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']

fig = plt.figure()
x = range(2)
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.boxenplot(merged_order_all['s2calc_x'], orient='v', ax=axes[0]).set(xlabel = 'Bound', ylabel = '')
p2 = sns.boxenplot(merged_order_all['s2calc_y'], orient='v', ax=axes[1]).set(xlabel = 'Unbound', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/apovholo_boxen.png')

fig = plt.figure()
x = range(2)
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.boxenplot(merged_order_5['s2calc_x'], orient='v', ax=axes[0]).set(xlabel = 'Bound', ylabel = '')
p2 = sns.boxenplot(merged_order_5['s2calc_y'], orient='v', ax=axes[1]).set(xlabel = 'Unbound', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/apovholo_boxen_order5.png')

print('Difference of s2calc on Holo V. Apo (5A)')
ind_MannWhit(merged_order_5['s2calc_x'], merged_order_far['s2calc_x_y'])

fig = plt.figure()

merged_order_far['Difference'] = merged_order_far['s2calc_x_x'] - merged_order_far['s2calc_x_y']
merged_order_5['Difference'] = merged_order_5['s2calc_x'] - merged_order_5['s2calc_y']
sns.kdeplot(merged_order_all['Difference'], label='Random', bw=0.02)
sns.kdeplot(merged_order_far['Difference'], label='>10A', bw=0.02)
sns.kdeplot(merged_order_5['Difference'], label='5A', bw=0.02)
#sns.kdeplot(merged_order_all['Difference'],label='All', bw=0.02 )
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_far_difference.png')
print('merged_order_far')
print(merged_order_far['Difference'].describe())

fig = plt.figure()
merged_order_5['Difference'] = merged_order_5['s2calc_x'] - merged_order_5['s2calc_y']
sns.kdeplot(merged_order_5['Difference'], label='Difference', bw=0.02)
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_5_difference.png')
print((merged_order_5['Difference'].describe()))

fig.show()

#merged_order_far.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_far.csv')



# #_________________________________###

# #SUBSET TO ONLY SEQUENCE UNIQUE
# print('SUBSET TO ONLY UNIQUE SEQUENCE:')
# PDB_to_keep = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/PDB_to_keep.csv')
# subset_10 = merged_order_10[merged_order_10['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# subset_5 = merged_order_5[merged_order_5['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# subset_far = merged_order_far[merged_order_far['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# print(subset_far.head())

# print('Subset Difference of s2calc on Side Chains with 5A between Bound/Unbound')
# paired_wilcoxon(subset_5['s2calc_x'], subset_5['s2calc_y'])


# print('Subset Difference of s2calc on Side Chains with 10A between Bound/Unbound')
# paired_wilcoxon(subset_10['s2calc_x'], subset_10['s2calc_y'])

# print('Subset Difference of s2calc on Side Chains Far between Bound/Unbound')
# paired_wilcoxon(subset_far['s2calc_x_x'], subset_far['s2calc_x_y'])


# make_dist_plot_AH(merged_order_10['s2calc_x'], merged_order_10['s2calc_y'], '', 'Number of Residues', 'Bound/Unbound within 5A', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2calc_5A_subset')

# make_boxenplot_AH(merged_order_5['s2calc_x'], merged_order_5['s2calc_y'], '', 'Number of Residues', 'Bound/Unbound within 10A', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2calc_10A_boxen_subset')

