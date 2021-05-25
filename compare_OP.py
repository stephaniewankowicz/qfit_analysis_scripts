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
from scipy.stats import hypergeom as hg
from sklearn.metrics import r2_score

merged_all_apo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_merged_apo.csv')
merged_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv')
merged_5 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_5.csv')
merged_10 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_10.csv')
#fields = ['s2calc_x_x', 's2calc_x_y','Apo', 'Holo'] , usecols=fields
merged_far = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_far.csv')

merged_far['Difference'] = merged_far['s2calc_x_x'] - merged_far['s2calc_x_y']
merged_5['Difference'] = merged_5['s2calc_x'] - merged_5['s2calc_y']
merged_10['Difference'] = merged_10['s2calc_x'] - merged_10['s2calc_y']
merged_all['Difference'] = merged_all['s2calc_x'] - merged_all['s2calc_y']
merged_all_apo['Difference'] = merged_all_apo['s2calc_x'] - merged_all_apo['s2calc_y']


merged_5.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_diff_5.csv')
li = []

#create random solvently exposed
sasa = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/sasa_all.csv')
sasa_exposed = sasa[sasa['Solvent_Exposed']=='Solvently Exposed']
sasa_notexposed = sasa[sasa['Solvent_Exposed']=='Not Solved Exposed']

merged_all_solvent = merged_all[merged_all['resi'].isin(sasa_exposed['resnum']) & merged_all['chain'].isin(sasa_exposed['chain']) & merged_all['resn_x'].isin(sasa_exposed['aa'])]
merged_not_solvent = merged_all[(~merged_all['resi'].isin(sasa_notexposed['resnum'])) & (~merged_all['chain'].isin(sasa_notexposed['chain'])) & (~merged_all['resn_x'].isin(sasa_notexposed['aa']))]
print(len(merged_not_solvent.index))
print(merged_not_solvent.head())
print('merged far')
print(merged_far.head())
print(len(merged_far.index))
far_solvent = merged_far[merged_far['resi'].isin(sasa_exposed['resnum']) & merged_all['chain'].isin(sasa_exposed['chain']) & merged_all['resn_x'].isin(sasa_exposed['aa'])]
far_not_solvent = merged_far[~merged_far['resi'].isin(sasa_exposed['resnum']) & merged_all['chain'].isin(sasa_exposed['chain']) & merged_all['resn_x'].isin(sasa_notexposed['aa'])]
print(len(far_not_solvent.index))
#print(far_not_solvent.head())
# print('all:')
# print(merged_all.head())

# print('close:')
# print(merged_5.head())
AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 
'Y']

for i in merged_all['Holo'].unique():
    tmp = merged_all_solvent[merged_all_solvent['Holo'] == i]
    tmp_close = merged_5[merged_5['Holo'] == i]
    for a in AA: #look at each amino acid
        num = len(tmp_close[tmp_close['resn_x'] == a]) #determine the number of AA in 
        if num > 0: #only select if that type of AA is located in a close residue
            tmp = merged_all_solvent[merged_all_solvent['resn_x'] == a].sample(n=num) #choose random residues corresponding to the type and number of AA found in binding site
            li.append(tmp)
random_solvent = pd.concat(li, axis=0, ignore_index=True)

for i in merged_all['Holo'].unique():
    tmp = merged_all[merged_all['Holo'] == i]
    tmp_close = merged_5[merged_5['Holo'] == i]
    for a in AA: #look at each amino acid
        num = len(tmp_close[tmp_close['resn_x'] == a]) #determine the number of AA in 
        if num > 0: #only select if that type of AA is located in a close residue
            tmp = merged_all[merged_all['resn_x'] == a].sample(n=num) #choose random residues corresponding to the type and number of AA found in binding site
            li.append(tmp)

random = pd.concat(li, axis=0, ignore_index=True)

print('random:')
print(random.head())
print(len(random.index))
print(random['Difference'].describe())

#HYPER TEST
#pvalue = hg.cdf(int(len_intersection)-1, int(pop_size), int(len_genes1), int(len_genes2))
print('random solvent v. merged 5')
ind_MannWhit(random_solvent['Difference'], merged_5['Difference'])
make_boxenplot_chem(random_solvent['Difference'], merged_5['Difference'], 'Control Dataset', 'Binding Site Residues', 'Difference in OP (Holo-Apo)', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_solvent_s2calc_diff')

print('random v. merged 5')
ind_MannWhit(random['Difference'], merged_5['Difference'])
make_boxenplot_chem(random['Difference'], merged_5['Difference'], 'Control Dataset', 'Binding Site Residues', 'Difference in OP (Holo-Apo)', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff')

print('not solvent v. merged 5')
ind_MannWhit(merged_not_solvent['Difference'], merged_5['Difference'])
make_boxenplot_chem(merged_not_solvent['Difference'], merged_5['Difference'], 'Control Dataset', 'Binding Site Residues', 'Difference in OP (Holo-Apo)', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_notsolvent')



print('random v. far')
#ind_MannWhit(random['Difference'], far_solvent['Difference'])
#make_boxenplot_chem(random['Difference'], far_solvent['Difference'], 'Control Dataset', 'Far Residues', 'Difference in OP (Holo-Apo)', '/Users/stephaniewankowicz/Downloads/qfit_paper/far_control_s2calc_diff.png')

print('random v. far (not solvent)')
ind_MannWhit(random['Difference'], far_not_solvent['Difference'])
make_boxenplot_chem(random['Difference'], far_not_solvent['Difference'], 'Control Dataset', 'Far Residues (Not Solvently Exposed)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/far_notsolvent_control_s2calc_diff')

print('apo v. A/H')
ind_MannWhit(merged_all_apo['Difference'], merged_all['Difference'])
make_boxenplot_chem(merged_all_apo['Difference'], merged_all['Difference'], 'Apo-Apo', 'Holo-Apo', 'Difference in OP', '/Users/stephaniewankowicz/Downloads/qfit_paper/apo_AH_s2calc_diff')

#COMPARE PAIRED CLOSE_FAR
print('number of holo structures:')
print(len(merged_5['Holo'].unique()))
print(len(merged_all['Holo'].unique()))
print(len(merged_far['Holo'].unique()))


merged_all = merged_all[~merged_all['Holo'].isin(['2q2n', '5e0r', '6gpc', '5tya'])]
merged_5 = merged_5[~merged_5['Holo'].isin(['2q2n', '5e0r', '6gpc', '5tya'])]
merged_far = merged_far[~merged_far['Holo'].isin(['2q2n', '5e0r', '6gpc', '5tya'])]

order5_summary = pd.DataFrame()
n = 1
for i in merged_all['Holo'].unique():
    tmp = merged_5[merged_5['Holo'] == i]
    for a in tmp['Apo'].unique():
    	tmp2 = tmp[tmp['Apo'] == a]
    	order5_summary.loc[n, 'Holo'] = i
    	order5_summary.loc[n, 'Apo'] = a
    	order5_summary.loc[n, 'Average_5_diff'] = tmp2['Difference'].mean()
    	order5_summary.loc[n, 'Average_5_holo_OP'] = tmp2['s2calc_x'].mean()
    	n += 1

print(n)
print('order5_summary')
print(len(order5_summary.index))
orderall_summary = pd.DataFrame()
n = 1
for i in merged_all['Holo'].unique():
    tmp = merged_all[merged_all['Holo'] == i]
    for a in tmp['Apo'].unique():
    	tmp2 = tmp[tmp['Apo'] == a]
    	orderall_summary.loc[n, 'Holo'] = i
    	orderall_summary.loc[n, 'Apo'] = a
    	orderall_summary.loc[n, 'Average_all_diff'] = tmp2['Difference'].mean()
    	orderall_summary.loc[n, 'Average_all_holo_OP'] = tmp2['s2calc_x'].mean()
    	n += 1

print('orderall_summary')
print(len(orderall_summary.index))

orderfar_summary = pd.DataFrame()
n = 1
for i in merged_far['Holo'].unique():
    tmp = merged_far[merged_far['Holo'] == i]
    for a in tmp['Apo'].unique():
    	tmp2 = tmp[tmp['Apo'] == a]
    	orderfar_summary.loc[n, 'Holo'] = i
    	orderfar_summary.loc[n, 'Apo'] = a
    	orderfar_summary.loc[n, 'Average_far_diff'] = tmp2['Difference'].mean()
    	orderfar_summary.loc[n, 'Average_far_holo_OP'] = tmp2['s2calc_x_x'].mean()
    	n += 1

orderfar_solvent_sum = pd.DataFrame()
n = 1
for i in far_not_solvent['Holo'].unique():
    tmp = far_not_solvent[far_not_solvent['Holo'] == i]
    for a in tmp['Apo'].unique():
        tmp2 = tmp[tmp['Apo'] == a]
        orderfar_solvent_sum.loc[n, 'Holo'] = i
        orderfar_solvent_sum.loc[n, 'Apo'] = a
        orderfar_solvent_sum.loc[n, 'Average_far_diff'] = tmp2['Difference'].mean()
        orderfar_solvent_sum.loc[n, 'Average_far_holo_OP'] = tmp2['s2calc_x_x'].mean()
        n += 1
print(orderfar_solvent_sum.head())

orderall_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/orderall_summary.csv')
order5_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order5_summary.csv')
orderfar_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/orderfar_summary.csv')

order5_summary.merge(far_not_solvent, on=['Apo', 'Holo'])

#not_sol_far_5 = orderfar_solvent_sum.merge(order5_summary, on=['Apo', 'Holo'])
#not_sol_far_5 = not_sol_far_5.drop_duplicates()
#not_sol_far_5['Difference_far_5'] = not_sol_far_5['Average_far_diff'] - not_sol_far_5['Average_5_diff']
#not_sol_far_5['Difference_Holo_far_5'] = not_sol_far_5['Average_far_holo_OP'] - not_sol_far_5['Average_5_holo_OP']

#print(not_sol_far_5.head())

combined_summary = orderall_summary.merge(order5_summary, on=['Apo', 'Holo'])
combined_summary1 = combined_summary.drop_duplicates()
combined_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/combined_summary1.csv')
combined_summary = combined_summary.drop_duplicates()
combined_summary = combined_summary.merge(orderfar_summary, on=['Apo', 'Holo'])
combined_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/combined_summary2.csv')

combined_summary = combined_summary.drop_duplicates()
combined_summary['Difference_all_5'] = combined_summary['Average_all_diff'] - combined_summary['Average_5_diff']
combined_summary['Difference_far_5'] = combined_summary['Average_far_diff'] - combined_summary['Average_5_diff']
combined_summary['Difference_Holo_far_5'] = combined_summary['Average_far_holo_OP'] - combined_summary['Average_5_holo_OP']
combined_summary['Difference_Holo_all_5'] = combined_summary['Average_all_holo_OP'] - combined_summary['Average_5_holo_OP']


combined_summary = combined_summary[combined_summary['Apo'] !='1uj4']
combined_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/combined_summary.csv', index=False)
# combined_summary['Difference_all_5'].hist()
# combined_summary['Average_all_diff'].hist()
# combined_summary['Average_5_diff'].hist()

#combined_summary.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/combined_summary.csv')
fig=plt.figure()
scatter_plot_with_linear_fit(combined_summary['Average_all_holo_OP'], combined_summary['Average_5_diff'], slope=None, y_intercept=None, label=None, color=None)
#scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Average Difference All Residues')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('Average Difference Binding Residues')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/all_5_diff.png')



# fig=plt.figure()
# scatter_plot_with_linear_fit(not_sol_far_5['Difference_Holo_far_5'], not_sol_far_5['Average_5_diff'], slope=None, y_intercept=None, label=None, color=None)
# print(r2_score(not_sol_far_5['Difference_Holo_far_5'], not_sol_far_5['Average_5_diff']))
# #scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
# plt.xlabel('Difference in Holo OP (Far-Binding Residues)')
# plt.legend(loc = 'upper left')
# #plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
# #plt.text(0.4, 0.2, 'Better Apo')
# plt.ylabel('Difference in OP Binding Residues (Holo-Apo)')
# #plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
# fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/far_5_diff_notsolvent.png')

print(combined_summary[combined_summary['Apo']=='1mb0'])
fig=plt.figure()
scatter_plot_with_linear_fit(combined_summary['Difference_Holo_far_5'], combined_summary['Average_5_diff'], slope=None, y_intercept=None, label=None, color=None)
print(r2_score(combined_summary['Difference_Holo_far_5'], combined_summary['Average_5_diff']))
#scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference in Holo OP (Far-Binding Residues)')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('Difference in OP Binding Residues (Holo-Apo)')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/far_5_diff.png')


# fig=plt.figure()
# #scatter_plot_with_linear_fit(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], slope=1, y_intercept=0, label='Rfree', color=None)
# scatter = plt.scatter(combined_summary['Average_far_diff'], combined_summary['Average_5_diff'], alpha=0.3, c='red')#['#A7E30E', '#226363'])
# scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='blue')#['#A7E30E', '#226363'])
# plt.xlabel('Average Difference Far Residues')
# plt.legend(loc = 'upper left')
# #plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
# #plt.text(0.4, 0.2, 'Better Apo')
# plt.ylabel('Average Difference Binding Residues')
# #plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
# fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/far_5_diff.png')


fig = plt.figure()
sns.distplot(combined_summary['Difference_all_5'], label='All-Binding', bins=20)
sns.distplot(combined_summary['Difference_far_5'], label='Far-Binding', bins=20)
#sns.distplot(merged_far['Difference'], label='Far')
#sns.distplot(merged_5['Difference'], label='Binding Site')
#sns.distplot(merged_10['Difference'], label='10A')
#sns.distplot(merged_all_apo['Difference'], label='Apo-Apo', bins=50, kde=False)
plt.legend()
#plt.xlim(-0.5,0.5)

plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/all_v_5.png')



print('all v. 5')
paired_wilcoxon(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'])


fig = plt.figure()
sns.kdeplot(random['Difference'], label='Random', bw=0.05)
#sns.kdeplot(random_solvent_exposed['Difference'], label='Random Solvently Exposed', bw=0.05)
sns.kdeplot(merged_all['Difference'], label='All', bw=0.05)
#sns.kdeplot(merged_far['Difference'], label='Far', bw=0.05)
#sns.kdeplot(merged_5['Difference'], label='Binding Site', bw=0.05)
#sns.kdeplot(merged_10['Difference'], label='10A', bw=0.05)
sns.kdeplot(merged_all_apo['Difference'], label='Apo', bw=0.1)
#plt.xlim(-0.15,0.15)
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/OP_compare.png')

fig = plt.figure()
sns.distplot(merged_all['Difference'], label='Holo-Apo', bins=50, kde=False)
#sns.distplot(merged_far['Difference'], label='Far')
#sns.distplot(merged_5['Difference'], label='Binding Site')
#sns.distplot(merged_10['Difference'], label='10A')
sns.distplot(merged_all_apo['Difference'], label='Apo-Apo', bins=50, kde=False)
plt.legend()
#plt.xlim(-0.5,0.5)

plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/OP_compare_dist.png')

