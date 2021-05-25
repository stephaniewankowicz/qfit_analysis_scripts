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


pd.set_option('display.max_columns', None)
PDB_bind = pd.read_excel("/Users/stephaniewankowicz/Downloads/qfit_paper/PDB_bind_refined.xlsx", skiprows=[0])
#print(PDB_bind.head())


#REFERENCE FILE
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ')
#pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()
#print(AH_pairs.head())
AH_key = create_AH_key(AH_pairs)

PDB_bind2 = PDB_bind[PDB_bind['PDB code'].isin(AH_pairs['Holo'])].reset_index()

PDB_bind2['Ki_Kd'] = PDB_bind2['Affinity Data'].str[:2]
PDB_bind2['value'] = PDB_bind2['Affinity Data'].str[3:7]
PDB_bind2['value'] = PDB_bind2['value'].str.rstrip("u")
PDB_bind2['value'] = PDB_bind2['value'].str.rstrip("n")
PDB_bind2['value'] = PDB_bind2['value'].str.rstrip("nM")
PDB_bind2['value'] = PDB_bind2['value'].str.rstrip("uM")
PDB_bind2['value'] = pd.to_numeric(PDB_bind2['value'], errors='coerce')
PDB_bind2['units'] = PDB_bind2['Affinity Data'].str[4:]
PDB_bind2['units'] = PDB_bind2['units'].str.lstrip(".")
PDB_bind2['units'] = PDB_bind2['units'].str.replace('\d+', '')
#strip all numerical values 


#print(PDB_bind2.head(20))
PDB_bind2_ki = PDB_bind2[PDB_bind2['Ki_Kd']=='Ki']
#bring in close residue summary
merged_5 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_5.csv')
merged_5['Difference'] = merged_5['s2calc_x'] - merged_5['s2calc_y']

order5_summary = pd.DataFrame()
n = 1
for i in merged_5['Holo'].unique():
    tmp = merged_5[merged_5['Holo'] == i]
    for a in tmp['Apo'].unique():
    	tmp2 = tmp[tmp['Apo'] == a]
    	order5_summary.loc[n, 'Holo'] = i
    	order5_summary.loc[n, 'Apo'] = a
    	order5_summary.loc[n, 'Average_5_diff'] = tmp2['Difference'].mean()
    	order5_summary.loc[n, 'Average_5_holo_OP'] = tmp2['s2calc_x'].mean()
    	n += 1

PDB_bind = PDB_bind2.merge(order5_summary, left_on='PDB code', right_on='Holo')
print(PDB_bind.head(20))
#PDB_bind_ki = PDB_bind[PDB_bind['Ki_Kd']=='Ki']
PDB_bind_ki_uM = PDB_bind[PDB_bind['units']=='uM'] #MAKE THIS PRETTIER
PDB_bind_ki_nM = PDB_bind[PDB_bind['units']=='nM']
print(PDB_bind_ki_nM.head())
PDB_bind_ki_nM['value'] = 1000 * PDB_bind_ki_nM['value']
PDB_bind_ki = PDB_bind_ki_nM.append(PDB_bind_ki_uM, ignore_index=True)
#PDB_bind_ki = pd.concat(PDB_bind_ki_nM, PDB_bind_ki_uM)
PDB_bind_ki['value_log']= np.log(PDB_bind_ki['value'])
#PDB_bind_ki = PDB_bind_ki[PDB_bind_ki['value']<=500]
print(len(PDB_bind_ki.index))

fig=plt.figure()
#scatter_plot_with_linear_fit(PDB_bind_ki['Average_5_holo_OP'], PDB_bind_ki['value'], slope=None, y_intercept=None, label=None, color=None)
scatter = plt.scatter(PDB_bind_ki['Average_5_holo_OP'], PDB_bind_ki['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Holo OP')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('log(ki)')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_close.png')


fig=plt.figure()
#scatter_plot_with_linear_fit(PDB_bind_ki['Average_5_holo_OP'], PDB_bind_ki['value'], slope=None, y_intercept=None, label=None, color=None)
scatter = plt.scatter(PDB_bind_ki['Average_5_diff'], PDB_bind_ki['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('log(ki)')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_difference.png')

compare = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/combined_summary.csv')
print(compare.head())
compare_bind = compare.merge(PDB_bind_ki, on='Holo')

fig=plt.figure()
#scatter_plot_with_linear_fit(PDB_bind_ki['Average_5_holo_OP'], PDB_bind_ki['value'], slope=None, y_intercept=None, label=None, color=None)
scatter = plt.scatter(compare_bind['Difference_Holo_far_5'], compare_bind['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference Holo Far 5')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('log(ki)')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_difference_far.png')

fig=plt.figure()
#scatter_plot_with_linear_fit(PDB_bind_ki['Average_5_holo_OP'], PDB_bind_ki['value'], slope=None, y_intercept=None, label=None, color=None)
scatter = plt.scatter(compare_bind['Difference_far_5'], compare_bind['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference Holo Far 5')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('log(ki)')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_far.png')

high_rot = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order5_chem_highRotBonds_PerHeavyAtom.csv')
print(compare_bind.head())

print(len(high_rot['Holo'].unique()))

high_rot_bind = compare_bind[compare_bind['Holo'].isin(high_rot['Holo'])]
print('high_rot_bind')
print(high_rot_bind.head())

fig=plt.figure()
scatter = plt.scatter(high_rot_bind['Average_5_diff_x'], high_rot_bind['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference Holo/Apo 5')
plt.legend(loc = 'upper left')
plt.ylabel('log(ki)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_rot_diff.png')

fig=plt.figure()
scatter = plt.scatter(high_rot_bind['Average_5_holo_OP_x'], high_rot_bind['value_log'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Difference Holo/Apo 5')
plt.legend(loc = 'upper left')
plt.ylabel('log(ki)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/binding_OP_rot_holo.png')

# merged_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_all.csv')
# merged_5 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_5.csv')
# merged_far = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_far.csv')

# merged_far['Difference'] = merged_far['s2calc_x_x'] - merged_far['s2calc_x_y']
# merged_5['Difference'] = merged_5['s2calc_x'] - merged_5['s2calc_y']
# merged_10['Difference'] = merged_10['s2calc_x'] - merged_10['s2calc_y']
# merged_all['Difference'] = merged_all['s2calc_x'] - merged_all['s2calc_y']
# merged_all_apo['Difference'] = merged_all_apo['s2calc_x'] - merged_all_apo['s2calc_y']

