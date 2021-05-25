#packages
import pandas as pd
import numpy as np
import os
import seaborn as sns
import sys
import matplotlib.pyplot as plt
from scipy import stats
from figure_functions import *

pd.set_option('display.max_columns', None)

#import DF
merged_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_all.csv')
order5_summary = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_summary_5.csv')
merged_5 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_5.csv')
chem = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/Chemical_Descriptors_PDB.csv')
remove = ['2q2n', '5e0r', '5tya']

order5_summary['Difference'] = order5_summary['Average_Order5_Calc_x'] - order5_summary['Average_Order5_Calc_y']

merged_5['Difference'] = merged_5['s2calc_x'] - merged_5['s2calc_y']
merged_all['Difference'] = merged_all['s2calc_x'] - merged_all['s2calc_y']

merged_all_chem = merged_all.merge(chem, on='Holo')
merged_all_chem = merged_all_chem.drop_duplicates()

merged_5_chem = merged_5.merge(chem, on='Holo')
merged_5_chem = merged_5_chem.drop_duplicates()

order5_chem = order5_summary.merge(chem, on='Holo')
order5_chem = merged_5_chem.drop_duplicates()

order5_chem = order5_chem[~order5_chem['PDB_x'].isin(remove)]
print(order5_chem.head())
#order5_chem['Average_Order5_Calc_x'] = order5_chem['Average_Order5_Calc_x'].clip(lower=0)

order5_chem['NumHAcceptors_PerHeavyAtom'] = order5_chem['NumHAcceptors']/order5_chem['HeavyAtomCount']
order5_chem['NumHDonors_PerHeavyAtom'] = order5_chem['NumHDonors']/order5_chem['HeavyAtomCount']
order5_chem['NumHTotal'] = order5_chem['NumHAcceptors'] + order5_chem['NumHDonors']
order5_chem['NumHTotal_PerHeavyAtom'] = order5_chem['NumHTotal']/order5_chem['HeavyAtomCount']
order5_chem['NumRotatableBonds_PerHeavyAtom'] = order5_chem['NumRotatableBonds']/order5_chem['HeavyAtomCount']

order5_chem_lowMolWeight = order5_chem[order5_chem['MolWeight'] <= order5_chem['MolWeight'].quantile(0.25)]
order5_chem_highMolWeight = order5_chem[order5_chem['MolWeight'] >= order5_chem['MolWeight'].quantile(0.75)]

merged_5_chem_lowRotBonds = merged_5_chem[merged_5_chem['NumRotatableBonds'] <= order5_chem['NumRotatableBonds'].quantile(0.25)]
merged_5_chem_highRotBonds = merged_5_chem[merged_5_chem['NumRotatableBonds'] >= order5_chem['NumRotatableBonds'].quantile(0.75)]

order5_chem_lowRotBonds = order5_chem[order5_chem['NumRotatableBonds'] <= order5_chem['NumRotatableBonds'].quantile(0.25)]
order5_chem_highRotBonds = order5_chem[order5_chem['NumRotatableBonds'] >= order5_chem['NumRotatableBonds'].quantile(0.75)]


random = build_random(merged_5_chem_highRotBonds, merged_all)

print('High Rotatable Bonds: Random, High')
ind_MannWhit(random['Difference'], merged_5_chem_highRotBonds['Difference'])
make_boxenplot_chem(random['Difference'], merged_5_chem_highRotBonds['Difference'], 'Control Dataset', 'Binding Site Residues (High Rotatable Bonds)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_highrot.png')

random = build_random(merged_5_chem_lowRotBonds, merged_all)
print('Low Rotatable Bonds: Random, Low')
ind_MannWhit(random['Difference'], merged_5_chem_lowRotBonds['Difference'])
make_boxenplot_chem(random['Difference'], merged_5_chem_lowRotBonds['Difference'], 'Control Dataset', 'Binding Site Residues (Low Rotatable Bonds)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_lowrot.png')


order5_chem_lowlogP = order5_chem[order5_chem['MolLogP'] <= order5_chem['MolLogP'].quantile(0.25)]
order5_chem_highlogP = order5_chem[order5_chem['MolLogP'] >= order5_chem['MolLogP'].quantile(0.75)]

merged_5_chem_lowlogP = merged_5_chem[merged_5_chem['MolLogP'] <= order5_chem['MolLogP'].quantile(0.25)]
merged_5_chem_highlogP = merged_5_chem[merged_5_chem['MolLogP'] >= order5_chem['MolLogP'].quantile(0.75)]

# random = build_random(merged_5_chem_lowlogP, merged_all)
# print('Low LogP: Random, Low')
# ind_MannWhit(random['Difference'], merged_5_chem_lowlogP['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_lowlogP['Difference'], 'Control Dataset', 'Binding Site Residues (Low LogP)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_lowlog.png')

# random = build_random(merged_5_chem_highlogP, merged_all)
# print('High LogP: Random, Low')
# ind_MannWhit(random['Difference'], merged_5_chem_highlogP['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_highlogP['Difference'], 'Control Dataset', 'Binding Site Residues (Low LogP)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_highlog.png')


order5_chem_lowRotBonds_PerHeavyAtom = order5_chem[order5_chem['NumRotatableBonds_PerHeavyAtom'] <= order5_chem['NumRotatableBonds_PerHeavyAtom'].quantile(0.25)]
order5_chem_highRotBonds_PerHeavyAtom = order5_chem[order5_chem['NumRotatableBonds_PerHeavyAtom'] >= order5_chem['NumRotatableBonds_PerHeavyAtom'].quantile(0.75)]

order5_chem_lowHdonor_PerHeavyAtom= order5_chem[order5_chem['NumHDonors_PerHeavyAtom'] <= order5_chem['NumHDonors_PerHeavyAtom'].quantile(0.25)]
order5_chem_highHdonor_PerHeavyAtom= order5_chem[order5_chem['NumHDonors_PerHeavyAtom'] >= order5_chem['NumHDonors_PerHeavyAtom'].quantile(0.75)]

order5_chem_lowHaccept_PerHeavyAtom = order5_chem[order5_chem['NumHAcceptors_PerHeavyAtom'] <= order5_chem['NumHAcceptors_PerHeavyAtom'].quantile(0.25)]
order5_chem_highHaccept_PerHeavyAtom = order5_chem[order5_chem['NumHAcceptors_PerHeavyAtom'] >= order5_chem['NumHAcceptors_PerHeavyAtom'].quantile(0.75)]

order5_chem_lowHtotal_PerHeavyAtom = order5_chem[order5_chem['NumHTotal_PerHeavyAtom'] <= order5_chem['NumHTotal_PerHeavyAtom'].quantile(0.25)]
order5_chem_highHtotal_PerHeavyAtom = order5_chem[order5_chem['NumHTotal_PerHeavyAtom'] >= order5_chem['NumHTotal_PerHeavyAtom'].quantile(0.75)]

order5_chem_lowHtotal = order5_chem[order5_chem['NumHTotal'] <= order5_chem['NumHTotal'].quantile(0.25)]
order5_chem_highHtotal= order5_chem[order5_chem['NumHTotal'] >= order5_chem['NumHTotal'].quantile(0.75)]

merged_5_chem['NumHTotal'] = merged_5_chem['NumHAcceptors'] + merged_5_chem['NumHDonors']

merged_5_chem['NumHTotal_PerHeavyAtom'] = merged_5_chem['NumHTotal']/merged_5_chem['HeavyAtomCount']
merged_5_chem_lowHtotal_PerHeavyAtom = merged_5_chem[merged_5_chem['NumHTotal_PerHeavyAtom'] <= merged_5_chem['NumHTotal_PerHeavyAtom'].quantile(0.25)]
merged_5_chem_highHtotal_PerHeavyAtom = merged_5_chem[merged_5_chem['NumHTotal_PerHeavyAtom'] >= merged_5_chem['NumHTotal_PerHeavyAtom'].quantile(0.75)]

# random = build_random(merged_5_chem_highHtotal_PerHeavyAtom, merged_all)
# print('High H per Heavy Atom: Random, High')
# ind_MannWhit(random['Difference'], merged_5_chem_highHtotal_PerHeavyAtom['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_highHtotal_PerHeavyAtom['Difference'], 'Control Dataset', 'Binding Site Residues (High H bond)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_highhbond.png')


fig=plt.figure()
plt.scatter(merged_5_chem['NumHTotal_PerHeavyAtom'], merged_5_chem['Difference'])
plt.xlabel('Number of H/Heavy Atoms')
plt.legend(loc = 'upper left')
plt.ylabel('Difference 5(Holo-Apo)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/Hbond_difference.png')


# random = build_random(merged_5_chem_lowHtotal_PerHeavyAtom, merged_all)
# print('Low H per Heavy Atom: Random, Low')
# ind_MannWhit(random['Difference'], merged_5_chem_lowHtotal_PerHeavyAtom['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_lowHtotal_PerHeavyAtom['Difference'], 'Control Dataset', 'Binding Site Residues (High H bond)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_lowhbond.png')


merged_5_chem_lowHtotal = merged_5_chem[merged_5_chem['NumHTotal'] <= merged_5_chem['NumHTotal'].quantile(0.25)]
merged_5_chem_highHtotal = merged_5_chem[merged_5_chem['NumHTotal'] >= merged_5_chem['NumHTotal'].quantile(0.75)]

# random = build_random(merged_5_chem_highHtotal, merged_all)
# print('High H: Random, High')
# ind_MannWhit(random['Difference'], merged_5_chem_highHtotal['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_highHtotal['Difference'], 'Control Dataset', 'Binding Site Residues (High H bond)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_highhbond.png')


# random = build_random(merged_5_chem_lowHtotal, merged_all)
# print('Low H: Random, Low')
# ind_MannWhit(random['Difference'], merged_5_chem_lowHtotal['Difference'])
# make_boxenplot_chem(random['Difference'], merged_5_chem_lowHtotal['Difference'], 'Control Dataset', 'Binding Site Residues (High H bond)', 'Difference in s2calc', '/Users/stephaniewankowicz/Downloads/qfit_paper/binding_control_s2calc_diff_lowhbond.png')



# fig = plt.figure()
# p1 = sns.scatterplot(order5_chem_highRotBonds_PerHeavyAtom['Average_Order5_Calc'], order5_chem_highRotBonds_PerHeavyAtom['Difference'], color='r')
# p2 = sns.scatterplot(order5_chem_lowRotBonds_PerHeavyAtom['Average_Order5_Calc'], order5_chem_lowRotBonds_PerHeavyAtom['Difference'], color='b')
# plt.xlabel('Bound OP')
# plt.ylabel('Bound-Unbound OP Difference')
# #plt.title(title)
# fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/Rotatable_Bonds_Difference_Order5.png')


# make_boxenplot_chem(order5_chem_lowlogP['Average_Order5_Calc_x'], order5_chem_highlogP['Average_Order5_Calc_x'], 
# 	'Low LogP', 'High LogP', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/logp_orderparameters')
# print('Ind ttest order5-logP, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowlogP['Average_Order5_Calc_x'], order5_chem_highlogP['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowlogP['Difference'], order5_chem_highlogP['Difference'], 
	'Low LogP', 'High LogP', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/logp_orderparameters_diff')
print('Ind ttest order5-logP Difference, column1=low, column2=high')
ind_MannWhit(order5_chem_lowlogP['Difference'], order5_chem_highlogP['Difference'])

print('ROTATABLE BONDS:')

print(order5_chem['NumRotatableBonds_PerHeavyAtom'].quantile(0.25))
print(order5_chem['NumRotatableBonds_PerHeavyAtom'].quantile(0.75))
print(order5_chem_lowRotBonds_PerHeavyAtom['NumRotatableBonds'].describe())
print(order5_chem_highRotBonds_PerHeavyAtom['NumRotatableBonds'].describe())
# make_boxenplot_chem(order5_chem_lowRotBonds_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highRotBonds_PerHeavyAtom['Average_Order5_Calc_x'], 
# 	'Low Rotable Bonds per Heavy Atom', 'High Rotable Bonds per Heavy Atom', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/rotatablebonds_peratom_orderparameters')
# print('Ind ttest order5-Rotatable Bonds (per Heavy Atom), column1=low, column2=high')
# ind_MannWhit(order5_chem_lowRotBonds_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highRotBonds_PerHeavyAtom['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowRotBonds_PerHeavyAtom['Difference'], order5_chem_highRotBonds_PerHeavyAtom['Difference'], 
	'Low Rotable Bonds per Heavy Atom', 'High Rotable Bonds per Heavy Atom', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/rotbonds_peratom_orderparameters_diff')
print('Ind ttest order5-Rotatable Bonds (Difference) per Heavy Atom, column1=low, column2=high')
ind_MannWhit(order5_chem_lowRotBonds_PerHeavyAtom['Difference'], order5_chem_highRotBonds_PerHeavyAtom['Difference'])
order5_chem_highRotBonds_PerHeavyAtom.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order5_chem_highRotBonds_PerHeavyAtom.csv')

print(order5_chem['NumRotatableBonds'].quantile(0.25))
print(order5_chem['NumRotatableBonds'].quantile(0.75))
order5_chem_lowRotBonds 
order5_chem_highRotBonds

# make_boxenplot_chem(order5_chem_lowRotBonds['Average_Order5_Calc_x'], order5_chem_highRotBonds['Average_Order5_Calc_x'], 
# 	'Low Rotable Bonds', 'High Rotable Bonds', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/rotatablebonds_orderparameters')
# print('Ind ttest order5-Rotatable Bonds (per Heavy Atom), column1=low, column2=high')
# ind_MannWhit(order5_chem_lowRotBonds['Average_Order5_Calc_x'], order5_chem_highRotBonds['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowRotBonds['Difference'], order5_chem_highRotBonds['Difference'], 
	'Low Rotable Bonds', 'High Rotable Bonds', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/highrotbonds_orderparameters_diff')
print('Ind ttest order5-Rotatable Bonds (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowRotBonds['Difference'], order5_chem_highRotBonds['Difference'])




print('MOLECULAR WEIGHT:')
print(order5_chem['MolWeight'].quantile(0.75))
print(order5_chem['MolWeight'].quantile(0.25))

# make_boxenplot_chem(order5_chem_lowMolWeight['Average_Order5_Calc_x'], order5_chem_highMolWeight['Average_Order5_Calc_x'], 
# 	'Low Molecular Weight', 'High Molecular Weight', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/MW_orderparameters')
# print('Ind ttest order5-MW, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowMolWeight['Average_Order5_Calc_x'], order5_chem_highMolWeight['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowMolWeight['Difference'], order5_chem_highMolWeight['Difference'], 
	'Low Molecular Weight', 'High Molecular Weight', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/MW_orderparameters_diff')
print('Ind ttest order5-Molecular Weight (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowMolWeight['Difference'], order5_chem_highMolWeight['Difference'])


# make_boxenplot_chem(order5_chem_lowHdonor_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHdonor_PerHeavyAtom['Average_Order5_Calc_x'], 
# 	'Low H-Bond Donors per Heavy Atom', 'High H-Bond Donors per Heavy Atom', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/Hdonor_orderparameters')
# print('Ind ttest order5-Hbond Donor, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowHdonor_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHdonor_PerHeavyAtom['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowHdonor_PerHeavyAtom['Difference'], order5_chem_highHdonor_PerHeavyAtom['Difference'], 
	'Low H-Bond Donors per Heavy Atom', 'High H-Bond Donors per Heavy Atom', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/Hdonor_orderparameters_diff')
print('Ind ttest order5-Hbonddonor (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowHdonor_PerHeavyAtom['Difference'], order5_chem_highHdonor_PerHeavyAtom['Difference'])


# make_boxenplot_chem(order5_chem_lowHaccept_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHaccept_PerHeavyAtom['Average_Order5_Calc_x'], 
# 	'Low H-Bond Accept per Heavy Atom', 'High H-Bond Accept per Heavy Atom', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/Haccept_orderparameters')
# print('Ind ttest order5-Hbond acceptor, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowHaccept_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHaccept_PerHeavyAtom['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowHaccept_PerHeavyAtom['Difference'], order5_chem_highHaccept_PerHeavyAtom['Difference'], 
	'Low H-Bond Accept per Heavy Atom', 'High H-Bond Accept per Heavy Atom', 'Scalc Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/Haccept_orderparameters_diff')
print('Ind ttest order5-Hbondaccept (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowHaccept_PerHeavyAtom['Difference'], order5_chem_highHaccept_PerHeavyAtom['Difference'])


print('Hbonds:')
print(order5_chem['NumHTotal'].quantile(0.75))
print(order5_chem['NumHTotal'].quantile(0.25))

print(order5_chem_lowHtotal_PerHeavyAtom['NumHTotal'].quantile(0.75))
print(order5_chem_lowHtotal_PerHeavyAtom['NumHTotal'].quantile(0.25))

# make_boxenplot_chem(order5_chem_lowHtotal_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHtotal_PerHeavyAtom['Average_Order5_Calc_x'], 
# 	'Low H-Bond Total per Heavy Atom', 'High H-Bond Total per Heavy Atom', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/Htotal_orderparameters_heavyatom')
# print('Ind ttest order5-Hbond total (heavy atom), column1=low, column2=high')
# ind_MannWhit(order5_chem_lowHtotal_PerHeavyAtom['Average_Order5_Calc_x'], order5_chem_highHtotal_PerHeavyAtom['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowHtotal_PerHeavyAtom['Difference'], order5_chem_highHtotal_PerHeavyAtom['Difference'], 
	'Low H-Bond Total per Heavy Atom', 'High H-Bond Total per Heavy Atom', 'Scalc Order Parameters Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/Htotal_orderparameters_diff_heavyatom')
print('Ind ttest order5-Hbond Total per heavy atom (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowHtotal_PerHeavyAtom['Difference'], order5_chem_highHtotal_PerHeavyAtom['Difference'])

order5_chem_highHtotal_PerHeavyAtom.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order5_chem_highHtotal_PerHeavyAtom.csv')
# make_boxenplot_chem(order5_chem_lowHtotal['Average_Order5_Calc_x'], order5_chem_highHtotal['Average_Order5_Calc_x'], 
# 	'Low H-Bond Total', 'High H-Bond Total', 'Scalc Order Parameters', '/Users/stephaniewankowicz/Downloads/qfit_paper/Htotal_orderparameters')
# print('Ind ttest order5-Hbond total, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowHtotal['Average_Order5_Calc_x'], order5_chem_highHtotal['Average_Order5_Calc_x'])

make_boxenplot_chem(order5_chem_lowHtotal['Difference'], order5_chem_highHtotal['Difference'], 
	'Low H-Bond Total per Heavy Atom', 'High H-Bond Total', 'Scalc Order Parameters Difference', '/Users/stephaniewankowicz/Downloads/qfit_paper/Htotal_orderparameters_diff')
print('Ind ttest order5-Hbond Total (Difference), column1=low, column2=high')
ind_MannWhit(order5_chem_lowHtotal['Difference'], order5_chem_highHtotal['Difference'])

merged_all_chem['Difference'] = merged_all_chem['scalc_x'] - merged_all_chem['scalc_y']
merged_all_chem_lowH = merged_all_chem[merged_all_chem['PDB'].isin(order5_chem_lowHtotal_PerHeavyAtom['PDB'])]
merged_all_chem_highH = merged_all_chem[merged_all_chem['PDB'].isin(order5_chem_lowHtotal_PerHeavyAtom['PDB'])]
ind_MannWhit(merged_all_chem_lowH['Difference'], merged_all_chem_highH['Difference'])


#Rotatable Bonds with Ligand B-Factor
# bfactor_ligand_merged = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/Bfactor_Ligand_Summary.csv')
# print(bfactor_ligand_merged.head())
# order5_chem_ligand_b = order5_chem.merge(bfactor_ligand_merged, left_on='PDB_x', right_on='PDB')

# order5_chem_lowRotBonds_PerHeavyAtom_ligand_b = order5_chem_ligand_b[order5_chem_ligand_b['NumRotatableBonds_PerHeavyAtom'] <= order5_chem_ligand_b['NumRotatableBonds_PerHeavyAtom'].quantile(0.25)]
# order5_chem_highRotBonds_PerHeavyAtom_ligand_b = order5_chem_ligand_b[order5_chem_ligand_b['NumRotatableBonds_PerHeavyAtom'] >= order5_chem_ligand_b['NumRotatableBonds_PerHeavyAtom'].quantile(0.75)]

# make_boxenplot_chem(order5_chem_lowRotBonds_PerHeavyAtom_ligand_b['Average_Bfactor_y'], order5_chem_highRotBonds_PerHeavyAtom_ligand_b['Average_Bfactor_y'], 
# 	'Low Rotable Bonds per Heavy Atom', 'High Rotable Bonds per Heavy Atom', 'Ligand B-factors', '/Users/stephaniewankowicz/Downloads/qfit_paper/rotatablebonds_ligandbfactors')
# print('Ind ttest ligand b factors-Rotatable Bonds, column1=low, column2=high')
# ind_MannWhit(order5_chem_lowRotBonds_PerHeavyAtom_ligand_b['Average_Bfactor_y'], order5_chem_highRotBonds_PerHeavyAtom_ligand_b['Average_Bfactor_y'])


fig=plt.figure()
plt.scatter(merged_5_chem['NumRotatableBonds'], merged_5_chem['NumHTotal'])
#scatter_plot_with_linear_fit(merged_5_chem['NumRotatableBonds'], merged_5_chem['NumHTotal'], slope=None, y_intercept=None, label=None, color=None)
#scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Num Rotatable Bonds')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('Num H Bonds')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/rot_v_hbonds.png')

fig=plt.figure()
#plt.scatter(merged_5_chem['NumRotatableBonds_PerHeavyAtom'], merged_5_chem['NumHTotal_PerHeavyAtom'])
plt.scatter(order5_chem['NumRotatableBonds_PerHeavyAtom'], order5_chem['NumHTotal_PerHeavyAtom'])
#scatter_plot_with_linear_fit(merged_5_chem['NumRotatableBonds'], merged_5_chem['NumHTotal'], slope=None, y_intercept=None, label=None, color=None)
#scatter = plt.scatter(combined_summary['Average_all_diff'], combined_summary['Average_5_diff'], alpha=0.3, cmap='tab10')#['#A7E30E', '#226363'])
plt.xlabel('Num Rotatable Bonds')
plt.legend(loc = 'upper left')
#plt.text(0.12, 0.4, 'Better Holo') #transform=ax.transAxes, 
#plt.text(0.4, 0.2, 'Better Apo')
plt.ylabel('Num H Bonds')
#plt.title('R Free Differences (Holo and Apo) (qFit Structures)')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/rot_v_hbonds_peratom.png')


#Ouput
order5_chem.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_summary_5_chem.csv')

