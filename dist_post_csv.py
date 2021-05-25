import pandas as pd
from figure_functions import *	
import numpy as np
import glob
#from scipy.interpolate import spline

os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200209.txt', sep=' ')
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
path=os.getcwd()

# apo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_apo.csv')
# holo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_holo.csv')
# dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')
# print(len(dist_all_collapse.index))
# dist_all_collapse = dist_all_collapse.merge(AH_key, on = ['PDB'])
# dist_all_collapse = dist_all_collapse.drop_duplicates()
# print(dist_all_collapse.head())
# print(len(dist_all_collapse.index))

#dist_all_collapse.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')
dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')
dist_all_collapse = dist_all_collapse.drop(['Unnamed: 0', 'Unnamed: 0.1', 'Unnamed: 0.1.1'], axis=1)
print(dist_all_collapse.head())
dist_all_collapse = dist_all_collapse.drop_duplicates()
dist_all_collapse_holo = dist_all_collapse[dist_all_collapse['Apo_Holo_x'] == 'Holo']
dist_all_collapse_apo = dist_all_collapse[dist_all_collapse['Apo_Holo_x'] == 'Apo']
print(len(dist_all_collapse_holo.index))
print(len(dist_all_collapse_apo.index))
test = dist_all_collapse_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
test = test.drop_duplicates()
print(test.head())
dist_all_collapse_m = test.merge(dist_all_collapse_apo, left_on='Apo', right_on='PDB')
dist_all_collapse_m = dist_all_collapse_m.drop_duplicates()
print('dist_all_collapse_m')
print(dist_all_collapse_m.head())
dist_all_collapse_m.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_merged.csv')
merged_order_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv')
merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']

#merge OP with distance

print(len(apo.index))
print(len(holo.index))
apo = apo.dropna(subset=['resi'])
holo = holo.dropna(subset=['resi'])


holo['resi'] = holo['resi'].astype(int)
holo['res'] = holo['resi'].astype(str) + holo['chain']

apo['resi'] = apo['resi'].astype(int)
apo['res'] = apo['resi'].astype(str) + apo['chain']

print(len(apo.index))
print(apo.head())
merged_order_all['resi'] = merged_order_all['resi'].astype(int)
merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

dist_holo_order = dist_all_collapse_holo.merge(holo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])
dist_apo_order = dist_all_collapse_apo.merge(apo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])
dist_holo_order = dist_holo_order.drop_duplicates()
dist_apo_order = dist_apo_order.drop_duplicates()
dist_apo_order.to_csv('dist_apo_order.csv')
dist_holo_order.to_csv('dist_holo_order.csv')
#test = pd.concat(dist_apo_order, dist_holo_order)
print('test')
#print(test.head())
print('dist_apo_order')
print(dist_apo_order.head())

fig = plt.figure()
sns.relplot(data=dist_apo_order, kind="line",
    x="Distance", y="s2calc", col="PDB",
    hue="PDB", style="Apo_Holo_x",
)
plt.show()
