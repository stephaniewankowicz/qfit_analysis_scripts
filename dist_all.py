import pandas as pd
from figure_functions import *	
import numpy as np
import glob
from scipy.interpolate import spline

CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
path=os.getcwd()


all_files = glob.glob(path + "/*_*_residue_dist.txt")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=None)
    df['PDB'] = filename[47:51]
    df['Ligand'] = filename[52:55]
    li.append(df)

dist_all = pd.concat(li, axis=0, ignore_index=True)
print(len(all_files))

dist_all = dist_all.rename(columns={0:'Resn', 1:'Residue', 2:'AltLoc', 3:'Distance'})
print(dist_all.head())
dist_all = dist_all[dist_all['PDB'].isin(CDK2)]


dist_all_collapse = pd.DataFrame()
n=1
for i in dist_all['PDB'].unique():
    #print(i)
    for lig in dist_all[dist_all['PDB'] == i]['Ligand'].unique():
        #print(lig)
        tmp = dist_all[(dist_all['PDB'] == i) & (dist_all['Ligand'] == lig)] 
        for res in tmp['Residue'].unique():
            #print(res)
            dist_all_collapse.loc[n, 'PDB'] = i
            dist_all_collapse.loc[n, 'res'] = res
            dist_all_collapse.loc[n, 'Ligand'] = lig
            dist_all_collapse.loc[n, 'Distance'] = np.min(dist_all[(dist_all['PDB'] == i) & (dist_all['Residue'] == res) & (dist_all['Ligand'] == lig)]['Distance'])
            n += 1
        
#dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Dropbox/Fraser_Rotation/qfit_output/200225/dist_all_collapse.csv')
dist_all_collapse.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_CDK2.csv')
print(len(dist_all_collapse.index))
dist_all_collapse = dist_all_collapse.merge(AH_key, on = ['PDB'])
dist_all_collapse = dist_all_collapse.drop_duplicates()
print(dist_all_collapse.head())
print(len(dist_all_collapse.index))

dist_all_collapse_holo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Holo']
dist_all_collapse_apo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Apo']
print(len(dist_all_collapse_holo.index))
print(len(dist_all_collapse_apo.index))
#test = dist_all_collapse_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
#dist_all_collapse_m = test.merge(dist_all_collapse_apo, left_on='Apo', right_on='PDB')
#dist_all_collapse_m = dist_all_collapse_m.drop_duplicates()

merged_order_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all_CDK2.csv')
merged_order_all = merged_order_all[merged_order_all['PDB_x'].isin(CDK2)]

merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']

apo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_apo.csv')
holo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_holo.csv')

print(len(apo.index))
print(len(holo.index))
apo = apo.dropna(subset=['resi'])
holo = holo.dropna(subset=['resi'])
#print(apo[apo['resi'].isnull()])


holo['resi'] = holo['resi'].astype(int)
holo['res'] = holo['resi'].astype(str) + holo['chain']

apo['resi'] = apo['resi'].astype(int)
apo['res'] = apo['resi'].astype(str) + apo['chain']

print(len(apo.index))
print(apo.head())
merged_order_all['resi'] = merged_order_all['resi'].astype(int)
merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

print(dist_all_collapse_apo.head())
dist_holo_order = dist_all_collapse_holo.merge(holo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])
dist_apo_order = dist_all_collapse_apo.merge(apo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])

dist_apo_order = dist_apo_order.drop_duplicates()
print(dist_apo_order.head(1000))
CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']
print(len(dist_apo_order.index))

CDK2_holo = dist_holo_order[dist_holo_order['PDB'].isin(CDK2)]
CDK2_apo = dist_apo_order[dist_apo_order['PDB']=='1pw2']
CDK2_apo = CDK2_apo.drop_duplicates()
print(CDK2_apo['PDB'].unique())
print(len(CDK2_apo.index))
print(CDK2_holo['PDB'].unique())

print(len(CDK2_holo.index))

CDK2_apo['PDB_lig'] = CDK2_apo['PDB'] + CDK2_apo['Ligand']
CDK2_holo['PDB_lig'] = CDK2_holo['PDB'] + CDK2_holo['Ligand']
CDK2_apo = CDK2_apo.drop_duplicates()


CDK2_merge = CDK2_apo.merge(CDK2_holo, on=['Ligand', 'res'])
CDK2_merge['Difference'] = CDK2_merge['s2calc_x'] - CDK2_merge['s2calc_y']
#CDK2_holo.to_csv('CDK2_holo.csv')
#CDK2_apo.to_csv('CDK2_apo.csv')
print(len(CDK2_apo.index))
#CDK2_apo = CDK2_apo.dropna(subset=['resi', 's2calc_y'])
#CDK2_holo = CDK2_holo.dropna(subset=['resi', 's2calc_x'])

fig = plt.figure()
sns.lineplot(x=CDK2_holo['Distance'], y=CDK2_holo['s2calc'], hue=CDK2_holo['PDB_lig'])
plt.show()

fig = plt.figure()
sns.lineplot(x=CDK2_apo['Distance'], y=CDK2_apo['s2calc'], hue=CDK2_apo['PDB_lig'])
plt.show()

fig = plt.figure()
sns.lineplot(x=CDK2_merge['resi_x'], y=CDK2_merge['Difference'], hue=CDK2_merge['Ligand'])
plt.show()

CDK2_merge.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/CDK2_merge.csv', index=False)

f, ax = plt.subplots(figsize=(8, 8))
piv = pd.pivot_table(CDK2_merge, values=["Difference"], index=['Ligand'], columns=['resi_x'])
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.heatmap(piv, square=True, cbar=False)
plt.show()

#CDK2_holo = CDK2_holo.dropna(subset=['Distance', 's2calc_x'])
#xnew = np.linspace(CDK2_holo['s2calc_x'].min(), CDK2_holo['s2calc_x'].max(), 100)
#power_smooth = spline(CDK2_holo['Distance'], CDK2_holo['s2calc_x'], xnew)

#plt.plot(xnew,power_smooth)
#plt.show()


