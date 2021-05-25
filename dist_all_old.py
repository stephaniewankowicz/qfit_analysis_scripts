import pandas as pd
from figure_functions import *	
import numpy as np
import glob
import scipy
from scipy.interpolate import interp1d


os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/')
path=os.getcwd()


all_files = glob.glob(path + "/*_*_residue_dist.txt")

li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=None)
    df['PDB'] = filename[54:58]
    df['Ligand'] = filename[59:62]
    li.append(df)

dist_all = pd.concat(li, axis=0, ignore_index=True)
print(dist_all.head())
dist_all = dist_all.rename(columns={0:'Resn', 1:'Residue', 2:'AltLoc', 3:'Distance'})
dist_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all.csv')

#dist_all_collapse = pd.DataFrame()
#n=1
#dist = []
#for i in dist_all['PDB'].unique():
#    print(i)
#    for lig in dist_all[dist_all['PDB'] == i]['Ligand'].unique():
#        #print(lig)
#        tmp = dist_all[(dist_all['PDB'] == i) & (dist_all['Ligand'] == lig)] 
        #print(tmp.head())
#        for res in tmp['Residue'].unique():
#            dist.append(tuple((i, res, lig, np.min(tmp[tmp['Residue']==res]['Distance']))))

#df = pd.DataFrame(dist, columns =['PDB', 'Residue', 'Ligand', 'Distance']) 
#print(df.head())


#df.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')
dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')
#print(len(dist_all_collapse.index))
dist_all_collapse = dist_all_collapse.merge(AH_key, on = ['PDB'])
dist_all_collapse = dist_all_collapse.drop_duplicates()
print(dist_all_collapse.head())

dist_all_collapse_holo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Holo']
dist_all_collapse_apo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Apo']

#test = dist_all_collapse_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
#dist_all_collapse_m = test.merge(dist_all_collapse_apo, left_on='Apo', right_on='PDB')
#dist_all_collapse_m = dist_all_collapse_m.drop_duplicates()

merged_order_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_all.csv')
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
dist_holo_order = dist_all_collapse_holo.merge(holo, left_on=['PDB', 'Residue'], right_on=['PDB', 'res'])
dist_apo_order = dist_all_collapse_apo.merge(apo, left_on=['PDB', 'Residue'], right_on=['PDB', 'res'])

dist_apo_order = dist_apo_order.drop_duplicates()
test = dist_holo_order.merge(AH_pairs, left_on='PDB', right_on='Holo')
df_merged = test.merge(dist_apo_order, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
df_merged = df_merged.drop_duplicates()
pd.set_option('display.max_columns', None)
#print(df_merged.head())
#df_merged.head().to_csv('df_merged.csv')
#print(dist_holo_order.head())

df_merged['Difference'] = df_merged['s2calc_x'] - df_merged['s2calc_y']

dist_holo_order_subset = dist_holo_order[['s2calc', 'Distance']]

#print(dist_holo_order_subset.head())
dist_holo_order_subset_melt = pd.melt(dist_holo_order_subset, id_vars=['Distance'], value_vars=['s2calc'])
dist_difference = pd.melt(df_merged, id_vars=['Distance'], value_vars=['Difference'])
print(dist_holo_order_subset_melt.head())

print('dist_difference:')
print(dist_difference.head())

fig = plt.figure()
sns.scatterplot(dist_difference['Distance'], dist_difference['value'])
#plt.xlabel(x_label)
#plt.legend()
#plt.ylabel(y_label)
plt.title('Distance by Difference Order Parameter')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_by_differenceOP.png')



fig = plt.figure()
sns.scatterplot(dist_holo_order_subset_melt['Distance'], dist_holo_order_subset_melt['value'])
#plt.xlabel(x_label)
#plt.legend()
#plt.ylabel(y_label)
plt.title('Distance by Order Parameter')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_by_OP.png')

#dist_holo_order_subset = dist_holo_order_subset.interpolate(method='cubic')
#dist_holo_order_subset.plot()
#dist_holo_order = dist_holo_order.dropna()
#f1 = interp1d(dist_holo_order.index, dist_holo_order['s2calc'],kind='cubic')

#df2 = pd.DataFrame()
#df2['s2calc'] = f1

#dist_holo_order['Distance_cubic']
#fig = plt.figure()
#sns.lineplot(x=dist_holo_order['Distance_cubic'], y=dist_holo_order['s2calc'])
plt.show()
#dist_holo_order.plot.line(Distance_cubic, s2calc)
#df_smooth = dist_holo_order.reindex(index='Distance').interpolate('cubic')
#df_smooth.plot(ax=axs, alpha=0.7)
#f.plot(ax=axs, alpha=0.7)
#fig.show()

#fig = plt.figure()
#scipy.interpolate.BSpline(dist_holo_order['Distance'], dist_holo_order['s2calc'], k=3)
#fig.show()

#f, ax = plt.subplots(figsize=(8, 8))
#piv = pd.pivot_table(CDK2_merge, values=["Difference"], index=['Ligand'], columns=['resi_x'])
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
#ax = sns.heatmap(piv, square=True, cbar=False)
#plt.show()

#print(dist_apo_order.head(1000))
CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']
#print(len(dist_apo_order.index))

#CDK2_holo = dist_holo_order[dist_holo_order['PDB'].isin(CDK2)]
#CDK2_apo = dist_apo_order[dist_apo_order['PDB']=='1pw2']
#CDK2_apo = CDK2_apo.drop_duplicates()
#print(CDK2_apo['PDB'].unique())
#print(len(CDK2_apo.index))
#print(CDK2_holo['PDB'].unique())

#print(len(CDK2_holo.index))

#CDK2_apo['PDB_lig'] = CDK2_apo['PDB'] + CDK2_apo['Ligand']
#CDK2_holo['PDB_lig'] = CDK2_holo['PDB'] + CDK2_holo['Ligand']
#CDK2_apo = CDK2_apo.drop_duplicates()


#CDK2_merge = CDK2_apo.merge(CDK2_holo, on=['Ligand', 'res'])
#CDK2_merge['Difference'] = CDK2_merge['s2calc_x'] - CDK2_merge['s2calc_y']
#CDK2_holo.to_csv('CDK2_holo.csv')
#CDK2_apo.to_csv('CDK2_apo.csv')
#print(len(CDK2_apo.index))
#CDK2_apo = CDK2_apo.dropna(subset=['resi', 's2calc_y'])
#CDK2_holo = CDK2_holo.dropna(subset=['resi', 's2calc_x'])

#fig = plt.figure()
#sns.lineplot(x=CDK2_holo['Distance'], y=CDK2_holo['s2calc'], hue=CDK2_holo['PDB_lig'])
#plt.show()

#fig = plt.figure()
#sns.lineplot(x=CDK2_apo['Distance'], y=CDK2_apo['s2calc'], hue=CDK2_apo['PDB_lig'])
#plt.show()

#fig = plt.figure()
#sns.lineplot(x=CDK2_merge['resi_x'], y=CDK2_merge['Difference'], hue=CDK2_merge['Ligand'])
#plt.show()

#CDK2_merge.to_csv('CDK2_merge.csv', index=False)

#f, ax = plt.subplots(figsize=(8, 8))
#piv = pd.pivot_table(CDK2_merge, values=["Difference"], index=['Ligand'], columns=['resi_x'])
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
#ax = sns.heatmap(piv, square=True, cbar=False)
#plt.show()

#CDK2_holo = CDK2_holo.dropna(subset=['Distance', 's2calc_x'])
#xnew = np.linspace(CDK2_holo['s2calc_x'].min(), CDK2_holo['s2calc_x'].max(), 100)
#power_smooth = spline(CDK2_holo['Distance'], CDK2_holo['s2calc_x'], xnew)

#plt.plot(xnew,power_smooth)
#plt.show()


