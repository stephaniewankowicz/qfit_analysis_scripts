import pandas as pd
from figure_functions import *  
import numpy as np
import glob
import scipy.cluster.hierarchy as sch
from scipy.spatial import distance_matrix
from scipy.cluster import hierarchy
from sklearn import manifold
from sklearn.decomposition import PCA
#from scipy.interpolate import spline

pd.set_option('display.max_columns', None)

CDK2 = ['3k7k','4q7s', '3mnu', '3sap', '5nee', '3bl0', '3rz8', '3hkn', '3ml2', '3t5u', '5n0e', '2wd3', '4qjm', '4pyy', '3rz7']
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()


all_files = glob.glob(path + "/*_*_residue_dist.txt")
li = []

for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',', header=None)
    df['PDB'] = filename[54:58]
    df['Ligand'] = filename[59:62]
    li.append(df)

dist_all = pd.concat(li, axis=0, ignore_index=True)

dist_all = dist_all.rename(columns={0:'Resn', 1:'Residue', 2:'AltLoc', 3:'Distance'})
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
        
#dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_Ca2.csv')
dist_all_collapse.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_Ca2.csv')
dist_all_collapse = dist_all_collapse.merge(AH_key, on = ['PDB'])
dist_all_collapse = dist_all_collapse.drop_duplicates()
print('dist_all_collapse:')
print(dist_all_collapse.head())
print(len(dist_all_collapse.index))

dist_all_collapse_holo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Holo']
dist_all_collapse_apo = dist_all_collapse[dist_all_collapse['Apo_Holo'] == 'Apo']
#test = dist_all_collapse_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
#dist_all_collapse_m = test.merge(dist_all_collapse_apo, left_on='Apo', right_on='PDB')
#dist_all_collapse_m = dist_all_collapse_m.drop_duplicates()

merged_order_all = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_all.csv')
print('merged_order_all')
print(merged_order_all.head())
#merged_order_all = merged_order_all[merged_order_all['PDB_x'].isin(CDK2)]

merged_order_all = merged_order_all[merged_order_all['Apo'] =='3atk']
merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']
print(merged_order_all.head())

# for i in CDK2:
#     print(i)
#     merged_order_all[merged_order_all['Holo']==i].to_csv(f'{i}_difference.csv', index=False)

merged_order_all = merged_order_all.merge(dist_all_collapse_apo, left_on=['PDB_y', 'Ligand'], right_on=['PDB', 'Ligand'])
print('merged_order_all:')
print(merged_order_all.head())

merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

#apo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_apo.csv')
#holo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_holo.csv')

#apo = apo.dropna(subset=['resi'])
#holo = holo.dropna(subset=['resi'])


# holo['resi'] = holo['resi'].astype(int)
# holo['res'] = holo['resi'].astype(str) + holo['chain']

# apo['resi'] = apo['resi'].astype(int)
# apo['res'] = apo['resi'].astype(str) + apo['chain']

#merged_order_all['resi'] = merged_order_all['resi'].astype(int)
#merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

# dist_holo_order = dist_all_collapse_holo.merge(holo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])
# dist_apo_order = dist_all_collapse_apo.merge(apo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])

# dist_apo_order = dist_apo_order.drop_duplicates()

#CDK2_holo = dist_holo_order[dist_holo_order['PDB'].isin(CDK2)]
#CDK2_apo = dist_apo_order[dist_apo_order['PDB']=='4q7s']
#CDK2_apo = CDK2_apo.drop_duplicates()

#print(CDK2_apo.head())
merged_order_all['PDB_lig'] = merged_order_all['Holo'] + merged_order_all['Ligand']
#CDK2_apo['PDB_lig'] = CDK2_apo['PDB'] + CDK2_apo['Ligand']
#CDK2_holo['PDB_lig'] = CDK2_holo['PDB'] + CDK2_holo['Ligand']
#CDK2_apo = CDK2_apo.drop_duplicates()


#CDK2_merge = CDK2_apo.merge(CDK2_holo, on=['Ligand', 'res'])
merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']
print(merged_order_all.head())


#fig = plt.figure()
#sns.lineplot(x=merged_order_all['Distance'], y=merged_order_all['s2calc_x'], hue=merged_order_all['PDB_lig'])
#plt.show()

#fig = plt.figure()
#sns.lineplot(x=merged_order_all['Distance'], y=merged_order_all['s2calc_y'], hue=merged_order_all['PDB_lig'])
#plt.show()

fig = plt.figure()
sns.lineplot(x=merged_order_all['resi'], y=merged_order_all['Difference'], hue=merged_order_all['Ligand'])
#plt.show()
merged_order_all = merged_order_all.drop(['Distance', 'Unnamed: 0_y', 'Unnamed: 0_x'], axis=1)
merged_order_all = merged_order_all.drop_duplicates()

merged_order_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/trypsin_merge.csv', index=False)

fig = plt.figure()
f, ax = plt.subplots(figsize=(20, 15))
#piv = pd.pivot_table(merged_order_all, values=["Difference"], index=['Ligand'], columns=['res'])
piv =  merged_order_all.pivot("PDB_lig", "res", "Difference")
#print(piv.head())
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.heatmap(piv, cbar=False, xticklabels=2, yticklabels=2) #, ax=ax
plt.xticks(rotation=45) 
#ax.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/cdk2_difference_heatmap.png')
#plt.show()
print(piv.head())


fig=plt.figure()
f, ax = plt.subplots(figsize=(20, 15))
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.clustermap(piv, yticklabels=True, method="ward") #, xticklabels=2, yticklabels=2#row_cluster=False,
plt.xticks(rotation=45)
ax.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/trypsin_difference_clustermap.png')
#plt.show()

distance_matrix = pd.DataFrame(distance_matrix(piv.values, piv.values), index=piv.index, columns=piv.index)
#distance_matrix.to_csv('distance_matrix.CSV',sep=',')
print(distance_matrix.head())


pca = PCA(n_components=3)
principalComponents = pca.fit_transform(piv)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['PCA1', 'PCA2', 'PCA3'])

print(piv.index)
finalDf = pd.concat([principalDf, merged_order_all[['res']]], axis = 1)

fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot(1,1,1) 
ax.set_xlabel('Principal Component 1', fontsize = 15)
ax.set_ylabel('Principal Component 2', fontsize = 15)
ax.set_title('2 component PCA', fontsize = 20)
ax.scatter(finalDf['PCA1']
               , finalDf['PCA2']
               , s = 50)
#ax.legend(targets)
ax.grid()
#plt.show()

mds = manifold.MDS(n_components=2, dissimilarity="precomputed", n_init=50, max_iter=1000, random_state=1)
results = mds.fit(distance_matrix.values)

Geo_bln = distance_matrix.columns
coords = results.embedding_

fig = plt.figure(figsize=(12,10))

plt.subplots_adjust(bottom = 0.1)
plt.scatter(coords[:, 0], coords[:, 1])

for label, x, y in zip(Geo_bln, coords[:, 0], coords[:, 1]):
    plt.annotate(
        label,
        xy = (x, y), 
        xytext = (-20, 20),
        textcoords = 'offset points'
    )
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/ca2_cluster_eucledian.png')




print('cluster:')
d = sch.distance.pdist(piv)
L = sch.linkage(d, method='complete')

clusters = sch.fcluster(L,0.5*d.max(), 'distance')

for i, cluster in enumerate(clusters):
    print(piv.index[i], cluster)





