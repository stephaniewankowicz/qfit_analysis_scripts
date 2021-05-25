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

CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']
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


# dist_all_collapse = pd.DataFrame()
# n=1
# for i in dist_all['PDB'].unique():
#     #print(i)
#     for lig in dist_all[dist_all['PDB'] == i]['Ligand'].unique():
#         #print(lig)
#         tmp = dist_all[(dist_all['PDB'] == i) & (dist_all['Ligand'] == lig)] 
#         for res in tmp['Residue'].unique():
#             #print(res)
#             dist_all_collapse.loc[n, 'PDB'] = i
#             dist_all_collapse.loc[n, 'res'] = res
#             dist_all_collapse.loc[n, 'Ligand'] = lig
#             dist_all_collapse.loc[n, 'Distance'] = np.min(dist_all[(dist_all['PDB'] == i) & (dist_all['Residue'] == res) & (dist_all['Ligand'] == lig)]['Distance'])
#             n += 1
        
dist_all_collapse = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_CDK2.csv')
#dist_all_collapse.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse_CDK2.csv')
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
merged_order_all = merged_order_all[merged_order_all['PDB_x'].isin(CDK2)]


merged_order_all['Difference'] = merged_order_all['s2calc_x'] - merged_order_all['s2calc_y']

CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']

for i in CDK2:
    print(i)
    merged_order_all[merged_order_all['Holo']==i].to_csv(f"{i}_difference.csv", index=False)

#merged_order_all = merged_order_all.merge(dist_all_collapse)
print('merged_order_all:')
print(merged_order_all.head())

merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

apo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_apo.csv')
holo = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/order_all_holo.csv')

apo = apo.dropna(subset=['resi'])
holo = holo.dropna(subset=['resi'])


holo['resi'] = holo['resi'].astype(int)
holo['res'] = holo['resi'].astype(str) + holo['chain']

apo['resi'] = apo['resi'].astype(int)
apo['res'] = apo['resi'].astype(str) + apo['chain']

merged_order_all['resi'] = merged_order_all['resi'].astype(int)
merged_order_all['res'] = merged_order_all['resi'].astype(str) + merged_order_all['chain']

dist_holo_order = dist_all_collapse_holo.merge(holo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])
dist_apo_order = dist_all_collapse_apo.merge(apo, left_on=['PDB', 'res'], right_on=['PDB', 'res'])

dist_apo_order = dist_apo_order.drop_duplicates()
CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']

CDK2_holo = dist_holo_order[dist_holo_order['PDB'].isin(CDK2)]
CDK2_apo = dist_apo_order[dist_apo_order['PDB']=='1pw2']
CDK2_apo = CDK2_apo.drop_duplicates()

CDK2_apo['PDB_lig'] = CDK2_apo['PDB'] + CDK2_apo['Ligand']
CDK2_holo['PDB_lig'] = CDK2_holo['PDB'] + CDK2_holo['Ligand']
CDK2_apo = CDK2_apo.drop_duplicates()


CDK2_merge = CDK2_apo.merge(CDK2_holo, on=['Ligand', 'res'])
CDK2_merge['Difference'] = CDK2_merge['s2calc_x'] - CDK2_merge['s2calc_y']

fig = plt.figure()
sns.lineplot(x=CDK2_holo['Distance'], y=CDK2_holo['s2calc'], hue=CDK2_holo['PDB_lig'])
#plt.show()

fig = plt.figure()
sns.lineplot(x=CDK2_apo['Distance'], y=CDK2_apo['s2calc'], hue=CDK2_apo['PDB_lig'])
#plt.show()

fig = plt.figure()
sns.lineplot(x=CDK2_merge['resi_x'], y=CDK2_merge['Difference'], hue=CDK2_merge['Ligand'])
#plt.show()

CDK2_merge.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/CDK2_merge.csv', index=False)

fig = plt.figure()
f, ax = plt.subplots(figsize=(20, 15))
#piv = pd.pivot_table(merged_order_all, values=["Difference"], index=['Ligand'], columns=['res'])
piv =  merged_order_all.pivot("Ligand", "res", "Difference")
#print(piv.head())
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.heatmap(piv, cbar=False, xticklabels=2, yticklabels=2) #, ax=ax
plt.xticks(rotation=45) 
#ax.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/cdk2_difference_heatmap.png')
#plt.show()
print(piv.head())

cmap = sns.palplot(sns.diverging_palette(240, 10, n=9))

fig=plt.figure()
f, ax = plt.subplots(figsize=(20, 15))
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.clustermap(piv, yticklabels=True, cmap="RdBu") #, xticklabels=2, yticklabels=2#row_cluster=False,
plt.xticks(rotation=45)
ax.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/cdk2_difference_clustermap.png')
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
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/cdk2_cluster_eucledian.png')




print('cluster:')
d = sch.distance.pdist(piv)
L = sch.linkage(d, method='complete')

clusters = sch.fcluster(L,0.5*d.max(), 'distance')

for i, cluster in enumerate(clusters):
    print(piv.index[i], cluster)









CDK2 = ['1pw2','3qu0', '3ygk', '3r1q', '3qtw', '3qwk', '3qx4', '3r7y', '3qql', '3rni', '3rm7', '3r7e', '3rk9', '3r7v', '3rjc', '3qtx', '3qzh', '3gl8', '3r9d', '1pxi', '2a0c', '1pw2', '3unj']

CDK2_holo = dist_holo_order[dist_holo_order['PDB'].isin(CDK2)]
CDK2_apo = dist_apo_order[dist_apo_order['PDB']=='1pw2']
CDK2_apo = CDK2_apo.drop_duplicates()


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

# fig = plt.figure()
# sns.lineplot(x=CDK2_holo['Distance'], y=CDK2_holo['s2calc'], hue=CDK2_holo['PDB_lig'])
# plt.show()

# fig = plt.figure()
# sns.lineplot(x=CDK2_apo['Distance'], y=CDK2_apo['s2calc'], hue=CDK2_apo['PDB_lig'])
# plt.show()

# fig = plt.figure()
# sns.lineplot(x=CDK2_merge['resi_x'], y=CDK2_merge['Difference'], hue=CDK2_merge['Ligand'])
# plt.show()

CDK2_merge.to_csv('CDK2_merge.csv', index=False)

f, ax = plt.subplots(figsize=(8, 8))
piv = pd.pivot_table(CDK2_merge, values=["Difference"], index=['Ligand'], columns=['res'])
#piv_AA = pd.pivot_table(order_CDK2, values=["resn"], index='PDB', columns='resi')
ax = sns.heatmap(piv, square=True, cbar=False)
#plt.show()


merged_order_5 = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_5.csv')
merged_order_far = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/merged_order_far.csv')
CDK2_merged_order_far = merged_order_far[merged_order_far['PDB_x'].isin(CDK2)]
CDK2_merged_order_5 = merged_order_5[merged_order_5['PDB_x'].isin(CDK2)]
CDK2_merged_order_5['Difference'] = CDK2_merged_order_5['s2calc_x'] - CDK2_merged_order_5['s2calc_y']
CDK2_merged_order_far['Difference'] = CDK2_merged_order_far['s2calc_x_x'] - CDK2_merged_order_far['s2calc_x_y']

print('CDK2_merged_order_far')
print(CDK2_merged_order_far.head())

make_dist_plot_AH(CDK2_merged_order_far['s2calc_x_x'], CDK2_merged_order_far['s2calc_x_y'], 's2ortho', 'Number of Residues', 'Bound v. Unbound s2calc (Further than 10A)', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2ortho_>10A_CDK2')

print('Difference of s2calc on Side Chains >10 A between Bound/Unbound [Entire Protein]')
paired_wilcoxon(CDK2_merged_order_far['s2calc_x_x'], CDK2_merged_order_far['s2calc_x_y'])


make_dist_plot_AH(CDK2_merged_order_5['s2calc_x'], CDK2_merged_order_5['s2calc_y'], 's2calc', 'Number of Residues', 'Bound v. Unbound s2calc Binding site', '/Users/stephaniewankowicz/Downloads/qfit_paper/AH_s2ortho_>10A_CDK2')

print('Difference of s2calc on Side Chains binding site')
paired_wilcoxon(CDK2_merged_order_5['s2calc_x'], CDK2_merged_order_5['s2calc_y'])
#CDK2_holo = CDK2_holo.dropna(subset=['Distance', 's2calc_x'])
#xnew = np.linspace(CDK2_holo['s2calc_x'].min(), CDK2_holo['s2calc_x'].max(), 100)
#power_smooth = spline(CDK2_holo['Distance'], CDK2_holo['s2calc_x'], xnew)

#plt.plot(xnew,power_smooth)
#plt.show()


