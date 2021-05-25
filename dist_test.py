import pandas as pd
from figure_functions import *	
import numpy as np
import glob
#from scipy.interpolate import spline

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

dist_all_collapse = pd.DataFrame()
n=1
dist = []
for i in dist_all['PDB'].unique():
    print(i)
    for lig in dist_all[dist_all['PDB'] == i]['Ligand'].unique():
        tmp = dist_all[(dist_all['PDB'] == i) & (dist_all['Ligand'] == lig)] 
        for res in tmp['Residue'].unique():
            dist.append(tuple((i, res, lig, np.min(tmp[tmp['Residue']==res]['Distance']))))

df = pd.DataFrame(dist, columns =['PDB', 'Residue', 'Ligand', 'Distance']) 
print(df.head())


df.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/dist_all_collapse.csv')