#packages
from __future__ import division
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


#REFERENCE FILE
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
AH_key=pd.read_csv('qfit_AH_key_191218.csv')


os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/water_partial')
path=os.getcwd()

#plot numbers
all_files = glob.glob(path + "/*qFit_3.0_water_numbers.txt")

li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=None)
    df['PDB'] = filename[61:65]
    li.append(df)
water_num = pd.concat(li, axis=0, ignore_index=True)
water_num = water_num.rename(columns={0: "PDB1", 1: "Num_Water", 2: "Num_Close_Lig", 3: "Num_Close_Protein", 4: "Num_Partial_Lig", 5: "Num_Partial_Protein"})

water_num['Percent_3A_Protein'] = water_num['Num_Close_Protein']/water_num['Num_Water']
water_num['Percent_3A_Ligand'] = water_num['Num_Close_Lig']/water_num['Num_Water']
water_num['Percent_3A_Partial_Protein'] = water_num['Num_Partial_Protein']/water_num['Num_Water']
water_num['Percent_3A_Partial_Ligand'] = water_num['Num_Partial_Lig']/water_num['Num_Water']
print(water_num.head())
print(len(water_num.index))
#plot hist of distance of partial versus full occupancy water molecules
fig = plt.figure()
x = range(3)
f, axes = plt.subplots(1, 3, sharey=True, sharex=True)
p1 = sns.boxenplot(water_num['Num_Water'], orient='v', ax=axes[0]).set(xlabel = 'All Water Molecules', ylabel = '')
p2 = sns.boxenplot(water_num['Num_Close_Protein'], orient='v', ax=axes[1]).set(xlabel = 'Within 3A of Protein', ylabel = '')
p3 = sns.boxenplot(water_num['Num_Close_Lig'], orient='v', ax=axes[2]).set(xlabel = 'Within 3A of Ligands', ylabel = '')

plt.savefig('water_numbers.png')

fig = plt.figure()
x = range(2)
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.boxenplot(water_num['Num_Partial_Protein'], orient='v', ax=axes[0]).set(xlabel = 'Within 3A of Protein (Partial)', ylabel = '')
p2 = sns.boxenplot(water_num['Num_Partial_Lig'], orient='v', ax=axes[1]).set(xlabel = 'Within 3A of Ligands (Partial)', ylabel = '')
plt.savefig('water_numbers2.png')

all_files = glob.glob(path + "/*qFit_3.0_partialwatercloseresidue.txt")

li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=None)
    df['PDB'] = filename[61:65]
    li.append(df)

water_partial = pd.concat(li, axis=0, ignore_index=True)
water_partial = water_partial.rename(columns={0: "Water_Resi", 1: "Water_Chain", 2: "Distance", 3: "Resi"})
water_partial['Full_Partial'] ='Partial'

print(water_partial.head())

all_files = glob.glob(path + "/*_qFit_3.0_waterclosestresidue.txt")

li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=None)
    df['PDB'] = filename[61:65]
    li.append(df)

water_closest = pd.concat(li, axis=0, ignore_index=True)
water_closest= water_closest.rename(columns={0: "Water_Resi", 1: "Water_Chain", 2:"PDB_Resi",3:"Distance", 4:"Resi_Bfactor", 5:"Water_Bfactor"})
water_closest['Water_Bfactor'] = water_closest.Water_Bfactor.str.replace('[','')
water_closest['Water_Bfactor'] = water_closest.Water_Bfactor.str.replace(']','')
water_closest['Water_Bfactor'] = pd.to_numeric(water_closest['Water_Bfactor'])

print(water_closest.head())










water_merged = water_closest.merge(water_partial, on=['PDB','Water_Resi', 'Water_Chain'], how='left')  
water_merged = water_merged.drop_duplicates()
water_closest.to_csv('water_closest.csv', index=False)
water_partial.to_csv('water_partial.csv',index=False)



water_partial_dist = water_merged[water_merged['Full_Partial']=='Partial']
print(water_partial_dist.head())

#plot hist of distance of partial versus full occupancy water molecules
fig = plt.figure()
x = range(2)
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.boxenplot(water_partial_dist['Distance_x'], orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
p2 = sns.boxenplot(water_merged['Distance_x'], orient='v', ax=axes[1]).set(xlabel = 'All', ylabel = '')
plt.savefig('water.png')

fig = plt.figure()
p1 = sns.scatterplot(water_partial_dist['Distance_x'], water_partial_dist['Water_Bfactor'])#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_scatter_partial.png')

fig = plt.figure()
p1 = sns.scatterplot(water_partial_dist['Resi_Bfactor'], water_partial_dist['Water_Bfactor'])#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_scatter_b_factor_partial.png')

fig = plt.figure()
p1 = sns.distplot(water_partial_dist['Water_Bfactor'], kde=False)#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_b_factor_partial_dist.png')

fig = plt.figure()
p1 = sns.distplot(water_merged['Water_Bfactor'], kde=False)#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_b_factor_all_dist.png')

print(water_merged['Water_Bfactor'].quantile(0.75))
print(water_merged['Water_Bfactor'].quantile(0.5))
print(water_merged['Water_Bfactor'].quantile(0.25))
high_bfactor_water = water_merged[water_merged['Water_Bfactor'] >= water_merged['Water_Bfactor'].quantile(0.75)]
high_bfactor_water.to_csv('high_bfactor_water.csv', index=False)

fig = plt.figure()
p1 = sns.scatterplot(water_merged['Distance_x'], water_merged['Water_Bfactor'])#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_scatter_all.png')


fig = plt.figure()
p1 = sns.scatterplot(water_merged['Resi_Bfactor'], water_merged['Water_Bfactor'])#orient='v', ax=axes[0]).set(xlabel = 'Partial', ylabel = 'Distance')
plt.savefig('water_scatter_b_factor_all.png')



 


