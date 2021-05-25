import pandas as pd
import numpy as np
import glob
import os
import matplotlib.pyplot as plt
from figure_functions import *	

#os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/water_analysis/ultra_high_res/')
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/water_analysis/')
path=os.getcwd()


# all_files = glob.glob(path + "/*_water_df.csv")

# li = []

# for filename in all_files:
#     df = pd.read_csv(filename, index_col=None, sep=',')
#     df['PDB'] = filename[77:81]
#     li.append(df)

# water = pd.concat(li, axis=0, ignore_index=True)
# print(water.head())


#plt.hist(water['bfactor'], bins = 200)
#plt.show()

#plt.hist(water['occ'], bins = 200)
#plt.show()
#ultra_high_res/
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/water_analysis/')
path=os.getcwd()


all_files = glob.glob(path + "/*_waterclosestresidue.txt")

li = []

for filename in all_files:
    df = pd.read_csv(filename, sep=',', header=None)
    df['PDB'] = filename[62:66]
    #df['PDB'] = filename[77:81]
    li.append(df)

water_close = pd.concat(li, axis=0, ignore_index=True)
water_close= water_close.rename(columns={0: "Protein_Resi", 1: "Protein_Chain", 2: "Water_Resi", 3: "Water_Chain", 4: "Water_Occ", 5: "Water_B"})

os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)
water_close = water_close[water_close['PDB'].isin(AH_key['PDB'])]
water_close = water_close[~water_close['PDB'].isin(['2tga'])]
print(water_close.head())
water_close[['Protein_Chain','Distance']] = water_close['Protein_Chain'].str.split(']',expand=True)
water_close['Protein_Resi'] = water_close.Protein_Resi.str.replace('[','')
water_close['Protein_Resi'] = water_close.Protein_Resi.str.replace(']','')
water_close['Protein_Chain'] = water_close.Protein_Chain.str.replace('[','')
water_close['Protein_Chain'] = water_close.Protein_Chain.str.replace(']','')
water_close['Protein_Chain'] = water_close.Protein_Chain.str.replace("'",'')
water_close['Protein_Chain'] = water_close.Protein_Chain.str.replace(']','')
water_close['Water_Occ'] = water_close.Water_Occ.str.replace('[','')
water_close['Water_Occ'] = water_close.Water_Occ.str.replace(']','')
water_close['Water_B'] = water_close.Water_B.str.replace('[','')
water_close['Water_B'] = water_close.Water_B.str.replace(']','')
water_close.Water_Occ.str.rstrip('.')
water_close['Water_Occ'] = water_close['Water_Occ'].astype('float64')
water_close[['Protein_Resi','Protein_Resi2']] = water_close["Protein_Resi"].str.split(" ", 1, expand=True)
water_close = water_close.drop(['Protein_Resi2'], axis=1)
water_close['Protein_Resi'] = water_close['Protein_Resi'].astype('int64')
water_close['Water_B'] = water_close['Water_B'].astype('float64')
water_close['Distance'] = water_close['Distance'].astype('float64')


print(water_close[water_close['Distance']>10])
print(water_close[water_close['Water_B']>200])
water_close = water_close[water_close['Water_B']<200]
print(water_close.head())
# plt.scatter(water_close['Distance'], water_close['Water_Occ'])
# plt.show()

# plt.scatter(water_close['Distance'], water_close['Water_B'])
# plt.show()

# plt.scatter(water_close['Water_B'], water_close['Water_Occ'])
# plt.show()

RMSF = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_summary.csv')

#print(RMSF.head())

RMSF_water = RMSF.merge(water_close, left_on=['PDB', 'resi'], right_on=['PDB', 'Protein_Resi'])
RMSF_water = RMSF_water[RMSF_water['Distance']<=2]
RMSF_water = RMSF_water[RMSF_water['Water_Occ']<1]
print(RMSF_water[RMSF_water['RMSF']>1])

plt.scatter(RMSF_water['Water_B'], RMSF_water['RMSF'])
plt.show()

plt.scatter(RMSF_water['Water_Occ'], RMSF_water['RMSF'])
plt.show()

plt.scatter(RMSF_water['Distance'], RMSF_water['RMSF'])
plt.show()