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


os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
path=os.getcwd()

all_files = glob.glob(path + "/*partialwatercloseresidue.txt")

li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=None)
    df['PDB'] = filename[47:51]
    li.append(df)

water_partial = pd.concat(li, axis=0, ignore_index=True)
water_partial = water_partial.rename(columns={0: "Water_Resi", 1: "Water_Chain", 2: "Chain", 3: "Resi"})

all_files = glob.glob(path + "/*_qFit_rotamer_output.txt")
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, header=0, sep=':')
    df['PDB'] = filename[47:51]
    li.append(df)

rotamer = pd.concat(li, axis=0, ignore_index=True)
print(len(rotamer.index))
#rotamer = rotamer[rotamer['residue'] != ['SUMMARY']]
rotamer = rotamer[rotamer['residue'] != ['SUMMARY']] #.drop(rotamer.filter(regex='SUMMARY').columns, axis=1)
split = rotamer['residue'].str.split(" ")
print(len(rotamer.index))
print(len(rotamer.index) - 1)
for i in range(0, len(rotamer.index) - 1):
	print(i)
	print(split[i])
	rotamer.loc[i,'chain'] = split[i][1]
	STUPID = str(rotamer.loc[i,'residue'])
	for s in STUPID.split():
		if s.isdigit():
			rotamer.loc[i,'resi'] = s
     #[int(s) for s in STUPID.split() if s.isdigit()]

print(rotamer.head())

rotamer_partial_water = rotamer.merge(water_partial, left_on=['PDB', 'chain', 'resi'], right_on=['PDB', 'Chain', 'Resi'], how='left')
rotamer_partial_water_yes = rotamer_partial_water[rotamer_partial_water.Water_Resi.notnull()]
rotamer_partial_water_no = rotamer_partial_water[rotamer_partial_water.Water_Resi.isnull()]


print(rotamer_partial_water_yes.head())












