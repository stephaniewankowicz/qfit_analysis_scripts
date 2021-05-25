#packages
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from figure_functions import *	

os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/') #_updated_200211
pairs = pd.read_csv('ligand_supplementary_table1_QCed.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key = create_AH_key(AH_pairs)

os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/201116/clashscore')
path=os.getcwd()

all_files = glob.glob(path + "/*_clashscore.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, sep=',', index_col=0)
    df['PDB'] = filename[65:69]
    li.append(df)

clash_all = pd.concat(li, axis=0, ignore_index=True)
print(clash_all.head())
clash_all = clash_all[clash_all['PDB'].isin(AH_key['PDB'])]

clash_all = clash_all.dropna()
clash_all.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/clash_all.csv')
clash_all[clash_all['Clashscore']>15].to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/clash_to_remove.csv')
print(clash_all['Clashscore'].describe())
print(len(clash_all[clash_all['Clashscore']>15].index))
fig = plt.figure()
sns.distplot(clash_all['Clashscore'], kde=False)
plt.xlabel('Clashscore')
plt.legend()
plt.ylabel('Number of Structures')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/clashscore.png')
