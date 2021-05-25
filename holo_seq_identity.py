#packages
from __future__ import division
import pandas as pd
import numpy as np
import os

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

holo_seq=pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/seq_holo_test.tab', sep='\t', header=None)
holo_seq.columns=['PDB1', 'PDB2', 'SeqID', 'Len', 'Blank1', 'B2', 'B3', 'B4', 'B5', 'B6', 'B7','B8']

to_keep = pd.DataFrame()
n = 1
for i in holo_seq['PDB1'].unique():
	to_remove=[]
	print(i)
	to_remove.append(i)
	tmp = holo_seq[holo_seq['PDB1']==i]
	for t in tmp['PDB2'].unique():
		if (tmp[tmp['PDB2']==t]['SeqID'] >= 90).any():
			to_remove.append(t)
	print(len(to_remove))
	subset = AH_pairs[AH_pairs['Holo'].isin(to_remove)]
	#print(subset['Holo_Res'].min())
	i2 = subset['Holo_Res'].min()
	try:
		to_keep.loc[n, 'PDB_to_keep'] = (subset[subset['Holo_Res'] == i2]['Holo'].unique()[0])
	except:
		continue
	n += 1
to_keep = to_keep.drop_duplicates()
to_keep.to_csv('PDB_to_keep.csv', index=False)



