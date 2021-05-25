import os
import pandas as pd
from figure_functions import *	

#reference files
pairs_ref = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ')
pairs_ref.drop_duplicates(inplace=True)

print(pairs_ref.head())
ah_key = create_AH_key(pairs_ref)
print(ah_key)
pairs = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/pre_refine_done_id.txt', sep=' ', header=None)
pairs.drop_duplicates(inplace=True)
pairs = pairs.rename(columns={0: "PDB"})

print(pairs.head())


test = pairs.merge(ah_key, on='PDB')
test.drop_duplicates(inplace=True)
print(test.head())

print(len(test[test['Apo_Holo']=='Apo'].index))
print(len(test[test['Apo_Holo']=='Holo'].index))