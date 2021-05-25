import pandas as pd
import os

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed.txt', sep=' ')
#pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()
print(len(AH_pairs.index))

rvalues = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/PDBs_to_keep_rvalue.csv') 
AH_pairs = AH_pairs[AH_pairs['Apo'].isin(rvalues['Apo']) & AH_pairs['Holo'].isin(rvalues['Holo'])]

clash_remove = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/clash_to_remove.csv') 
print(len(clash_remove['PDB'].index))
print(clash_remove.head())
AH_pairs = AH_pairs[~AH_pairs['Holo'].isin(clash_remove['PDB'])]
AH_pairs = AH_pairs[~AH_pairs['Apo'].isin(clash_remove['PDB'])]

to_remove = ['6fji']
AH_pairs = AH_pairs[~AH_pairs['Holo'].isin(to_remove)]
AH_pairs = AH_pairs[~AH_pairs['Apo'].isin(to_remove)]

print(len(AH_pairs.index))
AH_pairs.to_csv('ligand_supplementary_table1_QCed_updated_200211.txt', sep=' ', index=False)