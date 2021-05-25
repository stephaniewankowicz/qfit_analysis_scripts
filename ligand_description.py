#packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from rdkit import Chem
from figure_functions import *	


os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()
AH_key = create_AH_key(AH_pairs)
print(AH_key.head())

colors = ["#1b9e77","#d95f02", "#7570b3","#e7298a","#66a61e", "#e6ab02", "#666666"]
sns.set_palette(sns.color_palette(colors))

pd.set_option('display.max_columns', None)
chem = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/Chemical_Descriptors_PDB.csv')

print(chem.head())

chem = chem[chem['Holo'].isin(AH_key['PDB'])]

for i in range(len(chem['Smile_y'].index)):
	try:
		mol = Chem.MolFromSmiles(chem.loc[i,'Smile_y'])
		chem.loc[i,'Molecular_Size'] = mol.GetNumAtoms()
	except:
		continue

print('Median Chemical Molecular Size:')
print(np.nanmedian(chem['Molecular_Size']))

print('Median Chem Name Count:')
print(np.median(chem['ChemName'].value_counts()))

print('Median Molecular Weight:')
print(np.nanmedian(chem['MolWeight']))

fig = plt.figure()
sns.distplot(chem['Molecular_Size'], kde=False, label='')
plt.xlabel('Molecular_Size')
plt.ylabel('Number of Ligands')
plt.title('Ligand Size Distribution')
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/ligand_size_distribution.png')

print('Number of Unique Ligands:')
print(len(chem['ChemName'].unique()))

chem.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/chem_check.csv')

fig = plt.figure()
chemname_counts = chem['ChemName'].value_counts()[:30]
sns.barplot(chemname_counts.index, chemname_counts.values, alpha=0.8, palette=colors)
plt.title('Top 30 Frequent Ligands in Dataset')
plt.ylabel('Number of Ooccurrences')
plt.xlabel('Ligand PDB ID')
plt.xticks(rotation=45)
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/chem_name_dist.png')



# os.chdir('/Users/stephaniewankowicz/Downloads/apo_sdf_files/')
# path=os.getcwd()
# all_files = glob.glob(path + "/*sdf")
# li = []
# for filename in all_files:
#     print(filename)
#     suppl = Chem.SDMolSupplier(filename)
#     for mol in suppl:
#        print(mol)
#        print(mol.GetNumAtoms())



