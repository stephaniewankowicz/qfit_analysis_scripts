import pandas as pd
from figure_functions import *	
import numpy as np

#SUMMARY TABLES
bfactor_ligand_merged = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/Bfactor_Ligand_Summary.csv')
AH_rotamer_summary = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/AH_rotamer_summary.csv')

chem_properties = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/merged_order_summary_5.csv')
#water_5A = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/water_5A_merged.csv')
#water_5A_chem = water_5A.merge(chem_properties, left_on = 'Holo', right_on = 'PDB_x')

#water_5A_chem.head()

#water_5A_chem_lowMolWeight = water_5A_chem[water_5A_chem['MolWeight'] <= water_5A_chem['MolWeight'].quantile(0.25)]
#water_5A_chem_highMolWeight = water_5A_chem[water_5A_chem['MolWeight'] >= water_5A_chem['MolWeight'].quantile(0.75)]

#make_boxenplot_chem(water_5A_chem_highMolWeight['Water_Delta'], water_5A_chem_highMolWeight['Water_Delta'], 
#	'Low MolWeight', 'High MolWeight', 'Water Delta', '/Users/stephaniewankowicz/Downloads/qfit_paper/water_5A_Molweight')
#print('Ind ttest water moleweight, column1=low, column2=high')
#ind_ttest(water_5A_chem_lowMolWeight['Water_Delta'], water_5A_chem_highMolWeight['Water_Delta'])

#ind_ttest(water_5A_chem_lowMolWeight['num_water_holo'], water_5A_chem_highMolWeight['num_water_holo'])
#ind_ttest(water_5A_chem_lowMolWeight['num_water_apo'], water_5A_chem_highMolWeight['num_water_apo'])

AH_rotamer_summary_ligandb = AH_rotamer_summary.merge(bfactor_ligand_merged, left_on = 'Holo', right_on = 'PDB')
AH_rotamer_summary_chem = AH_rotamer_summary.merge(chem_properties, left_on = 'Holo', right_on = 'PDB_x')

#print(AH_rotamer_summary_ligandb.head())
#print(AH_rotamer_summary_chem.head())


AH_rotamer_summary_chem_lowMolWeight = AH_rotamer_summary_chem[AH_rotamer_summary_chem['MolWeight'] <= AH_rotamer_summary_chem['MolWeight'].quantile(0.25)]
AH_rotamer_summary_chem_highMolWeight = AH_rotamer_summary_chem[AH_rotamer_summary_chem['MolWeight'] >= AH_rotamer_summary_chem['MolWeight'].quantile(0.75)]

AH_rotamer_summary_chem_lowRotBonds = AH_rotamer_summary_chem[AH_rotamer_summary_chem['NumRotatableBonds'] <= AH_rotamer_summary_chem['NumRotatableBonds'].quantile(0.25)]
AH_rotamer_summary_chem_highRotBonds = AH_rotamer_summary_chem[AH_rotamer_summary_chem['NumRotatableBonds'] >= AH_rotamer_summary_chem['NumRotatableBonds'].quantile(0.75)]


different_low = len(AH_rotamer_summary_chem_lowMolWeight[AH_rotamer_summary_chem_lowMolWeight['Rotamer'] == 'Different'].index)
same_low = len(AH_rotamer_summary_chem_lowMolWeight[AH_rotamer_summary_chem_lowMolWeight['Rotamer'] == 'Same'].index)
both_low = len(AH_rotamer_summary_chem_lowMolWeight[AH_rotamer_summary_chem_lowMolWeight['Rotamer'] == 'Same and Different'].index)

different_high = len(AH_rotamer_summary_chem_highMolWeight[AH_rotamer_summary_chem_highMolWeight['Rotamer'] == 'Different'].index)
same_high = len(AH_rotamer_summary_chem_highMolWeight[AH_rotamer_summary_chem_highMolWeight['Rotamer'] == 'Same'].index)
both_high = len(AH_rotamer_summary_chem_highMolWeight[AH_rotamer_summary_chem_highMolWeight['Rotamer'] == 'Same and Different'].index)


fig, ax = plt.subplots()
barwidth = 0.4

high = [different_high, same_high, both_low]
low = [different_low, same_low, both_low]
print(high, low)

r1 = np.arange(len(high))
r2 = [x + barwidth for x in r1]

print(barwidth/2)
ax.bar(r1 - barwidth/2, high, width=barwidth, edgecolor='white', label='High')
ax.bar(r1 + barwidth/2, low, width=barwidth, edgecolor='white', label='Low')

ax.set_title('Rotamer Status by Rotatable Bonds', fontsize=18)
ax.set_ylabel('Number of Residues', fontsize=18)
ax.legend()
ax.set_xticks(r1)
ax.set_xticklabels(('All Different','All Same', 'Same & Different'), fontsize=12)
#plt.xticks(r1, ())
fig.tight_layout()
plt.savefig('RotamerStatus_by_Rotable Bonds'+ '.png')
plt.show()
