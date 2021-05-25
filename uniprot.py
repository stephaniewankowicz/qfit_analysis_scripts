import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns

colors = ["#1b9e77","#d95f02", "#7570b3","#e7298a","#66a61e", "#e6ab02", "#666666"]
sns.set_palette(sns.color_palette(colors))
pd.set_option('display.max_columns', None)

uniprot = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/uniprot_final_dataset.tab', sep='\t')

uniprot['EC_First'] = uniprot['EC number'].astype(str).str[0]

uniprot_enzyme = uniprot[uniprot['EC_First']!='n']
uniprot_enzyme['EC_Name'] = uniprot_enzyme['EC_First'].replace(['1','2','3','4','5','6','7'],['Oxidoreductase', 'Transferases', 'Hydrolases', 'Lyases', 'Isomerases', 'Ligases', 'Translocases'])

fig = plt.figure()
sns.countplot(x='EC_Name', data=uniprot_enzyme)
plt.xticks(rotation = 45)
plt.tight_layout()
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/enzyme_classification.png')


#Original Dataset
uniprot = pd.read_csv('/Users/stephaniewankowicz/Downloads/uniprot-yourlist_M20210325A94466D2655679D1FD8953E075198DA81965FAG.tab', sep='\t')
uniprot['Num_Apo']=uniprot['yourlist:M20210325A94466D2655679D1FD8953E075198DA81965FAG'].str.findall(',').str.len() + 1

print(uniprot['Num_Apo'].describe())
print(len(uniprot[uniprot['Num_Apo']==1].index))


uniprot_high = uniprot[uniprot['Num_Apo']>6]
uniprot_high = uniprot_high.sort_values(['Num_Apo'], ascending=False).reset_index(drop=True)

uniprot_high['Protein'] = uniprot_high['Protein names'].str.split("(").str[0]#[0]
uniprot_high['Protein'] = uniprot_high['Protein'].str.split("[").str[0]
fig = plt.figure()
sns.barplot(uniprot_high['Protein'], uniprot_high['Num_Apo'], palette=colors)
plt.title('Top 30 Frequent Ligands in Dataset')
plt.ylabel('Number of Ooccurrences')
plt.xlabel('Protein Name')
plt.xticks(rotation=20, ha='right')
plt.tight_layout()
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/uniprot_dist.png')
