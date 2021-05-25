#Alt Loc Summary
from __future__ import division
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from figure_functions import *

pd.set_option('display.max_columns', None)

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ')
pairs.drop_duplicates(inplace=True)

ah_key = create_AH_key(pairs)


#__________________POST QFIT_____________________________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()
all_files = glob.glob(path + "/*qfit_RMSF.csv")
li = []
for filename in all_files:
    df = pd.read_csv(filename, sep=',')
    df['PDB'] = filename[54:58]
    li.append(df)
altloc_qFit = pd.concat(li, axis=0, ignore_index=True)


altloc_qFit_sum= pd.DataFrame()
n = 1
for i in altloc_qFit['PDB'].unique():
    tmp = altloc_qFit[altloc_qFit['PDB'] == i]
    altloc_qFit_sum.loc[n, 'PDB'] = i
    altloc_qFit_sum.loc[n, 'Num_Residues'] = len(tmp.index)
    altloc_qFit_sum.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['Num_Alt']>0].index)
    altloc_qFit_sum.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
    altloc_qFit_sum.loc[n, 'test'] = 1
    n += 1

altloc_qFit_sum['per_altloc'] = altloc_qFit_sum['Num_Alt_Loc']/altloc_qFit_sum['Num_Residues']
altloc_qFit_sum['Num_Single'] = altloc_qFit_sum['Num_Residues'] - altloc_qFit_sum['Num_Alt_Loc']

altloc_qFit_sum = altloc_qFit_sum.merge(AH_key, on = ['PDB'])
altloc_qFit_sum_holo = altloc_qFit_sum[altloc_qFit_sum['Apo_Holo'] == 'Holo']
altloc_qFit_sum_apo = altloc_qFit_sum[altloc_qFit_sum['Apo_Holo'] == 'Apo']

test = altloc_qFit_sum_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
altloc_qFit_sum_m = test.merge(altloc_qFit_sum_apo, left_on=['Apo'], right_on=['PDB'])  
altloc_qFit_sum_m= altloc_qFit_sum_m.drop_duplicates()

altloc_qFit_sum_m['Apo_Holo_Diff'] = altloc_qFit_sum_m['per_altloc_x'] - altloc_qFit_sum_m['per_altloc_y']
altloc_qFit_sum_m['Apo_Holo_Diff_Num'] = altloc_qFit_sum_m['Num_Alt_Loc_x'] - altloc_qFit_sum_m['Num_Alt_Loc_y']


##___________________SUMMARY STATS_____________________

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_qFit_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
p2 = sns.boxenplot(altloc_qFit_sum[altloc_qFit_sum['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p3 = sns.boxenplot(altloc_qFit_sum[altloc_qFit_sum['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_fullprotein_postqFit.png')

print(altloc_qFit_sum.head())
fig = plt.figure()
#f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_qFit_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
#p2 = sns.violinplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
sns.violinplot(y='per_altloc', x='test', hue='Apo_Holo', data=altloc_qFit_sum, split=True, palette="muted", orient='v')#.set(xlabel='Holo', ylabel='')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_postqFit_violin.png')


print('Median Number of Alt Conf Post-qFit:')
print(altloc_qFit_sum['per_altloc'].median())

print('Median Number of Alt Conf Post-qFit (Holo):')
print(altloc_qFit_sum[altloc_qFit_sum['Apo_Holo']=='Holo']['per_altloc'].median())

print('Median Number of Alt Conf Post-qFit (Apo):')
print(altloc_qFit_sum[altloc_qFit_sum['Apo_Holo']=='Apo']['per_altloc'].median())


plt.figure()
sns.distplot(altloc_qFit_sum_m['Apo_Holo_Diff'], kde=False, label='Post qFit')
plt.xlabel('Difference in Number of Alternative Conformations (Holo-Apo)')
plt.text(-200, 200, 'More Alt Conf in Apo') 
plt.text(0, 200, 'More Alt Conf in Holo')
plt.ylabel('Number of Pairs')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/altloc_diff_dist.png')

fig = plt.figure()
sns.violinplot(x=altloc_qFit_sum_m['Apo_Holo_Diff'], orient='v').set(xlabel = 'Number of Alt Conf', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/altloc_diff_violin.png')


fig = plt.figure()
sns.swarmplot(altloc_qFit_sum_m['Apo_Holo_Diff'], orient='v').set(xlabel = 'Number of Alt Conf', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/altloc_diff_swarm.png')


#__________________POSTQFIT CLOSE_______________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210426/')
path=os.getcwd()

all_files = glob.glob(path + "/*5.0_rmsf_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, header=0)
    li.append(df)

close_RMSF = pd.concat(li, axis=0, ignore_index=True)

close_RMSF['PDB'] = close_RMSF['PDB_name'].astype(str).str[0:4]
close_RMSF = close_RMSF.rename(columns={"Chain": "chain", "resseq":"resi"})
merged_close_RMSF = merge_apo_holo_df(close_RMSF)

close_RMSF_summary = pd.DataFrame()
n = 1
for i in close_RMSF['PDB'].unique():
    tmp = close_RMSF[close_RMSF['PDB'] == i]
    close_RMSF_summary.loc[n, 'PDB'] = i
    close_RMSF_summary.loc[n, 'Num_Residues'] = len(tmp.index)
    close_RMSF_summary.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['Num_Alt']>0].index)
    close_RMSF_summary.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
    close_RMSF_summary.loc[n, 'test'] = 1
    n += 1

close_RMSF_summary['per_altloc'] = close_RMSF_summary['Num_Alt_Loc'] / close_RMSF_summary['Num_Residues']

close_RMSF_summary = close_RMSF_summary.merge(AH_key)
close_RMSF_sum_holo = close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Holo']
close_RMSF_sum_apo = close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']
test = close_RMSF_sum_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
merged_close_sum_RMSF = test.merge(close_RMSF_sum_apo, left_on='Apo', right_on='PDB')

merged_close_sum_RMSF['Percent_Holo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_x']/ merged_close_sum_RMSF['Num_Residues_x']
merged_close_sum_RMSF['Percent_Apo_Close'] = merged_close_sum_RMSF['Num_Alt_Loc_y']/ merged_close_sum_RMSF['Num_Residues_y']
merged_close_sum_RMSF['Apo_Holo_Diff'] = merged_close_sum_RMSF['Percent_Holo_Close'] - merged_close_sum_RMSF['Percent_Apo_Close']
merged_close_sum_RMSF['Apo_Holo_Diff_Num'] = merged_close_sum_RMSF['Percent_Holo_Close'] - merged_close_sum_RMSF['Percent_Apo_Close']

merged_close_sum_RMSF.drop_duplicates(inplace=True)

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_qFit_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
p2 = sns.boxenplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p3 = sns.boxenplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_5A_postqFit.png')

print(close_RMSF_summary['Apo_Holo'].head())

fig = plt.figure()
#f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_qFit_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
#p2 = sns.violinplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
sns.violinplot(y='per_altloc', x='test', hue='Apo_Holo', data=close_RMSF_summary, split=True, palette="muted", orient='v').set(xlabel='', ylabel='Fraction of Residues with Alt Locs')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_5A_postqFit_violin.png')

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.swarmplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p2 = sns.swarmplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_5A_postqFit_swarm.png')


print('Median Number of Alt Conf Post-qFit (5A):')
print(close_RMSF_summary['per_altloc'].median())

print('Median Number of Alt Conf Post-qFit (Holo):')
print(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Holo']['per_altloc'].median())

print('Median Number of Alt Conf Post-qFit (Apo):')
print(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'].median())

same = (len((altloc_qFit_sum_m[altloc_qFit_sum_m['Apo_Holo_Diff_Num']==0]).index))
gain = (len((altloc_qFit_sum_m[altloc_qFit_sum_m['Apo_Holo_Diff_Num']>0]).index))
loss = (len((altloc_qFit_sum_m[altloc_qFit_sum_m['Apo_Holo_Diff_Num']<0]).index))
whole = [same, gain, loss]

same = (len((merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Diff_Num']==0]).index))
gain = (len((merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Diff_Num']>0]).index))
loss = (len((merged_close_sum_RMSF[merged_close_sum_RMSF['Apo_Holo_Diff_Num']<0]).index))
close = [same, gain, loss]

print(f'whole: {whole}')
print(f'Binging Site: {close}')

fig = plt.figure()
a = range(3)
barWidth = 0.4
r1 = np.arange(len(close))
r2 = [x + barWidth for x in r1]
p1_close = plt.bar(r1, whole, color='#9400D3', width=barWidth, label='Entire Protein')
p2 = plt.bar(r2, close, color='#8B0000', width=barWidth, label='Binding Site Residues')
plt.ylabel('Number of Pairs')
plt.xticks(a, ('Same','Increase', 'Decrease'))
plt.legend()
fig.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/test_preqfit_samegainloss.png')

