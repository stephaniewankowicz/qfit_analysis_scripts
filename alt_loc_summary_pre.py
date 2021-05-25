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
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200209.txt', sep=' ')
pairs.drop_duplicates(inplace=True)

ah_key = create_AH_key(pairs)

#_____________________PRE QFIT_________________________________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
path=os.getcwd()
all_files = glob.glob(path + "/*_RMSF.csv")
li = []
for filename in all_files:
    df = pd.read_csv(filename, sep=',')
    df['PDB'] = df['PDB_name']
    li.append(df)
altloc = pd.concat(li, axis=0, ignore_index=True)


altloc_sum = pd.DataFrame()
n = 1
for i in altloc['PDB'].unique():
    tmp = altloc[altloc['PDB'] == i]
    altloc_sum.loc[n, 'PDB'] = i
    altloc_sum.loc[n, 'Num_Residues'] = len(tmp.index)
    altloc_sum.loc[n, 'Num_Alt_Loc'] =  len(tmp[tmp['Num_Alt']>0].index)#len(tmp[tmp['RMSF']>0].index)
    altloc_sum.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
    altloc_sum.loc[n, 'Test'] = 1
    n += 1

altloc_sum['per_altloc'] = altloc_sum['Num_Alt_Loc'] / altloc_sum['Num_Residues']
altloc_sum['Num_Single'] = altloc_sum['Num_Residues'] - altloc_sum['Num_Alt_Loc']

altloc_sum = altloc_sum.merge(AH_key, on = ['PDB'])
altloc_sum_holo = altloc_sum[altloc_sum['Apo_Holo'] == 'Holo']
altloc_sum_apo = altloc_sum[altloc_sum['Apo_Holo'] == 'Apo']

test = altloc_sum_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
altloc_sum_m = test.merge(altloc_sum_apo, left_on=['Apo'], right_on=['PDB'])  
altloc_sum_m= altloc_sum_m.drop_duplicates()

altloc_sum_m['Apo_Holo_Multi_Diff'] = altloc_sum_m['per_altloc_x'] - altloc_sum_m['per_altloc_y']
altloc_sum_m['Apo_Holo_Multi_Diff_Num'] = altloc_sum_m['Num_Alt_Loc_x'] - altloc_sum_m['Num_Alt_Loc_y']

altloc_sum = altloc_sum.drop_duplicates()


##___________________SUMMARY STATS_____________________

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel=' Fraction of Residues with Alt Locs')
p1 = sns.boxenplot(altloc_sum[altloc_sum['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p2 = sns.boxenplot(altloc_sum[altloc_sum['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_fullprotein_preqFit.png')

fig = plt.figure()
#f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(altloc_qFit_sum['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
#p2 = sns.violinplot(close_RMSF_summary[close_RMSF_summary['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
sns.violinplot(y='per_altloc', x='Test', hue='Apo_Holo', data=altloc_sum, split=True, orient='v').set(xlabel='', ylabel='Fraction of Residues with Alt Locs')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_preqFit_violin.png')

fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
p1 = sns.swarmplot(altloc_sum[altloc_sum['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p2 = sns.swarmplot(altloc_sum[altloc_sum['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.title('Percent of Alternative Conformers in qFit Structures')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_preqFit_swarm.png')



print('Median Number of Alt Conf Pre-qFit:')
print(altloc_sum['per_altloc'].median())
print('Median Number of Alt Conf Pre-qFit (Holo):')
print(altloc_sum[altloc_sum['Apo_Holo']=='Holo']['per_altloc'].median())
print(altloc_sum[altloc_sum['Apo_Holo']=='Holo']['per_altloc'].mean())

print('Median Number of Alt Conf Pre-qFit (Apo):')
print(altloc_sum[altloc_sum['Apo_Holo']=='Apo']['per_altloc'].median())
print(altloc_sum[altloc_sum['Apo_Holo']=='Apo']['per_altloc'].mean())


#_________________CLOSE TO LIGAND_________________
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/210127_preqfit/')
path=os.getcwd()

all_files = glob.glob(path + "/*5.0_rmsf_subset.csv")

li = []

for filename in all_files:
    df = pd.read_csv(filename, header=0)
    li.append(df)

close_RMSF = pd.concat(li, axis=0, ignore_index=True)

close_RMSF['PDB'] = close_RMSF['PDB_name'].astype(str).str[0:4]
close_RMSF = close_RMSF.rename(columns={"Chain": "chain", "resseq":"resi"})
merged_close_RMSF_pre = merge_apo_holo_df(close_RMSF)


close_RMSF_summary_pre = pd.DataFrame()
n = 1
for i in close_RMSF['PDB'].unique():
    tmp = close_RMSF[close_RMSF['PDB'] == i]
    close_RMSF_summary_pre.loc[n, 'PDB'] = i
    close_RMSF_summary_pre.loc[n, 'Num_Residues'] = len(tmp.index)
    close_RMSF_summary_pre.loc[n, 'Num_Alt_Loc'] = len(tmp[tmp['Num_Alt']>0].index)
    close_RMSF_summary_pre.loc[n, 'Average_RMSF'] = tmp['RMSF'].mean()
    n += 1

close_RMSF_summary_pre['per_altloc'] = close_RMSF_summary_pre['Num_Alt_Loc']/close_RMSF_summary_pre['Num_Residues']

close_RMSF_summary_pre = close_RMSF_summary_pre.merge(AH_key)
close_RMSF_sum_holo_pre = close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Holo']
close_RMSF_sum_apo_pre = close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Apo']
test_pre = close_RMSF_sum_holo_pre.merge(AH_pairs, left_on='PDB', right_on='Holo')
merged_close_sum_RMSF_pre = test_pre.merge(close_RMSF_sum_apo_pre, left_on='Apo', right_on='PDB')

merged_close_sum_RMSF_pre['Percent_Holo_Close'] = merged_close_sum_RMSF_pre['Num_Alt_Loc_x']/ merged_close_sum_RMSF_pre['Num_Residues_x']
merged_close_sum_RMSF_pre['Percent_Apo_Close'] = merged_close_sum_RMSF_pre['Num_Alt_Loc_y']/ merged_close_sum_RMSF_pre['Num_Residues_y']
merged_close_sum_RMSF_pre['Apo_Holo_Multi_Diff'] = merged_close_sum_RMSF_pre['Percent_Holo_Close'] - merged_close_sum_RMSF_pre['Percent_Apo_Close']
merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num'] = merged_close_sum_RMSF_pre['Num_Alt_Loc_x'] - merged_close_sum_RMSF_pre['Num_Alt_Loc_y']
merged_close_sum_RMSF_pre.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/close_RMSF_summary_pre.csv')
merged_close_sum_RMSF_pre = merged_close_sum_RMSF_pre.drop_duplicates()

##___________________SUMMARY STATS_____________________
close_RMSF_summary_pre = close_RMSF_summary_pre[close_RMSF_summary_pre['PDB']!='6es8']
close_RMSF_summary_pre.to_csv('close_RMSF_summary_pre.csv')
fig = plt.figure()
f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
#p1 = sns.boxenplot(close_RMSF_summary_pre['per_altloc'], orient='v', ax=axes[0]).set(xlabel='All', ylabel='% Residues with Alt Loc')
p2 = sns.boxenplot(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel='Apo', ylabel='Fraction of Residues with Alt Locs')
p3 = sns.boxenplot(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel='Holo', ylabel='')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/NumberAltConfBoundvUnbound_5A_preqFit.png')

print('Median Number of Alt Conf Pre-qFit (5A):')
print(close_RMSF_summary_pre['per_altloc'].median())

print('Median Number of Alt Conf Pre-qFit (Holo):')
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Holo']['per_altloc'].median())
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Holo']['per_altloc'].mean())
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Holo']['per_altloc'].describe())

print('Median Number of Alt Conf Pre-qFit (Apo):')
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Apo']['per_altloc'].median())
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Apo']['per_altloc'].mean())
print(close_RMSF_summary_pre[close_RMSF_summary_pre['Apo_Holo']=='Apo']['per_altloc'].describe())


#__________PRE QFIT COMPARE BINDING VERSUS ALL PROTEIN
plt.figure()
x = range(2)
f, axes = plt.subplots(1, 4, sharey=True, sharex=False)
p1 = sns.boxenplot(altloc_sum[altloc_sum['Apo_Holo']=='Holo']['per_altloc'], orient='v', ax=axes[0]).set(xlabel = 'Holo', ylabel = '')
p2 = sns.boxenplot(altloc_sum[altloc_sum['Apo_Holo']=='Apo']['per_altloc'], orient='v', ax=axes[1]).set(xlabel = 'Apo', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/test_preqfit.png')

altloc_sum_m = altloc_sum_m[altloc_sum_m['Holo']!='5tyk']
merged_close_sum_RMSF_pre = merged_close_sum_RMSF_pre[merged_close_sum_RMSF_pre['Holo']!='5tyk']


plt.figure()
x = range(2)
f, axes = plt.subplots(1, 2, sharey=True, sharex=False)
p1 = sns.boxenplot(altloc_sum_m['Apo_Holo_Multi_Diff_Num'], orient='v', ax=axes[0]).set(xlabel = 'Whole Protein', ylabel = '')
p2 = sns.boxenplot(merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num'], orient='v', ax=axes[1]).set(xlabel = 'Binding Sites', ylabel = '')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/test_preqfit.png')


plt.figure()
p1 = sns.distplot(altloc_sum_m['Apo_Holo_Multi_Diff_Num'], label= 'Whole Protein')
p2 = sns.distplot(merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num'], label="Binding Sites")
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/test_preqfit_distplot.png')


same = (len((altloc_sum_m[altloc_sum_m['Apo_Holo_Multi_Diff_Num']==0]).index))
gain = (len((altloc_sum_m[altloc_sum_m['Apo_Holo_Multi_Diff_Num']>0]).index))
loss = (len((altloc_sum_m[altloc_sum_m['Apo_Holo_Multi_Diff_Num']<0]).index))
whole = [same, gain, loss]

same = (len((merged_close_sum_RMSF_pre[merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num']==0]).index))
gain = (len((merged_close_sum_RMSF_pre[merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num']>0]).index))
loss = (len((merged_close_sum_RMSF_pre[merged_close_sum_RMSF_pre['Apo_Holo_Diff_Num']<0]).index))
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

#_________________________________###

#SUBSET TO ONLY SEQUENCE UNIQUE
# print('SUBSET TO ONLY UNIQUE SEQUENCE:')
# PDB_to_keep = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/PDB_to_keep.csv')

# altloc_qFit_sub = altloc_qFit_sum_m[altloc_qFit_sum_m['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# altloc_sum_sub = altloc_sum_m[altloc_sum_m['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# merged_close_sub = merged_close_sum_RMSF[merged_close_sum_RMSF['Holo'].isin(PDB_to_keep['PDB_to_keep'])]
# merged_close_pre_sub = merged_close_sum_RMSF_pre[merged_close_sum_RMSF_pre['Holo'].isin(PDB_to_keep['PDB_to_keep'])]

# print('All Protein Subset:')
# alt_loc_figure(altloc_qFit_sub, altloc_sum_sub,'/Users/stephaniewankowicz/Downloads/qfit_paper/GainSameLoss_subset.png')

# print('5A Subset:')
# alt_loc_figure(merged_close_sub, merged_close_pre_sub,'/Users/stephaniewankowicz/Downloads/qfit_paper/GainSameLoss_within5A_subset.png')
