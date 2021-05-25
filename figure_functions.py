#Figure Functions
from __future__ import division
from scipy import stats
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os
import pandas as pd
import numpy as np
from matplotlib import cm

colors = ["#1b9e77","#d95f02", "#7570b3","#e7298a","#66a61e", "#e6ab02", "#666666"]
sns.set_palette(sns.color_palette(colors))
plt.style.use('seaborn-dark-palette')

#reference files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit/pair_docs/')
pairs = pd.read_csv('ligand_supplementary_table1_QCed_updated_200422.txt', sep=' ', header=None)
pairs = pairs.rename(columns={0: "Apo", 1: "Apo_Res", 2: "Holo", 3: "Holo_Res", 5:"Ligand"})
AH_pairs = pairs.drop_duplicates()

AH_key1 = AH_pairs[['Apo']]
AH_key2 = AH_pairs[['Holo']]
AH_key2.columns = ['PDB']
AH_key1.columns = ['PDB']
AH_key1['Apo_Holo'] = 'Apo'
AH_key2['Apo_Holo'] = 'Holo'
AH_key = pd.concat([AH_key1, AH_key2])

def paired_wilcoxon(holo_col, apo_col): 
	print(stats.wilcoxon(holo_col, apo_col))

	print('Holo Mean:') 
	print(holo_col.mean())

	print('Apo Mean:')
	print(apo_col.mean())

	print('Holo Median:')
	print(holo_col.median())

	print('Apo Median:')
	print(apo_col.median())



def ind_MannWhit(col_1, col_2):
	print(stats.mannwhitneyu(col_1, col_2))

	print('Column 1 Mean:')
	print(col_1.mean())

	print('Column 2 Mean:')
	print(col_2.mean())

	print('Column 1 Median:')
	print(col_1.median())

	print('Column 2 Median:')
	print(col_2.median())



def merge_apo_holo_df(df):
	df = df.merge(AH_key, on=['PDB'])
	print(df.head())
	print(len(df['PDB'].unique()))
	df_holo = df[df['Apo_Holo'] == 'Holo']
	print(len(df_holo['PDB'].unique()))
	df_apo = df[df['Apo_Holo'] == 'Apo']
	test = df_holo.merge(AH_pairs, left_on='PDB', right_on='Holo')
	print('test:')
	print(len(test.index))
	df_merged = test.merge(df_apo, left_on=['Apo', 'chain', 'resi'], right_on=['PDB', 'chain', 'resi'])  
	print(len(df_merged.index))
	df_merged = df_merged.drop_duplicates()
	print(len(df_merged['Holo'].unique()))
	return df_merged



def make_dist_plot_AH(holo_col, apo_col, x_label, y_label, title, out_name):
    fig = plt.figure()
    sns.distplot(holo_col, kde=False, label='Bound')
    sns.distplot(apo_col, kde=False, label='Unbound')
    plt.xlabel(x_label)
    plt.legend()
    plt.ylabel(y_label)
    plt.title(title)
    fig.savefig(out_name + '.png')

def make_boxenplot_AH(holo_col, apo_col, xlabel, ylabel, title, out_name):
	fig = plt.figure()
	x = range(2)
	f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
	p1 = sns.boxenplot(holo_col, orient='v', ax=axes[0]).set(xlabel = 'Holo', ylabel = ylabel)
	p2 = sns.boxenplot(apo_col, orient='v', ax=axes[1]).set(xlabel = 'Apo', ylabel = '')
	plt.savefig(out_name + '.png')


def make_boxenplot_chem(low_col, high_col, xlabel_low, xlabel_high, ylabel, out_name):
	fig = plt.figure()
	x = range(2)
	f, axes = plt.subplots(1, 2, sharey=True, sharex=True)
	p1 = sns.boxenplot(low_col, orient='v', ax=axes[0]).set(xlabel = xlabel_low, ylabel = ylabel)
	p2 = sns.boxenplot(high_col, orient='v', ax=axes[1]).set(xlabel = xlabel_high, ylabel = '')
	plt.savefig(out_name + '.png')

def rotamer_compare(subset_rotamer):
	rotamer_ = []
	print('rotamer_compare')
	for i in AH_pairs['Holo'].unique():
		apo_list = AH_pairs[AH_pairs['Holo'] == i]['Apo'].unique()
		tmp1 = subset_rotamer[subset_rotamer['PDB'] == i]
		if tmp1.empty == True:
			continue
		else:
			for a in apo_list:
				apo = subset_rotamer[subset_rotamer['PDB'] == a]
				if apo.empty == True:
					continue
				else:
					for r in apo['chain_resi'].unique():
						num_alt_loc_a= len(apo[apo['chain_resi'] == r].index)
						if r in tmp1['chain_resi'].unique():
							num_alt_loc_h= len(tmp1[tmp1['chain_resi'] == r].index)
							if num_alt_loc_h > 1 or num_alt_loc_a > 1:
								for alt in tmp1[tmp1['chain_resi'] == r]['altloc'].unique():
									for alt2 in apo[apo['chain_resi'] == r]['altloc'].unique():
										tmp_alt = tmp1[(tmp1['chain_resi'] == r) & (tmp1['altloc'] == alt)]
										apo_alt = apo[(apo['chain_resi'] == r) & (apo['altloc'] == alt2)]
										if (apo_alt['rotamer'].unique() == tmp_alt['rotamer'].unique()).all() == True:
											rotamer = 'Same'
										else:                                              
											rotamer = 'Different'
										rotamer_.append(tuple((i, a, r, num_alt_loc_a, num_alt_loc_h, apo_alt[apo_alt['chain_resi'] == r]['rotamer'].unique(), tmp_alt[tmp_alt['chain_resi'] == r]['rotamer'].unique(), rotamer, tmp_alt[tmp_alt['chain_resi'] == r]['AA'].astype(str))))
							else:
								if (apo[apo['chain_resi'] == r]['rotamer'].unique() == tmp1[tmp1['chain_resi'] == r]['rotamer'].unique()).all() == True:
									rotamer = 'Same'
								else:
									rotamer = 'Different'
								rotamer_.append(tuple((i, a, r, num_alt_loc_a, num_alt_loc_h, apo[apo['chain_resi'] == r]['rotamer'].unique(), tmp1[tmp1['chain_resi'] == r]['rotamer'].unique(), rotamer, tmp1[tmp1['chain_resi'] == r]['AA'].astype(str))))
	AH_rotamer = pd.DataFrame(rotamer_, columns =['Holo', 'Apo', 'chain_resi', 'Apo_altloc', 'Holo_altloc', 'Apo_Rotamer', 'Holo_Rotamer', 'Rotamer', 'AA']) #'Apo_altloc_name', 'Holo_altloc_name',
	return AH_rotamer

def rotamer_AH_summary(multi, AH_rotamer_singlealt):
	AH_rotamer_summary_multi = pd.DataFrame()
	n=0
	for i in multi['Holo'].unique():
		tmp = multi[multi['Holo']==i]
		for a in tmp['Apo'].unique():
			tmp2 = tmp[tmp['Apo']== a]
			for res in tmp2['chain_resi'].unique():
				AH_rotamer_summary_multi.loc[n, 'chain_resi'] = res
				AH_rotamer_summary_multi.loc[n, 'Holo'] = i
				AH_rotamer_summary_multi.loc[n, 'Apo'] = a
				if len(tmp2[tmp2['chain_resi']==res]['Rotamer'].unique())==1:
					#print(type(tmp2[tmp2['chain_resi']==res]['rotamer'].unique()))
					AH_rotamer_summary_multi.loc[n, 'Rotamer'] = str(tmp2[tmp2['chain_resi']==res]['Rotamer'].unique())
				else:
					AH_rotamer_summary_multi.loc[n, 'Rotamer'] = 'Same and Different'
				n += 1
	
	AH_rotamer_summary_single = pd.DataFrame()
	n = 1
	for i in AH_rotamer_singlealt['Holo'].unique():
		tmp = AH_rotamer_singlealt[AH_rotamer_singlealt['Holo'] == i]
		for a in tmp['Apo'].unique():
			tmp2 = tmp[tmp['Apo'] == a]
			for res in tmp2['chain_resi'].unique():
				AH_rotamer_summary_single.loc[n, 'Holo'] = i
				AH_rotamer_summary_single.loc[n, 'Apo'] = a
				AH_rotamer_summary_single.loc[n, 'chain_resi'] = res
				AH_rotamer_summary_single.loc[n, 'Rotamer'] = str(tmp2[tmp2['chain_resi']==res]['Rotamer'].unique())
				n += 1

	AH_rotamer_summary_single['Rotamer'] = AH_rotamer_summary_single.Rotamer.str.replace('[','')
	AH_rotamer_summary_single['Rotamer'] = AH_rotamer_summary_single.Rotamer.str.replace(']','')
	AH_rotamer_summary_single['Rotamer'] = AH_rotamer_summary_single.Rotamer.str.replace("\'", '')
	
	AH_rotamer_summary_multi['Rotamer'] = AH_rotamer_summary_multi.Rotamer.str.replace('[','')
	AH_rotamer_summary_multi['Rotamer'] = AH_rotamer_summary_multi.Rotamer.str.replace("\'", '')
	AH_rotamer_summary_multi['Rotamer'] = AH_rotamer_summary_multi.Rotamer.str.replace(']','')

	return AH_rotamer_summary_multi, AH_rotamer_summary_single

def rotamer_summary(subset_rotamer_holo):
	summary = pd.DataFrame()
	n = 1 
	for i in subset_rotamer_holo['PDB'].unique():
		for res in subset_rotamer_holo[subset_rotamer_holo['PDB'] == i]['chain_resi'].unique():
			summary.loc[n,'PDB'] = i
			summary.loc[n,'Residue'] = res
			tmp = subset_rotamer_holo[(subset_rotamer_holo['chain_resi'] == res) & (subset_rotamer_holo['PDB'] == i)]
			summary.loc[n,'num_altloc'] = len(tmp.index)
			if len(tmp[tmp['chain_resi'] == res]['rotamer'].unique()) == 1:
				summary.loc[n,'Rotamer_Status'] = 'same'
			else:
				summary.loc[n,'Rotamer_Status'] = 'different'
			n += 1
	return summary

def create_AH_key(AH_pairs):
    AH_key1 = AH_pairs[['Apo']]
    AH_key2 = AH_pairs[['Holo']]
    AH_key2.columns = ['PDB']
    AH_key1.columns = ['PDB']
    AH_key1['Apo_Holo'] = 'Apo'
    AH_key2['Apo_Holo'] = 'Holo'
    AH_key = pd.concat([AH_key1, AH_key2])
    return AH_key


def scatter_plot_with_linear_fit(x, y, slope=None, y_intercept=None, label=None, color=None):
    """
    :param x: an array
    :param y: an array
    :param slope: slope of the fitted line
    :param y_intercept: y-intercept of the fitted line
    If slope or y_intercept is not specified, these parameters will be generated
    by linear fit.
    :return: Pearson correlation coefficient and p-value
    """
    c_norm = matplotlib.colors.BoundaryNorm(boundaries=[0,1,2], ncolors=3)
    scatter = plt.scatter(x, y, alpha=0.3, label=label, c=color, cmap='tab10', norm=c_norm)#['#A7E30E', '#226363'])
    # fitted line
    if slope is None or y_intercept is None:
        slope, y_intercept = np.polyfit(x, y, 1)
        print(f'slope: {slope}')
    x_fit = np.linspace(np.min(x), np.max(x), 100)
    y_fit = slope * x_fit + y_intercept
    plt.plot(x_fit, y_fit, color='black', alpha=0.5)

def alt_loc_figure(pre_qfit, post_qfit, path):
	same = (len((pre_qfit[pre_qfit['Apo_Holo_Multi_Diff']==0]).index))
	gain = (len((pre_qfit[pre_qfit['Apo_Holo_Multi_Diff']>0]).index))
	loss = (len((pre_qfit[pre_qfit['Apo_Holo_Multi_Diff']<0]).index))
	total = same + gain + loss
	print('Pre Qfit:')
	print('Percent pre qFit Same:')
	print(same/total)
	print('Percent pre qFit Gain:')
	print(gain/total)
	print('Percent pre qFit Loss:')
	print(loss/total)


	same_qFit = (len((post_qfit[post_qfit['Apo_Holo_Multi_Diff']==0]).index))
	gain_qFit = (len((post_qfit[post_qfit['Apo_Holo_Multi_Diff']>0]).index))
	loss_qFit = (len((post_qfit[post_qfit['Apo_Holo_Multi_Diff']<0]).index))
	total_qFit = same_qFit + gain_qFit + loss_qFit
	print('Post Qfit:')
	print('Percent qFit Same:')
	print(same_qFit/total_qFit)
	print('Percent qFit Gain:')
	print(gain_qFit/total_qFit)
	print('Percent qFit Loss:')
	print(loss_qFit/total_qFit)

	#___________FIGURE___________
	fig = plt.figure()
	a = range(3)

	org = [same, gain, loss]
	qFit = [same_qFit, gain_qFit, loss_qFit]

	barWidth = 0.4
	r1 = np.arange(len(org))
	r2 = [x + barWidth for x in r1]

	p1_close = plt.bar(r1, org, color='#9400D3', width=barWidth, label='Refined PDB')
	p2 = plt.bar(r2, qFit, color='#8B0000', width=barWidth, label='qFit')

	plt.ylabel('Number of Pairs')
	plt.xticks(a, ('Same','Increase', 'Decrease'))
	plt.legend()
	fig.savefig(path)

def build_random(data, all_data):
	AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 
'Y']
	li = []
	for i in data['Holo'].unique():
		tmp = all_data[all_data['Holo'] == i]
		tmp_close = data[data['Holo'] == i]
		for a in AA: #look at each amino acid
			num = len(tmp_close[tmp_close['resn_x'] == a]) #determine the number of AA in 
			if num > 0: #only select if that type of AA is located in a close residue
				tmp = all_data[all_data['resn_x'] == a].sample(n=num) #choose random residues corresponding to the type and number of AA found in binding site
				li.append(tmp)

	random = pd.concat(li, axis=0, ignore_index=True)
	return random

