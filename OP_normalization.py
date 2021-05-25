#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import scipy
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from scipy import stats
import matplotlib.patches as mpatches
from matplotlib import lines as mpl_lines


os.chdir('/Users/stephaniewankowicz/Downloads/OP_normalization/')
path=os.getcwd()

all_files = glob.glob(path + "/Wilson*")
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',')
    li.append(df)

Wilson = pd.concat(li, axis=0, ignore_index=True)
Wilson.to_csv('Wilson.csv')

all_files = glob.glob(path + "/*_methyl.out")
li = []
for filename in all_files:
    df = pd.read_csv(filename, index_col=None, sep=',')
    df['PDB'] = filename[53:57]
    li.append(df)

order = pd.concat(li, axis=0, ignore_index=True)


summary = pd.DataFrame()
n = 1
for i in order['PDB'].unique():
    summary.loc[n, 'PDB'] = i
    subset = order[order['PDB']==i]
    summary.loc[n, 'Res']= Wilson[Wilson['PDB']==i]['Resolution'].unique()
    Res = (Wilson[Wilson['PDB']==i]['Resolution'].unique())*0.1
    #summary['Res'] = subset_w['Resolution'].unique()
    summary.loc[n, 's2calc_median'] = subset['s2calc'].median()
    summary.loc[n, 's2calc_median_normal'] = (subset['s2calc']/Res).median()
    summary.loc[n, 's2calc_mean'] = subset['s2calc'].mean()
    n+=1

#print(summary)
order_Wilson=Wilson.merge(order, on='PDB')
order_Wilson = order_Wilson[order_Wilson['resi']<=200]

print(order_Wilson.head())


sns.scatterplot(order_Wilson[order_Wilson['PDB']=='5lan']['s2ortho'], order_Wilson[order_Wilson['PDB']!='5lan']['s2ortho'])

# for norm in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.9, 1, 1.1, 1.2, 1.8, 2]:
#      order_Wilson['s2ortho_norm'] = order_Wilson['s2ortho']*order_Wilson['Resolution']*norm   
#      print(norm)
#      s2ortho_all = []
#      for i in order_Wilson['resi']:
#          #range_s2ortho.append(order_Wilson[order_Wilson['resi']==i]['s2ortho_norm'].quantile(0.75) - order_Wilson[order_Wilson['resi']==i]['s2ortho_norm'].quantile(0.25))
#          sd = np.std(order_Wilson[order_Wilson['resi']==i]['s2ortho_norm']) #sd_s2ortho.append
#          mean = np.mean(order_Wilson[order_Wilson['resi']==i]['s2ortho_norm'])
# #         #print('sd')
# #         #print(sd)
# #         #print(mean)
#          z_s2ortho = []
#          for p in order_Wilson[order_Wilson['resi']==i]['PDB']:
# #         	#print('zscore:') sd_s2ortho.append(
# #         	#print(order_Wilson[(order_Wilson['resi']==i)&(order_Wilson['PDB']==p)])
#          	s2ortho = ((order_Wilson[(order_Wilson['resi']==i)&(order_Wilson['PDB']==p)]['s2ortho_norm'].values[0]-mean)/sd)
# #         	#print(s2ortho)
#          	z_s2ortho.append(s2ortho)
#      	 s2ortho_all = np.median(z_s2ortho) 
#      print(np.median(s2ortho_all))

from sklearn.preprocessing import MinMaxScaler
scaler = MinMaxScaler(feature_range=(0,1))
scaler.fit(OW_melt_res['s2calc_x'])
OW_melt_res['normalized_s2calc_x'] = scaler.transform(OW_melt_res['s2calc_x'])
