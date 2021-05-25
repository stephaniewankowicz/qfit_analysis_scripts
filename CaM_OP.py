#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import seaborn as sns
import glob
import sys
import matplotlib.pyplot as plt
import matplotlib.pylab as plb
from scipy import stats
import matplotlib.patches as mpatches
from figure_functions import *	

#read in files
os.chdir('/Users/stephaniewankowicz/Downloads/qfit_paper/')
path=os.getcwd()


CaM = pd.read_csv('CaM_OP.csv')
#print(CaM.head())

#CaM_NMR = CaM['s2calc_NMR', 'residue']
#CaM_Xray = CaM['residue', 's2calc-xray', 's2calc-xray_2', 's2calc-xray_1xer']


CaM_melt = pd.melt(CaM, id_vars=['residue'], value_vars=['s2calc-xray', 's2calc-xray_2', 's2calc-xray_1xer'])
#print(CaM_melt)
CaM_melt = CaM_melt.dropna()

print(CaM)

#figure
fig = plt.figure()
ax = sns.boxplot(x=CaM_melt['residue'], y=CaM_melt['value']) #hue=CaM_melt['variable']
#ax = sns.lineplot(x=CaM['residue'], y=CaM['s2calc_NMR'], linewidth=2)
plt.xlabel('Residue')
plt.legend()
plt.ylabel('Order Parameter')
labels = ax.axes.get_xticklabels()
ax.axes.set_xticklabels(labels, rotation=45)
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/CaM_OP.png')
