#packages
import pandas as pd
import numpy as np
import os
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import sys
from scipy import stats
import matplotlib.patches as mpatches
from figure_functions import *	

bfactor = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/bfactor.csv')
RMSF = pd.read_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_summary.csv')

RMSF_b = bfactor.merge(RMSF, on=['PDB', 'resi', 'chain'])
RMSF_b.to_csv('/Users/stephaniewankowicz/Downloads/qfit_paper/RMSF_bfactor.csv')

fig = plt.figure()
scatter_plot_with_linear_fit(RMSF_b['RMSF'], RMSF_b['Average_Bfactor'], label='')
plt.xlabel('Residue RMSF')
plt.ylabel('Residue B-Factors')
plt.savefig('/Users/stephaniewankowicz/Downloads/qfit_paper/bfactor_v_RMSF.png')
