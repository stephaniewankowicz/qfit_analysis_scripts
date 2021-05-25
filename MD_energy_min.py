import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns

pot = pd.read_csv('/Users/stephaniewankowicz/Downloads/potential_7kqo.csv', sep=' ',header=None)
pot = pot.rename(columns={0: "Step", 1: "Potential"})
print(pot.head())

fig = plt.figure()
ax = sns.lineplot(y=pot['Potential'], x=pot['Step']) #hue=CaM_melt['variable']
# #ax = sns.lineplot(x=CaM['residue'], y=CaM['s2calc_NMR'], linewidth=2)
plt.xlabel('Step')
# plt.legend()
plt.ylabel('Potential Energy')
# labels = ax.axes.get_xticklabels()
# ax.axes.set_xticklabels(labels, rotation=45)
#plt.show()
plt.savefig('/Users/stephaniewankowicz/Downloads/7kqo_pot_energy.png')
