# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 10:26:20 2024

@author: mre283
"""

import pandas as pd
import cbmpy
import matplotlib.pyplot as plt
import os

#%%

path = os.getcwd() + r'\Results\Ecoli_carbonsources\glucose.xlsx'
macromolecules = ['protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']

protein_df = pd.read_excel(path, sheet_name='protein', header=0, index_col=0)
lipid_df = pd.read_excel(path, sheet_name='lipid', header=0, index_col=0)
lps_df = pd.read_excel(path, sheet_name='lps', header=0, index_col=0)
rna_df = pd.read_excel(path, sheet_name='rna', header=0, index_col=0)
dna_df = pd.read_excel(path, sheet_name='dna', header=0, index_col=0)
murein_df = pd.read_excel(path, sheet_name='murein', header=0, index_col=0)
ion_df = pd.read_excel(path, sheet_name='ion', header=0, index_col=0)
cofactor_df = pd.read_excel(path, sheet_name='cofactor', header=0, index_col=0)
GAM_df = pd.read_excel(path, sheet_name='GAM', header=0, index_col=0)
dfs = [protein_df, lipid_df, lps_df, rna_df, dna_df, murein_df, ion_df, cofactor_df, GAM_df]

#%%

ATP_dem_df = pd.DataFrame(columns=['ATP_demand_tot', 'ATP_demand_cat'])

for i, df in enumerate(dfs):
    ATP_dem_df.loc[macromolecules[i], 'ATP_demand_tot'] = df.loc['R_ATP_demand', 'catabolic']+df.loc['R_ATP_demand', 'anabolic']
    ATP_dem_df.loc[macromolecules[i], 'ATP_demand_cat'] = df.loc['R_ATP_demand', 'catabolic']
    
    
#%%

path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism'
stouthamer_df = pd.read_excel(path+r'\Stouthamer_table_data_pythonreadable.xlsx', index_col=0)



# # Define the data for the DataFrame
# data = {
#     'E. coli': [42.2, 1.55, 9.58, 1.59, 1.26, 0.047, 0.14, 16.12],
#     'E. coli (Stouthamer method)': [10.15, -2.09, 4.27, 0.033, -0.55, 0.05, 0.03, 16.12],
#     'Stouthamer': [20.50, 0.14, 4.37, 1.05, 2.05, 5.21, 0, 1.39]
# }

# index = ['protein', 'lipid', 'rna', 'dna', 'polysaccharide', 'transport', 'production of cofactors', 'GAM']

# # Create the DataFrame
# stouthamer_df = pd.DataFrame(data, index=index)

# Plotting the stacked bar plot with bars for each column
ax = stouthamer_df.T.plot(kind='barh', stacked=True, zorder=3)

# Adding title and labels
plt.xlabel('$(Y_{X/ATP})^{-1}$')


# Adding grid and x=0 line
plt.grid(True, zorder=0)
plt.axvline(x=0, color='black', linewidth=1)
plt.legend(bbox_to_anchor=(1.05, 1))
ax.set_yticklabels([r'$\it{S.\ cerevisiae}$', r'$\it{E.\ coli}$', r'$\it{E.\ coli}$ (Stouthamer method)', r'Stouthamer'])
# Display the plot

path = r'C:\Users\mre283\Dropbox\maaike remijer paper separating catabolism from analbolism\Figures'
plt.savefig(r'Stouthamer_barplot.png', dpi=300, bbox_inches='tight')
plt.show()
