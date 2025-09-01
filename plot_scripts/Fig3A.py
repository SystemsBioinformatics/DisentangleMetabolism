# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 13:30:36 2025

@author: mre283
"""


import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


plt.rcParams.update({'text.usetex': True, 'font.size': 24})

def format_label(index_label):
    # Define replacements and special cases
    replacements = {
        'Pputida': r'P.~putida',
        'Ecoli_carbonsources': r'E.~coli',
        'Ecoli_nitrogensources': r'E.~coli', 
        'Scerevisiae': r'S.~cerevisiae',
        'Gmetallireducens': r'G.~metallireducens',
        'Cljungdahlii': r'C.~ljungdahlii',
        'Mbarkeri': r'M.~barkeri',
        'Synechocystis': r'Synechocystis',
        'akg': r'alpha-ketoglutarate',
        'carbonmonoxidehydrogen':r'CO+H2',
        'carbonmonoxide':'CO',
        'carbondioxidehydrogen':r'CO2+H2'
    }
    
    # Replace special cases and remove additional parts
    for key, value in replacements.items():
        if key in index_label:
            index_label = index_label.replace(key, value)
    
    # General formatting for the rest of the string
    parts = index_label.split('\\')
    organism = parts[0]
    substrate = parts[1] if len(parts) > 1 else ''
    
    # Format organism in italics and concatenate with substrate
    formatted_label = f'${organism}$ {substrate}'
    formatted_label = formatted_label.replace('_', ' ')
    
    return formatted_label


df_dict = {}


directories = ['Results\Cljungdahlii', 'Results\Gmetallireducens', 'Results\Mbarkeri', 'Results\Scerevisiae', 'Results\Pputida','Results\Ecoli', 'Results\Synechocystis']

for directory in directories:
    for filename in os.scandir(directory):
        if filename.is_file():
            if 'nadh' in filename.path and not 'Synechocystis' in filename.path:
                continue
            if 'atp' in filename.path and 'Synechocystis' in filename.path:
                continue
            print(filename.path)
            df_dict[filename.path[8:-5]] = pd.read_excel(filename.path, sheet_name='atp_overall', index_col=0)
df_fraction = pd.DataFrame(columns=['ATP prod ana','ATP transferred'], index=df_dict.keys())
# df_fraction_new = pd.DataFrame(columns=['Y_ATP/E'], index=df_dict.keys())

    
    

for k, df in df_dict.items():
    df_fraction.loc[k, 'ATP prod ana'] = df.loc['ana', 'prod']
    df_fraction.loc[k, 'ATP transferred'] = df.loc['cat', 'net']
    df_fraction.loc[k, 'ATP fraction'] = df.loc['ana', 'prod']/df.loc['ana', 'cons']



df_fraction = df_fraction.sort_values(by=['ATP fraction'], ascending=False)


# Assuming `df_fraction` and `df_fraction_new` are your dataframes
cat_atp = df_fraction['ATP transferred']
ana_atp = df_fraction['ATP prod ana']

# Align the indices of both series
common_index = cat_atp.index.union(ana_atp.index)
cat_atp_aligned = cat_atp.reindex(common_index).fillna(0)
ana_atp_aligned = ana_atp.reindex(common_index).fillna(0)

# Sort by the sum of the two series
sum_values = cat_atp_aligned + ana_atp_aligned
sorted_index = sum_values.sort_values().index
cat_atp_sorted = cat_atp_aligned[sorted_index]
ana_atp_sorted = ana_atp_aligned[sorted_index]

# Apply the format_label function to the indices
cat_atp_sorted.index = [format_label(idx) for idx in cat_atp_sorted.index]
ana_atp_sorted.index = [format_label(idx) for idx in ana_atp_sorted.index]


# Set the font sizes
plt.rcParams.update({
    'font.size': 28,          # Default text size
    'axes.titlesize': 36,     # Title font size
    'axes.labelsize': 28,     # Axis labels font size
    'xtick.labelsize': 28,    # X-tick labels font size
    'ytick.labelsize': 32,    # Y-tick labels font size
    'legend.fontsize': 28,    # Legend font size
    'figure.titlesize': 32    # Figure title font size
})

fig = plt.figure(figsize=(44, 24), constrained_layout=True) 

ax3 = plt.subplot2grid((12, 12), (0, 8), colspan=3, rowspan=12)
ax4 = plt.subplot2grid((12, 12), (0, 11), colspan=1, rowspan=12)

n = len(cat_atp_sorted)
y_pos = np.arange(n) * 3  # Increase spacing multiplier (default ~1.0)
bar_height = 2.8

# Plot on ax3 (first part of the broken x-axis)
ax3.set_xlim(0, 160)
ax3.barh(y_pos, cat_atp_sorted, label='eCAT', color='tab:orange', height=bar_height, zorder=3)
ax3.barh(y_pos, ana_atp_sorted, left=cat_atp_sorted, label='pCAT and ANA', color='tab:blue', height=bar_height, zorder=3)

# Set custom y-ticks and labels
ax3.set_yticks(y_pos)
ax3.set_yticklabels(cat_atp_sorted.index)
ax3.grid(True, axis='x', zorder=0)
# ax3.axvline(x=34.7, color='black', linestyle='--', linewidth=3, label='$(Y_{X/ATP})^{-1}$ from Stouthamer',zorder=4)
ax3.set_xlabel('ATP requirement  [mmolATP/gDW]')

# Plot on ax4 (second part of the broken x-axis)
ax4.set_xlim(235, sum_values.max() + 20)
ax4.barh(y_pos, cat_atp_sorted, color='tab:orange', height=bar_height, zorder=3)
ax4.barh(y_pos, ana_atp_sorted, left=cat_atp_sorted, color='tab:blue', height=bar_height, zorder=3)
ax4.set_yticks(y_pos)
ax4.set_yticklabels([])  # Hide labels to avoid duplication
ax4.grid(True, axis='x', zorder=0)

ax3.legend(loc='upper center', bbox_to_anchor=(0.4, -0.05))

for spine in ax3.spines.values():
    spine.set_linewidth(2.5)

for spine in ax4.spines.values():
    spine.set_linewidth(2.5)
figpath = r'/Plots'
# plt.savefig(figpath + r'\Bars_allorgs_new.png', bbox_inches='tight', dpi=300)

plt.show()