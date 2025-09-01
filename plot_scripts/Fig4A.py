# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 13:54:30 2025

@author: mre283
"""


import numpy as np
import pandas as pd
import re
import sympy
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.patches as mpatches
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'


df_mceqs_aCP = pd.read_excel('Results/all_mceqs_Ecoli_anabolism.xlsx', header=0, index_col=0)

directory = 'Results\Ecoli'
df_dict = {}


for filename in os.scandir(directory):
    if filename.is_file():
        if 'nadh' in filename.path:
            continue
        print(filename.path)
        df_dict[filename.path[14:-5]] = pd.read_excel(filename.path, sheet_name='mceqs', index_col=0)

def extract_metabolites(expr):
    metabolites = {}
    for term in expr.as_ordered_terms():
        coeff, metabolite = term.as_coeff_Mul()
        metabolite_str = str(metabolite)
        coeff = np.float64(coeff)
        if metabolite_str in metabolites:
            metabolites[metabolite_str] += coeff
        else:
            metabolites[metabolite_str] = coeff
    return metabolites



df_mceqs_cat = pd.DataFrame()
df_mceqs_tot = pd.DataFrame()
df_mceqs_aS = pd.DataFrame()

for sim, df in df_dict.items():
    mceq_cat = df.loc['catabolic', 'overall']
    mceq_tot = df.loc['biomass', 'overall']
    mceq_aS = df.loc['anabolic', 'overall']
    
    if mceq_cat[0] == '[':
        mceq_cat = mceq_cat[1:-1]
    df_mceqs_cat.loc[sim, :] = np.zeros(len(df_mceqs_cat.columns))
    symp = sympy.parsing.sympy_parser.parse_expr(mceq_cat)
    metabolite_dict = extract_metabolites(symp)
    for metab, coeff in metabolite_dict.items():
        df_mceqs_cat.loc[sim, metab] = coeff
        
    if mceq_tot[0] == '[':
        mceq_tot = mceq_tot[1:-1]
    df_mceqs_tot.loc[sim, :] = np.zeros(len(df_mceqs_tot.columns))
    symp = sympy.parsing.sympy_parser.parse_expr(mceq_tot)
    metabolite_dict = extract_metabolites(symp)
    for metab, coeff in metabolite_dict.items():
        df_mceqs_tot.loc[sim, metab] = coeff
        
    if mceq_aS[0] == '[':
        mceq_aS = mceq_aS[1:-1]
    df_mceqs_aS.loc[sim, :] = np.zeros(len(df_mceqs_aS.columns))
    symp = sympy.parsing.sympy_parser.parse_expr(mceq_aS)
    metabolite_dict = extract_metabolites(symp)
    for metab, coeff in metabolite_dict.items():
        df_mceqs_aS.loc[sim, metab] = coeff
    
    
df_mceqs_cat.fillna(0, inplace=True)
df_mceqs_tot.fillna(0, inplace=True)
df_mceqs_aS.fillna(0, inplace=True)


new_i = []
for i in df_mceqs_aCP.index:
    
    if i[:7] == 'aerobic':
        new_i.append(i[8:])
    elif i == 'alphaketoglutarate':
        new_i.append('akg')
    elif i =='anaerobic_glucose':
        new_i.append('glucose_anaerobic')
    else:
        new_i.append(i)


df_mceqs_aCP.index = new_i




all_cols = df_mceqs_tot.columns.union(df_mceqs_aCP.columns).union(df_mceqs_cat.columns)

df_mceqs_tot = df_mceqs_tot.reindex(columns=all_cols, fill_value=0)
df_mceqs_aCP = df_mceqs_aCP.reindex(columns=all_cols, fill_value=0)
df_mceqs_cat = df_mceqs_cat.reindex(columns=all_cols, fill_value=0)

df_mceqs_PP = df_mceqs_tot - df_mceqs_aCP - df_mceqs_cat


import matplotlib.pyplot as plt

# Define a common index for alignment based on the union of all indices
common_index = df_mceqs_cat.index.union(
    df_mceqs_PP.index).union(
    df_mceqs_aCP.index
)

# Reindex all DataFrames to the common index, filling missing values with 0
df_mceqs_cat_reindexed = df_mceqs_cat.reindex(common_index).fillna(0)
df_mceqs_PP_reindexed = df_mceqs_PP.reindex(common_index).fillna(0)
df_mceqs_aCP_reindexed = df_mceqs_aCP.reindex(common_index).fillna(0)



# Define bar colors for ATP, NADPH, NADH
bar_colors = {'atp frac': 'tab:blue', 'nadph frac': 'tab:green', 'nadh frac': 'tab:orange'}
names = {'atp frac': 'ATP', 'nadph frac': 'NADPH', 'nadh frac': 'NADH'}

# Define bar colors for ATP, NADPH, NADH
bar_colors2 = {'M_atp_c': 'tab:blue', 'M_nadph_c': 'tab:green', 'M_nadh_c': 'tab:orange'}
names2 = {'M_atp_c': 'ATP', 'M_nadph_c': 'NADPH', 'M_nadh_c': 'NADH'}

# Subplot titles
titles = ['Energy catabolism', 'Precursor catabolism', 'Anabolism']

# DataFrames to plot
dataframes = [
    df_mceqs_cat_reindexed,
    df_mceqs_PP_reindexed,
    df_mceqs_aCP_reindexed,

]


# Define the sorting dictionary and corresponding colors
sort_dict = {
    'Sugars': ['xylose', 'sorbitol', 'glucose', 'glucose_anaerobic', 'galactose', 'ribose', 'fructose', 'gluconate'],
    'TCA Cycle Intermediates': ['succinate', 'malate', 'acetate', 'fumarate', 'citrate', 'akg'], 
    'lower glycolytic intermediates': ['pyruvate', 'lactate', 'glycerol'],
    'amino acids': ['tryptophan', 'glutamine', 'glutamate', 'aspartate', 'arginine', 'alanine']
}

category_colors = {
    'amino acids': 'tab:red',
    'lower glycolytic intermediates': 'tab:purple',
    'TCA Cycle Intermediates': 'tab:green',
    'Sugars': 'tab:orange'
    }

# Flatten the dictionary into a sorted list of substrates
sorted_substrates = [item for category in sort_dict.values() for item in category]

# Filter `common_index` based on `sorted_substrates`, keeping only those present in `common_index`
filtered_sorted_index = [substrate for substrate in sorted_substrates if substrate in common_index]

# Map substrates to their categories
substrate_to_category = {substrate: category 
                         for category, substrates in sort_dict.items() 
                         for substrate in substrates}

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(18, 11.5))  # Adjust figure size for 4 subplots

axes[1].tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

axes[2].tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

# axes[3].tick_params(
#     axis='y',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     left=False,      # ticks along the bottom edge are off
#     right=False,         # ticks along the top edge are off
#     labelleft=False) # labels along the bottom edge are off

twin_axes = []
# Loop through each DataFrame and plot
for ax, df, title in zip(axes, dataframes, titles):
    # Calculate bar positions
    indices = np.arange(len(filtered_sorted_index))  # Position of substrates on y-axis
    bar_width = 0.25  # Width of each bar
    offsets = [-bar_width, 0, bar_width]  # Offsets for grouped bars

    # Plot ATP on bottom x-axis
    atp_col = 'M_atp_c'
    ax.barh(indices + offsets[0], df[atp_col], height=bar_width, label=atp_col, color=bar_colors2[atp_col])
    
    # Create a twin x-axis for NADPH and NADH
    twin_ax = ax.twiny()
    twin_ax.barh(indices + offsets[1], df['M_nadph_c'], height=bar_width, label='M_nadph_c', color=bar_colors2['M_nadph_c'])
    twin_ax.barh(indices + offsets[2], df['M_nadh_c'], height=bar_width, label='M_nadh_c', color=bar_colors2['M_nadh_c'])

    # Set tick parameters for twin axes
    ax.tick_params(axis='x', labelsize=28, labelcolor=bar_colors2[atp_col])
    twin_ax.tick_params(axis='x', labelsize=28, labelcolor='black')
    
    if title == 'Precursor catabolism':  # Update as needed to match the second subplot
        x_min = min(df[atp_col].min(), df[['M_nadph_c', 'M_nadh_c']].min().min())  # Get shared min value
        x_max = max(df[atp_col].max(), df[['M_nadph_c', 'M_nadh_c']].max().max())  # Get shared max value
        
        ax.set_xlim(x_min -5, 100)  # Set the primary x-axis limits
        twin_ax.set_xlim((x_min-5)/6, 100/6)  # Set the secondary x-axis limits to match
        
        # ax.set_xlim(x_min -5, x_max + 5)  # Set the primary x-axis limits
        # twin_ax.set_xlim((x_min-5)/4, (x_max+5)/4)  # Set the secondary x-axis limits to match
    twin_axes.append(twin_ax)
    # Add category shading
    current_y = 0
    for category, substrates in sort_dict.items():
        count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
        ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
        current_y += count


    ax.set_title(title, fontsize=28, pad=20)  # Increase subplot title font size
    ax.axvline(x=0, color='black', linewidth=1)  # Add a vertical line at x=0 for reference
    
    lines = np.arange(-0.5, len(filtered_sorted_index)-0.5, 1)
    for lin in lines:
        ax.axhline(y=lin, color='black', linewidth=1)


    # Remove whitespace above and below the bars
    ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)

# Set y-axis labels
axes[0].set_yticks(indices)
axes[0].set_yticklabels(filtered_sorted_index, fontsize=28)  # Increase y-axis tick font size

axes[0].set_xlabel('$(Y_{X/ATP}^{eCAT})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])
axes[1].set_xlabel('$(Y_{X/ATP}^{pCAT})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])
axes[2].set_xlabel('$(Y_{X/ATP}^{ANA})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])

twin_axes[0].set_xlabel('$(Y_{X/NAD(P)H}^{eCAT})^{-1}$ [mmol/gDW]', fontsize=28, labelpad=15)
twin_axes[1].set_xlabel('$(Y_{X/NAD(P)H}^{pCAT})^{-1}$ [mmol/gDW]', fontsize=28, labelpad=15)
twin_axes[2].set_xlabel('$(Y_{X/NAD(P)H}^{ANA})^{-1}$ [mmol/gDW]', fontsize=28, labelpad=15)

axes[0].set_xlim(0, 100)
twin_axes[0].set_xlim(0, 100/6)


# Add a common legend for bar colors (energy carriers)
bar_legend = fig.legend(
    handles=[mpatches.Patch(color=color, label=names[col]) for col, color in bar_colors.items()],
    loc='lower center', bbox_to_anchor=(0.4, -0.1),  # Center below the plot
    frameon=False, fontsize=28, title_fontsize=28, title='Energy Carriers'
)

# Add a common legend for category shading
shading_legend_patches = [
    mpatches.Patch(color=color, label=category, alpha=0.2)
    for category, color in category_colors.items()
]
shading_legend = fig.legend(
    handles=shading_legend_patches,
    loc='lower center', bbox_to_anchor=(0.65, -0.15),  # Slightly lower than the first legend
    frameon=False, fontsize=28, title_fontsize=28, title="Carbon source pathway"
)


# Adjust layout to make space for the legends
plt.tight_layout(rect=[0, 0.1, 1, 1])  # Leave extra space at the bottom
plt.savefig(r'Plots/all_energycarriers_threeparts.png', dpi=300, bbox_inches='tight')
plt.show()
