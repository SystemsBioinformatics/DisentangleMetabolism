# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 09:35:52 2024

@author: mre283
"""

import cbmpy
import separate_cat_ana as sca
import numpy as np
import pandas as pd
import re
import sympy
import os
import matplotlib.pyplot as plt

#%%


result_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\Precursor_simulations_Ecoli'
df_mceqs_aCP = pd.read_excel(result_path+'\\all_mceqs.xlsx', header=0, index_col=0)



directory = 'Results\Ecoli_carbonsources'
df_dict = {}

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        if 'NRC' in filename.path:
            continue
        if 'respirofermentative' in filename.path:
            continue
        if 'nadh' in filename.path:
            continue
        df_dict[filename.path[28:-5]] = pd.read_excel(filename.path, sheet_name='mceqs', index_col=0)

directory = 'Results\Ecoli_nitrogensources'

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.path[30:-5] == 'glucose_nitrate' or filename.path[30:-5] == 'glucose_nitrate_nadh':
        continue
    if filename.is_file():
        df_dict[filename.path[30:-5]] = pd.read_excel(filename.path, sheet_name='mceqs', index_col=0)
        
#%%

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

#%%
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



#%%

all_cols = df_mceqs_tot.columns.union(df_mceqs_aCP.columns).union(df_mceqs_cat.columns)

df_mceqs_tot = df_mceqs_tot.reindex(columns=all_cols, fill_value=0)
df_mceqs_aCP = df_mceqs_aCP.reindex(columns=all_cols, fill_value=0)
df_mceqs_cat = df_mceqs_cat.reindex(columns=all_cols, fill_value=0)

df_mceqs_PP = df_mceqs_tot - df_mceqs_aCP - df_mceqs_cat

#%%

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

# Create subplots
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Plot each reindexed DataFrame's column
df_mceqs_cat_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0])
df_mceqs_PP_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1])
df_mceqs_aCP_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2])

# Add titles to each subplot
ax[0].set_title("Energy conservation")
ax[1].set_title("Precursor production")
ax[2].set_title("Anabolism from carbon precursors")

ax[0].set_xlabel('$(Y_{X/ATP}^{EC})^{-1}$')
ax[1].set_xlabel('$(Y_{X/ATP}^{PP})^{-1}$')
ax[2].set_xlabel('$(Y_{X/ATP}^{ana})^{-1}$')

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(common_index)))  # Set y-ticks to match the length of the common index
    axis.set_yticklabels(common_index)        # Set y-tick labels to the common index
    axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Adjust layout and show plot
plt.suptitle('ATP')
plt.tight_layout()
plt.show()


#%%
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

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

# Reindex the DataFrames to the filtered and sorted index
df_mceqs_cat_reindexed = df_mceqs_cat.reindex(filtered_sorted_index).fillna(0)
df_mceqs_PP_reindexed = df_mceqs_PP.reindex(filtered_sorted_index).fillna(0)
df_mceqs_aCP_reindexed = df_mceqs_aCP.reindex(filtered_sorted_index).fillna(0)

# Create subplots
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Calculate shading ranges
current_y = 0
shading_ranges = []
for category, substrates in sort_dict.items():
    count = sum(1 for substrate in substrates if substrate in common_index)
    shading_ranges.append((current_y, current_y + count, category_colors[category]))
    current_y += count

# Add background shading
for axis in ax:
    for start, end, color in shading_ranges:
        axis.axhspan(start - 0.5, end - 0.5, color=color, alpha=0.2)  # Add shaded band

# Plot each reindexed DataFrame's column
df_mceqs_cat_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0], color='tab:blue')
df_mceqs_PP_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1], color='tab:blue')
df_mceqs_aCP_reindexed['M_atp_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2], color='tab:blue')

# Add titles to each subplot
ax[0].set_title("Energy Conservation")
ax[1].set_title("Precursor Production")
ax[2].set_title("Anabolism from Carbon Precursors")

ax[0].set_xlabel('$(Y_{X/ATP}^{EC})^{-1}$ [mmol/gDW]')
ax[1].set_xlabel('$(Y_{X/ATP}^{PP})^{-1}$ [mmol/gDW]')
ax[2].set_xlabel('$(Y_{X/ATP}^{ana})^{-1}$ [mmol/gDW]')

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(filtered_sorted_index)))  # Set y-ticks to match the length of the filtered and sorted index
    axis.set_yticklabels(filtered_sorted_index)         # Set y-tick labels to the custom labels
    # axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Add legend for shading
legend_patches = [mpatches.Patch(color=color, label=category, alpha=0.2) 
                  for category, color in category_colors.items()]
fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(0.85, 0.5), title="Categories", frameon=False)

# Adjust layout and show plot
plt.suptitle('ATP')
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for the legend
savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\_02 Presentaties posters abstracts\241209_generallab_figs'
# plt.savefig(savepath+r'\ATP_threeparts_Ecoli.png', dpi=300, bbox_inches='tight')
plt.show()


# according to chatGPT: Substrate bars are grouped and aligned with the shaded regions, making the visualization informative and aesthetically pleasing.



#%%

# Create subplots for NADPH
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Add background shading
for axis in ax:
    for start, end, color in shading_ranges:
        axis.axhspan(start - 0.5, end - 0.5, color=color, alpha=0.2)  # Add shaded band

# Plot each reindexed DataFrame's column
df_mceqs_cat_reindexed['M_nadph_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0], color='tab:blue')
df_mceqs_PP_reindexed['M_nadph_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1], color='tab:blue')
df_mceqs_aCP_reindexed['M_nadph_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2], color='tab:blue')

# Add titles to each subplot
ax[0].set_title("Energy Conservation")
ax[1].set_title("Precursor Production")
ax[2].set_title("Anabolism from Carbon Precursors")

ax[0].set_xlabel('$(Y_{X/NADPH}^{EC})^{-1}$')
ax[1].set_xlabel('$(Y_{X/NADPH}^{PP})^{-1}$')
ax[2].set_xlabel('$(Y_{X/NADPH}^{ana})^{-1}$')

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(filtered_sorted_index)))  # Set y-ticks to match the length of the filtered and sorted index
    axis.set_yticklabels(filtered_sorted_index)         # Set y-tick labels to the custom labels
    # axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Add legend for shading
legend_patches = [mpatches.Patch(color=color, label=category, alpha=0.2) 
                  for category, color in category_colors.items()]
fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(0.85, 0.5), title="Categories", frameon=False)

# Adjust layout and show plot
plt.suptitle('NADPH')
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for the legend
plt.show()


#%%

# Create subplots for NADH
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Add background shading
for axis in ax:
    for start, end, color in shading_ranges:
        axis.axhspan(start - 0.5, end - 0.5, color=color, alpha=0.2)  # Add shaded band

# Plot each reindexed DataFrame's column
df_mceqs_cat_reindexed['M_nadh_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0], color='tab:blue')
df_mceqs_PP_reindexed['M_nadh_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1], color='tab:blue')
df_mceqs_aCP_reindexed['M_nadh_c'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2], color='tab:blue')

# Add titles to each subplot
ax[0].set_title("Energy Conservation")
ax[1].set_title("Precursor Production")
ax[2].set_title("Anabolism from Carbon Precursors")

ax[0].set_xlabel('$(Y_{X/NADH}^{EC})^{-1}$')
ax[1].set_xlabel('$(Y_{X/NADH}^{PP})^{-1}$')
ax[2].set_xlabel('$(Y_{X/NADH}^{ana})^{-1}$')

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(filtered_sorted_index)))  # Set y-ticks to match the length of the filtered and sorted index
    axis.set_yticklabels(filtered_sorted_index)         # Set y-tick labels to the custom labels
    # axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Add legend for shading
legend_patches = [mpatches.Patch(color=color, label=category, alpha=0.2) 
                  for category, color in category_colors.items()]
fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(0.85, 0.5), title="Categories", frameon=False)

# Adjust layout and show plot
plt.suptitle('NADH')
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for the legend
plt.show()

#%%
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Define bar colors for ATP, NADPH, NADH
bar_colors = {'M_atp_c': 'tab:blue', 'M_nadph_c': 'tab:green', 'M_nadh_c': 'tab:orange'}
names = {'M_atp_c': 'ATP', 'M_nadph_c': 'NADPH', 'M_nadh_c': 'NADH'}

# Subplot titles
titles = ['Energy Conservation', 'Precursor Production', 'Anabolism from Carbon Precursors']

# DataFrames to plot
dataframes = [df_mceqs_cat_reindexed, df_mceqs_PP_reindexed, df_mceqs_aCP_reindexed]

# Create subplots
fig, axes = plt.subplots(1, 3, figsize=(30, 11.5))

axes[0].set_xlabel('$(Y_{X/ec}^{EC})^{-1}$', fontsize=20)
axes[1].set_xlabel('$(Y_{X/ec}^{PP})^{-1}$', fontsize=20)
axes[2].set_xlabel('$(Y_{X/ec}^{ana})^{-1}$', fontsize=20)

# Loop through each DataFrame and plot
for ax, df, title in zip(axes, dataframes, titles):
    # Calculate bar positions
    indices = np.arange(len(filtered_sorted_index))  # Position of substrates on y-axis
    bar_width = 0.25  # Width of each bar
    offsets = [-bar_width, 0, bar_width]  # Offsets for grouped bars

    # Plot bars for each energy carrier
    for i, (col, color) in enumerate(bar_colors.items()):
        ax.barh(indices + offsets[i], df[col], height=bar_width, label=col, color=color)

    # Add category shading
    current_y = 0
    for category, substrates in sort_dict.items():
        count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
        ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
        current_y += count

    # Set y-axis labels
    ax.set_yticks(indices)
    ax.set_yticklabels(filtered_sorted_index, fontsize=20)  # Increase y-axis tick font size
    ax.set_title(title, fontsize=20)  # Increase subplot title font size
    ax.axvline(x=0, color='black', linewidth=1)  # Add a vertical line at x=0 for reference

    ax.tick_params(axis='x', labelsize=20)
    # Add xlabel for each subplot
    # ax.set_xlabel('Inverse Yield [mmol/gDW]', fontsize=20)  # Increase x-axis label font size

    # Remove whitespace above and below the bars
    ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)

# Add a common legend for bar colors (energy carriers)
bar_legend = fig.legend(
    handles=[mpatches.Patch(color=color, label=names[col]) for col, color in bar_colors.items()],
    loc='upper right', title='Energy Carriers', bbox_to_anchor=(0.95, 0.75), frameon=False, fontsize=20, title_fontsize=20
)

# Add a common legend for category shading
shading_legend_patches = [
    mpatches.Patch(color=color, label=category, alpha=0.2)
    for category, color in category_colors.items()
]
fig.legend(
    handles=shading_legend_patches,
    loc='center right', title="Categories", bbox_to_anchor=(1.05, 0.5), frameon=False, fontsize=20, title_fontsize=20
)

# Adjust layout and show the plot
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legends
savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\_02 Presentaties posters abstracts\241209_generallab_figs'
# plt.savefig(savepath+r'\all_energycarriers_threeparts_Ecoli.png', dpi=300, bbox_inches='tight')
plt.show()




#%%

directory = 'Results\Precursor_simulations_Ecoli'
df_dict = {}

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        if 'mceqs' in filename.path:
            continue
        print(filename)
        df_dict[filename.path[36:-5]] = pd.read_excel(filename.path, sheet_name='info', index_col=0)


#%%


df_atpprod = pd.DataFrame(columns=['atp prod', 'nadph prod', 'nadh prod', 'atp frac', 'nadph frac', 'nadh frac'], index=list(df_dict.keys()))

for i in list(df_dict.keys()):
    df_atpprod.loc[i, 'atp prod'] = df_dict[i].loc['atp production', 0]
    df_atpprod.loc[i, 'nadph prod'] = df_dict[i].loc['nadph production', 0]
    df_atpprod.loc[i, 'nadh prod'] = df_dict[i].loc['nadh production', 0]
    
    df_atpprod.loc[i, 'atp frac'] = df_dict[i].loc['atp production', 0]/df_dict[i].loc['atp consumption', 0]
    df_atpprod.loc[i, 'nadph frac'] = df_dict[i].loc['nadph production', 0]/df_dict[i].loc['nadph consumption', 0]
    df_atpprod.loc[i, 'nadh frac'] = df_dict[i].loc['nadh production', 0]/df_dict[i].loc['nadh consumption', 0]


#%%

# Create subplots
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Plot each reindexed DataFrame's column
df_atpprod['atp frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0])
df_atpprod['nadph frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1])
df_atpprod['nadh frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2])

# Add titles to each subplot
ax[0].set_title("ATP")
ax[1].set_title("NADPH")
ax[2].set_title("NADH")

ax[0].set_xlabel('ATP fraction')
ax[1].set_xlabel('NADPH fraction')
ax[2].set_xlabel('NADH fraction')

index = list(df_atpprod.index)

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(index)))  # Set y-ticks to match the length of the common index
    axis.set_yticklabels(index)        # Set y-tick labels to the common index
    
    axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Adjust layout and show plot
plt.suptitle('X produced in aCP/X consumed in aCP')
plt.tight_layout()
plt.show()

#%%

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

new_i = []
for i in df_atpprod.index:
    
    if i[:7] == 'aerobic':
        new_i.append(i[8:])
    elif i == 'alphaketoglutarate':
        new_i.append('akg')
    elif i =='anaerobic_glucose':
        new_i.append('glucose_anaerobic')
    else:
        new_i.append(i)


df_atpprod.index = new_i


# Ensure df_atpprod uses the same sorted index
df_atpprod_reindexed = df_atpprod.reindex(filtered_sorted_index).fillna(0)

# Create subplots
fig, ax = plt.subplots(1, 3, figsize=(15, 5))  # Adjust figsize as needed

# Calculate shading ranges
current_y = 0
shading_ranges = []
for category, substrates in sort_dict.items():
    count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
    shading_ranges.append((current_y, current_y + count, category_colors[category]))
    current_y += count

# Add background shading
for axis in ax:
    for start, end, color in shading_ranges:
        axis.axhspan(start - 0.5, end - 0.5, color=color, alpha=0.2)  # Add shaded band

# Plot each reindexed DataFrame's column
df_atpprod_reindexed['atp frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[0])
df_atpprod_reindexed['nadph frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[1])
df_atpprod_reindexed['nadh frac'].plot(kind='barh', stacked=False, zorder=3, ax=ax[2])

# Add titles to each subplot
ax[0].set_title("ATP")
ax[1].set_title("NADPH")
ax[2].set_title("NADH")

ax[0].set_xlabel('ATP fraction')
ax[1].set_xlabel('NADPH fraction')
ax[2].set_xlabel('NADH fraction')

# Set consistent y-axis labels
for axis in ax:
    axis.set_yticks(range(len(filtered_sorted_index)))  # Set y-ticks to match the length of the filtered and sorted index
    axis.set_yticklabels(filtered_sorted_index)         # Set y-tick labels to the custom labels
    axis.grid(True, zorder=0)
    axis.axvline(x=0, color='black', linewidth=1)

# Add legend for shading
legend_patches = [mpatches.Patch(color=color, label=category, alpha=0.2) 
                  for category, color in category_colors.items()]
fig.legend(handles=legend_patches, loc='center left', bbox_to_anchor=(0.85, 0.5), title="Categories", frameon=False)

# Adjust layout and show plot
plt.suptitle('X produced in aCP/X consumed in aCP')
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Adjust layout to make space for the legend
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

# Define bar colors for ATP, NADPH, NADH
bar_colors = {'atp frac': 'tab:blue', 'nadph frac': 'tab:green', 'nadh frac': 'tab:orange'}
names = {'atp frac': 'ATP', 'nadph frac': 'NADPH', 'nadh frac': 'NADH'}

# Define bar colors for ATP, NADPH, NADH
bar_colors2 = {'M_atp_c': 'tab:blue', 'M_nadph_c': 'tab:green', 'M_nadh_c': 'tab:orange'}
names2 = {'M_atp_c': 'ATP', 'M_nadph_c': 'NADPH', 'M_nadh_c': 'NADH'}

# Subplot titles
titles = ['Energy Conservation', 'Precursor Production', 'Anabolism from Carbon Precursors', 'X produced/X consumed']

# DataFrames to plot
dataframes = [
    df_mceqs_cat_reindexed,
    df_mceqs_PP_reindexed,
    df_mceqs_aCP_reindexed,
    df_atpprod_reindexed,
]

# Create subplots
fig, axes = plt.subplots(1, 4, figsize=(40, 11.5))  # Adjust figure size for 4 subplots

# Loop through each DataFrame and plot
for ax, df, title in zip(axes[:-1], dataframes[:-1], titles[:-1]):
    # Calculate bar positions
    indices = np.arange(len(filtered_sorted_index))  # Position of substrates on y-axis
    bar_width = 0.25  # Width of each bar
    offsets = [-bar_width, 0, bar_width]  # Offsets for grouped bars

    # Plot bars for each energy carrier
    for i, (col, color) in enumerate(bar_colors2.items()):
        ax.barh(indices + offsets[i], df[col], height=bar_width, label=col, color=color)

    # Add category shading
    current_y = 0
    for category, substrates in sort_dict.items():
        count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
        ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
        current_y += count

    # Set y-axis labels
    ax.set_yticks(indices)
    ax.set_yticklabels(filtered_sorted_index, fontsize=20)  # Increase y-axis tick font size
    ax.set_title(title, fontsize=20)  # Increase subplot title font size
    ax.axvline(x=0, color='black', linewidth=1)  # Add a vertical line at x=0 for reference

    ax.tick_params(axis='x', labelsize=20)
    # Add xlabel for each subplot
    # ax.set_xlabel('Inverse Yield [mmol/gDW]', fontsize=20)  # Increase x-axis label font size

    # Remove whitespace above and below the bars
    ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)


# Loop through each DataFrame and plot
ax, df, title = axes[-1], dataframes[-1], titles[-1]
# Calculate bar positions
indices = np.arange(len(filtered_sorted_index))  # Position of substrates on y-axis
bar_width = 0.2  # Width of each bar
offsets = [-bar_width, 0, bar_width]  # Offsets for grouped bars

# Plot bars for each energy carrier
for i, (col, color) in enumerate(bar_colors.items()):
    if col in df.columns:  # Ensure the column exists in the DataFrame
        ax.barh(indices + offsets[i], df[col], height=bar_width, label=col, color=color)
        
        

# Add category shading
current_y = 0
for category, substrates in sort_dict.items():
    count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
    ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
    current_y += count

# Set y-axis labels
ax.set_yticks(indices)
ax.set_yticklabels(filtered_sorted_index, fontsize=20)  # Increase y-axis tick font size
ax.set_title(title, fontsize=20)  # Increase subplot title font size
ax.axvline(x=0, color='black', linewidth=1)  # Add a vertical line at x=0 for reference

ax.tick_params(axis='x', labelsize=20)
ax.set_xlabel('Fraction' if title == 'X produced/X consumed' else 'Inverse Yield [mmol/gDW]', fontsize=20)

# Remove whitespace above and below the bars
ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)

# Add a common legend for bar colors (energy carriers)
bar_legend = fig.legend(
    handles=[mpatches.Patch(color=color, label=names[col]) for col, color in bar_colors.items()],
    loc='upper right', title='Energy Carriers', bbox_to_anchor=(0.95, 0.75), frameon=False, fontsize=20, title_fontsize=20
)

# Add a common legend for category shading
shading_legend_patches = [
    mpatches.Patch(color=color, label=category, alpha=0.2)
    for category, color in category_colors.items()
]
fig.legend(
    handles=shading_legend_patches,
    loc='center right', title="Categories", bbox_to_anchor=(1.05, 0.5), frameon=False, fontsize=20, title_fontsize=20
)

# Adjust layout and show the plot
plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space for legends
savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\_02 Presentaties posters abstracts\241209_generallab_figs'
# plt.savefig(savepath + r'\all_energycarriers_fourparts_Ecoli.png', dpi=300, bbox_inches='tight')
plt.show()

#%%

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches



# Define bar colors for ATP, NADPH, NADH
bar_colors = {'atp frac': 'tab:blue', 'nadph frac': 'tab:green', 'nadh frac': 'tab:orange'}
names = {'atp frac': 'ATP', 'nadph frac': 'NADPH', 'nadh frac': 'NADH'}

# Define bar colors for ATP, NADPH, NADH
bar_colors2 = {'M_atp_c': 'tab:blue', 'M_nadph_c': 'tab:green', 'M_nadh_c': 'tab:orange'}
names2 = {'M_atp_c': 'ATP', 'M_nadph_c': 'NADPH', 'M_nadh_c': 'NADH'}

# Subplot titles
titles = ['Energy Conservation', 'Precursor Production', 'Anabolism from Carbon Precursors', 'Fraction EC produced in anaCP']

# DataFrames to plot
dataframes = [
    df_mceqs_cat_reindexed,
    df_mceqs_PP_reindexed,
    df_mceqs_aCP_reindexed,
    df_atpprod_reindexed,
]

# Create subplots
fig, axes = plt.subplots(1, 4, figsize=(30, 11.5))  # Adjust figure size for 4 subplots

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

axes[3].tick_params(
    axis='y',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    left=False,      # ticks along the bottom edge are off
    right=False,         # ticks along the top edge are off
    labelleft=False) # labels along the bottom edge are off

twin_axes = []
# Loop through each DataFrame and plot
for ax, df, title in zip(axes[:-1], dataframes[:-1], titles[:-1]):
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
    
    if title == 'Precursor Production':  # Update as needed to match the second subplot
        x_min = min(df[atp_col].min(), df[['M_nadph_c', 'M_nadh_c']].min().min())  # Get shared min value
        x_max = max(df[atp_col].max(), df[['M_nadph_c', 'M_nadh_c']].max().max())  # Get shared max value
        
        
        
        ax.set_xlim(x_min -5, x_max + 5)  # Set the primary x-axis limits
        twin_ax.set_xlim((x_min-5)/4, (x_max+5)/4)  # Set the secondary x-axis limits to match
    twin_axes.append(twin_ax)
    # Add category shading
    current_y = 0
    for category, substrates in sort_dict.items():
        count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
        ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
        current_y += count


    ax.set_title(title, fontsize=28)  # Increase subplot title font size
    ax.axvline(x=0, color='black', linewidth=1)  # Add a vertical line at x=0 for reference
    
    lines = np.arange(-0.5, len(filtered_sorted_index)-0.5, 1)
    for lin in lines:
        ax.axhline(y=lin, color='black', linewidth=1)


    # Remove whitespace above and below the bars
    ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)

# Set y-axis labels
axes[0].set_yticks(indices)
axes[0].set_yticklabels(filtered_sorted_index, fontsize=28)  # Increase y-axis tick font size

axes[0].set_xlabel('$(Y_{X/ATP}^{EC})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])
axes[1].set_xlabel('$(Y_{X/ATP}^{PP})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])
axes[2].set_xlabel('$(Y_{X/ATP}^{ana})^{-1}$ [mmol/gDW]', fontsize=28, color=bar_colors['atp frac'])

twin_axes[0].set_xlabel('$(Y_{X/NAD(P)H}^{EC})^{-1}$ [mmol/gDW]', fontsize=28)
twin_axes[1].set_xlabel('$(Y_{X/NAD(P)H}^{PP})^{-1}$ [mmol/gDW]', fontsize=28)
twin_axes[2].set_xlabel('$(Y_{X/NAD(P)H}^{ana})^{-1}$ [mmol/gDW]', fontsize=28)


# Plot for the last subplot (X produced/X consumed)
ax, df, title = axes[-1], dataframes[-1], titles[-1]
indices = np.arange(len(filtered_sorted_index))
bar_width = 0.25
offsets = [-bar_width, 0, bar_width]

# Plot ATP on bottom x-axis
ax.barh(indices + offsets[0], df['atp frac'], height=bar_width, label='ATP', color=bar_colors['atp frac'])

# Create twin x-axis for NADPH and NADH
twin_ax = ax.twiny()
twin_ax.barh(indices + offsets[1], df['nadph frac'], height=bar_width, label='NADPH', color=bar_colors['nadph frac'])
twin_ax.barh(indices + offsets[2], df['nadh frac'], height=bar_width, label='NADH', color=bar_colors['nadh frac'])

ax.tick_params(axis='x', labelsize=28, labelcolor=bar_colors2[atp_col])
twin_ax.tick_params(axis='x', labelsize=28, labelcolor='black')

# Add category shading
current_y = 0
for category, substrates in sort_dict.items():
    count = sum(1 for substrate in substrates if substrate in filtered_sorted_index)
    ax.axhspan(current_y - 0.5, current_y + count - 0.5, color=category_colors[category], alpha=0.2)
    current_y += count


ax.set_title(title, fontsize=28)
ax.axvline(x=0, color='black', linewidth=1)

lines = np.arange(-0.5, len(filtered_sorted_index)-0.5, 1)
for lin in lines:
    ax.axhline(y=lin, color='black', linewidth=1)


# Add xlabel for ATP (bottom x-axis) and NADPH/NADH (top x-axis)
ax.set_xlabel('ATP Fraction', fontsize=28, color=bar_colors['atp frac'])
twin_ax.set_xlabel('NADPH/NADH Fraction', fontsize=28, color='black')

# Remove whitespace above and below the bars
ax.set_ylim(-0.5, len(filtered_sorted_index) - 0.5)

# Add a common legend for bar colors (energy carriers)
bar_legend = fig.legend(
    handles=[mpatches.Patch(color=color, label=names[col]) for col, color in bar_colors.items()],
    loc='lower center', bbox_to_anchor=(0.4, -0.2),  # Center below the plot
    frameon=False, fontsize=28, title_fontsize=28, title='Energy Carriers'
)

# Add a common legend for category shading
shading_legend_patches = [
    mpatches.Patch(color=color, label=category, alpha=0.2)
    for category, color in category_colors.items()
]
shading_legend = fig.legend(
    handles=shading_legend_patches,
    loc='lower center', bbox_to_anchor=(0.6, -0.2),  # Slightly lower than the first legend
    frameon=False, fontsize=28, title_fontsize=28, title="Carbon source pathway"
)


# Adjust layout to make space for the legends
plt.tight_layout(rect=[0, 0.1, 1, 1])  # Leave extra space at the bottom
plt.savefig(savepath+r'\all_energycarriers_threeparts_fractions.png', dpi=300, bbox_inches='tight')
plt.show()

