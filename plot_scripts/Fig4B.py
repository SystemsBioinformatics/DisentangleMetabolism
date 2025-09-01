# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 14:05:37 2025

@author: mre283
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

df_mceqs = pd.read_excel('Results/all_mceqs_Ecoli_anabolism.xlsx', header=0, index_col=0)

df_mceqs_plot = df_mceqs.drop(['M_biomass_e', 'M_btn_e', 'M_ca2_e', 'M_cl_e', 'M_cobalt2_e', 'M_cu2_e', 'M_fe2_e', 'M_h2o_e',
                               'M_h_e', 'M_k_e', 'M_kdo2lipid4_e', 'M_mg2_e', 'M_mn2_e', 'M_mobd_e', 'M_nh4_e', 'M_ni2_e',
                               'M_so4_e', 'M_zn2_e', 'M_adp_c', 'M_fad_c', 'M_nad_c', 'M_nadp_c', 'M_q8_c', 'M_gdp_e', 'M_udp_e',
                               'M_fe3_e', 'M_gtp_e', 'M_utp_e', 'M_3mob_e', 'M_fgam_e', 'M_ala__D_e'], axis=1)

# Dictionary with groups and corresponding metabolites
groups = {
    'upper glycolysis': ['fructose', 'glucose', 'galactose', 'gluconate', 'xylose', 'ribose', 'sorbitol', 'anaerobic_glucose'],
    'TCA-cycle intermediates': ['acetate', 'succinate', 'alphaketoglutarate', 'citrate', 'fumarate', 'malate'],
    'lower glycolysis': ['pyruvate', 'glycerol', 'lactate'],
    'amino acids': ['alanine', 'arginine', 'aspartate', 'glutamate', 'glutamine', 'tryptophan']
}

# Function to categorize each index based on the groups dictionary
def categorize_index(index):
    for group, metabolites in groups.items():
        if any(metabolite in index for metabolite in metabolites):
            return group
    return 'Other'  # Default category if no match is found

# Assuming df_mceqs_plot is your pre-defined DataFrame
# Columns to move to the right-hand subplot
columns_to_move = ['M_atp_c', 'M_nadh_c', 'M_nadph_c', 'M_fadh2_c', 'M_q8h2_c']

# Check if all columns exist in the DataFrame
existing_columns = [col for col in columns_to_move if col in df_mceqs_plot.columns]

# Create two DataFrames: one for the main heatmap and one for the specified columns
df_main = df_mceqs_plot.drop(columns=existing_columns)  # Remove specified columns

# Categorize the indices of df_main
df_main['Category'] = df_main.index.map(categorize_index)

# Sort the DataFrame by Category
df_main = df_main.sort_values(by='Category')

# Extract the sorted index and categories for labeling
sorted_index = df_main.index
sorted_index_main = df_main.index
sorted_categories = df_main['Category'].values

# Remove the 'Category' column before plotting
df_main = df_main.drop(columns='Category')



groups = {
    'upper glycolysis': {
        'metabolites': ['fructose', 'glucose', 'galactose', 'gluconate', 'xylose', 'ribose', 'sorbitol', 'anaerobic_glucose'],
        'color': 'tab:orange'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['acetate', 'succinate', r'$\alpha$-ketoglutarate', 'citrate', 'fumarate', 'malate'],
        'color': 'tab:green'
    },
    'lower glycolysis': {
        'metabolites': ['pyruvate', 'glycerol', 'lactate'],
        'color': 'tab:purple'
    },
    'amino acids': {
        'metabolites': ['alanine', 'arginine', 'aspartate', 'glutamate', 'glutamine', 'tryptophan'],
        'color': 'tab:red'
    }
}


# BiGG ID to full names/formulas mapping
bigg_to_name = {
    'ac': 'acetate',
    'accoa': 'acetyl-CoA',
    'akg': r'$\alpha$-ketoglutarate',
    'co2': r'CO$_2$',  # Latex formatted
    'coa': 'CoA',
    'fum': 'fumarate',
    'o2': r'O$_2$',
    'oaa': 'oxaloacetate',
    'pi': 'phosphate',
    'succ': 'succinate',
    'succoa': 'succinyl-CoA',
    '3pg': '3P-glycerate',
    'dhap': 'dihydroxyacetone-P',
    'e4p': 'erythrose 4P',
    'f6p': 'fructose 6P',
    'g3p': 'glyceraldehyde 3P',
    'gtp': 'Guanosine triphosphate',
    'pep': 'PEP',
    'pyr': 'pyruvate',
    'r5p': 'ribose 5P',
    'ru5p__D': 'ribulose 5P',
    'dha': 'dihydroxyacetone',
    'glyc': 'glycerol',
    'icit': 'isocitrate',
    'utp': 'Uridine triphosphate',
    '3mob': '3-Methyl-2-oxobutanoate',
    'ala__L': 'alanine',
    'glu__L': 'glutamate',
    'gly': 'glycine',
    'val__L': 'valine',
    'for': 'formate',
    'arg__L': 'arginine',
    'asp__L': 'aspartate',
    '10fthf': 'formyl-THF',
    'fgam': 'Formylglycinamide ribonucleotide',
    'mlthf': 'methylene-THF',
    'ser__L': 'serine',
    'thf': 'THF',
    'ala__D': 'D-alanine',
    'gln__L': 'glutamine',
    'trp__L': 'tryptophan'
}

groups_y = {
    'upper glycolysis': {
        'metabolites': ['f6p'],
        'color': 'tab:orange'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['ac', 'accoa', 'akg', 'fum', 'oaa', 'succ', 'succoa', 'icit'],
        'color': 'tab:green'
    },
    'lower glycolysis': {
        'metabolites': ['3pg', 'dhap', 'g3p', 'pep', 'pyr', 'dha', 'glyc', 'for'],
        'color': 'tab:purple'
    },
    'amino acids': {
        'metabolites': ['ala__L', 'glu__L', 'gly', 'val__L', 'arg__L', 'asp__L', 'ser__L', 'ala__D', 'gln__L', 'trp__L'],
        'color': 'tab:red'
    },
    
    'conserved moieties': {
        'metabolites': ['thf', 'coa', 'pi'],
        'color': 'blue'
    },
    'pentose phosphate pathway': {
        'metabolites': ['e4p', 'r5p', 'ru5p__D'],
        'color': 'brown'
    }, 
    'inorganic': {
        'metabolites': ['co2', 'o2'],
        'color': 'pink'
    },
    'other': {
        'metabolites': ['10fthf', 'fgam', 'mlthf'],
        'color': 'grey'
    }
}


# Function to categorize each index based on the groups dictionary
def categorize_index(index):
    for group, data in groups.items():
        if any(metabolite in index for metabolite in data['metabolites']):
            return group
    return 'Other'



# Define grid-drawing function
def add_grid(ax, df):
    n_rows, n_cols = df.shape
    ax.hlines(np.arange(-0.5, n_rows), xmin=-0.5, xmax=n_cols - 0.5, color='black', linewidth=0.5)
    ax.vlines(np.arange(-0.5, n_cols), ymin=-0.5, ymax=n_rows - 0.5, color='black', linewidth=0.5)

# Function to remove 'aerobic_' prefix from y-tick labels
def clean_label(label):
    nw_lb = label
    if label.startswith('aerobic_'):
        nw_lb = label[8:]
    
    if 'alphaketoglutarate' in label:
        nw_lb = r'$\alpha$-ketoglutarate'
    
    return nw_lb


# Customize y-axis labels with cleaned metabolite indices and their colors
yticks_main = [clean_label(index) for index in sorted_index_main]  # Clean the labels
xticks_main = [idd[2:-2] for idd in df_main.columns]


# Define a function to map x-axis labels to their categories
def categorize_x_label(label):
    for group, data in groups.items():
        if any(metabolite in label for metabolite in data['metabolites']):
            return group
    return 'Other'

# Define a function to get the color for the x-axis labels
def get_x_label_color(label):
    category = categorize_x_label(label)
    return groups.get(category, {}).get('color', 'black')  # Default to black


# Define a function to map x-axis labels to their categories
def categorize_y_label(label):
    for group, data in groups_y.items():
        if any(metabolite in label for metabolite in data['metabolites']):
            return group
    return 'Other'

# Define a function to get the color for the x-axis labels
def get_y_label_color(label):
    category = categorize_y_label(label[2:-2])
    return groups_y.get(category, {}).get('color', 'black')  # Default to black

# Create a mapping of x-axis labels to their categories
x_label_categories = {label: categorize_x_label(label[2:-2]) for label in df_main.columns}

# Sort the columns based on the category order in groups_y
category_order = list(groups_y.keys()) + ['Other']  # Ensure 'Other' is sorted last
sorted_columns = sorted(df_main.columns, key=lambda x: category_order.index(x_label_categories[x]))



# --- Preprocessing ---
# Sort columns (metabolites) by category
df_main.columns.name = None
sorted_columns_main = sorted(df_main.columns, key=lambda x: categorize_y_label(x[2:-2]))[5:]
df_main = df_main[sorted_columns_main]
sorted_categories_main = [categorize_y_label(col[2:-2]) for col in sorted_columns_main]

# Get tick labels
xticks_main = [bigg_to_name.get(col[2:-2], col[2:-2]) for col in sorted_columns_main]
# Desired substrate order (reversed)
desired_order = ['xylose', 'sorbitol', 'glucose', 'anaerobic_glucose', 'galactose',
                 'ribose', 'fructose', 'gluconate', 'succinate', 'malate', 'acetate',
                 'fumarate', 'citrate', 'alphaketoglutarate', 'pyruvate', 'lactate', 'glycerol',
                 'tryptophan', 'glutamine', 'glutamate', 'aspartate', 'arginine',
                 'alanine'][::-1]  # Reverse order

# Function to extract substrate name from index (assuming index like "aerobic_glucose")
def extract_substrate(index):
    if index[:8] == 'aerobic_':
        return index.replace('aerobic_', '')  
    else: 
        return index

# Map indices to substrate names
df_main['substrate'] = df_main.index.map(extract_substrate)

# Drop rows not in desired order
df_main = df_main[df_main['substrate'].isin(desired_order)]

# Sort by reversed desired order
df_main['substrate'] = pd.Categorical(df_main['substrate'], categories=desired_order, ordered=True)
df_main = df_main.sort_values('substrate')

# Drop helper column
df_main = df_main.drop(columns='substrate')

# Update y-tick labels
yticks_main = [clean_label(idx) for idx in df_main.index]

# --- Plotting ---
fig, ax1 = plt.subplots(figsize=(20, 12))

# Heatmap
cmap_main = plt.get_cmap('seismic')
norm_main = mcolors.TwoSlopeNorm(vmin=df_main.values.min(), vcenter=0, vmax=df_main.values.max())
cax1 = ax1.matshow(df_main.values, cmap=cmap_main, norm=norm_main)

# Colorbar
colorbar = plt.colorbar(cax1, ax=ax1)
colorbar.ax.tick_params(labelsize=18) 
colorbar.set_label('yield [mmol/gDW]', fontsize=22, labelpad=10, loc='center')

# Clear default ticks (we'll add custom labels)
ax1.set_xticks([])
ax1.set_yticks([])

# X-axis labels (metabolites)
for i, label in enumerate(xticks_main):
    color = get_y_label_color(sorted_columns_main[i])
    ax1.text(i, -1.2, label, color=color, fontsize=18, ha='center', va='bottom', rotation=90)

# for i, col in enumerate(sorted_columns_main):
#     bigg_id = col[2:-2]
#     label = bigg_to_name.get(bigg_id, bigg_id)
#     color = get_y_label_color(bigg_id)
#     ax1.text(i, -1.2, label, color=color, fontsize=18, ha='center', va='bottom', rotation=90)

# Y-axis labels (carbon sources)
for i, label in enumerate(yticks_main):
    color = get_x_label_color(label)
    ax1.text(-1.2, i, label, color=color, fontsize=18, ha='right', va='center')

# Labels
ax1.set_xlabel('Precursor metabolite', fontsize=18, labelpad=160)
ax1.xaxis.set_label_position('top')
ax1.set_ylabel('Carbon and energy Source', fontsize=18, labelpad=170)
ax1.yaxis.set_label_position('left')

# Grid
add_grid(ax1, df_main)

# Legend
legend_handles = [
    plt.Line2D([0], [0], color=data['color'], marker='o', markersize=14,
               linestyle='None', label=group)
    for group, data in groups_y.items()
]
ax1.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.25),
           ncol=2, frameon=False, fontsize=18)


# Add a cross for each cell where the value is exactly zero
zero_locs = np.argwhere(df_main.values == 0)
for y, x in zero_locs:
    ax1.plot(x, y, marker='x', color='black', markersize=10, markeredgewidth=2)

# Layout and display
plt.tight_layout()

plt.savefig(r'Plots\heatmap_Ecoli_carbonsubs_horizontal.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()