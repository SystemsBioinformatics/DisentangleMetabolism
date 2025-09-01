# -*- coding: utf-8 -*-
"""
Created on Fri Dec 13 09:10:19 2024

@author: mre283
"""

import matplotlib.pyplot as plt
import pandas as pd

result_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\Precursor_simulations_Ecoli'


df_mceqs = pd.read_excel(result_path+'\\all_mceqs.xlsx', index_col=0)

df_mceqs_plot = df_mceqs.drop(['M_biomass_e', 'M_btn_e', 'M_ca2_e', 'M_cl_e', 'M_cobalt2_e', 'M_cu2_e', 'M_fe2_e', 'M_h2o_e',
                               'M_h_e', 'M_k_e', 'M_kdo2lipid4_e', 'M_mg2_e', 'M_mn2_e', 'M_mobd_e', 'M_nh4_e', 'M_ni2_e',
                               'M_so4_e', 'M_zn2_e', 'M_adp_c', 'M_fad_c', 'M_nad_c', 'M_nadp_c', 'M_q8_c', 'M_gdp_e', 'M_udp_e',
                               'M_fe3_e', 'M_gtp_e', 'M_utp_e', 'M_3mob_e', 'M_fgam_e', 'M_ala__D_e'], axis=1)


#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

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
        'color': 'blue'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['acetate', 'succinate', 'alphaketoglutarate', 'citrate', 'fumarate', 'malate'],
        'color': 'green'
    },
    'lower glycolysis': {
        'metabolites': ['pyruvate', 'glycerol', 'lactate'],
        'color': 'red'
    },
    'amino acids': {
        'metabolites': ['alanine', 'arginine', 'aspartate', 'glutamate', 'glutamine', 'tryptophan'],
        'color': 'purple'
    }
}


# BiGG ID to full names/formulas mapping
bigg_to_name = {
    'ac': 'acetate',
    'accoa': 'acetyl-CoA',
    'akg': r'alphaketoglutarate',
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
        'color': 'blue'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['ac', 'accoa', 'akg', 'fum', 'oaa', 'succ', 'succoa', 'icit'],
        'color': 'green'
    },
    'lower glycolysis': {
        'metabolites': ['3pg', 'dhap', 'g3p', 'pep', 'pyr', 'dha', 'glyc'],
        'color': 'red'
    },
    'amino acids': {
        'metabolites': ['ala__L', 'glu__L', 'gly', 'val__L', 'arg__L', 'asp__L', 'ser__L', 'ala__D', 'gln__L', 'trp__L'],
        'color': 'purple'
    },
    
    'conserved moieties': {
        'metabolites': ['thf', 'coa', 'pi'],
        'color': 'orange'
    },
    'ppp': {
        'metabolites': ['e4p', 'r5p', 'ru5p__D'],
        'color': 'brown'
    }, 
    'external': {
        'metabolites': ['co2', 'o2', 'for'],
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



#%%
# Define grid-drawing function
def add_grid(ax, df):
    n_rows, n_cols = df.shape
    ax.hlines(np.arange(-0.5, n_rows), xmin=-0.5, xmax=n_cols - 0.5, color='black', linewidth=0.5)
    ax.vlines(np.arange(-0.5, n_cols), ymin=-0.5, ymax=n_rows - 0.5, color='black', linewidth=0.5)

# Function to remove 'aerobic_' prefix from y-tick labels
def clean_label(label):
    return label[8:] if label.startswith('aerobic_') else label


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

#%%



# Transpose the DataFrame for the main heatmap
df_main = df_main.transpose()

# Reassign categories and sort the new rows (previously columns)
df_main['Category'] = df_main.index.map(lambda x: categorize_y_label(x[2:-2]))
df_main = df_main.sort_values(by='Category')
sorted_index_main = df_main.index
sorted_categories_main = df_main['Category'].values
df_main = df_main.drop(columns='Category')

# Update y-axis ticks and labels (was x-axis before transposition)
yticks_main = [bigg_to_name.get(idx[2:-2], idx[2:-2]) for idx in sorted_index_main]
xticks_main = [clean_label(col) for col in df_main.columns]

# Update color mapping for y-axis labels
# def get_y_label_color(label):
#     category = categorize_x_label(label[2:-2])
#     return groups_y.get(category, {}).get('color', 'black')

# Plot the transposed heatmap
fig, ax1 = plt.subplots(figsize=(20, 12))

ax1.set_xticklabels([])  # Start with empty labels to add colored text manually
ax1.set_yticklabels([])  # Start with empty labels to add colored text manually
cmap_main = plt.get_cmap('seismic')
norm_main = mcolors.TwoSlopeNorm(vmin=df_main.values.min(), vcenter=0, vmax=df_main.values.max())
cax1 = ax1.matshow(df_main, cmap=cmap_main, norm=norm_main)
ax1.set_title("Anabolic carbon precursor requirement", fontsize=16)
colorbar = plt.colorbar(cax1, ax=ax1)
colorbar.set_label('yield [mmol/gDW]', fontsize=14, labelpad=10, loc='center')  # Add label above



# Add x-axis labels (formerly y-axis labels)
for i, (tick, category) in enumerate(zip(xticks_main, sorted_categories_main)):
    ax1.text(i, -1.2, tick, color=get_x_label_color(tick), fontsize=12, ha='center', va='bottom', rotation=90)

# Add y-axis labels (formerly x-axis labels)
for i, label in enumerate(yticks_main):
    color = get_y_label_color(sorted_index_main[i])
    ax1.text(-1.2, i, label, color=color, fontsize=12, ha='right', va='center')

# Add grid to the transposed heatmap
add_grid(ax1, df_main)

# Adjust legend for y-axis colors
legend_handles = [plt.Line2D([0], [0], color=data['color'], marker='o', markersize=10, linestyle='None', label=group)
                  for group, data in groups_y.items()]
ax1.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=12)

# Move labels
ax1.set_xlabel('Carbon and Energy Source', fontsize=14, labelpad=140)
ax1.xaxis.set_label_position('top')
ax1.set_ylabel('Metabolite', fontsize=14, labelpad=125)
ax1.yaxis.set_label_position('left')

plt.tight_layout()
plt.savefig('heatmap_Ecoli_carbonsubs_transposed.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()


# #%%

# # Sort x-axis columns by their categories
# # Create a mapping of x-axis labels to their categories
# x_label_categories = {label: categorize_x_label(label) for label in df_main.columns}

# # Define the order of categories
# category_order = list(groups.keys()) + ['Other']  # Ensure 'Other' is sorted last

# # Sort the columns based on the category order
# sorted_columns = sorted(df_main.columns, key=lambda x: category_order.index(x_label_categories[x]))

# # Reorder the DataFrame's columns
# df_main = df_main[sorted_columns]

# # Update x-axis ticks and labels based on the sorted order
# xticks_main = [clean_label(col) for col in sorted_columns]
# sorted_categories_main = [x_label_categories[col] for col in sorted_columns]

# # Update the plot with sorted columns and labels
# fig, ax1 = plt.subplots(figsize=(20, 12))

# ax1.set_xticklabels([])  # Start with empty labels to add colored text manually
# ax1.set_yticklabels([])  # Start with empty labels to add colored text manually
# cmap_main = plt.get_cmap('seismic')
# norm_main = mcolors.TwoSlopeNorm(vmin=df_main.values.min(), vcenter=0, vmax=df_main.values.max())
# cax1 = ax1.matshow(df_main, cmap=cmap_main, norm=norm_main)
# ax1.set_title("Anabolic carbon precursor requirement", fontsize=16)
# colorbar = plt.colorbar(cax1, ax=ax1)
# colorbar.set_label('yield [mmol/gDW]', fontsize=14, labelpad=10, loc='center')  # Add label above

# # Add x-axis labels (formerly y-axis labels)
# for i, (tick, category) in enumerate(zip(xticks_main, sorted_categories_main)):
#     ax1.text(i, -1.2, tick, color=get_x_label_color(tick), fontsize=12, ha='center', va='bottom', rotation=90)

# # Add y-axis labels (formerly x-axis labels)
# for i, label in enumerate(yticks_main):
#     color = get_y_label_color(sorted_index_main[i])
#     ax1.text(-1.2, i, label, color=color, fontsize=12, ha='right', va='center')

# # Add grid to the transposed heatmap
# add_grid(ax1, df_main)

# # Adjust legend for y-axis colors
# legend_handles = [plt.Line2D([0], [0], color=data['color'], marker='o', markersize=10, linestyle='None', label=group)
#                   for group, data in groups.items()]
# ax1.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=12)

# # Move labels
# ax1.set_xlabel('Carbon and Energy Source', fontsize=14, labelpad=150)
# ax1.xaxis.set_label_position('top')
# ax1.set_ylabel('Metabolite', fontsize=14, labelpad=125)
# ax1.yaxis.set_label_position('left')

# plt.tight_layout()
# plt.savefig('heatmap_Ecoli_carbonsubs_sorted_transposed.png', dpi=300, bbox_inches='tight')
# plt.show()

#%%

df_comparison = pd.read_excel(r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\Ecolicore_Neidhardt_table_comparison.xlsx', sheet_name='pythonreadable', index_col=0)

df_comparison_ec = df_comparison.iloc[-5:-2, :]

df_comparison = df_comparison.drop(['NADH', 'NADPH', 'ATP'], axis=0)



#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

groups_y = {
    'upper glycolysis': {
        'metabolites': ['fructose-6P', 'glucose-6P'],
        'color': 'blue'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['acetate', 'acetyl-CoA', 'alphaketoglutarate', 'fumarate', 'oxaloacetate', 'succinate', 'succinyl-CoA'],
        'color': 'green'
    },
    
    'lower glycolysis': {
        'metabolites': ['3P-glycerate', 'dihydroxyacetone-P', 'glyceraldehyde-3P', 'PEP', 'pyruvate'],
        'color': 'red'
    },
    'amino acids': {
        'metabolites': ['glutamine', 'glutamate'],
        'color': 'purple'
    }, 
    
    'ppp': {
        'metabolites': ['erythrose-4P', 'ribose-5P', 'ribulose-5P'],
        'color': 'brown'
    }, 
    'external': {
        'metabolites': ['CO2', 'O2'],
        'color': 'pink'
    }
}

def categorize_x_label(label):
    for group, data in groups_y.items():
        if any(metabolite in label for metabolite in data['metabolites']):
            return group
    return 'Other'


def get_x_label_color(label):
    category = categorize_x_label(label)
    # print(category)
    # print(label)
    return groups_y.get(category, {}).get('color', 'black')  # Default to black


# Ensure groups and color definitions remain the same
# Assuming the groups and bigg_to_name dictionaries are already defined as in your code

def categorize_index(index):
    for group, data in groups_y.items():
        if any(metabolite in index for metabolite in data['metabolites']):
            return group
    return 'Other'

# Prepare the df_comparison DataFrame (assuming it is already loaded)
# Ensure the index represents the metabolites and the columns match the specified ones

# Categorize the indices of df_comparison
df_comparison['Category'] = df_comparison.index.map(categorize_index)
df_comparison = df_comparison.sort_values(by='Category')
sorted_index = df_comparison.index
sorted_categories = df_comparison['Category'].values
df_comparison = df_comparison.drop(columns='Category')

# Transpose the DataFrame for better visualization
# df_comparison = df_comparison.transpose()
df_comparison['Category'] = df_comparison.index.map(lambda x: categorize_y_label(x))
df_comparison = df_comparison.sort_values(by='Category')

sorted_index_main = df_comparison.index
sorted_categories_main = df_comparison['Category'].values
df_comparison = df_comparison.drop(columns='Category')

# Prepare labels
yticks_main = [bigg_to_name.get(idx, idx) for idx in sorted_index_main]
xticks_main = [clean_label(col) for col in df_comparison.columns]

# Plot setup
fig, ax1 = plt.subplots(figsize=(20, 12))
ax1.set_xticklabels([])
ax1.set_yticklabels([])

cmap_main = plt.get_cmap('seismic')
norm_main = mcolors.TwoSlopeNorm(vmin=df_comparison.values.min(), vcenter=0, vmax=df_comparison.values.max())
cax1 = ax1.matshow(df_comparison, cmap=cmap_main, norm=norm_main)

# ax1.set_title("Comparison Heatmap", fontsize=16)
colorbar = plt.colorbar(cax1, ax=ax1)
colorbar.set_label('yield [mmol/gDW]', fontsize=14, labelpad=10)

# Add labels
for i, (tick, category) in enumerate(zip(xticks_main, sorted_categories_main)):
    ax1.text(i, -1.2, tick, color=get_x_label_color(tick), fontsize=12, ha='center', va='bottom', rotation=90)

for i, label in enumerate(yticks_main):
    color = get_x_label_color(sorted_index_main[i])
    ax1.text(-1.2, i, label, color=color, fontsize=12, ha='right', va='center')

# Add grid
add_grid(ax1, df_comparison)

# Adjust legend
legend_handles = [plt.Line2D([0], [0], color=data['color'], marker='o', markersize=10, linestyle='None', label=group)
                  for group, data in groups_y.items()]
ax1.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=12)

# Move labels
# ax1.set_xlabel('Source', fontsize=14, labelpad=140)
ax1.xaxis.set_label_position('top')
ax1.set_ylabel('Metabolite', fontsize=14, labelpad=140)
ax1.yaxis.set_label_position('left')

plt.tight_layout()
plt.savefig('comparison_heatmap.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()

#%%

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

groups_y = {
    'upper glycolysis': {
        'metabolites': ['fructose-6P', 'glucose-6P'],
        'color': 'blue'
    },
    'TCA-cycle intermediates': {
        'metabolites': ['acetate', 'acetyl-CoA', 'alphaketoglutarate', 'fumarate', 'oxaloacetate', 'succinate', 'succinyl-CoA'],
        'color': 'green'
    },
    'lower glycolysis': {
        'metabolites': ['3P-glycerate', 'dihydroxyacetone-P', 'glyceraldehyde-3P', 'PEP', 'pyruvate'],
        'color': 'red'
    },
    'amino acids': {
        'metabolites': ['glutamine', 'glutamate'],
        'color': 'purple'
    },
    'ppp': {
        'metabolites': ['erythrose-4P', 'ribose-5P', 'ribulose-5P'],
        'color': 'brown'
    },
    'external': {
        'metabolites': ['CO2', 'O2'],
        'color': 'pink'
    }
}

def categorize_x_label(label):
    for group, data in groups_y.items():
        if any(metabolite in label for metabolite in data['metabolites']):
            return group
    return 'Other'

def get_x_label_color(label):
    category = categorize_x_label(label)
    return groups_y.get(category, {}).get('color', 'black')

def categorize_index(index):
    for group, data in groups_y.items():
        if any(metabolite in index for metabolite in data['metabolites']):
            return group
    return 'Other'

df_comparison_ec['Category'] = df_comparison_ec.index.map(categorize_index)
df_comparison_ec = df_comparison_ec.sort_values(by='Category')
sorted_index = df_comparison_ec.index
sorted_categories = df_comparison_ec['Category'].values
df_comparison_ec = df_comparison_ec.drop(columns='Category')
df_comparison_ec['Category'] = df_comparison_ec.index.map(lambda x: categorize_index(x))
df_comparison_ec = df_comparison_ec.sort_values(by='Category')
sorted_index_main = df_comparison_ec.index
sorted_categories_main = df_comparison_ec['Category'].values
df_comparison_ec = df_comparison_ec.drop(columns='Category')
yticks_main = [idx for idx in sorted_index_main]
xticks_main = [col for col in df_comparison_ec.columns]

fig, ax1 = plt.subplots(figsize=(10, 3))
ax1.set_xticklabels([])
ax1.set_yticklabels([])
cmap_main = plt.get_cmap('seismic')
norm_main = mcolors.TwoSlopeNorm(vmin=-70, vcenter=0, vmax=10)
cax1 = ax1.matshow(df_comparison_ec, cmap=cmap_main, norm=norm_main)
colorbar = plt.colorbar(cax1, ax=ax1, ticks=[-70, -35, 0, 5, 10], aspect=7.5)
colorbar.set_label('yield [mmol/gDW]', fontsize=12, labelpad=10)
add_grid(ax1, df_comparison_ec)

for i, (tick, category) in enumerate(zip(xticks_main, sorted_categories_main)):
    ax1.text(i, -1.2, tick, color=get_x_label_color(tick), fontsize=12, ha='center', va='bottom', rotation=90)

for i, label in enumerate(yticks_main):
    color = get_x_label_color(sorted_index_main[i])
    ax1.text(-1.2, i, label, color=color, fontsize=12, ha='right', va='center')

plt.tight_layout()
plt.savefig('comparison_heatmap_ec.png', dpi=300, bbox_inches='tight', transparent=True)
plt.show()
