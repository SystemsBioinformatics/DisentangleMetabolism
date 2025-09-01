# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:49:18 2024

@author: mre283
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

#%%

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
        df_dict[filename.path[28:-5]] = pd.read_excel(filename.path, sheet_name='atp_overall', index_col=0)

directory = 'Results\Ecoli_nitrogensources'

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.path[30:-5] == 'glucose_nitrate' or filename.path[30:-5] == 'glucose_nitrate_nadh':
        continue
    if filename.is_file():
        df_dict[filename.path[30:-5]] = pd.read_excel(filename.path, sheet_name='atp_overall', index_col=0)

#%%

for k, df in df_dict.items():
    print(k)
    print(df)

#%%
plt.figure()
for k, df in df_dict.items():
    plt.bar(k, df.loc['ana', 'prod']/df.loc['ana']['cons'],zorder=3)

plt.xticks(rotation=90)
plt.grid(zorder=0)
plt.ylabel('fraction of ATP used produced in anabolism [-]')

#%% ordered by value
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Set the font sizes
plt.rcParams.update({
    'font.size': 28,          # Default text size
    'axes.titlesize': 28,     # Title font size
    'axes.labelsize': 28,     # Axis labels font size
    'xtick.labelsize': 28,    # X-tick labels font size
    'ytick.labelsize': 32,    # Y-tick labels font size
    'legend.fontsize': 28,    # Legend font size
    'figure.titlesize': 24    # Figure title font size
})

# Assuming df_dict and NADPH_exch are defined as in the user's provided code
NADPH_exch = {'glutamine': False, 'glutamate': False, 'sorbitol': True, 'arginine': False, 'akg': False, 'fructose': True,
              'glucose': True, 'galactose': True, 'xylose': True, 'ribose': True, 'citrate': True, 'fumarate': True,
              'succinate': True, 'alanine': False, 'malate': True, 'glycerol': False, 'lactate': False, 'gluconate': True,
              'glucose_anaerobic': True, 'aspartate': False, 'acetate': False, 'pyruvate': False, 'tryptophan': False}

# Create df_fraction
df_fraction = pd.DataFrame(columns=['ATP fraction'], index=df_dict.keys())
df_fraction = df_fraction.drop('aer_glc_nadh')

for k, df in df_dict.items():
    df_fraction['ATP fraction'].loc[k] = df.loc['ana', 'prod'] / df.loc['ana', 'cons']

# Sort df_fraction
df_fraction = df_fraction.sort_values(by=['ATP fraction'], ascending=False)

# Customize labels for amino acids and akg
df_fraction.index = [f"{label}*" if NADPH_exch.get(label) == False and label in ['glutamine', 'glutamate', 'arginine', 'alanine', 'aspartate', 'tryptophan'] else label for label in df_fraction.index]
df_fraction = df_fraction.rename(index={'akg': r'$\alpha$-ketoglutarate'})

# Plotting with the correct figsize
fig, ax = plt.subplots(figsize=(14, 10))  # Adjust the figure size here

# Horizontal bar plot
ax = df_fraction.plot.barh(y='ATP fraction', zorder=3, legend=False, ax=ax)  # Pass ax to the plot method

# Customize colors based on NADPH_exch dictionary
blue_rgba = (0.12, 0.47, 0.71, 1.0)
orange_rgba = (1.0, 0.50, 0.05, 1.0)
colors = [blue_rgba if NADPH_exch.get(idx.rstrip('*')) else orange_rgba for idx in df_fraction.index]

# Applying the colors to the bars
bars = ax.patches
for bar, color in zip(bars, colors):
    bar.set_color(color)

# Annotating the 'glucose_anaerobic' bar
for bar, idx in zip(bars, df_fraction.index):
    if idx == 'glucose_anaerobic':
        ax.text(bar.get_width() + 0.01, bar.get_y() + bar.get_height() / 2, 'anaerobic', ha='left', va='center', fontsize=28)

labels = df_fraction.index.tolist()
labels[labels.index('glucose_anaerobic')] = 'glucose'  # Change 'glucose_anaerobic' to 'glucose'
ax.set_yticklabels(labels)

# Adding legend with larger title font size
legend_elements = [
    Line2D([0], [0], color=blue_rgba, lw=4, label='ATP and NAD(P)H'),
    Line2D([0], [0], color=orange_rgba, lw=4, label='ATP')
]

ax.legend(handles=legend_elements, loc='upper right', title='Energy carriers exchanged \n between catabolism and anabolism', title_fontsize=28)

plt.yticks(rotation=0)  # Rotate the y-axis labels to horizontal
plt.grid(zorder=0)
plt.xlabel('anabolic ATP fraction [-]')
plt.tight_layout()

fig_path = 'C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\Paper latex files\Figures'
plt.savefig(fig_path+'\\ATPfrac_Ecoli.pdf')
plt.savefig(fig_path+'\\ATPfrac_Ecoli.png')
# Show plot
plt.show()




#%%

nC_dict = {'glucose': 6, 
           'glutamate': 5, 
           'sorbitol': 6,
           'fructose': 6,
           'xylose': 5,
           'galactose': 6, 
           'ribose': 5,
           'alanine': 3,
           'succinate': 4,
           'glycerol': 3,
           'gluconate': 6,
           'aspartate': 4,
           'acetate': 2,
           'tryptophan': 11,
           'pyruvate': 3, 
           'citrate': 6,
           'lactate': 3, 
           'malate': 4, 
           'fumarate': 4,
           'akg': 5,
           'arginine': 6,
           'glutamine': 5}

dor_dict = {'glucose': 24, 
           'glutamate': 24, 
           'sorbitol': 26,
           'fructose': 24,
           'xylose': 20,
           'galactose': 24, 
           'ribose': 20,
           'alanine': 18,
           'succinate': 12,
           'glycerol': 14,
           'gluconate': 22,
           'aspartate': 18,
           'acetate': 8,
           'tryptophan': 58,
           'pyruvate': 10, 
           'citrate': 18,
           'lactate': 12, 
           'malate': 12, 
           'fumarate': 12,
           'akg': 16,
           'arginine': 22,
           'glutamine': 18}

exch_reac_dict = {'glucose': 'R_EX_glc__D_e_rev', 
           'glutamate': 'R_EX_glu__L_e_rev', 
           'sorbitol': 'R_EX_sbt__D_e_rev',
           'fructose': 'R_EX_fru_e_rev',
           'xylose': 'R_EX_xyl__D_e_rev',
           'galactose': 'R_EX_gal_e_rev', 
           'ribose': 'R_EX_rib__D_e_rev',
           'alanine': 'R_EX_ala__L_e_rev',
           'succinate': 'R_EX_succ_e_rev',
           'glycerol': 'R_EX_glyc_e_rev',
           'gluconate': 'R_EX_glcn_e_rev',
           'aspartate': 'R_EX_asp__L_e_rev',
           'acetate': 'R_EX_ac_e_rev',
           'tryptophan': 'R_EX_trp__L_e_rev',
           'pyruvate': 'R_EX_pyr_e_rev', 
           'citrate': 'R_EX_cit_e_rev',
           'lactate': 'R_EX_lac__D_e_rev', 
           'malate': 'R_EX_mal__L_e_rev', 
           'fumarate': 'R_EX_fum_e_rev',
           'akg': 'R_EX_akg_e_rev',
           'arginine': 'R_EX_arg__L_e_rev',
           'glutamine': 'R_EX_gln__L_e_rev'}

#%%

directory = 'Results\Ecoli_carbonsources'
df_dict_allfluxes = {}


# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        if 'NRC' in filename.path:
            continue
        df_dict_allfluxes[filename.path[28:-5]] = pd.read_excel(filename.path, sheet_name='overall', index_col=0)

directory = 'Results\Ecoli_nitrogensources'

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        if filename.path[30:-5] == 'glucose_nitrate' or filename.path[30:-5] == 'glucose_nitrate_nadh':
            continue
        df_dict_allfluxes[filename.path[30:-5]] = pd.read_excel(filename.path, sheet_name='overall', index_col=0)

subsyieldC_dict = {}
subsyielde_dict = {}
subsyield_dict = {}

for name, df in df_dict_allfluxes.items():
    for substrate in list(exch_reac_dict.keys()):
        if substrate in name:
            exch = exch_reac_dict[substrate]
            subsname = substrate
            # print(exch)
    

    print(name)
    print(exch)
    subsyield = 1000/df_dict_allfluxes[name].loc[exch, 'biomass'] # [gDW/mol substrate]
    subsyieldC = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*nC_dict[subsname]) # [gDW/Cmol substrate]
    subsyielde = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*dor_dict[subsname]) # [gDW/emol substrate]
    subsyieldeC = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*dor_dict[subsname]/nC_dict[subsname])
    
    subsyield_dict[name] = subsyield
    subsyieldC_dict[name] = subsyieldC
    subsyielde_dict[name] = subsyielde
    
    df_fraction.loc[name, 'yield on substrate'] = subsyield
    df_fraction.loc[name, 'yield on C'] = subsyieldC
    df_fraction.loc[name, 'yield on e'] = subsyielde
    df_fraction.loc[name, 'yield on e per C'] = subsyieldeC
    df_fraction.loc[name, 'e/C'] = dor_dict[subsname]/nC_dict[subsname]
    
#%%

fig, ax = plt.subplots(3, 4, figsize=(20,10))

df_fraction.plot.bar(y='ATP fraction', ax=ax[0,0], zorder=3)

plt.xticks(rotation=90)
ax[0,0].grid(zorder=0)
ax[0,0].set_ylabel('fraction of ATP used produced in anabolism [-]')


df_fraction.plot.bar(y='yield on substrate', ax=ax[0,1], zorder=3)

plt.xticks(rotation=90)
ax[0,1].grid(zorder=0)
ax[0,1].set_ylabel('biomass yield on substrate [gDW/mol substrate]')


df_fraction.plot.bar(y='yield on C', ax=ax[0,2], zorder=3)

plt.xticks(rotation=90)
ax[0,2].grid(zorder=0)
ax[0,2].set_ylabel('biomass yield on carbon [gDW/Cmol substrate]')


df_fraction.plot.bar(y='yield on e', ax=ax[0,3], zorder=3)

plt.xticks(rotation=90)
ax[0,3].grid(zorder=0)
ax[0,3].set_ylabel('biomass yield on electrons [gDW/emol substrate]')


df_fraction.plot.bar(y='yield on e per C', ax=ax[1,0], zorder=3)

plt.xticks(rotation=90)
ax[1,0].grid(zorder=0)
ax[1,0].set_ylabel('biomass yield on electrons per carbon in the substrate')


df_fraction.plot.bar(y='e/C', ax=ax[1,1], zorder=3)

plt.xticks(rotation=90)
ax[1,1].grid(zorder=0)
ax[1,1].set_ylabel('emol/Cmol')

for name, df in df_dict_allfluxes.items():
    
    if 'R_EX_co2_e_fwd' in df_dict_allfluxes[name].index:
        byproduct = df_dict_allfluxes[name].loc['R_EX_co2_e_fwd', 'biomass']
    else:
        byproduct = 0

    df_fraction.loc[name, 'byproduct'] = byproduct



df_fraction.plot.bar(y='byproduct', ax=ax[1,2], zorder=3)

plt.xticks(rotation=90)
ax[1,2].grid(zorder=0)
ax[1,2].set_ylabel('byproduct CO2 (mmol/gDW/h)')


for k, df in df_dict.items():
    df_fraction.loc[k, 'Yatp'] = 1000/df.loc['tot', 'cons']


df_fraction.plot.bar(y='Yatp', ax=ax[1,3], zorder=3)

plt.xticks(rotation=90)
ax[1,3].grid(zorder=0)
ax[1,3].set_ylabel('Yatp [gDW/molATP]')
plt.tight_layout()

for k, df in df_dict.items():
    df_fraction.loc[k, 'ATP transferred'] = df.loc['cat', 'net']


df_fraction.plot.bar(y='ATP transferred', ax=ax[2,0], zorder=3)

plt.xticks(rotation=90)
ax[2,0].grid(zorder=0)
ax[2,0].set_ylabel('')
plt.tight_layout()


#%% making a matrix of plots against each other
# each point should be a different color/thing, corresponding to the growth condition

# we want to plot all the things in df_fraction against each other, so we need a 9x9 matrix, that is filled half
# ax[1,1].set_visible(False) --> to clear the plots that we don't need


from matplotlib.lines import Line2D

fig, ax = plt.subplots(9, 9, figsize=(40,20))
 # --> to clear the plots that we don't need
# added 7 things

colors = np.array(["red","green","blue","yellow","pink","black","orange","purple","beige","brown","gray","cyan"])
not_done_list = list(df_fraction.columns)
for i, col1 in enumerate(df_fraction.columns):
    for j, col2 in enumerate(not_done_list):
        # if col1 != col2:
        if i < j:
            ax[i,j].scatter(df_fraction[col2][:12], df_fraction[col1][:12], color=colors, marker='o')
    
colors = np.array(["red","green","blue","yellow","pink","black","orange","purple","beige","brown","gray","cyan"])
not_done_list = list(df_fraction.columns)
for i, col1 in enumerate(df_fraction.columns):
    for j, col2 in enumerate(not_done_list):
        # if col1 != col2:
        if i < j:
            ax[i,j].scatter(df_fraction[col2][12:], df_fraction[col1][12:], color=colors[:-1], marker='^')

cols = df_fraction.columns
rows = df_fraction.columns

for a, col in zip(ax[0], cols):
    a.set_title(col)

i=1
j=0
for row in rows[:-1]:
    ax[j,i].set_ylabel(row, size='large')
    i+=1
    j+=1

for i in range(9):
    for j in range(9):
        if i >= j:  # Only plot if it's the first column or row number is greater than or equal to column number minus 1
            ax[i, j].axis('off')

plt.suptitle('x axis -->', size='x-large')
plt.text(-0.1, 0.5, 'y axis -->', rotation=90, size = 'x-large', va='center', ha='center', transform=ax[4, 0].transAxes)

handles = list(df_fraction.index)
# Create custom legend
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) for color, label in zip(colors, handles[:12])] + [Line2D([0], [0], marker='^', color='w', markerfacecolor=color, markersize=10, label=label) for color, label in zip(colors[:-1], handles[12:])]



fig.legend(handles=legend_elements, loc='lower center')
plt.tight_layout()

plt.savefig('matrixfig.pdf')














