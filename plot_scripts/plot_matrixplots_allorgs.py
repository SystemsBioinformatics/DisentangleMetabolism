# -*- coding: utf-8 -*-
"""
Created on Tue May  7 15:50:42 2024

@author: mre283
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import seaborn as sns


#%%

df_dict = {}


directories = ['Results\Cljungdahlii', 'Results\Gmetallireducens', 'Results\Mbarkeri', 'Results\Scerevisiae', 'Results\Pputida', 'Results\Ecoli_carbonsources', 'Results\Ecoli_nitrogensources', 'Results\Synechocystis'] # , 'Results\Synechocystis'

for directory in directories:
    for filename in os.scandir(directory):
        if filename.is_file():
            if 'mba' in filename.path:
                continue
            if 'nadh' in filename.path and not 'Synechocystis' in filename.path:
                continue
            if 'atp' in filename.path and 'Synechocystis' in filename.path:
                continue
            if 'NRC' in filename.path or 'nitrate' in filename.path:
                continue
            if 'respirofermentative' in filename.path or 'overflow' in filename.path:
                continue
            print(filename.path)
            df_dict[filename.path[8:-5]] = pd.read_excel(filename.path, sheet_name='atp_overall', index_col=0)

colordict = {}
colors = [
    'red',
    'blue',
    'green',
    'yellow',
    'cyan',
    'magenta',
    'orange',
    'pink',
    'purple',
    'brown',
    'lime',
    'navy',
    'gold',
    'teal',
    'indigo',
    'violet',
    'maroon',
    'olive',
    'coral',
    'turquoise']

for i, directory in enumerate(directories):
    colordict[directory[8:]] = colors[i]
    # if directory[8:] == 'Ecoli_nitrogensources':
    #     colordict[directory[8:]] = colordict['Ecoli_carbonsources']

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
           'glutamine': 5,
           'carbondioxide': 1,
           'carbonmonoxide': 1,
           'ethanol': 2,
           'formate': 1, 
           'photon':0,
           'bicarbonate': 1, 
           'hydrogen': 0,
           'butyrate': 4,
           'hexadecanoate161': 16,
           'methanol': 1,
           'octadecanoate': 18,
           'putrescine': 4,
           'maltotriose': 18}

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
           'photon':1, # what is the amount of electrons in a photon???
           'malate': 12, 
           'fumarate': 12,
           'akg': 16,
           'arginine': 22,
           'glutamine': 18,
           'carbonmonoxide': 2,
           'hydrogen': 2,
           'ethanol': 12,
           'formate': 2,
           'butyrate': 20,
           'hexadecanoate161': 90,
           'methanol': 6,
           'octadecanoate': 104,
           'putrescine': 22, 
           'maltotriose': 72}

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
           'glutamine': 'R_EX_gln__L_e_rev',
           'carbondioxide': 'R_EX_co2_e_rev',
           'hydrogen': 'R_EX_h2_e_rev',
           'carbonmonoxide': 'R_EX_co_e_rev',
           'ethanol': 'R_EX_etoh_e_rev',
           'formate': 'R_EX_for_e_rev',
           'butyrate': 'R_EX_but_e_rev',
           'hexadecanoate161': 'R_EX_hdcea_e_rev',
           'methanol': 'R_EX_meoh_e_rev',
           'octadecanoate': 'R_EX_ocdca_e_rev',
           'putrescine': 'R_EX_ptrc_e_rev',
           'photon':'R_EX_photon_e_rev',
           'maltotriose': 'R_EX_malttr_e_rev'}#,
           # 'acetateMbarkeri': 'R_EX_ac_LPAREN_e_RPAREN__rev',
           # 'hydrogenMbarkeri': 'R_EX_h2_LPAREN_e_RPAREN__rev'}


# M barkeri acetate exchange: R_EX_ac_LPAREN_e_RPAREN_
# Synechocystis: what is the electron donor?????

#%%

df_dict_allfluxes = {}
df_fraction = pd.DataFrame(columns=['yield on substrate'], index=df_dict.keys())
df_fraction_new = pd.DataFrame(columns=['Y_ATP/E'], index=df_dict.keys())

# iterate over files in 
# that directory
for directory in directories:
    for filename in os.scandir(directory):
        if filename.is_file():
            if 'mba' in filename.path:
                continue
            if 'nadh' in filename.path and not 'Synechocystis' in filename.path:
                continue
            if 'atp' in filename.path and 'Synechocystis' in filename.path:
                continue
            if 'NRC' in filename.path or 'nitrate' in filename.path:
                continue
            if 'respirofermentative' in filename.path or 'overflow' in filename.path:
                continue
            print(filename.path)
            df_dict_allfluxes[filename.path[8:-5]] = pd.read_excel(filename.path, sheet_name='overall', index_col=0)


subsyieldC_dict = {}
subsyielde_dict = {}
subsyield_dict = {}

for name, df in df_dict_allfluxes.items():
    if 'Mbarkeri' in name:
        if 'acetate' in name:
    
            exch = 'R_EX_ac_LPAREN_e_RPAREN__rev'
            
            subsname = 'acetate'
        else:
            exch = 'R_EX_h2_LPAREN_e_RPAREN__rev'
        
            subsname = 'hydrogen'
    elif 'Scerevisiae' in name and 'lactate' in name:
        exch = 'R_EX_lac__L_e_rev'
        subsname='lactate'
            
    else:
        for substrate in list(exch_reac_dict.keys()):
            
            if substrate in name:
                exch = exch_reac_dict[substrate]
                subsname = substrate
            
            # print(exch)
    
    
    
    if subsname == 'hydrogen' or subsname == 'photon':
        for subs in ['carbondioxide', 'carbonmonoxide', 'bicarbonate']:
            if subs in name:
                Csubsname = subs
    else:
        Csubsname = subsname
    
    print(exch)
    print(name)
    subsyield = 1000/df_dict_allfluxes[name].loc[exch, 'biomass'] # [gDW/mol substrate]
    subsyieldC = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*nC_dict[Csubsname]) # [gDW/Cmol substrate]
    subsyielde = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*dor_dict[subsname]) # [gDW/emol substrate]
    subsyieldeC = 1000/(df_dict_allfluxes[name].loc[exch, 'biomass']*dor_dict[subsname]/nC_dict[Csubsname])
    
    if 'R_ATP_demand' in df_dict_allfluxes[name].columns:
        ATPyieldE = df_dict_allfluxes[name].loc['R_ATP_demand', 'catabolic']/df_dict_allfluxes[name].loc[exch, 'catabolic']
        ATPyieldE_atp = df_dict_allfluxes[name].loc['R_ATP_demand', 'R_ATP_demand']/df_dict_allfluxes[name].loc[exch, 'R_ATP_demand']
    else:
        ATPyieldE = np.nan
        ATPyieldE_atp = np.nan
    if 'R_NADPH_demand' in df_dict_allfluxes[name].columns:
        NADPHyieldE = df_dict_allfluxes[name].loc['R_NADPH_demand', 'catabolic']/df_dict_allfluxes[name].loc[exch, 'catabolic']
        print(name)
        NADPHyieldE_nadph = df_dict_allfluxes[name].loc['R_NADPH_demand', 'R_NADPH_demand']/df_dict_allfluxes[name].loc[exch, 'R_NADPH_demand']
    else:
        NADPHyieldE = np.nan
        NADPHyieldE_nadph = np.nan
    if 'R_NADH_demand' in df_dict_allfluxes[name].columns:
        NADHyieldE = df_dict_allfluxes[name].loc['R_NADH_demand', 'catabolic']/df_dict_allfluxes[name].loc[exch, 'catabolic']
        NADHyieldE_nadh = df_dict_allfluxes[name].loc['R_NADH_demand', 'R_NADH_demand']/df_dict_allfluxes[name].loc[exch, 'R_NADH_demand']
    else:
        NADHyieldE = np.nan
        NADHyieldE_nadh=np.nan
    if 'R_FDX_demand' in df_dict_allfluxes[name].columns:
        FDyieldE = df_dict_allfluxes[name].loc['R_FDX_demand', 'catabolic']/df_dict_allfluxes[name].loc[exch, 'catabolic']
        FDyieldE_fd = df_dict_allfluxes[name].loc['R_FDX_demand', 'R_FDX_demand']/df_dict_allfluxes[name].loc[exch, 'R_FDX_demand']
    else:
        FDyieldE = np.nan
        FDyieldE_fd = np.nan
    
    subsyield_dict[name] = subsyield
    subsyieldC_dict[name] = subsyieldC
    subsyielde_dict[name] = subsyielde
    
    df_fraction.loc[name, 'yield on substrate'] = subsyield
    df_fraction.loc[name, 'yield on C'] = subsyieldC
    df_fraction.loc[name, 'yield on e'] = subsyielde
    df_fraction.loc[name, 'yield on e per C'] = subsyieldeC
    df_fraction.loc[name, 'e/C'] = dor_dict[subsname]/nC_dict[Csubsname]
    df_fraction_new.loc[name, 'Y_ATP/E[cat]'] = ATPyieldE
    df_fraction_new.loc[name, 'Y_NADPH/E[cat]'] = NADPHyieldE
    df_fraction_new.loc[name, 'Y_NADH/E[cat]'] = NADHyieldE
    df_fraction_new.loc[name, 'Y_FD/E[cat]'] = FDyieldE
    df_fraction_new.loc[name, 'Y_ATP/E[atp]'] = ATPyieldE_atp
    df_fraction_new.loc[name, 'Y_NADPH/E[nadph]'] = NADPHyieldE_nadph
    df_fraction_new.loc[name, 'Y_NADH/E[nadh]'] = NADHyieldE_nadh
    df_fraction_new.loc[name, 'Y_FD/E[fd]'] = FDyieldE_fd
    
    
    
for name, df in df_dict_allfluxes.items():
    
    if 'R_EX_co2_e_fwd' in df_dict_allfluxes[name].index:
        byproduct = df_dict_allfluxes[name].loc['R_EX_co2_e_fwd', 'biomass']
    else:
        byproduct = 0

    df_fraction.loc[name, 'byproduct'] = byproduct

for k, df in df_dict.items():
    df_fraction.loc[k, 'Yatp'] = 1000/df.loc['tot', 'cons']
    df_fraction_new.loc[k, 'ATP prod ana'] = df.loc['ana', 'prod']

for k, df in df_dict.items():
    df_fraction.loc[k, 'ATP transferred'] = df.loc['cat', 'net']


for k, df in df_dict.items():
    df_fraction.loc[k, 'ATP fraction'] = df.loc['ana', 'prod']/df.loc['ana', 'cons']

df_fraction = df_fraction.sort_values(by=['ATP fraction'], ascending=False)


# df_fraction_new = df_fraction_new.fillna(0)

#%%

df_fraction_plot = df_fraction.drop('byproduct', axis=1)


#%%

import matplotlib.pyplot as plt

name_dict = {
    'yield on substrate': '$Y_{X/E}$ [gDW/mol E]',
    'yield on C': '$Y_{X/C}$ [gDW/Cmol E]',
    'yield on e': '$Y_{X/e}$ [gDW/emol E]',
    'yield on e per C': '$Y_{X/eC}$ [gDW/(emol/Cmol)]',
    'e/C': '$\gamma _E$ [mol e/mol C]',
    'Yatp': '$Y_{X/ATP}$ [gDW/molATP]', 
    'ATP transferred': '$(Y_{X/ATP}^{cat})^{-1}$ [mmolATP/gDW]',
    'ATP fraction': '$(Y_{X/ATP}^{ana})^{-1}/(Y_{X/ATP})^{-1}$ [-]'
}

df_fraction_plot = df_fraction_plot.rename(name_dict, axis=1)
sns.set_theme(style="whitegrid")

corr_mat = df_fraction_plot.corr()

corr_mat = corr_mat.stack().reset_index(name="correlation")

corr_mat2 = df_fraction_plot.corr()**2
corr_mat2 = corr_mat2.stack().reset_index(name="r squared")
corr_mat['r squared'] = corr_mat2['r squared']


g = sns.relplot(
    data=corr_mat,
    x="level_0", y="level_1", hue="correlation", size="r squared",
    palette="vlag", hue_norm=(-1, 1), edgecolor=".7",
    height=12, sizes=(200, 1200), size_norm=(-.2, .8),
)

# Tweak the figure to finalize
g.set(xlabel="", ylabel="", aspect="equal")
g.despine(left=True, bottom=True)
g.ax.margins(.05)

# Increase tick label size
g.ax.tick_params(axis='x', labelsize=24)
g.ax.tick_params(axis='y', labelsize=24)

# Rotate x labels
for label in g.ax.get_xticklabels():
    label.set_rotation(90)


# Increase legend title and text size
plt.setp(g.legend.get_texts(), fontsize=24)  # for legend text

# Make dots in the legend larger
for legend_handle in g.legend.legendHandles[:6]:
    legend_handle._sizes = [700]  # Adjust the size of the legend dots

# Move the legend to the right and up, and increase the space between entries
g.legend.set_bbox_to_anchor((1.05, 0.7))  # Shift legend to the right and up
g.legend.set_loc('center left')  # Align legend to the center left
g.legend.set_frame_on(False)  # Optional: Remove the legend border

plt.savefig(r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\corellogram.svg', bbox_inches='tight')

plt.show()
    


#%%


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib as mpl

# Set up LaTeX text rendering and font properties
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

# Define colordict and display_names
colordict = {'Cljungdahlii': '#1f77b4',  # Blue
             'Gmetallireducens': '#ff7f0e',  # Orange
             'Mbarkeri': '#2ca02c',  # Green
             'Scerevisiae': '#d62728',  # Red
             'Pputida': '#9467bd',  # Purple
             'Ecoli': '#8c564b',  # Brown
             'Synechocystis': '#e377c2'}  # Pink

display_names = {
    'Cljungdahlii': r'\textit{C. ljungdahlii}',
    'Gmetallireducens': r'\textit{G. metallireducens}',
    'Mbarkeri': r'\textit{M. barkeri}',
    'Scerevisiae': r'\textit{S. cerevisiae}',
    'Pputida': r'\textit{P. putida}',
    'Ecoli': r'\textit{E. coli}',
    'Synechocystis': r'\textit{Synechocystis}'
}

# Assuming df_fraction_plot is defined
fig, ax = plt.subplots(8, 8, figsize=(20, 25))

# Assign colors based on colordict
colors = []
for sim in df_fraction_plot.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict[key])

not_done_list = list(df_fraction_plot.columns)
for i, col1 in enumerate(df_fraction_plot.columns):
    for j, col2 in enumerate(not_done_list):
        if i < j:
            ax[i, j].scatter(df_fraction_plot[col2], df_fraction_plot[col1], color=colors)

cols = df_fraction_plot.columns
rows = df_fraction_plot.columns

for a, col in zip(ax[0][1:], cols[1:]):
    a.set_title(col, fontsize=28)  # Set font size for column titles

i = 1
j = 0
for row in rows[:-1]:
    ax[j, i].set_ylabel(row, size=28)  # Set font size for row labels
    i += 1
    j += 1

for i in range(8):
    for j in range(8):
        if i >= j:
            ax[i, j].axis('off')
        else:
            ax[i, j].tick_params(axis='both', which='major', labelsize=20)  # Set font size for tick labels

# Set font size for x-axis labels
for i in range(1, 8):
    ax[7, i].set_xlabel(cols[i], fontsize=28)

# Create custom legend handles with formatted labels
handles = []
for key, color in colordict.items():
    display_key = display_names[key]
    handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
                          markerfacecolor=color, markersize=15))  # Set font size for legend markers

# Add legend to the figure
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.2, 0.5), title='Organisms', title_fontsize=28, fontsize=28)  # Set font size for legend

plt.tight_layout()

# plt.savefig('matrixfig_allorgs.pdf', bbox_inches='tight')
plt.show()


#%%

import matplotlib as mpl
mpl.rcParams['text.usetex'] = True
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = 'DejaVu Sans'
mpl.rcParams['text.latex.preamble'] = r'\usepackage{amsmath}'

colordict = {'Cljungdahlii': '#1f77b4',  # Blue
             'Gmetallireducens': '#ff7f0e',  # Orange
             'Mbarkeri': '#2ca02c',  # Green
             'Scerevisiae': '#d62728',  # Red
             'Pputida': '#9467bd',  # Purple
             'Ecoli': '#8c564b',  # Brown
             'Synechocystis': '#e377c2'}  # Pink

# Create display name mapping for the legend with LaTeX formatting
display_names = {
    'Cljungdahlii': r'\textit{C. ljungdahlii}',
    'Gmetallireducens': r'\textit{G. metallireducens}',
    'Mbarkeri': r'\textit{M. barkeri}',
    'Scerevisiae': r'\textit{S. cerevisiae}',
    'Pputida': r'\textit{P. putida}',
    'Ecoli': r'\textit{E. coli}',
    'Synechocystis': r'\textit{Synechocystis}'
}

fig, ax = plt.subplots(1, 2)

colors = []
for sim in df_fraction.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict[key])

# Plotting the data
ax[0].scatter(df_fraction['yield on e'], df_fraction['ATP fraction'], color=colors)
ax[0].set_xlabel(r'Y\textsubscript{X/e}\, [gDW mmol\textsubscript{e}\textsuperscript{-1}]')
ax[0].set_ylabel('ATP fraction')

ax[1].scatter(df_fraction['Yatp'], df_fraction['ATP fraction'], color=colors)
ax[1].set_ylabel('ATP fraction')
ax[1].set_xlabel(r'Y\textsubscript{ATP}\, [gDW mol\textsubscript{ATP}\textsuperscript{-1}]')

# Create custom legend handles with formatted labels
handles = []
for key, color in colordict.items():
    display_key = display_names[key]
    handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
                          markerfacecolor=color, markersize=10))

# Add legend to the figure
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), title='Organisms')

plt.tight_layout()

# plt.savefig('yields_ATPfraction.pdf', bbox_inches='tight')
plt.show()



#%%

df_dict_anabolic_mceq = {}

# iterate over files in 
# that directory
for directory in directories:
    for filename in os.scandir(directory):
        if filename.is_file():
            if 'mba' in filename.path:
                continue
            if 'nadh' in filename.path and not 'Synechocystis' in filename.path:
                continue
            if 'atp' in filename.path and 'Synechocystis' in filename.path:
                continue
            if 'NRC' in filename.path or 'nitrate' in filename.path:
                continue
            if 'respirofermentative' in filename.path or 'overflow' in filename.path:
                continue
            print(filename.path)
            
            df_dict_anabolic_mceq[filename.path[8:-5]] = pd.read_excel(filename.path, sheet_name='mceqs', index_col=0).loc['anabolic', 'overall']

#%%
import sympy

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

# # Example chemical equation
# mceq = "2*A + 3*B - C + 0.5*D - 4.2*E"

# # Parse and extract
# expr = sympy.parsing.sympy_parser.parse_expr(mceq)
# metabolite_dict = extract_metabolites(expr)

# Print the resulting dictionary
# print(metabolite_dict)

df_mceqs = pd.DataFrame()
for sim, mceq in df_dict_anabolic_mceq.items():
    if mceq[0] == '[':
        mceq = mceq[1:-1]
    
    symp = sympy.parsing.sympy_parser.parse_expr(mceq)
    metabolite_dict = extract_metabolites(symp)
    for metab, coeff in metabolite_dict.items():
        df_mceqs.loc[sim, metab] = coeff
    
    
df_mceqs.fillna(0, inplace=True)
    
#%%

dor_dict = {'M_glc__D_e': 24, 
           'M_glu__L_e': 24, 
           'M_sbt__D_e': 26,
           'M_fru_e': 24,
           'M_xyl__D_e': 20,
           'M_gal_e': 24, 
           'M_rib__D_e': 20,
           'M_ala__L_e': 18,
           'M_succ_e': 12,
           'M_glyc_e': 14,
           'M_glcn_e': 22,
           'M_asp__L_e': 18,
           'M_ac_e': 8,
           'M_trp__L_e': 58,
           'M_pyr_e': 10, 
           'M_cit_e': 18,
           'M_lac__D_e': 12,
           'M_lac__L_e':12,
           # 'M_photon_e':1, # TODO what is the amount of electrons in a photon??? 
           'M_mal__L_e': 12, 
           'M_fum_e': 12,
           'M_akg_e': 16,
           'M_arg__L_e': 22,
           'M_gln__L_e': 18,
           'M_co_e': 2,
           'M_etoh_e': 12,
           'M_for_e': 2,
           'M_but_e': 20,
           'M_hdcea_e': 90,
           'M_meoh_e': 6,
           'M_ocdca_e': 104,
           'M_ptrc_e': 22, 
           'M_malttr_e': 72,
           'M_indole_e': 36}

df_deltae = pd.DataFrame(np.zeros(len(df_mceqs.index)),index=df_mceqs.index, columns=['delta e'])
for metab, dor in dor_dict.items():
    for sim in df_mceqs.index:
        df_deltae.loc[sim, 'delta e'] += df_mceqs.loc[sim, metab]*dor

for sim in df_deltae.index:
    df_deltae.loc[sim, 'delta e']+= 42*4.2

#%%
df_deltae['ATP fraction'] = np.zeros(len(df_deltae))
for sim in df_deltae.index:
    df_deltae.loc[sim, 'ATP fraction'] = df_fraction.loc[sim, 'ATP fraction']

plt.figure()
plt.plot(df_deltae['delta e'], df_deltae['ATP fraction'], 'o', label='simulations')
# plt.plot([-42*4.2, -42*4.2], [-0.2, 1], '--', label='dor of biomass')

plt.legend()
plt.xlabel('delta e')
plt.ylabel('ATP fraction in ana')

#%%

# dor_dict = {'M_glc__D_e': 24, 
#            'M_glu__L_e': 24, 
#            'M_sbt__D_e': 26,
#            'M_fru_e': 24,
#            'M_xyl__D_e': 20,
#            'M_gal_e': 24, 
#            'M_rib__D_e': 20,
#            'M_ala__L_e': 18,
#            'M_succ_e': 12,
#            'M_glyc_e': 14,
#            'M_glcn_e': 22,
#            'M_asp__L_e': 18,
#            'M_ac_e': 8,
#            'M_trp__L_e': 58,
#            'M_pyr_e': 10, 
#            'M_cit_e': 18,
#            'M_lac__D_e': 12,
#            'M_lac__L_e':12,
#            # 'M_photon_e':1, # TODO what is the amount of electrons in a photon??? 
#            'M_mal__L_e': 12, 
#            'M_fum_e': 12,
#            'M_akg_e': 16,
#            'M_arg__L_e': 22,
#            'M_gln__L_e': 18,
#            'M_co_e': 2,
#            'M_etoh_e': 12,
#            'M_for_e': 2,
#            'M_but_e': 20,
#            'M_hdcea_e': 90,
#            'M_meoh_e': 6,
#            'M_ocdca_e': 104,
#            'M_ptrc_e': 22, 
#            'M_malttr_e': 72,
#            'M_indole_e': 36,
#            'M_o2_e': -2}

# df_deltae = pd.DataFrame(np.zeros(len(df_mceqs.index)),index=df_mceqs.index, columns=['delta e'])
# for metab, dor in dor_dict.items():
#     for sim in df_mceqs.index:
#         df_deltae.loc[sim, 'delta e'] += df_mceqs.loc[sim, metab]*dor

# for sim in df_deltae.index:
#     df_deltae.loc[sim, 'delta e']+= 42*4.2

# #%%

# plt.figure()
# plt.plot(df_deltae['delta e'], df_deltae['ATP fraction'], 'o', label='simulations')
# # plt.plot([-42*4.2, -42*4.2], [-0.2, 1], '--', label='dor of biomass')

# plt.legend()
# plt.xlabel('delta e')
# plt.ylabel('ATP fraction in ana')

#%%

# df_excel = pd.read_excel('C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\python files\generalized_scripts\\results_categories.xlsx', header=0)

# for i, sim in enumerate(list(df_fraction.index)):
#     df_excel.loc[i, 'Simulation'] = sim

# df_excel.to_excel('C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\python files\generalized_scripts\\results_categories.xlsx')


#%%

df_cate = pd.read_excel('C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\python files\generalized_scripts\\results_categories.xlsx', header=0)
df_cate = df_cate.drop(index=44)
cate_dict = dict(zip(df_cate['Simulation'], df_cate['Category']))

colors = {'1 external e.a. used': 'tab:blue', '2 internal e.d. used': 'tab:orange',  '3 external e.d. used':'tab:purple', '4 internal e.d. used + external e.a. used':'tab:green'}
color_dict = {}

for sim, cat in cate_dict.items():
    color_dict[sim] = colors[cat]

fig, ax = plt.subplots()

ax.scatter(df_deltae['delta e'], df_deltae['ATP fraction'], color=list(color_dict.values()))

# handles = ['1 external e.a. used', '2 internal e.d. used',  '3 external e.d. used','4 internal e.d. used + external e.a. used']
# Create custom legend
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) for label, color in colors.items()]

fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.9, 0.5))
ax.set_ylabel('Anabolic ATP fraction')
ax.set_xlabel('$\Delta$e in carbon-containing molecules')


#%%
colors = []
for sim in df_deltae.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict[key])



fig, ax = plt.subplots()

ax.scatter(df_deltae['delta e'], df_deltae['ATP fraction'], color=colors)

# handles = ['1 external e.a. used', '2 internal e.d. used',  '3 external e.d. used','4 internal e.d. used + external e.a. used']
# Create custom legend
legend_elements = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label) for label, color in colordict.items()]

fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.9, 0.5))
ax.set_ylabel('Anabolic ATP fraction')
ax.set_xlabel('$\Delta$e in carbon-containing molecules')

#%%
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

colors = []
for sim in df_deltae.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict[key])

# Separate the data points into two groups based on the marker type
data_o = {'x': [], 'y': [], 'color': []}
data_triangle = {'x': [], 'y': [], 'color': []}

for i, sim in enumerate(df_deltae.index):
    if '_test' in sim:
        sim = 'Gmetallireducens\\acetateFe3'
    if sim not in list(cate_dict.keys()):
        data_o['x'].append(df_deltae['delta e'][i])
        data_o['y'].append(df_deltae['ATP fraction'][i])
        data_o['color'].append(colors[i])
    elif cate_dict[sim] == '4 internal e.d. used + external e.a. used':
        data_triangle['x'].append(df_deltae['delta e'][i])
        data_triangle['y'].append(df_deltae['ATP fraction'][i])
        data_triangle['color'].append(colors[i])
    else:
        data_o['x'].append(df_deltae['delta e'][i])
        data_o['y'].append(df_deltae['ATP fraction'][i])
        data_o['color'].append(colors[i])

fig, ax = plt.subplots()

# Plot the data points with 'o' markers
ax.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')

# Plot the data points with '^' markers
ax.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# Create custom legend for markers
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=10, label='EA or ED'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label='both EA and ED')
]

# Add color legends
for label, color in colordict.items():
    legend_elements.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=label))

fig.legend(handles=legend_elements, loc='center left', bbox_to_anchor=(0.9, 0.5))
ax.set_ylabel('Anabolic ATP fraction')
ax.set_xlabel('$\Delta$e in carbon-containing molecules')

plt.show()


#%% Making the plot of electron acceptor used in anabolic MCEQ

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Assuming df_mceqs, df_fraction, and cate_dict are defined in your code

colordict = {
    'Cljungdahlii': '#1f77b4',  # Blue
    'Gmetallireducens': '#ff7f0e',  # Orange
    'Mbarkeri': '#2ca02c',  # Green
    'Scerevisiae': '#d62728',  # Red
    'Pputida': '#9467bd',  # Purple
    'Ecoli': '#8c564b'  # Brown
}

# Create display name mapping for the legend with LaTeX formatting
display_names = {
    'Cljungdahlii': r'\textit{C. ljungdahlii}',
    'Gmetallireducens': r'\textit{G. metallireducens}',
    'Mbarkeri': r'\textit{M. barkeri}',
    'Scerevisiae': r'\textit{S. cerevisiae}',
    'Pputida': r'\textit{P. putida}',
    'Ecoli': r'\textit{E. coli}'
}

electronacceptors = ['M_o2_e', 'M_fe3_e', 'M_no3_e']

df_elecacc = df_mceqs[electronacceptors]
df_elecacc.loc[:, 'sum'] = np.abs(df_elecacc['M_o2_e'] * 4 + df_elecacc['M_fe3_e'] + df_elecacc['M_no3_e'] * 2)

df_elecacc['ATP fraction'] = np.zeros(len(df_elecacc))
for sim in df_elecacc.index:
    df_elecacc.loc[sim, 'ATP fraction'] = df_fraction.loc[sim, 'ATP fraction']

colors = []
data_o = {'x': [], 'y': [], 'color': []}
data_triangle = {'x': [], 'y': [], 'color': []}

for i, sim in enumerate(df_elecacc.index):
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'
    if 'Synechocystis' in key:
        continue
    colors.append(colordict[key])
    if '_test' in sim:
        sim = 'Gmetallireducens\\acetateFe3'
    if sim not in cate_dict or cate_dict[sim] != '4 internal e.d. used + external e.a. used':
        data_o['x'].append(df_elecacc['sum'][i])
        data_o['y'].append(df_elecacc['ATP fraction'][i])
        data_o['color'].append(colors[i])
    else:
        data_triangle['x'].append(df_elecacc['sum'][i])
        data_triangle['y'].append(df_elecacc['ATP fraction'][i])
        data_triangle['color'].append(colors[i])

# Set font size for the entire plot
plt.rcParams.update({
    'font.size': 22,           # General font size
    'axes.titlesize': 26,      # Title font size
    'axes.labelsize': 24,      # Axis labels font size
    'xtick.labelsize': 24,     # X-tick labels font size
    'ytick.labelsize': 24,     # Y-tick labels font size
    'legend.fontsize': 24,     # Legend font size
    'legend.title_fontsize': 24,  # Legend title font size
})

fig, ax = plt.subplots(figsize=(15, 10))

# Plot the data points with 'o' markers
ax.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED', s=200)  # s controls marker size

# Plot the data points with '^' markers
ax.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED', s=200)

# Create custom legend for markers with larger markers
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=15, label='EA or ED'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=15, label='both EA and ED')
]

# Add color legends with larger markers
for key, color in colordict.items():
    display_key = display_names[key]
    handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
                          markerfacecolor=color, markersize=15))

# Add legend to the figure with larger title
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.95, 0.5), title='Organisms')

# Set axis labels with larger font size
ax.set_ylabel('Anabolic ATP fraction')
ax.set_xlabel(r'mmol e$^-$ accepted/gDW in ana')

plt.savefig('ATPfrac_eacc_ana.svg', bbox_inches='tight')
plt.show()


#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from brokenaxes import brokenaxes

# Adding Synechocystis to colordict and display_names
colordict['Synechocystis'] = '#e377c2'  # Pink color
display_names['Synechocystis'] = r'\textit{Synechocystis}'

# Assuming df_mceqs and df_fraction are already defined
# Electron acceptors
electronacceptors = ['M_o2_e', 'M_fe3_e', 'M_no3_e']

# Subset df_mceqs
df_elecacc = df_mceqs[electronacceptors]

# Calculate sum
df_elecacc['sum'] = np.abs(df_elecacc['M_o2_e'] * 4 + df_elecacc['M_fe3_e'] + df_elecacc['M_no3_e'] * 2)

# Initialize colors list
colors = []
for sim in df_elecacc.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict.get(key, '#000000'))  # Default to black if not found

# Add ATP fraction column
df_elecacc['ATP fraction'] = np.zeros(len(df_elecacc))
for sim in df_elecacc.index:
    df_elecacc.loc[sim, 'ATP fraction'] = df_fraction.loc[sim, 'ATP fraction']

# Data for plotting
data_o = {'x': [], 'y': [], 'color': []}
data_triangle = {'x': [], 'y': [], 'color': []}

# Prepare data for plotting
for sim in df_elecacc.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'
    if 'Synechocystis' in key:
        colors.append(colordict['Synechocystis'])
    else:
        colors.append(colordict.get(key, '#000000'))  # Default to black if not found
    if '_test' in sim:
        sim = 'Gmetallireducens\\acetateFe3'
    if sim not in cate_dict or cate_dict[sim] != '4 internal e.d. used + external e.a. used':
        data_o['x'].append(df_elecacc.loc[sim, 'sum'])
        data_o['y'].append(df_elecacc.loc[sim, 'ATP fraction'])
        data_o['color'].append(colors[-1])
    else:
        data_triangle['x'].append(df_elecacc.loc[sim, 'sum'])
        data_triangle['y'].append(df_elecacc.loc[sim, 'ATP fraction'])
        data_triangle['color'].append(colors[-1])

# Create brokenaxes object
fig = plt.figure(figsize=(8, 6))
bax = brokenaxes(xlims=((0, 110), (390, 410)), hspace=.05, fig=fig)

# Plot the data points with 'o' markers
bax.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')

# Plot the data points with '^' markers
bax.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# Create custom legend for markers
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=10, label='EA or ED'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label='both EA and ED')
]

# Add color legends
for key, color in colordict.items():
    display_key = display_names[key]
    handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
                          markerfacecolor=color, markersize=10))

# Add legend to the figure
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(1.0, 0.5), title='Organisms')

bax.set_ylabel('Anabolic ATP fraction')
bax.set_xlabel('mmol e$^-$ accepted/gDW in ana')

# plt.savefig('ATPfrac_eacc_ana.svg', bbox_inches='tight')
plt.show()


#%%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# Adding Synechocystis to colordict and display_names
colordict = {'Cljungdahlii': '#1f77b4',  # Blue
             'Gmetallireducens': '#ff7f0e',  # Orange
             'Mbarkeri': '#2ca02c',  # Green
             'Scerevisiae': '#d62728',  # Red
             'Pputida': '#9467bd',  # Purple
             'Ecoli': '#8c564b',  # Brown
             'Synechocystis': '#e377c2'}  # Pink

display_names = {
    'Cljungdahlii': r'\textit{C. ljungdahlii}',
    'Gmetallireducens': r'\textit{G. metallireducens}',
    'Mbarkeri': r'\textit{M. barkeri}',
    'Scerevisiae': r'\textit{S. cerevisiae}',
    'Pputida': r'\textit{P. putida}',
    'Ecoli': r'\textit{E. coli}',
    'Synechocystis': r'\textit{Synechocystis}'
}

# # Example DataFrames, replace these with your actual data
# df_mceqs = pd.DataFrame({
#     'M_o2_e': np.random.rand(10),
#     'M_fe3_e': np.random.rand(10),
#     'M_no3_e': np.random.rand(10)
# }, index=['Cljungdahlii\\test1', 'Gmetallireducens\\test2', 'Mbarkeri\\test3', 'Scerevisiae\\test4',
#           'Pputida\\test5', 'Ecoli\\test6', 'Synechocystis\\test7', 'Ecoli\\test8', 'Ecoli\\test9', 'Gmetallireducens\\test10'])

# df_fraction = pd.DataFrame({
#     'ATP fraction': np.random.rand(10)
# }, index=['Cljungdahlii\\test1', 'Gmetallireducens\\test2', 'Mbarkeri\\test3', 'Scerevisiae\\test4',
#           'Pputida\\test5', 'Ecoli\\test6', 'Synechocystis\\test7', 'Ecoli\\test8', 'Ecoli\\test9', 'Gmetallireducens\\test10'])

electronacceptors = ['M_o2_e', 'M_fe3_e', 'M_no3_e', 'M_h2_e']

# Subset df_mceqs to only include rows where M_h2_e is positive
df_elecacc = df_mceqs[df_mceqs['M_h2_e'] >= 0][electronacceptors]


# Calculate sum
df_elecacc['sum'] = np.abs(df_elecacc['M_o2_e'] * 4 + df_elecacc['M_fe3_e'] + df_elecacc['M_no3_e'] * 2  + df_elecacc['M_h2_e'] * 2)

# Initialize colors list
colors = []
for sim in df_elecacc.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
    colors.append(colordict.get(key, '#000000'))  # Default to black if not found

# Add ATP fraction column
df_elecacc['ATP fraction'] = np.zeros(len(df_elecacc))
for sim in df_elecacc.index:
    df_elecacc.loc[sim, 'ATP fraction'] = df_fraction.loc[sim, 'ATP fraction']

# Data for plotting
data_o = {'x': [], 'y': [], 'color': []}
data_triangle = {'x': [], 'y': [], 'color': []}

# Prepare data for plotting
# cate_dict = {}  # Assuming you have this dictionary defined
for sim in df_elecacc.index:
    key = sim.split('\\')[0]
    if 'Ecoli' in key:
        key = 'Ecoli'
    if 'Synechocystis' in key:
        colors.append(colordict['Synechocystis'])
    else:
        colors.append(colordict.get(key, '#000000'))  # Default to black if not found
    if '_test' in sim:
        sim = 'Gmetallireducens\\acetateFe3'
    if sim not in cate_dict or cate_dict[sim] != '4 internal e.d. used + external e.a. used':
        data_o['x'].append(df_elecacc.loc[sim, 'sum'])
        data_o['y'].append(df_elecacc.loc[sim, 'ATP fraction'])
        data_o['color'].append(colors[-1])
    else:
        data_triangle['x'].append(df_elecacc.loc[sim, 'sum'])
        data_triangle['y'].append(df_elecacc.loc[sim, 'ATP fraction'])
        data_triangle['color'].append(colors[-1])

# Create the figure and the subplots
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [8, 1]}, figsize=(8, 5))
fig.subplots_adjust(wspace=0.05)

# Plot the data points with 'o' markers on ax1
ax1.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')
ax1.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# Plot the data points with 'o' markers on ax2 (outliers)
ax2.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')
ax2.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# Setting the x-axis limits for the broken effect
ax1.set_xlim(-2, 155)
ax2.set_xlim(420, 425)

# Hiding the spines between ax1 and ax2
ax1.spines['right'].set_visible(False)
ax2.spines['left'].set_visible(False)
ax1.yaxis.tick_left()
ax1.tick_params(labelright=False)
ax2.yaxis.tick_right()

# Adding break markers
d = .015  # Size of break marker
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (-d, +d), **kwargs)
ax2.plot((-d, +d), (1-d, 1+d), **kwargs)

# Custom legend for markers
handles = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=10, label='EA or ED'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label='both EA and ED')
]

# Add color legends
for key, color in colordict.items():
    display_key = display_names[key]
    handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
                          markerfacecolor=color, markersize=10))

# Add legend to the figure
fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.95, 0.5), title='Organisms')

# Set labels
ax1.set_ylabel('Anabolic ATP fraction')
ax1.set_xlabel('mmol e$^-$ accepted/gDW in ana')
ax2.set_xlabel(' ')

# plt.savefig('ATPfrac_eacc_ana.svg', bbox_inches='tight')
plt.show()


#%%

# import numpy as np
# import pandas as pd
# import matplotlib.pyplot as plt
# from matplotlib.lines import Line2D

# # Adding Synechocystis to colordict and display_names
# colordict = {'Cljungdahlii': '#1f77b4',  # Blue
#              'Gmetallireducens': '#ff7f0e',  # Orange
#              'Mbarkeri': '#2ca02c',  # Green
#              'Scerevisiae': '#d62728',  # Red
#              'Pputida': '#9467bd',  # Purple
#              'Ecoli': '#8c564b',  # Brown
#              'Synechocystis': '#e377c2'}  # Pink

# display_names = {
#     'Cljungdahlii': r'\textit{C. ljungdahlii}',
#     'Gmetallireducens': r'\textit{G. metallireducens}',
#     'Mbarkeri': r'\textit{M. barkeri}',
#     'Scerevisiae': r'\textit{S. cerevisiae}',
#     'Pputida': r'\textit{P. putida}',
#     'Ecoli': r'\textit{E. coli}',
#     'Synechocystis': r'\textit{Synechocystis}'
# }


# electronacceptors = ['M_o2_e', 'M_fe3_e', 'M_no3_e', 'M_h2_e']

# # Subset df_mceqs to only include rows where M_h2_e is positive
# df_elecacc = df_mceqs[df_mceqs['M_h2_e'] >= 0][electronacceptors]


# # Calculate sum
# df_elecacc['sum'] = np.abs(df_elecacc['M_o2_e'] * 4 + df_elecacc['M_fe3_e'] + df_elecacc['M_no3_e'] * 2 + df_elecacc['M_h2_e'] * 2)

# # Initialize colors list
# colors = []
# for sim in df_elecacc.index:
#     key = sim.split('\\')[0]
#     if 'Ecoli' in key:
#         key = 'Ecoli'  # Treat Ecoli_carbonsources and Ecoli_nitrogensources the same
#     colors.append(colordict.get(key, '#000000'))  # Default to black if not found

# # Add ATP fraction column
# df_elecacc['ATP transferred'] = np.zeros(len(df_elecacc))
# for sim in df_elecacc.index:
#     df_elecacc.loc[sim, 'ATP transferred'] = df_fraction.loc[sim, 'ATP transferred']

# # Data for plotting
# data_o = {'x': [], 'y': [], 'color': []}
# data_triangle = {'x': [], 'y': [], 'color': []}

# # Prepare data for plotting
# # cate_dict = {}  # Assuming you have this dictionary defined
# for sim in df_elecacc.index:
#     key = sim.split('\\')[0]
#     if 'Ecoli' in key:
#         key = 'Ecoli'
#     if 'Synechocystis' in key:
#         colors.append(colordict['Synechocystis'])
#     else:
#         colors.append(colordict.get(key, '#000000'))  # Default to black if not found
#     if '_test' in sim:
#         sim = 'Gmetallireducens\\acetateFe3'
#     if sim not in cate_dict or cate_dict[sim] != '4 internal e.d. used + external e.a. used':
#         data_o['x'].append(df_elecacc.loc[sim, 'sum'])
#         data_o['y'].append(df_elecacc.loc[sim, 'ATP transferred'])
#         data_o['color'].append(colors[-1])
#     else:
#         data_triangle['x'].append(df_elecacc.loc[sim, 'sum'])
#         data_triangle['y'].append(df_elecacc.loc[sim, 'ATP transferred'])
#         data_triangle['color'].append(colors[-1])

# # Create the figure and the subplots
# fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, gridspec_kw={'width_ratios': [8, 1]}, figsize=(8, 5))
# fig.subplots_adjust(wspace=0.05)

# # Plot the data points with 'o' markers on ax1
# ax1.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')
# ax1.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# # Plot the data points with 'o' markers on ax2 (outliers)
# ax2.scatter(data_o['x'], data_o['y'], color=data_o['color'], marker='o', label='EA or ED')
# ax2.scatter(data_triangle['x'], data_triangle['y'], color=data_triangle['color'], marker='^', label='both EA and ED')

# # Setting the x-axis limits for the broken effect
# ax1.set_xlim(-2, 115)
# ax2.set_xlim(420, 425)

# # Hiding the spines between ax1 and ax2
# ax1.spines['right'].set_visible(False)
# ax2.spines['left'].set_visible(False)
# ax1.yaxis.tick_left()
# ax1.tick_params(labelright=False)
# ax2.yaxis.tick_right()

# # Adding break markers
# d = .015  # Size of break marker
# kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
# ax1.plot((1-d, 1+d), (-d, +d), **kwargs)
# ax1.plot((1-d, 1+d), (1-d, 1+d), **kwargs)
# kwargs.update(transform=ax2.transAxes)
# ax2.plot((-d, +d), (-d, +d), **kwargs)
# ax2.plot((-d, +d), (1-d, 1+d), **kwargs)

# # Custom legend for markers
# handles = [
#     Line2D([0], [0], marker='o', color='w', markerfacecolor='k', markersize=10, label='EA or ED'),
#     Line2D([0], [0], marker='^', color='w', markerfacecolor='k', markersize=10, label='both EA and ED')
# ]

# # Add color legends
# for key, color in colordict.items():
#     display_key = display_names[key]
#     handles.append(Line2D([0], [0], marker='o', color='w', label=display_key,
#                           markerfacecolor=color, markersize=10))

# # Add legend to the figure
# fig.legend(handles=handles, loc='center left', bbox_to_anchor=(0.95, 0.5), title='Organisms')

# # Set labels
# ax1.set_ylabel('ATP transferred')
# ax1.set_xlabel('mmol e$^-$ accepted/gDW in ana')
# ax2.set_xlabel(' ')

# plt.savefig('ATPfrac_eacc_ana.svg', bbox_inches='tight')
# plt.show()

#%%

table_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\yields_table_v2.xlsx'
# df_yields = pd.DataFrame(columns = ['Organism', 'Carbon source', 'Energy source', 'Electron acceptor', 'Growth mode', 'Y_X/ATP[cat]', 'Y_ATP/E', 'Y_X/ATP[tot]', 'Y_X/E', 'dor (e/C)'], index=df_fraction.index)
df_yields = pd.read_excel(table_path, index_col=0, header=0)
# df_yields2.index = list(df_yields.index).remove('Gmetallireducens\\acetateFe3_test')
# df_yields['Organism'] = df_yields2['Organism']
# df_yields['Carbon source'] = df_yields2['Carbon source']
# df_yields['Energy source'] = df_yields2['Energy source']
# df_yields['Growth mode'] = df_yields2['Growth mode']
# df_yields['Electron acceptor'] = df_yields2['Electron acceptor']

df_yields['Y_X/E'] = df_fraction['yield on substrate']
df_yields['dor (e/C)'] = df_fraction['e/C']
df_yields['Y_X/ATP[tot]'] = df_fraction['Yatp']
df_yields['Y_X/ATP[cat]'] = 1000/df_fraction['ATP transferred']
df_yields['Y_ATP/E[cat]'] = df_fraction_new['Y_ATP/E[cat]']
df_yields['Y_NADPH/E[cat]'] = df_fraction_new['Y_NADPH/E[cat]']
df_yields['Y_NADH/E[cat]'] = df_fraction_new['Y_NADH/E[cat]']
df_yields['Y_FD/E[cat]'] = df_fraction_new['Y_FD/E[cat]']


df_yields['(Y_X/ATP[tot])^-1'] = 1/df_yields['Y_X/ATP[tot]']
df_yields['(Y_X/ATP[cat])^-1'] = 1/df_yields['Y_X/ATP[cat]']
df_yields['Y_ATP/E[ATP]'] = df_fraction_new['Y_ATP/E[atp]']
df_yields['Y_NADPH/E[NADPH]'] = df_fraction_new['Y_NADPH/E[nadph]']
df_yields['Y_NADH/E[NADH]'] = df_fraction_new['Y_NADH/E[nadh]']
df_yields['Y_FD/E[FD]'] = df_fraction_new['Y_FD/E[fd]']

result_path = 'C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\python files\generalized_scripts\Results'

with pd.ExcelWriter(result_path+"\\yields_table.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    df_yields.to_excel(writer)

#%%

plt.figure()
plt.plot(df_yields['Y_X/ATP[tot]'], df_yields['Y_X/ATP[cat]'], 'o')
plt.plot(range(16), range(16), '-')
plt.ylim(0, 100)
plt.xlim(0, 15)
plt.xlabel('Y_X/ATP[tot] (method)')
plt.ylabel('Y_X/ATP[cat] (Stouthamer)')
plt.show()

#%%
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D  # Import Line2D for custom legend entries
import pandas as pd

# Enable LaTeX rendering and set general font size
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
        'akg': 'alpha-ketoglutarate'
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

# Assuming `df_fraction` and `df_fraction_new` are your dataframes
cat_atp = df_fraction['ATP transferred']
ana_atp = df_fraction_new['ATP prod ana']

# Drop specified indices
cat_atp = cat_atp.drop(['Synechocystis\\photonbicarbonate_nadh', 'Gmetallireducens\\acetateFe3_test'])
ana_atp = ana_atp.drop(['Synechocystis\\photonbicarbonate_nadh', 'Gmetallireducens\\acetateFe3_test'])

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

# Plot with a broken x-axis
fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(28, 20), gridspec_kw={'width_ratios': [3, 1]})

# Define the x-axis limits for the two plots
ax1.set_xlim(0, 160)  # Upper part
ax2.set_xlim(560, sum_values.max() + 20)  # Lower part

# Plot the data on both axes
ax1.barh(cat_atp_sorted.index, cat_atp_sorted, label='$(Y_{X/ATP}^{cat})^{-1}$', color='tab:blue')
ax1.barh(ana_atp_sorted.index, ana_atp_sorted, left=cat_atp_sorted, label='$(Y_{X/ATP}^{ana})^{-1}$', color='tab:orange')

ax2.barh(cat_atp_sorted.index, cat_atp_sorted, label='$(Y_{X/ATP}^{cat})^{-1}$', color='tab:blue')
ax2.barh(ana_atp_sorted.index, ana_atp_sorted, left=cat_atp_sorted, label='$(Y_{X/ATP}^{ana})^{-1}$', color='tab:orange')

# Add a red dotted vertical line at x=34.7 on the upper plot
line = ax1.axvline(x=34.7, color='tab:blue', linestyle='--', linewidth=1.5)

# Break marks
ax1.spines['left'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.tick_params(labelright=False)  # Don't put tick labels on the right
ax1.yaxis.tick_left()

# Diagonal lines for the break
d = .015  # how big to make the diagonal lines in axes coordinates
kwargs = dict(color='k', clip_on=False)
# Lines for the axis break
# Upper right corner of ax1 and left top corner of ax2
fig.add_artist(plt.Line2D((1 - d, 1 + d), (-d, +d), transform=ax1.transAxes, **kwargs))
fig.add_artist(plt.Line2D((1 - d, 1 + d), (1 - d, 1 + d), transform=ax1.transAxes, **kwargs))

# Lower left corner of ax2 and right bottom corner of ax1
fig.add_artist(plt.Line2D((-d, +d), (1 - d, 1 + d), transform=ax2.transAxes, **kwargs))
fig.add_artist(plt.Line2D((-d, +d), (- d, + d), transform=ax2.transAxes, **kwargs))

# Rotate the y-axis labels by 0 degrees (default is vertical, no rotation needed)
plt.yticks(rotation=0, fontsize=32)  # Increase y-tick labels font size

# Coloring y-axis labels
for label in ax1.get_yticklabels():
    text = label.get_text()
    if 'C.~ljungdahlii' in text or 'bicarbonate' in text:
        label.set_color('green')
    elif 'Synechocystis' in text:
        label.set_color('purple')
    else:
        label.set_color('black')

# Adding labels and title
ax1.set_xlabel('$(Y_{X/ATP})^{-1}$  [mmolATP/gDW/h]', fontsize=28)
ax1.tick_params(axis='y', labelsize=24)

# Set the x-tick labels font size
ax1.tick_params(axis='x', labelsize=24)  # Set x-tick labels font size for ax1
ax2.tick_params(axis='x', labelsize=24)  # Set x-tick labels font size for ax2

# Create custom legend entries for organism classifications
chemoautotroph = Line2D([0], [0], color='green', marker='', linestyle='', label='chemoautotroph')
photoautotroph = Line2D([0], [0], color='purple', marker='', linestyle='', label='photoautotroph')
chemoheterotroph = Line2D([0], [0], color='black', marker='', linestyle='', label='chemoheterotroph')

# Create a custom legend entry for the red dotted line
custom_line = Line2D([0], [0], color='tab:blue', linestyle='--', linewidth=1.5, label='Stouthamer estimation of $(Y_{X/ATP}^{cat})^{-1}$')

# Add the legend, including the custom line and classifications
legend = ax2.legend(handles=[ax1.get_legend_handles_labels()[0][0], ax1.get_legend_handles_labels()[0][1], custom_line, chemoautotroph, photoautotroph, chemoheterotroph], loc='lower right', fontsize=32)

# Change the color of the text in the legend
for text in legend.get_texts():
    if 'chemoautotroph' in text.get_text():
        text.set_color('green')
    elif 'photoautotroph' in text.get_text():
        text.set_color('purple')
    elif 'chemoheterotroph' in text.get_text():
        text.set_color('black')


# Save the figure
figpath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\Paper latex files\Figures'
plt.savefig(figpath + r'\AllOrganisms_YatpCatAna.pdf', bbox_inches='tight')
plt.savefig(figpath + r'\AllOrganisms_YatpCatAna.png', bbox_inches='tight')

# Display the plot
plt.show()

