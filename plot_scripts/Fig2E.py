# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 13:50:54 2025

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
        'Ecoli': r'E.~coli',
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


directories = ['Results\Cljungdahlii', 'Results\Gmetallireducens', 'Results\Mbarkeri', 'Results\Scerevisiae', 'Results\Pputida', 'Results\Ecoli', 'Results\Synechocystis']

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

anaer_cat = cat_atp.loc['Ecoli\glucose_anaerobic']
anaer_ana = ana_atp.loc['Ecoli\glucose_anaerobic']

plt.figure()
plt.barh('$E.~coli$ glucose anaerobic', anaer_cat, color='tab:orange', label='eCAT')
plt.barh('$E.~coli$ glucose anaerobic', anaer_ana, color='tab:blue', label='pCAT and ANA', left=anaer_cat)

plt.xlabel('ATP requirement  [mmolATP/gDW]')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
plt.tight_layout()

# Save the figure
plt.savefig(r'Plots/examplebar_ecoli_anaerobic.png', bbox_inches='tight', dpi=300)
plt.show()
