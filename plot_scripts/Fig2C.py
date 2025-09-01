# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 13:46:53 2025

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

stouthamer_df = pd.read_excel(r'Results\Stouthamer_table_data_pythonreadable.xlsx', index_col=0)

fig, ax = plt.subplots(figsize=(12, 6)) 
# Plotting the stacked bar plot with bars for each column
ax = stouthamer_df.T.plot(kind='barh', stacked=True, zorder=3, ax=ax)

# Adding title and labels
plt.xlabel('$(Y_{X/ATP})^{-1}$ [mmolATP/gDW]')


# Adding grid and x=0 line
plt.grid(True, zorder=0)
plt.axvline(x=0, color='black', linewidth=1)
plt.legend(bbox_to_anchor=(1.05, 1), ncol=2)
ax.set_yticklabels([r'$\it{S.\ cerevisiae}$', r'$\it{E.\ coli}$', r'$\it{E.\ coli}$ (Stouthamer method)', r'Stouthamer'])
# Display the plot
plt.title('Aerobic growth on glucose')

figpath = r'/Plots'
# plt.savefig(figpath+r'/Stouthamer_barplot.png', dpi=300, bbox_inches='tight')
plt.show()
