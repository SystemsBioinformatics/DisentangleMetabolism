# -*- coding: utf-8 -*-
"""
Created on Mon Aug 25 14:11:42 2025

@author: mre283
"""


import xml.etree.ElementTree as ET
import os
import re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as mcolors


def parse_metabolite_expression(expr: str) -> dict:
    expr = expr.strip().strip('[]')
    pattern = r'([+-]?\s*\d*\.?\d+(?:e[+-]?\d+)?)[*]\s*M_([a-zA-Z0-9]+?)_[a-z]+'
    matches = re.findall(pattern, expr)
    return {met_name: float(coeff_str.replace(' ', '')) for coeff_str, met_name in matches}


def is_dark_color(hex_color):
    rgb = mcolors.hex2color(hex_color)
    luminance = 0.299 * rgb[0] + 0.587 * rgb[1] + 0.114 * rgb[2]
    return luminance < 0.5


def map_values_to_colors(values_dict, cmap_name='Blues', vmin=None, vmax=None):
    values = list(values_dict.values())
    vmin = vmin if vmin is not None else min(values)
    vmax = vmax if vmax is not None else max(values)
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap_name)
    return {met: mcolors.to_hex(cmap(norm(val))) for met, val in values_dict.items()}


def update_style(style_str, new_fill):
    styles = dict(s.split(':') for s in style_str.split(';') if ':' in s)
    styles['fill'] = new_fill
    return ';'.join(f'{k}:{v}' for k, v in styles.items())


def apply_coloring_to_svg(svg_path, color_dict, output_path, text_color_decision=True):
    tree = ET.parse(svg_path)
    root = tree.getroot()
    ns = {'svg': 'http://www.w3.org/2000/svg'}

    for elem in root.iter():
        elem_id = elem.attrib.get('id')
        if elem_id and elem_id.startswith('r_'):
            met = elem_id[2:]
            if met in color_dict:
                new_fill = color_dict[met]

                if elem.tag.endswith('rect'):
                    original_style = elem.attrib.get('style', '')
                    updated_style = update_style(original_style, new_fill)
                    elem.set('style', updated_style)

                text_elem = root.find(f".//svg:text[@id='t_{met}']", ns)
                if text_elem is not None and text_color_decision:
                    text_color = '#ffffff' if is_dark_color(new_fill) else '#000000'
                    text_style = text_elem.attrib.get('style', '')
                    updated_text_style = update_style(text_style, text_color)
                    text_elem.set('style', updated_text_style)

    tree.write(output_path)


def save_colorbar(cmap_name, vmin, vmax, filename='colorbar_legend.svg', label='Stoichiometric coefficient'):
    fig, ax = plt.subplots(figsize=(1.5, 4))
    norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
    cb = plt.colorbar(cm.ScalarMappable(norm=norm, cmap=cm.get_cmap(cmap_name)), cax=ax)
    cb.set_label(label, rotation=90)
    plt.tight_layout()
    plt.savefig(filename, format='svg')
    plt.close()


# === CONFIGURATION ===
cwd = os.getcwd()
legend_folder = cwd+r'/Plots/Network_plots'
figpath = cwd+r'/Plots/Network_plots'
resultpath = cwd+r'/Results'


figure_configs = [
    {
        'svg': figpath+r'\networkdrawings_CO_Clungdahlii.svg',
        'data': resultpath + r'/Cljungdahlii/carbonmonoxide.xlsx',
        'carbon_mets': ['co2', '10fthf', '5mthf', 'accoa', 'ac'],
        'energy_mets': ['atp', 'nadh', 'nadph']
    },
    {
        'svg': figpath+r'\networkdrawings_acetate_Ecoli.svg',
        'data': resultpath+r'/Ecoli/acetate.xlsx',
        'carbon_mets': ['akg', 'accoa', 'oaa', 'succoa'],
        'energy_mets': ['atp', 'nadph', 'fadh2', 'nadh']
    },
    {
        'svg': figpath+r'\networkdrawings_acetate_Mbarkeri.svg',
        'data': resultpath + r'/Mbarkeri/acetate.xlsx',
        'carbon_mets': ['ac', 'accoa', 'co2', 'hco3', 'mh4spt'],
        'energy_mets': ['atp', 'fdred']
    },
    {
        'svg': figpath+r'\networkdrawings_acetate_Pputida.svg',
        'data': resultpath + r'/Pputida/acetate.xlsx',
        'carbon_mets': ['accoa', 'icit', 'akg', 'succoa', 'oaa'],
        'energy_mets': ['atp', 'nadh']
    },
    
    {
        'svg': figpath+r'\networkdrawings_glucose_Pputida.svg',
        'data': resultpath + r'/Pputida/glucose.xlsx',
        'carbon_mets': ['accoa', 'akg', 'succoa', 'oaa', '6pgc', 'g3p', '3pg', 'pep', 'pyr'],
        'energy_mets': ['atp', 'nadph', 'nadh']
    },
    {
        'svg': figpath+r'\networkdrawings_glucose_Ecoli.svg',
        'data': resultpath+r'/Ecoli/glucose.xlsx',
        'carbon_mets': ['f6p', 'ru5p__D', 'r5p', 'e4p', 'dhap', '3pg', 'pep', 'pyr', 'accoa', 'oaa', 'akg', 'succoa'],
        'energy_mets': ['atp', 'nadph', 'nadh']
    },
    {
        'svg': figpath+r'\networkdrawings_glutamine_Ecoli.svg',
        'data': resultpath+r'/Ecoli/glutamine.xlsx',
        'carbon_mets': ['glu__L', 'gln__L', '3pg', 'pep', 'pyr', 'accoa', 'oaa', 'succoa', 'ser__L', 'gly', 'mlthf', '10fthf', 'ala__D'],
        'energy_mets': ['atp', 'nadph', 'nadh', 'fadh2']
    },
]

carbon_cmap = 'Blues'
energy_cmap = 'Reds'
output_svg_suffix = '_colored.svg'

# === AGGREGATE VALUES ACROSS FIGURES ===
all_carbon_vals = {}
all_energy_vals = {}

for config in figure_configs:
    st_mceq = pd.read_excel(config['data'], sheet_name='info_BP', index_col=0).loc['mceq', 0]
    mceq_dict = parse_metabolite_expression(st_mceq)

    for met in config['carbon_mets']:
        if met in mceq_dict:
            all_carbon_vals[met+config['svg'][137:-4]] = mceq_dict[met]

    for met in config['energy_mets']:
        if met in mceq_dict:
            all_energy_vals[met+config['svg'][137:-4]] = mceq_dict[met]

# Generate color maps based on global range
carbon_colors = map_values_to_colors(all_carbon_vals, carbon_cmap)
energy_colors = map_values_to_colors(all_energy_vals, energy_cmap)
combined_colors_global = {**carbon_colors, **energy_colors}

# === APPLY COLORING TO EACH FIGURE ===
for config in figure_configs:
    st_mceq = pd.read_excel(config['data'], sheet_name='info_BP', index_col=0).loc['mceq', 0]
    mceq_dict = parse_metabolite_expression(st_mceq)

    # Filter relevant metabolites for this figure
    carbon_local = {met: mceq_dict[met] for met in config['carbon_mets'] if met in mceq_dict}
    energy_local = {met: mceq_dict[met] for met in config['energy_mets'] if met in mceq_dict}

    # Use globally normalized colors
    local_colors = {**{k: carbon_colors[k+config['svg'][137:-4]] for k in carbon_local},
                    **{k: energy_colors[k+config['svg'][137:-4]] for k in energy_local}}

    svg_path = config['svg']
    fig_dir, fig_name = os.path.split(svg_path)
    output_path = os.path.join(fig_dir, fig_name.replace('.svg', output_svg_suffix))

    apply_coloring_to_svg(svg_path, local_colors, output_path)

# === SAVE LEGENDS ===
os.makedirs(legend_folder, exist_ok=True)

if all_carbon_vals:
    save_colorbar(carbon_cmap, min(all_carbon_vals.values()), max(all_carbon_vals.values()),
                  filename=os.path.join(legend_folder, 'carbon_legend.svg'), label='Carbon metabolites')

if all_energy_vals:
    save_colorbar(energy_cmap, min(all_energy_vals.values()), max(all_energy_vals.values()),
                  filename=os.path.join(legend_folder, 'energy_legend.svg'), label='Energy cofactors')
