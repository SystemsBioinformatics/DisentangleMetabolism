# -*- coding: utf-8 -*-
"""
Created on Tue Nov 12 11:46:15 2024

@author: mre283
"""

import cbmpy
import separate_cat_ana as sca
import numpy as np
import pandas as pd
import re
import sympy



def process_inputs(model, result_df, possible_energy_carriers):
    """
    Process input reactions and species to identify catabolic intermediates and external metabolites.
    
    Args:
        model: The metabolic model
        result_df: DataFrame containing reaction results, including 'catabolic' values
        energy_carriers: List of energy carrier metabolites
    
    Returns:
        catabolic_intermediates: List of catabolic intermediate metabolites
        catabolic_external: List of catabolic external metabolites
    """
    catabolic_intermediates = []
    catabolic_external = []
    energy_carriers = []

    for r in result_df.index:
        if np.abs(result_df.loc[r, 'catabolic']) > 1e-5:
            if 'demand' in r:
                continue

            species = model.getReaction(r).getSpeciesIds()

            for s in species:
                if s[-1] == 'e' and s not in catabolic_external:
                    catabolic_external.append(s)

            for s in species:
                if s not in catabolic_intermediates and s[-1] == 'c' and s not in possible_energy_carriers:
                    catabolic_intermediates.append(s)
                elif s not in energy_carriers and s in possible_energy_carriers:
                    energy_carriers.append(s)
    
    # for m in catabolic_external:
    #     if 'R_EX_'+m[2:]+'_rev' in result_df.index:
    #         s = m[:-1] + 'c'
            
    #         if s in catabolic_intermediates:
    #             print('removed from catabolic intermediates:')
    #             print(s)
    #             catabolic_intermediates.remove(s)
    
    if 'M_h_e' not in catabolic_external:
        catabolic_external.append('M_h_e')
    if 'M_h_c' not in catabolic_intermediates:
        catabolic_intermediates.append('M_h_c')


    return catabolic_intermediates, catabolic_external, energy_carriers


def create_species_and_reactions(model, catabolic_intermediates, energy_carriers):
    """
    Create external species and corresponding transport and exchange reactions.
    
    Args:
        model: The metabolic model
        catabolic_intermediates: List of catabolic intermediate metabolites
        energy_carriers: List of energy carrier metabolites
    """
    for s in catabolic_intermediates:
        if s in energy_carriers:
            continue

        name = model.getSpecies(s).name
        idd = s[:-1] + 'e'
        print(idd)
        # Create external species if it does not exist
        if model.getSpecies(idd) is None:
            model.createSpecies(idd, name=name + ' ex', compartment='e')
    
    ex_added = []
    for s in catabolic_intermediates:
        # Create transport reaction
        create_transport_reaction(model, s)

        # Create exchange reaction
        r_name = create_exchange_reaction(model, s)
        if r_name:
            ex_added.append(r_name)
    
    return ex_added

def create_transport_reaction(model, s):
    """
    Create a transport reaction for a given species.
    
    Args:
        model: The metabolic model
        s: The metabolite species ID
    """
    model.createReaction(f'R_{s[2:-2]}_trans', name=f'{s[2:-2]} transport', reversible=True)
    model.createReactionReagent(f'R_{s[2:-2]}_trans', s, -1)
    model.createReactionReagent(f'R_{s[2:-2]}_trans', s[:-1] + 'e', 1)
    model.setReactionBounds(f'R_{s[2:-2]}_trans', -1000, 1000)
    print(model.getReaction(f'R_{s[2:-2]}_trans').getEquation())


def create_exchange_reaction(model, s):
    """
    Create an exchange reaction for a given species.
    
    Args:
        model: The metabolic model
        s: The metabolite species ID
    """
    # r_ex_st = f'R_EX_{s[2:-2]}_e'
    # if r_ex_st not in model.getReactionIds() and r_ex_st+'_fwd' not in model.getReactionIds() and r_ex_st+'_rev' not in model.getReactionIds():
    #     model.createReaction(r_ex_st, name=f'{s[2:-2]} exchange', reversible=True)
    #     model.createReactionReagent(r_ex_st, s[:-1] + 'e', -1)
    #     model.setReactionBounds(r_ex_st, -1000, 1000)
    #     print(model.getReaction(r_ex_st).getEquation())
    
    
    name = model.getSpecies(s).name
    chemform = model.getSpecies(s).chemFormula
    print(chemform)
    # Use a regular expression to find the number after 'P'
    match = re.search(r'P(\d*)', chemform)
    print(match)
    if match:
        # Check if there's a number after P, otherwise default to 1
        number_after_P = int(match.group(1)) if match.group(1) else 1
        
    r_ex_st = 'R_EX_' + s[2:-2] + '_e'
    if r_ex_st not in model.getReactionIds() and r_ex_st + '_rev' not in model.getReactionIds() and r_ex_st+'_fwd' not in model.getReactionIds():
        model.createReaction('R_EX_' + s[2:-2] + '_e', name=s[2:-2]+' exchange', reversible=True)
        model.createReactionReagent('R_EX_' + s[2:-2] + '_e', s[:-1]+'e', -1)
        if 'CoA' in name and s != 'M_coa_c':
            model.createReactionReagent('R_EX_' + s[2:-2] + '_e', 'M_coa_e', 1)
        elif 'thf' in s and s != 'M_thf_c':
            model.createReactionReagent('R_EX_' + s[2:-2] + '_e', 'M_thf_e', 1)
        elif ('phosph' in name or 'Phosph' in name) and s!= 'M_pi_c':
            model.createReactionReagent('R_EX_' + s[2:-2] + '_e', 'M_pi_e', number_after_P)
        model.setReactionBounds('R_EX_' + s[2:-2] + '_e', -1000, 1000)
        print(model.getReaction('R_EX_' + s[2:-2] + '_e').getEquation())
    
    
        return 'R_EX_' + s[2:-2] + '_e'


def make_Ncat(result_df, model):
    model_Ncat = model.clone()
    
    demand_dict = {'R_ATP_demand': [('M_atp_c', -1.0), ('M_h2o_c', -1.0), ('M_adp_c', 1.0), ('M_h_c', 1.0), ('M_pi_c', 1.0)],
    'R_NADH_demand': [('M_nadh_c', -1.0), ('M_nad_c', 1.0), ('M_h_c', 1.0)],
    'R_NADPH_demand': [('M_nadph_c', -1.0), ('M_nadp_c', 1.0), ('M_h_c', 1.0)],
    'R_FADH2_demand': [('M_fadh2_c', -1.0), ('M_fad_c', 1.0), ('M_h_c', 2.0)],
    'R_Q8H2_demand': [('M_q8h2_c', -1.0), ('M_q8_c', 1.0)]}
    
    
    for name, reagents in demand_dict.items():
        # Check if any of the metabolites for the current energy carrier are in the energy_carriers list
        if name in result_df.index:
            # Create the reaction for the current energy carrier
            model_Ncat.createReaction(name, reversible=True, name=name[2:])

            # Add each reagent to the reaction
            for metabolite, coefficient in reagents:
                # if metabolite in energy_carriers:  # Check if the metabolite is in the energy_carriers list
                model_Ncat.createReactionReagent(name, metabolite=metabolite, coefficient=coefficient)
            
            # Set reaction bounds
            model_Ncat.setReactionBounds(name, -1000, 1000)

    
    cat_fluxes = result_df['catabolic']

    N = np.array(model_Ncat.buildStoichMatrix(matrix_type='numpy', only_return=True).array)
    metabs = model_Ncat.buildStoichMatrix(matrix_type='numpy', only_return=True).row
    i_inactive = []

    for i, f in enumerate(cat_fluxes.index):
        if np.abs(cat_fluxes.loc[f]) < 1e-6:
            i_inactive.append(i)

    N_activeflux = np.delete(N, i_inactive, 1)

    i_inactive = []

    # Loop over metabolites and store inactive indices
    for i in range(N_activeflux.shape[0]):
        if np.all(np.abs(N_activeflux[i, :])<1e-10):
            i_inactive.append(i)

    # Delete inactive metabolites from matrix with only active fluxes
    N_cat = np.delete(N_activeflux, i_inactive, 0)
    metabs_cat = np.delete(metabs, i_inactive)
    
    return N_cat, metabs_cat

def gaussian_elimination_with_operations(matrix):
    """
    Perform Gaussian elimination on a given matrix while tracking operations.
    
    Parameters:
        matrix (np.ndarray): The augmented matrix to be transformed.
    
    Returns:
        tuple: The row echelon form of the matrix and the operation matrix.
    """
    # Number of rows and columns
    n, m = matrix.shape
    # Create an identity matrix for tracking row operations
    operation_matrix = np.eye(n)

    for i in range(min(n, m)):  # Loop only over the smaller dimension
        # Partial pivoting: find the maximum element in the current column
        max_row_index = np.argmax(np.abs(matrix[i:, i])) + i
        
        # Swap the current row with the row containing the maximum element
        if i != max_row_index:
            matrix[[i, max_row_index]] = matrix[[max_row_index, i]]
            operation_matrix[[i, max_row_index]] = operation_matrix[[max_row_index, i]]
        
        # Make the diagonal element 1 and eliminate below
        for j in range(i + 1, n):
            if matrix[j, i] != 0:  # Avoid unnecessary calculations
                factor = matrix[j, i] / matrix[i, i]
                matrix[j] -= factor * matrix[i]
                
                # Track the operation in the operation matrix
                operation_matrix[j] -= factor * operation_matrix[i]

    return matrix, operation_matrix

# def print_matrix(matrix):
#     """ Helper function to print the matrix nicely. """
#     for row in matrix:
#         print(" ".join(f"{num:.4f}" for num in row))

def find_conserved_moieties(N, metabs):
    row_echelon_form, operations_matrix = gaussian_elimination_with_operations(N.copy())
    
    # Find rows with all zeros in row echelon form
    i_zero = np.all(np.abs(row_echelon_form) < 1e-5, axis=1)
    dictt = {}
    for idx, is_zero in enumerate(i_zero):  # Iterate over indices and boolean values
        if is_zero:  # Check if the row is all zeros
            print(f"Row {idx} is zero.")
            listt = []
            
            for n, k in enumerate(operations_matrix[idx]):  # Access the corresponding row in operations_matrix
                # print(k)
                # print(operations_matrix[idx])
                if np.abs(k) >1e-6:  # Check if k is non-zero
                    print(f"{k} * {metabs[int(n)]}")  # Use int(k) only if k is a valid index
                    listt.append((k.round(4), str(metabs[int(n)])))
                    
            dictt[f'Row {idx}'] = listt
    
    return dictt

    
def check_conserved_moieties(cons_dict, supplies):
    
    metabs=[]
    for name, reac in supplies.items():
        for tup in reac:
            if tup[0] not in metabs:
                metabs.append(tup[0])
    
    
    for row, moi in cons_dict.items():
        if (1.0, 'M_pi_c') in moi or (1.0, 'M_coa_c') in moi or (1.0, 'M_thf_c') in moi:
            continue
        for m in moi:
            
            
            if m[1] not in metabs:
                raise Exception(f'There is a conserved moiety that is not added! {moi}')



def add_supply_reactions(model, energy_carriers, result_df):
    """
    Add supply reactions for energy carriers that are in the energy_carriers list.
    
    Args:
        model: The metabolic model
        energy_carriers: List of energy carrier metabolites
    """
    # Supply reactions defined with metabolite IDs and coefficients
    supplies = {
        'ATP': [('M_atp_c', 1.0), ('M_h2o_c', 1.0), ('M_adp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)],
        'FDXR': [('M_fdxo_42_c', -1.0), ('M_fdxr_42_c', 1.0)],
        'NADH': [('M_nadh_c', 1.0), ('M_nad_c', -1.0), ('M_h_c', -1.0)],
        'NADPH': [('M_nadph_c', 1.0), ('M_nadp_c', -1.0), ('M_h_c', -1.0)],
        'FADH': [('M_fadh2_c', 1.0), ('M_fad_c', -1.0), ('M_h_c', -2.0)],
        'Q8H2': [('M_q8h2_c', 1.0), ('M_q8_c', -1.0)],
        'ATP_m': [('M_atp_m', 1.0), ('M_h2o_m', 1.0), ('M_adp_m', -1.0), ('M_h_m', -1.0), ('M_pi_m', -1.0)],
        'NADH_m': [('M_nadh_m', 1.0), ('M_nad_m', -1.0), ('M_h_m', -1.0)],
        'NADPH_m': [('M_nadph_m', 1.0), ('M_nadp_m', -1.0), ('M_h_m', -1.0)],
        'F420R': [('M_f420_DASH_2h2_c', 1.0), ('M_f420_DASH_2_c', -1.0)],
        'FDRED':[('M_fdred_c', 1.0), ('M_fdox_c', -1.0)],
        'MPHEN': [('M_mphenh2_c', 1.0), ('M_mphen_c', -1.0)],
        'OMCHR': [('M_omchr_e', 1.0), ('M_omcho_e', -1.0)],
        'MQN8': [('M_mqn8_c', 1.0), ('M_mql8_c', -1.0)],
        'FLXR': [('M_flxr_c', 1.0), ('M_flxso_c', -1.0)],
        'FICYTC': [('M_ficytC_c', 1.0), ('M_focytC_c', -1.0)],
        'FICYTC_m':[('M_focytc_m', -1.0), ('M_ficytc_m', 1.0)],
        'Q6H2_m':[('M_q6h2_m', 1.0), ('M_q6_m', -1.0)],
        'FADH2_m':[('M_fadh2_m', 1.0), ('M_fad_m', -1.0)],
        'DHLAM_m':[('M_dhlam_m', 1.0), ('M_lpam_m', -1.0)], 
        'GTP': [('M_gtp_c', 1.0), ('M_h2o_c', 1.0), ('M_gdp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)],
        'UTP': [('M_utp_c', 1.0), ('M_h2o_c', 1.0), ('M_udp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)]
    }
    added_supplies = []
    
    N_cat, metabs_cat = make_Ncat(result_df, model)
    cons_dict = find_conserved_moieties(N_cat, metabs_cat)
    check_conserved_moieties(cons_dict, supplies)
    # Iterate over the supplies and add reactions only if the energy carrier is in the list
    for name, reagents in supplies.items():
        # Check if any of the metabolites for the current energy carrier are in the energy_carriers list
        if any(metabolite in energy_carriers for metabolite, _ in reagents):
            # Create the reaction for the current energy carrier
            model.createReaction(f'R_{name}_supply', reversible=True, name=f'{name} supply')

            # Add each reagent to the reaction
            for metabolite, coefficient in reagents:
                # if metabolite in energy_carriers:  # Check if the metabolite is in the energy_carriers list
                model.createReactionReagent(f'R_{name}_supply', metabolite=metabolite, coefficient=coefficient)
            
            # Set reaction bounds
            model.setReactionBounds(f'R_{name}_supply', -1000, 1000)
            added_supplies.append(f'R_{name}_supply')
    
    return added_supplies


def adjust_reaction_bounds(model, catabolic_external, catabolic_intermediates, energy_carriers):
    """
    Adjust reaction bounds based on given external metabolites.
    
    Args:
        model: The metabolic model
        catabolic_external: List of catabolic external metabolites
    """
    
    # catabolic_external.append('M_h_e')
    # catabolic_intermediates.remove('M_h_c')
    for m in catabolic_external:
        rid = 'R_EX' + m[1:]
        
        if rid in model.getReactionIds():
            model.setReactionBounds(rid, -1000, 1000)
            model.getReaction(rid).reversible=True
        elif rid+'_fwd' in model.getReactionIds():
            model.setReactionBounds(rid+'_fwd', -1000, 1000)
            model.getReaction(rid+'_fwd').reversible=True
        elif rid+'_rev' in model.getReactionIds():
            model.setReactionBounds(rid+'_rev', -1000, 1000)
            model.getReaction(rid+'_rev').reversible=True
        else:
            raise Exception('Missing exchange reaction?')

    kos = []
    for r in result_df.index:
        if 'demand' in r or 'ATPS' in r:
            continue
        catabolic_reac = False
        if np.abs(result_df.loc[r, 'catabolic']) > 1e-5:
            catabolic_reac = True
            species = model.getReaction(r).getSpeciesIds()
            for s in species:
                if s not in catabolic_intermediates and s not in energy_carriers:
                    catabolic_reac = False
                    # print(s)
                    # print(r)
            
        if catabolic_reac is True:
            print(r)
            model.setReactionBounds(r, 0, 0)
            kos.append(r)
        if 'R_Htex_fwd' in model.getReactionIds():
            model.setReactionBounds('R_Htex_fwd', 0, 0)
            kos.append('R_Htex_fwd')
        if 'R_ACtex_rev' in model.getReactionIds():
            model.setReactionBounds('R_ACtex_rev', 0, 0)
            kos.append('R_ACtex_rev')
            
    return kos
    

def get_MCEQ_precursors(model, obj, fluxes, added_supplies=['R_ATP_supply'], biomass_reac='R_EX_BIOMASS', energy_carriers=['M_atp_c', 'M_adp_c']):
    
    mceq_model = model.clone()
    fluxes_mceq = fluxes.copy()
    for r in added_supplies:
        print(r)
        mceq_model.deleteReactionAndBounds(r+'_rev')
        mceq_model.deleteReactionAndBounds(r+'_fwd')
        fluxes_mceq.pop(r+'_rev')
        fluxes_mceq.pop(r+'_fwd')
    
    mceq = obj.get_MCEQ(mceq_model, biomass_reac=biomass_reac, flux_dist= pd.Series(fluxes_mceq), 
                         additional_active_metabolites=energy_carriers)
    
    
    return mceq

# def merge_model(model):
#     merged_model = model.clone()
    
#     reactions = merged_model.getReactionIds()
    
#     done = []
#     for r in reactions:
#         if r[-4:] == '_fwd' and r not in done:
#             if r[:-4]+'_rev' in reactions:
#                 # merge reactions
#                 merged_model.setReactionBounds(r, -1000, 1000)
#                 merged_model.deleteReactionsAndBounds(r[:-4]+'_rev')
#                 done.append(r)
#                 done.append(r[:-4]+'_rev')
#         if r[-4:] == '_rev' and r not in done:
#             if r[:-4]+'_fwd' in reactions:
#                 # merge reactions
#                 merged_model.setReactionBounds(r, -1000, 1000)
#                 merged_model.deleteReactionsAndBounds(r[:-4]+'_fwd')
#                 done.append(r)
#                 done.append(r[:-4]+'_fwd')   
        
    
    

def main(model, result_df, biomass_reac='R_EX_BIOMASS'):
    """
    Main function to process the model and update it with reactions and species.
    
    Args:
        model: The metabolic model
        result_df: DataFrame containing reaction results, including 'catabolic' values
    """
    possible_energy_carriers = [
        'M_nad_c', 'M_nadh_c', 'M_nadph_c', 'M_nadp_c', 
        'M_atp_c', 'M_adp_c', 'M_fdxr_42_c', 'M_fdxo_42_c', 
        'M_fad_c', 'M_fadh2_c', 'M_q8_c', 'M_q8h2_c',
        'M_atp_m', 'M_adp_m', 'M_nad_m', 'M_nadh_m',
        'M_nadph_m', 'M_nadp_m', 'M_f420_DASH_2h2_c',
        'M_f420_DASH_2_c', 'M_gdp_c', 'M_gtp_c', 'M_ttp_c',
        'M_tdp_c', 'M_utp_c', 'M_udp_c'
    ]

    # Step 1: Process inputs to identify species
    catabolic_intermediates, catabolic_external, energy_carriers = process_inputs(model, result_df, possible_energy_carriers)
    print(catabolic_external)
    print(catabolic_intermediates)
    print(energy_carriers)

    # Step 2: Create species and reactions
    ex_added = create_species_and_reactions(model, catabolic_intermediates, energy_carriers)

    # Step 3: Add supply reactions
    added_supplies = add_supply_reactions(model, energy_carriers, result_df)
    print(added_supplies)
    # Step 4: Adjust reaction bounds
    kos = adjust_reaction_bounds(model, catabolic_external, catabolic_intermediates, energy_carriers)
    
    
    # Set the objective for optimization
    model.setObjectiveFlux(biomass_reac, osense='minimize')
    # model.setObjectiveFlux('R_THD2pp', osense='minimize')
    
    model_notsplit = model.clone()
    # model.setReactionBounds(substrate, 0, 0)
    # Analyze the model using CBGLPK
    model = cbmpy.CBTools.splitReversibleReactions(model)
    
    # cbmpy.CBWrite.writeSBML3FBCV2(model, 'biomassprod_model.xml')
    # model.setReactionBounds('R_EX_oaa_e_fwd', 0, 0)
    result = cbmpy.CBGLPK.glpk_analyzeModel(model, method='s')
    fluxes = model.getReactionValues()
    mceq = get_MCEQ_precursors(model, gluc, fluxes, added_supplies=added_supplies, biomass_reac = biomass_reac, energy_carriers=energy_carriers)
    # print(model_notsplit.getReactionIds())
    print(ex_added+added_supplies)
    
    result = cbmpy.CBGLPK.glpk_analyzeModel(model_notsplit, method='s')
    fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model_notsplit, selected_reactions=ex_added+added_supplies)

    df_fva = pd.DataFrame(fva[0], index=fva[1], 
                          columns=['Reaction', 'Reduced Costs', 'Variability Min', 'Variability Max', 'abs(Max-Min)', 'MinStatus', 'MaxStatus'])

    
    return fluxes, mceq, df_fva, kos


# file_list = ['readfile_aerobic_acetate_Ecoli.xlsx', 'readfile_aerobic_fructose_Ecoli.xlsx','readfile_aerobic_galactose_Ecoli.xlsx','readfile_aerobic_gluconate_Ecoli.xlsx',
#               'readfile_aerobic_glucose_Ecoli.xlsx','readfile_aerobic_glycerol_Ecoli.xlsx','readfile_aerobic_pyruvate_Ecoli.xlsx','readfile_aerobic_ribose_Ecoli.xlsx',
#               'readfile_aerobic_sorbitol_Ecoli.xlsx']
# file_list = ['readfile_aerobic_succinate_Ecoli.xlsx','readfile_aerobic_xylose_Ecoli.xlsx','readfile_alanine_Ecoli.xlsx',
#               'readfile_alphaketoglutarate_Ecoli.xlsx','readfile_anaerobic_glucose_Ecoli.xlsx','readfile_arginine_Ecoli.xlsx','readfile_aspartate_Ecoli.xlsx',
#               'readfile_citrate_Ecoli.xlsx','readfile_fumarate_Ecoli.xlsx','readfile_glutamate_Ecoli.xlsx','readfile_glutamine_Ecoli.xlsx','readfile_lactate_Ecoli.xlsx',
#               'readfile_malate_Ecoli.xlsx','readfile_tryptophan_Ecoli.xlsx']

# file_list=['readfile_aerobic_acetate_Ecoli.xlsx', 'readfile_aerobic_pyruvate_Ecoli.xlsx', 'readfile_aerobic_succinate_Ecoli.xlsx', 'readfile_alphaketoglutarate_Ecoli.xlsx',
#            'readfile_arginine_Ecoli.xlsx','readfile_aspartate_Ecoli.xlsx','readfile_citrate_Ecoli.xlsx','readfile_fumarate_Ecoli.xlsx','readfile_glutamate_Ecoli.xlsx','readfile_glutamine_Ecoli.xlsx','readfile_lactate_Ecoli.xlsx',
#            'readfile_malate_Ecoli.xlsx','readfile_tryptophan_Ecoli.xlsx']

file_list = ['readfile_aerobic_glucose_Ecoli.xlsx']

df_dict_anabolic_mceq = {}
result_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\Precursor_simulations_Ecoli'


# file_list = [file_list[0]]
for file in file_list:
    gluc =  sca.separate_cat_ana(file, all_cofactors=False)
    
    gluc.do_cat_ana_split()
    
    dfs = gluc.generate_dfs()
    
    result_df, atp_df, cofactor_df = dfs
    

    model = gluc.active_model.clone()
    N_cat, metabs_cat = make_Ncat(result_df, model)
    cons_dict = find_conserved_moieties(N_cat, metabs_cat)
    
    fluxes, mceq, df_fva, kos = main(model, result_df, biomass_reac='R_EX_BIOMASS')
    
    energy_obj = sca.energy_generation(model)
    atp_consumption, atp_production, atp_special_consumption = energy_obj.calculate_energy(fluxes, split_names=True)
    nadh_prod, nadh_cons, transhydrogenase_prod, transhydrogenase_cons = energy_obj.calculate_NADH_prod(fluxes, split_names=True)
    nadph_prod, nadph_cons, transhydrogenase_prod, transhydrogenase_cons = energy_obj.calculate_NADPH_prod(fluxes, split_names=True)
    
    info_dict = {}
    info_dict['mceq'] = mceq
    info_dict['atp production'] = atp_production
    info_dict['atp consumption'] = atp_consumption
    info_dict['nadh production'] = nadh_prod
    info_dict['nadh consumption'] = nadh_cons
    info_dict['nadph production'] = nadph_prod
    info_dict['nadph consumption'] = nadph_cons
    
    with pd.ExcelWriter(result_path+"\\"+file[9:-11]+".xlsx") as writer:
       
        # use to_excel function and specify the sheet_name and index 
        # to store the dataframe in specified sheet
        pd.Series(fluxes).to_excel(writer, sheet_name="fluxes")
        pd.Series(info_dict).to_excel(writer, sheet_name="info")
        df_fva.to_excel(writer, sheet_name="fva")

        
        
    df_dict_anabolic_mceq[file[9:-11]] = str(mceq)




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



df_mceqs = pd.read_excel(result_path+'\\all_mceqs.xlsx', header=0, index_col=0)
for sim, mceq in df_dict_anabolic_mceq.items():
    if mceq[0] == '[':
        mceq = mceq[1:-1]
    df_mceqs.loc[sim, :] = np.zeros(len(df_mceqs.columns))
    symp = sympy.parsing.sympy_parser.parse_expr(mceq)
    metabolite_dict = extract_metabolites(symp)
    for metab, coeff in metabolite_dict.items():
        df_mceqs.loc[sim, metab] = coeff
    
    
df_mceqs.fillna(0, inplace=True)

# df_mceqs.to_excel(result_path+'\\all_mceqs.xlsx')

#%%
# result_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\Precursor_simulations_Ecoli'


# df_mceqs = pd.read_excel(result_path+'\\all_mceqs.xlsx', index_col=0)

# df_mceqs_plot = df_mceqs.drop(['M_biomass_e', 'M_btn_e', 'M_ca2_e', 'M_cl_e', 'M_cobalt2_e', 'M_cu2_e', 'M_fe2_e', 'M_h2o_e',
#                                'M_h_e', 'M_k_e', 'M_kdo2lipid4_e', 'M_mg2_e', 'M_mn2_e', 'M_mobd_e', 'M_nh4_e', 'M_ni2_e',
#                                'M_so4_e', 'M_zn2_e', 'M_adp_c', 'M_fad_c', 'M_nad_c', 'M_nadp_c', 'M_q8_c', 'M_gdp_e', 'M_udp_e',
#                                'M_fe3_e', 'M_gtp_e', 'M_utp_e', 'M_3mob_e', 'M_fgam_e'], axis=1)


# #%%
# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors

# # Dictionary with groups and corresponding metabolites
# groups = {
#     'upper glycolysis': ['fructose', 'glucose', 'galactose', 'gluconate', 'xylose', 'ribose', 'sorbitol', 'anaerobic_glucose'],
#     'TCA-cycle intermediates': ['acetate', 'succinate', 'alphaketoglutarate', 'citrate', 'fumarate', 'malate'],
#     'lower glycolysis': ['pyruvate', 'glycerol', 'lactate'],
#     'amino acids': ['alanine', 'arginine', 'aspartate', 'glutamate', 'glutamine', 'tryptophan']
# }

# # Function to categorize each index based on the groups dictionary
# def categorize_index(index):
#     for group, metabolites in groups.items():
#         if any(metabolite in index for metabolite in metabolites):
#             return group
#     return 'Other'  # Default category if no match is found

# # Assuming df_mceqs_plot is your pre-defined DataFrame
# # Columns to move to the right-hand subplot
# columns_to_move = ['M_atp_c', 'M_nadh_c', 'M_nadph_c', 'M_fadh2_c', 'M_q8h2_c']

# # Check if all columns exist in the DataFrame
# existing_columns = [col for col in columns_to_move if col in df_mceqs_plot.columns]

# # Create two DataFrames: one for the main heatmap and one for the specified columns
# df_main = df_mceqs_plot.drop(columns=existing_columns)  # Remove specified columns
# df_right = df_mceqs_plot[existing_columns]  # Keep only the specified columns

# # Categorize the indices of df_main
# df_main['Category'] = df_main.index.map(categorize_index)

# # Sort the DataFrame by Category
# df_main = df_main.sort_values(by='Category')

# # Extract the sorted index and categories for labeling
# sorted_index = df_main.index
# sorted_categories = df_main['Category'].values

# # Remove the 'Category' column before plotting
# df_main = df_main.drop(columns='Category')

# # Categorize the indices of df_main
# df_right['Category'] = df_right.index.map(categorize_index)

# # Sort the DataFrame by Category
# df_right = df_right.sort_values(by='Category')

# # Extract the sorted index and categories for labeling
# sorted_index = df_right.index
# sorted_categories = df_right['Category'].values

# # Remove the 'Category' column before plotting
# df_right = df_right.drop(columns='Category')


# # Plotting the heatmaps
# fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 8))

# # Main heatmap (excluding the specified columns)
# cmap_main = plt.get_cmap('seismic')
# norm_main = mcolors.TwoSlopeNorm(vmin=df_main.values.min(), vcenter=0, vmax=df_main.values.max())
# cax1 = ax1.matshow(df_main, cmap=cmap_main, norm=norm_main)
# ax1.set_title("Carbon precursors")

# # Adding a color bar for the main heatmap
# plt.colorbar(cax1, ax=ax1)

# # Customize y-axis labels with grouped categories and their metabolites
# yticks = [f"{category}: {index}" for category, index in zip(sorted_categories, sorted_index)]
# xticks_right = [idd[2:-2] for idd in df_right.columns]
# xticks_main = [idd[2:-2] for idd in df_main.columns]

# # Set axis ticks for the main heatmap
# ax1.set_xticks(range(len(df_main.columns)))
# ax1.set_yticks(range(len(df_main.index)))
# ax1.set_xticklabels(xticks_main, rotation=90)
# ax1.set_yticklabels(yticks)


# # Subplot for the specified columns
# cmap_right = plt.get_cmap('seismic')  # Same colormap but separate normalization
# norm_right = mcolors.TwoSlopeNorm(vmin=df_right.values.min(), vcenter=0, vmax=df_right.values.max())
# cax2 = ax2.matshow(df_right, cmap=cmap_right, norm=norm_right)
# ax2.set_title("Energy carriers")

# # Adding a color bar for the specified columns heatmap
# plt.colorbar(cax2, ax=ax2)

# # Set axis ticks for the specified columns heatmap
# ax2.set_xticks(range(len(df_right.columns)))
# ax2.set_yticks(range(len(df_right.index)))
# ax2.set_xticklabels(xticks_right, rotation=90)
# ax2.set_yticklabels(yticks)

# # Show the plot
# plt.tight_layout()
# plt.show()

# #%%

# import pandas as pd
# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.colors as mcolors

# # Dictionary with groups, corresponding metabolites, and their colors
# groups = {
#     'upper glycolysis': {
#         'metabolites': ['fructose', 'glucose', 'galactose', 'gluconate', 'xylose', 'ribose', 'sorbitol', 'anaerobic_glucose'],
#         'color': 'blue'
#     },
#     'TCA-cycle intermediates': {
#         'metabolites': ['acetate', 'succinate', 'alphaketoglutarate', 'citrate', 'fumarate', 'malate'],
#         'color': 'green'
#     },
#     'lower glycolysis': {
#         'metabolites': ['pyruvate', 'glycerol', 'lactate'],
#         'color': 'red'
#     },
#     'amino acids': {
#         'metabolites': ['alanine', 'arginine', 'aspartate', 'glutamate', 'glutamine', 'tryptophan'],
#         'color': 'purple'
#     }
# }


# # BiGG ID to full names/formulas mapping
# bigg_to_name = {
#     'ac': 'Acetate',
#     'accoa': 'Acetyl-CoA',
#     'akg': 'Alpha-Ketoglutarate',
#     'co2': r'CO$_2$',  # Latex formatted
#     'coa': 'CoA',
#     'fum': 'Fumarate',
#     'o2': r'O$_2$',
#     'oaa': 'Oxaloacetate',
#     'pi': 'Phosphate',
#     'succ': 'Succinate',
#     'succoa': 'Succinyl-CoA',
#     '3pg': '3P-glycerate',
#     'dhap': 'Dihydroxyacetone-P',
#     'e4p': 'Erythrose 4P',
#     'f6p': 'Fructose 6P',
#     'g3p': 'Glyceraldehyde 3P',
#     'gtp': 'Guanosine triphosphate',
#     'pep': 'Phosphoenolpyruvate',
#     'pyr': 'Pyruvate',
#     'r5p': 'Ribose 5-phosphate',
#     'ru5p__D': 'D-Ribulose 5P',
#     'dha': 'Dihydroxyacetone',
#     'glyc': 'Glycerol',
#     'icit': 'Isocitrate',
#     'utp': 'Uridine triphosphate',
#     '3mob': '3-Methyl-2-oxobutanoate',
#     'ala__L': 'L-Alanine',
#     'glu__L': 'L-Glutamate',
#     'gly': 'Glycine',
#     'val__L': 'L-Valine',
#     'for': 'Formate',
#     'arg__L': 'L-Arginine',
#     'asp__L': 'L-Aspartate',
#     '10fthf': '10-Formyl-THF',
#     'fgam': 'Formylglycinamide ribonucleotide',
#     'mlthf': '5,10-Methylene-THF',
#     'ser__L': 'L-Serine',
#     'thf': 'THF',
#     'ala__D': 'D-Alanine',
#     'gln__L': 'L-Glutamine',
#     'trp__L': 'L-Tryptophan'
# }

# groups_y = {
#     'upper glycolysis': {
#         'metabolites': ['f6p'],
#         'color': 'blue'
#     },
#     'TCA-cycle intermediates': {
#         'metabolites': ['ac', 'accoa', 'akg', 'fum', 'oaa', 'succ', 'succoa', 'icit'],
#         'color': 'green'
#     },
#     'lower glycolysis': {
#         'metabolites': ['3pg', 'dhap', 'g3p', 'pep', 'pyr', 'dha', 'glyc'],
#         'color': 'red'
#     },
#     'amino acids': {
#         'metabolites': ['ala__L', 'glu__L', 'gly', 'val__L', 'arg__L', 'asp__L', 'ser__L', 'ala__D', 'gln__L', 'trp__L'],
#         'color': 'purple'
#     },
    
#     'conserved moieties': {
#         'metabolites': ['thf', 'coa', 'pi'],
#         'color': 'orange'
#     },
#     'ppp': {
#         'metabolites': ['e4p', 'r5p', 'ru5p__D'],
#         'color': 'brown'
#     }, 
#     'external': {
#         'metabolites': ['co2', 'o2', 'for'],
#         'color': 'pink'
#     },
#     'other': {
#         'metabolites': ['10fthf', 'fgam', 'mlthf'],
#         'color': 'grey'
#     }
# }


# # Function to categorize each index based on the groups dictionary
# def categorize_index(index):
#     for group, data in groups.items():
#         if any(metabolite in index for metabolite in data['metabolites']):
#             return group
#     return 'Other'

# # Function to get text color based on category
# def get_text_color(category):
#     return groups.get(category, {}).get('color', 'black')  # Default to black if no color is defined

# # Assuming df_mceqs_plot is your pre-defined DataFrame
# # Columns to move to the right-hand subplot
# columns_to_move = ['M_atp_c', 'M_nadh_c', 'M_nadph_c', 'M_fadh2_c', 'M_q8h2_c']

# # Check if all columns exist in the DataFrame
# existing_columns = [col for col in columns_to_move if col in df_mceqs_plot.columns]

# # Create two DataFrames: one for the main heatmap and one for the specified columns
# df_main = df_mceqs_plot.drop(columns=existing_columns)  # Remove specified columns
# df_right = df_mceqs_plot[existing_columns]  # Keep only the specified columns

# # Categorize the indices of df_main
# df_main['Category'] = df_main.index.map(categorize_index)
# df_main = df_main.sort_values(by='Category')
# sorted_index_main = df_main.index
# sorted_categories_main = df_main['Category'].values
# df_main = df_main.drop(columns='Category')

# # Define grid-drawing function
# def add_grid(ax, df):
#     n_rows, n_cols = df.shape
#     ax.hlines(np.arange(-0.5, n_rows), xmin=-0.5, xmax=n_cols - 0.5, color='black', linewidth=0.5)
#     ax.vlines(np.arange(-0.5, n_cols), ymin=-0.5, ymax=n_rows - 0.5, color='black', linewidth=0.5)

# # Function to remove 'aerobic_' prefix from y-tick labels
# def clean_label(label):
#     return label[8:] if label.startswith('aerobic_') else label


# # Customize y-axis labels with cleaned metabolite indices and their colors
# yticks_main = [clean_label(index) for index in sorted_index_main]  # Clean the labels
# xticks_main = [idd[2:-2] for idd in df_main.columns]


# # Define a function to map x-axis labels to their categories
# def categorize_x_label(label):
#     for group, data in groups_y.items():
#         if any(metabolite in label for metabolite in data['metabolites']):
#             return group
#     return 'Other'

# # Define a function to get the color for the x-axis labels
# def get_x_label_color(label):
#     category = categorize_x_label(label)
#     return groups_y.get(category, {}).get('color', 'black')  # Default to black

# # Create a mapping of x-axis labels to their categories
# x_label_categories = {label: categorize_x_label(label[2:-2]) for label in df_main.columns}

# # Sort the columns based on the category order in groups_y
# category_order = list(groups_y.keys()) + ['Other']  # Ensure 'Other' is sorted last
# sorted_columns = sorted(df_main.columns, key=lambda x: category_order.index(x_label_categories[x]))

# # Reorder the DataFrame and xticks_main accordingly
# df_main = df_main[sorted_columns]
# xticks_main = [bigg_to_name.get(col[2:-2], col[2:-2]) for col in sorted_columns]

# # Update the plot code to reflect the sorted columns
# fig, ax1 = plt.subplots(figsize=(20, 12))
# cmap_main = plt.get_cmap('seismic')
# norm_main = mcolors.TwoSlopeNorm(vmin=df_main.values.min(), vcenter=0, vmax=df_main.values.max())
# cax1 = ax1.matshow(df_main, cmap=cmap_main, norm=norm_main)
# ax1.set_title("Anabolic carbon precursor requirement", fontsize=16)
# colorbar = plt.colorbar(cax1, ax=ax1)
# colorbar.set_label('yield [mmol/gDW]', fontsize=14, labelpad=10, loc='center')  # Add label above


# # Update the x-axis ticks and labels
# ax1.set_xticks(range(len(df_main.columns)))
# ax1.set_xticklabels([])  # Start with empty labels to add colored text manually
# ax1.set_yticklabels([])  # Start with empty labels to add colored text manually

# # Add x-axis labels with their corresponding colors
# for i, label in enumerate(xticks_main):
#     color = get_x_label_color(sorted_columns[i][2:-2])  # Get the color based on the label
#     ax1.text(i, -1.2, label, color=color, fontsize=12, ha='center', va='bottom', rotation=90)

# # Add y-axis labels as before
# yticks_main = [clean_label(index) for index in sorted_index_main]
# for i, (tick, category) in enumerate(zip(yticks_main, sorted_categories_main)):
#     ax1.text(-1.2, i, tick, color=get_text_color(category), fontsize=12, ha='right', va='center')

# # Add grid to the main heatmap
# add_grid(ax1, df_main)

# # Add legend for text colors
# legend_handles = [plt.Line2D([0], [0], color=data['color'], marker='o', markersize=10, linestyle='None', label=group)
#                   for group, data in groups_y.items()]
# ax1.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.25), ncol=2, frameon=False, fontsize=12)

# # Move x-axis label to the top
# ax1.set_xlabel('Metabolite', fontsize=14, labelpad=150)
# ax1.xaxis.set_label_position('top')

# # Move y-axis label to the left
# ax1.set_ylabel('Carbon and Energy Source', fontsize=14, labelpad=125)
# ax1.yaxis.set_label_position('left')

# plt.tight_layout()

# plt.savefig('heatmap_Ecoli_carbonsubs.png', dpi=300, bbox_inches='tight')
# plt.show()


# #%%
# # Categorize the indices of df_right
# df_right['Category'] = df_right.index.map(categorize_index)
# df_right = df_right.sort_values(by='Category')  # Sort by category
# sorted_index_right = df_right.index
# sorted_categories_right = df_right['Category'].values
# df_right = df_right.drop(columns='Category')  # Remove the temporary 'Category' column

# # Function to remove 'aerobic_' prefix from y-tick labels
# def clean_label(label):
#     return label[8:] if label.startswith('aerobic_') else label

# # Right-hand heatmap plot
# fig2, ax2 = plt.subplots(figsize=(6, 8))
# cmap_right = plt.get_cmap('seismic')
# norm_right = mcolors.TwoSlopeNorm(vmin=df_right.values.min(), vcenter=0, vmax=df_right.values.max())
# cax2 = ax2.matshow(df_right, cmap=cmap_right, norm=norm_right)
# ax2.set_title("Energy carriers")
# plt.colorbar(cax2, ax=ax2)

# # Customize y-axis labels with cleaned metabolite indices and their colors
# yticks_right = [clean_label(index) for index in sorted_index_right]  # Clean the labels
# xticks_right = [idd[2:-2] for idd in df_right.columns]

# # Set x-axis labels
# ax2.set_xticks(range(len(df_right.columns)))
# ax2.set_xticklabels(xticks_right, rotation=90)

# # Disable default y-axis labels
# ax2.set_yticks(range(len(df_right.index)))
# ax2.set_yticklabels([])  # Remove default y-axis tick labels

# # Apply custom text color for y-tick labels, shifted to the left
# for i, (tick, category) in enumerate(zip(yticks_right, sorted_categories_right)):
#     ax2.text(-1.2, i, tick, color=get_text_color(category), fontsize=10, ha='right', va='center')  # Adjusted x-coordinate

# # Add grid to the right-hand heatmap
# add_grid(ax2, df_right)

# # Add legend for text colors
# legend_handles = [plt.Line2D([0], [0], color=data['color'], marker='o', markersize=10, linestyle='None', label=group)
#                   for group, data in groups.items()]
# ax2.legend(handles=legend_handles, loc='lower center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False)

# # Move x-axis label to the top
# ax2.set_xlabel('Consumed or produced in anabolism from precursor [mmol/gDW]')
# ax2.xaxis.set_label_position('top')  # Position the x-label at the top

# # Move y-axis label to the left
# ax2.set_ylabel('Carbon and energy source', labelpad=100)  # Increase the labelpad to move it further left
# ax2.yaxis.set_label_position('left')  # Position the y-label to the left

# # Show right-hand plot
# plt.tight_layout()
# plt.show()

# #%%

# # Transpose the DataFrame for the main heatmap
# df_main = df_main.transpose()

# # Reassign categories and sort the new rows (previously columns)
# df_main['Category'] = df_main.index.map(lambda x: categorize_x_label(x[2:-2]))
# df_main = df_main.sort_values(by='Category')
# sorted_index_main = df_main.index
# sorted_categories_main = df_main['Category'].values
# df_main = df_main.drop(columns='Category')

# # Update y-axis ticks and labels (was x-axis before transposition)
# yticks_main = [bigg_to_name.get(idx[2:-2], idx[2:-2]) for idx in sorted_index_main]
# xticks_main = [clean_label(col) for col in df_main.columns]

# # Update color mapping for y-axis labels
# def get_y_label_color(label):
#     category = categorize_x_label(label[2:-2])
#     return groups_y.get(category, {}).get('color', 'black')

# # Plot the transposed heatmap
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
#     ax1.text(i, -1.2, tick, color=get_text_color(category), fontsize=12, ha='center', va='bottom', rotation=90)

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
# plt.savefig('heatmap_Ecoli_carbonsubs_transposed.png', dpi=300, bbox_inches='tight')
# plt.show()

