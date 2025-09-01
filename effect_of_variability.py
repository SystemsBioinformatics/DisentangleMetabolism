# -*- coding: utf-8 -*-
"""
Created on Thu Jul  4 15:35:06 2024

@author: mre283
"""

import separate_cat_ana as sca
import pandas as pd
import numpy as np
import os
import cbmpy

#%%

modelPath = 'C:\\Users\\mre283\\OneDrive - Vrije Universiteit Amsterdam\\89 Genome scale models' + '\iML1515_macroBOF_core.xml'
model = cbmpy.loadModel(modelPath)

model.setReactionBounds('R_ATPM', 0, 1000)
model.setReactionBounds('R_EX_BIOMASS', 1.0, 1000)
model.setReactionBounds('R_EX_glc__D_e', -1000, 1000)

model.setObjectiveFlux('R_EX_glc__D_e', osense='maximize')

result = cbmpy.CBGLPK.glpk_analyzeModel(model)
minsum = cbmpy.CBCPLEX.cplx_MinimizeSumOfAbsFluxes(model)

fluxes = model.getReactionValues()

#%%

model.setReactionBounds('R_EX_glc__D_e', result, result)

fva = cbmpy.CBCPLEX.cplx_FluxVariabilityAnalysis(model)

df_fva = pd.DataFrame(fva[0], index=fva[1], columns = ['Reaction', 'Reduced Costs', 'Variability Min', 'Variability Max', 'abs(Max-Min)', 'MinStatus', 'MaxStatus'])

#%%

model2 = model.clone()

# model3 = model.clone()
# model3.setReactionBounds('R_EX_glc__D_e', -1000, 1000)

fluxlistlist = []


for r in df_fva.index:
    if 1e-6 < df_fva['abs(Max-Min)'].loc[r] < 100:
        if np.abs(df_fva['Variability Max'].loc[r] - df_fva['Reaction'].loc[r]) > np.abs(df_fva['Reaction'].loc[r] - df_fva['Variability Max'].loc[r]):
            sense = 'maximize'
        else:
            sense = 'minimize'
        
        print(r)
        print(sense)
        model2.setObjectiveFlux(r, osense=sense)
        result2 = cbmpy.CBGLPK.glpk_analyzeModel(model2)
        fluxes2 = model2.getReactionValues()
        fluxlist = []
        for flux, value in fluxes2.items():
            
            # print(flux)
            # print(value)
            
            if np.abs(value) <1e-6 and np.abs(df_fva['Reaction'].loc[flux]) > 1e-6:
                print(flux)
                fluxlist.append(flux)
            
        fluxlistlist.append(fluxlist)
        
        if len(fluxlistlist) > 99:
            break

        
unique_fluxlistlist = [list(x) for x in set(tuple(x) for x in fluxlistlist)]

#%%

for i in range(len(unique_fluxlistlist)):
    aer_glc_nadh =  sca.separate_cat_ana('readfile_aerobic_glucose_Ecoli_var.xlsx', all_cofactors=False)
    
    for r in unique_fluxlistlist[i]:
        aer_glc_nadh.model.setReactionBounds(r, 0.0, 0.0)
    
    aer_glc_nadh.do_cat_ana_split()
    
    dfs = aer_glc_nadh.generate_dfs()
    
    
    result_df, atp_df, cofactor_df = dfs
    
    mceq_dict_biomass = aer_glc_nadh.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadh_c', 'M_adp_c', 'M_nad_c', 'M_nadph_c', 'M_nadp_c'])
    
    df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T
    
    result_path = 'C:\\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\\01 Separating catabolism and anabolism\python files\generalized_scripts\Results\\Variability_check'
    
    
    with pd.ExcelWriter(result_path+"\\EFM_nr_{}.xlsx".format(i)) as writer:
       
        # use to_excel function and specify the sheet_name and index 
        # to store the dataframe in specified sheet
        result_df.to_excel(writer, sheet_name="overall")
        atp_df.to_excel(writer, sheet_name="atp_overall")
        cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
        df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    
 

# for a reaction with value 0 and a nonzero variability:
    # maximize or minimize this reaction and see reactions that are now zero in FVA
    # soft-KO these reactions (first check: does it now work that we get the target reaction active in the model if we do FBA?)
    # perform FBA and split cat and ana
    # Look at some characteristics of splitting: ATP transfer, ATP fraction, NADPH transfer

#%%


directory = 'Results\Variability_check'
df_dict_atp = {}

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        
        df_dict_atp[filename.path[26:-5]] = pd.read_excel(filename.path, sheet_name='atp_overall', index_col=0)



df_dict_cofactor = {}

# iterate over files in 
# that directory
for filename in os.scandir(directory):
    if filename.is_file():
        
        df_dict_cofactor[filename.path[26:-5]] = pd.read_excel(filename.path, sheet_name='cofactor_overall', index_col=0)


#%%

df_compare = pd.DataFrame(index = list(df_dict_atp.keys()), columns = ['ATP transferred', 'ATP fraction', 'NADPH transferred'])

for EFM, df in df_dict_atp.items():
    df_compare.loc[EFM, 'ATP transferred'] = df.loc['cat', 'net']
    df_compare.loc[EFM, 'ATP fraction'] = df.loc['ana', 'prod']/df.loc['ana', 'cons']

for EFM, df in df_dict_cofactor.items():
    df_compare.loc[EFM, 'NADPH transferred'] = df.loc['cat', 'NADPH net']
    
EFMlist = list(df_dict_atp.keys())
for i in range(len(unique_fluxlistlist[:-2])):
    
    df_compare.loc[EFMlist[i], 'Deactivated fluxes'] = str(unique_fluxlistlist[i])


with pd.ExcelWriter(result_path+"\\comparison.xlsx") as writer:
    
    df_compare.to_excel(writer, sheet_name="comparison")








