# -*- coding: utf-8 -*-
"""
Created on Mon Apr  8 11:41:10 2024

@author: mre283
"""


import pandas as pd
import separate_cat_ana as sca
import PrecursorSeparator as ps
import os


#%%

cwd = os.getcwd()
result_path = cwd+'\Results\\Cljungdahlii'
input_files = cwd+r'\input_files/'

#%%

carbonmonoxide =  sca.separate_cat_ana(input_files+'readfile_carbonmonoxide_Cljungdahlii.xlsx', all_cofactors=False)
carbonmonoxide.do_cat_ana_split()
dfs = carbonmonoxide.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = carbonmonoxide.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_fdxr_42_c', 'M_adp_c', 'M_fdxo_42_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T


model = carbonmonoxide.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=carbonmonoxide.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(carbonmonoxide)

energy_obj = sca.energy_generation(model)

split_names = True
atp_consumption, atp_production = energy_obj.calculate_energy(fluxes, split_names=split_names)
nadh_prod, nadh_cons = energy_obj.calculate_NADH_prod(fluxes, split_names=split_names)
nadph_prod, nadph_cons = energy_obj.calculate_NADPH_prod(fluxes, split_names=split_names)

info_dict = {}
info_dict['mceq'] = mceq
info_dict['atp production'] = atp_production
info_dict['atp consumption'] = atp_consumption
info_dict['nadh production'] = nadh_prod
info_dict['nadh consumption'] = nadh_cons
info_dict['nadph production'] = nadph_prod
info_dict['nadph consumption'] = nadph_cons

with pd.ExcelWriter(result_path+"\\carbonmonoxide.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")

#%%

carbonmonoxide_hydrogen =  sca.separate_cat_ana(input_files+'readfile_carbonmonoxidehydrogen_Cljungdahlii.xlsx', all_cofactors=False)
carbonmonoxide_hydrogen.do_cat_ana_split()
dfs = carbonmonoxide_hydrogen.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = carbonmonoxide_hydrogen.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_fdxr_42_c', 'M_adp_c', 'M_fdxo_42_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = carbonmonoxide_hydrogen.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=carbonmonoxide_hydrogen.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(carbonmonoxide_hydrogen)

energy_obj = sca.energy_generation(model)

split_names = True
atp_consumption, atp_production = energy_obj.calculate_energy(fluxes, split_names=split_names)
nadh_prod, nadh_cons = energy_obj.calculate_NADH_prod(fluxes, split_names=split_names)
nadph_prod, nadph_cons = energy_obj.calculate_NADPH_prod(fluxes, split_names=split_names)

info_dict = {}
info_dict['mceq'] = mceq
info_dict['atp production'] = atp_production
info_dict['atp consumption'] = atp_consumption
info_dict['nadh production'] = nadh_prod
info_dict['nadh consumption'] = nadh_cons
info_dict['nadph production'] = nadph_prod
info_dict['nadph consumption'] = nadph_cons

with pd.ExcelWriter(result_path+"\\carbonmonoxidehydrogen.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

carbondioxide_hydrogen =  sca.separate_cat_ana(input_files+'readfile_carbondioxidehydrogen_Cljungdahlii.xlsx', all_cofactors=False)
carbondioxide_hydrogen.do_cat_ana_split()
dfs = carbondioxide_hydrogen.generate_dfs()
result_df, atp_df, cofactor_df = dfs


mceq_dict_biomass = carbondioxide_hydrogen.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_fdxr_42_c', 'M_adp_c', 'M_fdxo_42_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = carbondioxide_hydrogen.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=carbondioxide_hydrogen.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(carbondioxide_hydrogen)

energy_obj = sca.energy_generation(model)

split_names = True
atp_consumption, atp_production = energy_obj.calculate_energy(fluxes, split_names=split_names)
nadh_prod, nadh_cons = energy_obj.calculate_NADH_prod(fluxes, split_names=split_names)
nadph_prod, nadph_cons = energy_obj.calculate_NADPH_prod(fluxes, split_names=split_names)

info_dict = {}
info_dict['mceq'] = mceq
info_dict['atp production'] = atp_production
info_dict['atp consumption'] = atp_consumption
info_dict['nadh production'] = nadh_prod
info_dict['nadh consumption'] = nadh_cons
info_dict['nadph production'] = nadph_prod
info_dict['nadph consumption'] = nadph_cons

with pd.ExcelWriter(result_path+"\\carbondioxidehydrogen.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")
