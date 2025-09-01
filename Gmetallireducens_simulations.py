# -*- coding: utf-8 -*-
"""
Created on Fri Apr 19 13:36:10 2024

@author: mre283
"""

import pandas as pd
import separate_cat_ana as sca
import PrecursorSeparator as ps
import os

#%%

cwd = os.getcwd()
result_path = cwd+'\Results\\Gmetallireducens'
input_files = cwd+r'\input_files/'


#%%

acetateFe3 =  sca.separate_cat_ana(input_files+'readfile_acetateFe3_Gmetallireducens.xlsx', all_cofactors=True)
acetateFe3.do_cat_ana_split()
dfs = acetateFe3.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = acetateFe3.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = acetateFe3.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=acetateFe3.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(acetateFe3)

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

with pd.ExcelWriter(result_path+"\\acetateFe3_test.xlsx") as writer:
   
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

acetateNO3 =  sca.separate_cat_ana(input_files+'readfile_acetateNO3_Gmetallireducens.xlsx', all_cofactors=True)
acetateNO3.do_cat_ana_split()
dfs = acetateNO3.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = acetateNO3.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = acetateNO3.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=acetateNO3.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(acetateNO3)

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

with pd.ExcelWriter(result_path+"\\acetateNO3.xlsx") as writer:
   
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

ethanolNO2 =  sca.separate_cat_ana(input_files+'readfile_ethanolNO2_Gmetallireducens.xlsx', all_cofactors=True)
ethanolNO2.do_cat_ana_split()
dfs = ethanolNO2.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = ethanolNO2.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = ethanolNO2.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=ethanolNO2.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(ethanolNO2)

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


with pd.ExcelWriter(result_path+"\\ethanolNO2.xlsx") as writer:
   
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

formateFe3 =  sca.separate_cat_ana(input_files+'readfile_formateFe3_Gmetallireducens.xlsx', all_cofactors=True)
formateFe3.do_cat_ana_split()
dfs = formateFe3.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = formateFe3.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nad_c', 'M_nadh_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = formateFe3.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=formateFe3.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(formateFe3)

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

with pd.ExcelWriter(result_path+"\\formateFe3.xlsx") as writer:
   
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


butyrateFe3 =  sca.separate_cat_ana(input_files+'readfile_butyrateFe3_Gmetallireducens.xlsx', all_cofactors=True)
butyrateFe3.do_cat_ana_split()
dfs = butyrateFe3.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = butyrateFe3.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nad_c', 'M_nadh_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = butyrateFe3.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=butyrateFe3.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(butyrateFe3)

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


with pd.ExcelWriter(result_path+"\\butyrateFe3.xlsx") as writer:
   
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


ethanolFe3 =  sca.separate_cat_ana(input_files+'readfile_ethanolFe3_Gmetallireducens.xlsx', all_cofactors=True)
ethanolFe3.do_cat_ana_split()
dfs = ethanolFe3.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = ethanolFe3.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = ethanolFe3.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=ethanolFe3.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(ethanolFe3)

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

with pd.ExcelWriter(result_path+"\\ethanolFe3.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")
    
