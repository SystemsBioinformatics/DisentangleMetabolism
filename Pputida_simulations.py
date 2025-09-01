# -*- coding: utf-8 -*-
"""
Created on Wed May 15 08:37:37 2024

@author: mre283
"""

import pandas as pd
import separate_cat_ana as sca
import PrecursorSeparator as ps
import os

#%%

cwd = os.getcwd()
result_path = cwd+'\Results\\Pputida'
input_files = cwd+r'\input_files/'

#%%

glucose =  sca.separate_cat_ana(input_files+'readfile_glucose_Pputida.xlsx', all_cofactors=False)
glucose.do_cat_ana_split()
dfs = glucose.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = glucose.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = glucose.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=glucose.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(glucose)

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

with pd.ExcelWriter(result_path+"\\glucose.xlsx") as writer:
   
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


acetate =  sca.separate_cat_ana(input_files+'readfile_acetate_Pputida.xlsx', all_cofactors=False)
acetate.do_cat_ana_split()
dfs = acetate.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = acetate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = acetate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=acetate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(acetate)

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

with pd.ExcelWriter(result_path+"\\acetate.xlsx") as writer:
   
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


ethanol =  sca.separate_cat_ana(input_files+'readfile_ethanol_Pputida.xlsx', all_cofactors=False)
ethanol.do_cat_ana_split()
dfs = ethanol.generate_dfs()
result_df, atp_df, cofactor_df = dfs


mceq_dict_biomass = ethanol.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = ethanol.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=ethanol.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(ethanol)

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

with pd.ExcelWriter(result_path+"\\ethanol.xlsx") as writer:
   
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


methanol =  sca.separate_cat_ana(input_files+'readfile_methanol_Pputida.xlsx', all_cofactors=False)
methanol.do_cat_ana_split()
dfs = methanol.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = methanol.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = methanol.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=methanol.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(methanol)

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

with pd.ExcelWriter(result_path+"\\methanol.xlsx") as writer:
   
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


lactate =  sca.separate_cat_ana(input_files+'readfile_lactate_Pputida.xlsx', all_cofactors=False)
lactate.do_cat_ana_split()
dfs = lactate.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = lactate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = lactate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=lactate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(lactate)

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

with pd.ExcelWriter(result_path+"\\lactate.xlsx") as writer:
   
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
glycerol =  sca.separate_cat_ana(input_files+'readfile_glycerol_Pputida.xlsx', all_cofactors=False)
glycerol.do_cat_ana_split()
dfs = glycerol.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = glycerol.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = glycerol.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=glycerol.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(glycerol)

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

with pd.ExcelWriter(result_path+"\\glycerol.xlsx") as writer:
   
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


octadecanoate =  sca.separate_cat_ana(input_files+'readfile_octadecanoate_Pputida.xlsx', all_cofactors=False)
octadecanoate.do_cat_ana_split()
dfs = octadecanoate.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = octadecanoate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = octadecanoate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=octadecanoate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(octadecanoate)

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

with pd.ExcelWriter(result_path+"\\octadecanoate.xlsx") as writer:
   
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

hexadecanoate161 =  sca.separate_cat_ana(input_files+'readfile_hexadecanoate161_Pputida.xlsx', all_cofactors=False)
hexadecanoate161.do_cat_ana_split()
dfs = hexadecanoate161.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = hexadecanoate161.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = hexadecanoate161.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=hexadecanoate161.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(hexadecanoate161)

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

with pd.ExcelWriter(result_path+"\\hexadecanoate161.xlsx") as writer:
   
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


putrescine =  sca.separate_cat_ana(input_files+'readfile_putrescine_Pputida.xlsx', all_cofactors=False)
putrescine.do_cat_ana_split()
dfs = putrescine.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = putrescine.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = putrescine.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=putrescine.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(putrescine)

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

with pd.ExcelWriter(result_path+"\\putrescine.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")