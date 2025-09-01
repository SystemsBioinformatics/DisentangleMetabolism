# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 16:01:45 2024

@author: mre283
"""

import pandas as pd
import separate_cat_ana as sca
import PrecursorSeparator as ps
import os

#%%

cwd = os.getcwd()
result_path = cwd+'\Results\\Mbarkeri'
input_files = cwd+r'\input_files/'

#%%

acetate =  sca.separate_cat_ana(input_files+'readfile_acetate_Mbarkeri.xlsx', all_cofactors=False)
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


bicarbonate =  sca.separate_cat_ana(input_files+'readfile_bicarbonateH2_Mbarkeri.xlsx', all_cofactors=False)
bicarbonate.do_cat_ana_split()
dfs = bicarbonate.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = bicarbonate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_f420_DASH_2h2_c', 'M_adp_c', 'M_f420_DASH_2_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = bicarbonate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=bicarbonate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(bicarbonate)

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

with pd.ExcelWriter(result_path+"\\bicarbonateH2.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")