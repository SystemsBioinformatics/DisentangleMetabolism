# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 08:27:47 2024

@author: mre283
"""

import pandas as pd
import separate_cat_ana as sca
import PrecursorSeparator as ps
import os

#%%

cwd = os.getcwd()
result_path = cwd+'\Results\\Scerevisiae'
input_files = cwd+r'\input_files/'

#%%

aer_glc =  sca.separate_cat_ana(input_files+'readfile_aerobic_glucose_yeast.xlsx', all_cofactors=False)
aer_glc.do_cat_ana_split()
dfs = aer_glc.generate_dfs()

result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_carbohydrate, df_rna, df_dna, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_carbohydrate_atp, df_rna_atp, df_dna_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_carbohydrate_cofactor, df_rna_cofactor, df_dna_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = aer_glc.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = aer_glc.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = aer_glc.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_carbohydrate = aer_glc.generate_MCEQs(df_carbohydrate, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = aer_glc.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = aer_glc.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = aer_glc.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = aer_glc.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = aer_glc.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_carbohydrate, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'carbohydrate', 'rna', 'dna', 'ion', 'cofactor', 'GAM']).T

model = aer_glc.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=aer_glc.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(aer_glc)

energy_obj = sca.energy_generation(model)

split_names = False
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

with pd.ExcelWriter(result_path+"\\aerobic_glucose.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_carbohydrate.to_excel(writer, sheet_name="carbohydrate")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_carbohydrate_atp.to_excel(writer, sheet_name="atp_carbohydrate")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_carbohydrate_cofactor.to_excel(writer, sheet_name="cofactor_carbohydrate")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")
    
#%%

lactate =  sca.separate_cat_ana(input_files+'readfile_lactate_yeast.xlsx', all_cofactors=False)
lactate.do_cat_ana_split()
dfs = lactate.generate_dfs()
result_df, atp_df, cofactor_df = dfs


mceq_dict_biomass = lactate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = lactate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=lactate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(lactate)

energy_obj = sca.energy_generation(model)

split_names = False
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

galactose =  sca.separate_cat_ana(input_files+'readfile_galactose_yeast.xlsx', all_cofactors=False)
galactose.do_cat_ana_split()
dfs = galactose.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = galactose.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = galactose.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=galactose.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(galactose)

energy_obj = sca.energy_generation(model)

split_names = False
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

with pd.ExcelWriter(result_path+"\\galactose.xlsx") as writer:
   
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

ribose =  sca.separate_cat_ana(input_files+'readfile_ribose_yeast.xlsx', all_cofactors=False)
ribose.do_cat_ana_split()
dfs = ribose.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = ribose.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = ribose.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=ribose.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(ribose)

energy_obj = sca.energy_generation(model)

split_names = False
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

with pd.ExcelWriter(result_path+"\\ribose.xlsx") as writer:
   
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

maltotriose =  sca.separate_cat_ana(input_files+'readfile_maltotriose_yeast.xlsx', all_cofactors=False)
maltotriose.do_cat_ana_split()
dfs = maltotriose.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = maltotriose.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = maltotriose.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=maltotriose.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(maltotriose)

energy_obj = sca.energy_generation(model)

split_names = False
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

with pd.ExcelWriter(result_path+"\\maltotriose.xlsx") as writer:
   
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

pyruvate =  sca.separate_cat_ana(input_files+'readfile_pyruvate_yeast.xlsx', all_cofactors=False)
pyruvate.do_cat_ana_split()
dfs = pyruvate.generate_dfs()
result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = pyruvate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = pyruvate.active_model.clone()
prec = ps.PrecursorSeparator(model, result_df, biomass_reac=pyruvate.biomass_reac)
fluxes, mceq, fva_df, kos = prec.run(pyruvate)

energy_obj = sca.energy_generation(model)

split_names = False
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

with pd.ExcelWriter(result_path+"\\pyruvate.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")