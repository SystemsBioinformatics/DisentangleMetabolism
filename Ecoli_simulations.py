# -*- coding: utf-8 -*-
"""
Created on Thu Mar 28 09:40:19 2024

@author: mre283
"""

import separate_cat_ana as sca
import pandas as pd
import os
import PrecursorSeparator as ps

#%%
cwd = os.getcwd()
result_path = cwd+'\Results\\Ecoli'
input_files = cwd+r'\input_files/'

#%%


aer_glc =  sca.separate_cat_ana(input_files+'readfile_aerobic_glucose_Ecoli.xlsx', all_cofactors=False)
aer_glc.do_cat_ana_split()
dfs = aer_glc.generate_dfs()


result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs # 


df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs


mceq_dict_biomass = aer_glc.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = aer_glc.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = aer_glc.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = aer_glc.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = aer_glc.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = aer_glc.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = aer_glc.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = aer_glc.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = aer_glc.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = aer_glc.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                      index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T



model = aer_glc.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(aer_glc)

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
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")



#%%
aer_glc_nadh =  sca.separate_cat_ana(input_files+'readfile_aerobic_glucose_Ecoli_nadh.xlsx', all_cofactors=False)
aer_glc_nadh.do_cat_ana_split()
dfs = aer_glc_nadh.generate_dfs()

result_df, atp_df, cofactor_df = dfs
mceq_dict_biomass = aer_glc_nadh.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadh_c', 'M_adp_c', 'M_nad_c', 'M_nadph_c', 'M_nadp_c'])
df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = aer_glc_nadh.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(aer_glc_nadh)

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


with pd.ExcelWriter(result_path+"\\aer_glc_nadh.xlsx") as writer:
   
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

anaer_glc = sca.separate_cat_ana(input_files+'readfile_anaerobic_glucose_Ecoli.xlsx', all_cofactors=False)

anaer_glc.do_cat_ana_split()

dfs = anaer_glc.generate_dfs()

result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs


mceq_dict_biomass = anaer_glc.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_protein = anaer_glc.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_lipid = anaer_glc.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_lps = anaer_glc.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_rna = anaer_glc.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_dna = anaer_glc.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_murein = anaer_glc.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_ion = anaer_glc.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_cofactor = anaer_glc.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])
mceq_dict_GAM = anaer_glc.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_nadh_c', 'M_nad_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T


model = anaer_glc.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(anaer_glc)

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


with pd.ExcelWriter(result_path+"\\glucose_anaerobic.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")

    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

aer_ac = sca.separate_cat_ana(input_files+'readfile_aerobic_acetate_Ecoli.xlsx', all_cofactors=False)

aer_ac.do_cat_ana_split()

dfs = aer_ac.generate_dfs()


result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = aer_ac.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = aer_ac.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = aer_ac.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = aer_ac.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = aer_ac.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = aer_ac.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = aer_ac.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = aer_ac.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = aer_ac.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = aer_ac.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T


model = aer_ac.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(aer_ac)

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
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")

#%%


aer_xyl = sca.separate_cat_ana(input_files+'readfile_aerobic_xylose_Ecoli.xlsx', all_cofactors=False)


aer_xyl.do_cat_ana_split()

dfs = aer_xyl.generate_dfs()


result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = aer_xyl.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = aer_xyl.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = aer_xyl.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = aer_xyl.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = aer_xyl.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = aer_xyl.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = aer_xyl.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = aer_xyl.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = aer_xyl.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = aer_xyl.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T



model = aer_xyl.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(aer_xyl)

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



with pd.ExcelWriter(result_path+"\\xylose.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")



#%%

fruct = sca.separate_cat_ana(input_files+'readfile_aerobic_fructose_Ecoli.xlsx', all_cofactors=False)


fruct.do_cat_ana_split()

dfs = fruct.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = fruct.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = fruct.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = fruct.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = fruct.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = fruct.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = fruct.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = fruct.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = fruct.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = fruct.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = fruct.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T



model = fruct.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(fruct)

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


with pd.ExcelWriter(result_path+"\\fructose.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

galact = sca.separate_cat_ana(input_files+'readfile_aerobic_galactose_Ecoli.xlsx', all_cofactors=False)


galact.do_cat_ana_split()

dfs = galact.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = galact.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = galact.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = galact.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = galact.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = galact.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = galact.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = galact.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = galact.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = galact.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = galact.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T


model = galact.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(galact)

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


with pd.ExcelWriter(result_path+"\\galactose.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

glucon = sca.separate_cat_ana(input_files+'readfile_aerobic_gluconate_Ecoli.xlsx', all_cofactors=False)


glucon.do_cat_ana_split()

dfs = glucon.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = glucon.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = glucon.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = glucon.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = glucon.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = glucon.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = glucon.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = glucon.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = glucon.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = glucon.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = glucon.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T

model = glucon.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(glucon)

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

with pd.ExcelWriter(result_path+"\\gluconate.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

glyc = sca.separate_cat_ana(input_files+'readfile_aerobic_glycerol_Ecoli.xlsx', all_cofactors=False)


glyc.do_cat_ana_split()

dfs = glyc.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = glyc.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_gtp_c', 'M_gdp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = glyc.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(glyc)

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

pyr = sca.separate_cat_ana(input_files+'readfile_aerobic_pyruvate_Ecoli.xlsx', all_cofactors=False)


pyr.do_cat_ana_split()

dfs = pyr.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = pyr.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = pyr.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = pyr.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = pyr.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = pyr.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = pyr.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = pyr.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = pyr.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = pyr.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = pyr.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T

model = pyr.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(pyr)

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

with pd.ExcelWriter(result_path+"\\pyruvate.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

ribose = sca.separate_cat_ana(input_files+'readfile_aerobic_ribose_Ecoli.xlsx', all_cofactors=False)


ribose.do_cat_ana_split()

dfs = ribose.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = ribose.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = ribose.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = ribose.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = ribose.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = ribose.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = ribose.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = ribose.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = ribose.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = ribose.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = ribose.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T

model = ribose.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(ribose)

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

with pd.ExcelWriter(result_path+"\\ribose.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

sorb = sca.separate_cat_ana(input_files+'readfile_aerobic_sorbitol_Ecoli.xlsx', all_cofactors=False)


sorb.do_cat_ana_split()

dfs = sorb.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = sorb.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = sorb.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = sorb.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = sorb.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = sorb.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = sorb.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = sorb.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = sorb.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = sorb.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = sorb.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T

model = sorb.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(sorb)

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

with pd.ExcelWriter(result_path+"\\sorbitol.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_GAM_cofactor.to_excel(writer, sheet_name="cofactor_GAM")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

succ = sca.separate_cat_ana(input_files+'readfile_aerobic_succinate_Ecoli.xlsx', all_cofactors=False)


succ.do_cat_ana_split()

dfs = succ.generate_dfs()



result_df, atp_df, cofactor_df, macromolecule_dfs, atp_macromolecule_dfs, cofactor_macromolecule_dfs = dfs
df_protein, df_lipid, df_lps, df_rna, df_dna, df_murein, df_ion, df_cofactor, df_GAM = macromolecule_dfs
df_protein_atp, df_lipid_atp, df_lps_atp, df_rna_atp, df_dna_atp, df_murein_atp, df_ion_atp, df_cofactor_atp, df_GAM_atp = atp_macromolecule_dfs
df_protein_cofactor, df_lipid_cofactor, df_lps_cofactor, df_rna_cofactor, df_dna_cofactor, df_murein_cofactor, df_ion_cofactor, df_cofactor_cofactor, df_GAM_cofactor = cofactor_macromolecule_dfs

mceq_dict_biomass = succ.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_protein = succ.generate_MCEQs(df_protein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lipid = succ.generate_MCEQs(df_lipid, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_lps = succ.generate_MCEQs(df_lps, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_rna = succ.generate_MCEQs(df_rna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_dna = succ.generate_MCEQs(df_dna, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_murein = succ.generate_MCEQs(df_murein, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_ion = succ.generate_MCEQs(df_ion, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_cofactor = succ.generate_MCEQs(df_cofactor, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])
mceq_dict_GAM = succ.generate_MCEQs(df_GAM, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records([mceq_dict_biomass, mceq_dict_protein, mceq_dict_lipid, mceq_dict_lps, mceq_dict_rna,
                                      mceq_dict_dna, mceq_dict_murein, mceq_dict_ion, mceq_dict_cofactor, mceq_dict_GAM], 
                                     index=['overall', 'protein', 'lipid', 'lps', 'rna', 'dna', 'murein', 'ion', 'cofactor', 'GAM']).T

model = succ.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(succ)

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

with pd.ExcelWriter(result_path+"\\succinate.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_protein.to_excel(writer, sheet_name="protein")
    df_lipid.to_excel(writer, sheet_name="lipid")
    df_lps.to_excel(writer, sheet_name="lps")
    df_rna.to_excel(writer, sheet_name="rna")
    df_dna.to_excel(writer, sheet_name="dna")
    df_murein.to_excel(writer, sheet_name="murein")
    df_ion.to_excel(writer, sheet_name="ion")
    df_cofactor.to_excel(writer, sheet_name="cofactor")
    df_GAM.to_excel(writer, sheet_name="GAM")
    df_protein_atp.to_excel(writer, sheet_name="atp_protein")
    df_lipid_atp.to_excel(writer, sheet_name="atp_lipid")
    df_lps_atp.to_excel(writer, sheet_name="atp_lps")
    df_rna_atp.to_excel(writer, sheet_name="atp_rna")
    df_dna_atp.to_excel(writer, sheet_name="atp_dna")
    df_murein_atp.to_excel(writer, sheet_name="atp_murein")
    df_ion_atp.to_excel(writer, sheet_name="atp_ion")
    df_cofactor_atp.to_excel(writer, sheet_name="atp_cofactor")
    df_GAM_atp.to_excel(writer, sheet_name="atp_GAM")
    df_protein_cofactor.to_excel(writer, sheet_name="cofactor_protein")
    df_lipid_cofactor.to_excel(writer, sheet_name="cofactor_lipid")
    df_lps_cofactor.to_excel(writer, sheet_name="cofactor_lps")
    df_rna_cofactor.to_excel(writer, sheet_name="cofactor_rna")
    df_dna_cofactor.to_excel(writer, sheet_name="cofactor_dna")
    df_murein_cofactor.to_excel(writer, sheet_name="cofactor_murein")
    df_ion_cofactor.to_excel(writer, sheet_name="cofactor_ion")
    df_cofactor_cofactor.to_excel(writer, sheet_name="cofactor_cofactor")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")


#%%

citrate = sca.separate_cat_ana(input_files+'readfile_citrate_Ecoli.xlsx', all_cofactors=False)

citrate.do_cat_ana_split()

dfs = citrate.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = citrate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_gtp_c', 'M_gdp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = citrate.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(citrate)

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

with pd.ExcelWriter(result_path+"\\citrate.xlsx") as writer:
   
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

akg = sca.separate_cat_ana(input_files+'readfile_alphaketoglutarate_Ecoli.xlsx', all_cofactors=False)

akg.do_cat_ana_split()

dfs = akg.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = akg.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = akg.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(akg)

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

with pd.ExcelWriter(result_path+"\\akg.xlsx") as writer:
   
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

arginine = sca.separate_cat_ana(input_files+'readfile_arginine_Ecoli.xlsx', all_cofactors=False)

arginine.do_cat_ana_split()

dfs = arginine.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = arginine.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = arginine.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(arginine)

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

with pd.ExcelWriter(result_path+"\\arginine.xlsx") as writer:
   
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

glutamine = sca.separate_cat_ana(input_files+'readfile_glutamine_Ecoli.xlsx', all_cofactors=False)

glutamine.do_cat_ana_split()

dfs = glutamine.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = glutamine.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T


model = glutamine.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(glutamine)

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


with pd.ExcelWriter(result_path+"\\glutamine.xlsx") as writer:
   
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

lactate = sca.separate_cat_ana(input_files+'readfile_lactate_Ecoli.xlsx', all_cofactors=False)

lactate.do_cat_ana_split()

dfs = lactate.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = lactate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T


model = lactate.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(lactate)

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

malate = sca.separate_cat_ana(input_files+'readfile_malate_Ecoli.xlsx', all_cofactors=False)

malate.do_cat_ana_split()

dfs = malate.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = malate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_utp_c', 'M_udp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = malate.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(malate)

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




with pd.ExcelWriter(result_path+"\\malate.xlsx") as writer:
   
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

fumarate = sca.separate_cat_ana(input_files+'readfile_fumarate_Ecoli.xlsx', all_cofactors=False)

fumarate.do_cat_ana_split()

dfs = fumarate.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = fumarate.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadph_c', 'M_adp_c', 'M_nadp_c', 'M_gtp_c', 'M_gdp_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T

model = fumarate.active_model.clone()
prec_glc = ps.PrecursorSeparator(model, result_df)
fluxes, mceq, fva_df, kos = prec_glc.run(fumarate)

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

with pd.ExcelWriter(result_path+"\\fumarate.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")