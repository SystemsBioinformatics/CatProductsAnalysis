# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 15:32:04 2025

@author: mre283
"""


import time
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cbmpy
import matplotlib.patches as patches
import separate_cat_ana as sca
import PrecursorSeparator as ps

#%%

path = r"C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application\readfile_ecoli_NMN_GAPNx_anaerobic.xlsx"
obj =  sca.separate_cat_ana(path, all_cofactors=False)

obj.do_cat_ana_split()

dfs = obj.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = obj.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadh_c', 'M_adp_c', 'M_nad_c', 'M_nadph_c', 'M_nadp_c', 'M_pi_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T



specific_model = obj.active_model.clone()
prec_glc = ps.PrecursorSeparator(specific_model, result_df, biomass_reac = 'R_EX_BIOMASS')
fluxes, mceq, fva_df, kos = prec_glc.run(obj)

energy_obj = sca.energy_generation(specific_model)

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

result_path = r"C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application\Result files"

with pd.ExcelWriter(result_path+f"/ecoli_NMN_GAPNx_lactate.xlsx") as writer:
   
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

path = r"C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application\readfile_ecoli_anaerobic_GAPNx.xlsx"
obj =  sca.separate_cat_ana(path, all_cofactors=False)

obj.do_cat_ana_split()

dfs = obj.generate_dfs()

result_df, atp_df, cofactor_df = dfs

mceq_dict_biomass = obj.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadh_c', 'M_adp_c', 'M_nad_c', 'M_nadph_c', 'M_nadp_c', 'M_pi_c'])

df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T



specific_model = obj.active_model.clone()
prec_glc = ps.PrecursorSeparator(specific_model, result_df, biomass_reac = 'R_EX_BIOMASS')
fluxes, mceq, fva_df, kos = prec_glc.run(obj)

energy_obj = sca.energy_generation(specific_model)

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

result_path = r"C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application\Result files"

with pd.ExcelWriter(result_path+f"/ecoli_anaerobic_GAPNx.xlsx") as writer:
   
    # use to_excel function and specify the sheet_name and index 
    # to store the dataframe in specified sheet
    result_df.to_excel(writer, sheet_name="overall")
    atp_df.to_excel(writer, sheet_name="atp_overall")
    cofactor_df.to_excel(writer, sheet_name="cofactor_overall")
    df_mceqs.to_excel(writer, sheet_name="mceqs")
    
    pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
    pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
    fva_df.to_excel(writer, sheet_name="fva_BP")