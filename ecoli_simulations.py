# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 11:16:18 2025

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
                if s[-1] == 'e' and s not in catabolic_external and s not in possible_energy_carriers:
                    catabolic_external.append(s)

            for s in species:
                if s not in catabolic_intermediates and (s[-1] == 'c' or s[-1] == 'm') and s not in possible_energy_carriers:
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
    model.createReaction(f'R_{s[2:]}_trans', name=f'{s[2:-2]} transport', reversible=True)
    model.createReactionReagent(f'R_{s[2:]}_trans', s, -1)
    model.createReactionReagent(f'R_{s[2:]}_trans', s[:-1] + 'e', 1)
    model.setReactionBounds(f'R_{s[2:]}_trans', -1000, 1000)
    print(model.getReaction(f'R_{s[2:]}_trans').getEquation())


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
            if 'M_thf_e' not in model.getSpeciesIds():
                model.createSpecies('M_thf_e')
            model.createReactionReagent('R_EX_' + s[2:-2] + '_e', 'M_thf_e', 1)
        elif ('phosph' in name or 'Phosph' in name) and s!= 'M_pi_c':
            model.createReactionReagent('R_EX_' + s[2:-2] + '_e', 'M_pi_e', number_after_P)
        model.setReactionBounds('R_EX_' + s[2:-2] + '_e', -1000, 1000)
        print(model.getReaction('R_EX_' + s[2:-2] + '_e').getEquation())
    
    
        return 'R_EX_' + s[2:-2] + '_e'


def add_supply_reactions(model, energy_carriers):
    """
    Add supply reactions for energy carriers that are in the energy_carriers list.
    
    Args:
        model: The metabolic model
        energy_carriers: List of energy carrier metabolites
    """
    # Supply reactions defined with metabolite IDs and coefficients
    supplies = {
        'NADPH': [('M_nadph_c', 1.0), ('M_nadp_c', -1.0), ('M_h_c', -1.0)],
        'NADH': [('M_nadh_c', 1.0), ('M_nad_c', -1.0), ('M_h_c', -1.0)],
        'ATP': [('M_atp_c', 1.0), ('M_h2o_c', 1.0), ('M_adp_c', -1.0), ('M_h_c', -1.0), ('M_pi_c', -1.0)],
        'FDXR': [('M_fdxo_42_c', -1.0), ('M_fdxr_42_c', 1.0)],
        'FADH': [('M_fadh2_c', 1.0), ('M_fad_c', -1.0), ('M_h_c', -2.0)],
        'Q8H2': [('M_q8h2_c', 1.0), ('M_q8_c', -1.0)],
        'ATP_m': [('M_atp_m', 1.0), ('M_h2o_m', 1.0), ('M_adp_m', -1.0), ('M_h_m', -1.0), ('M_pi_m', -1.0)],
        'NADH_m': [('M_nadh_m', 1.0), ('M_nad_m', -1.0), ('M_h_m', -1.0)],
        'NADPH_m': [('M_nadph_m', 1.0), ('M_nadp_m', -1.0), ('M_h_m', -1.0)],
        'F420R': [('M_f420_DASH_2h2_c', 1.0), ('M_f420_DASH_2_c', -1.0)],
        'FDRED':[('M_fdred_c', 1.0), ('M_fdox_c', -1.0)],
        'FDXR2': [('M_fdxrd_c', 1.0), ('M_fdxo_2_2_c', -1.0)],
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
    
    # 'M_focytc_m', 'M_ficytc_m', 'M_q6h2_m', 'M_q6_m',
    # 'M_fad_m', 'M_fadh2_m', 'M_lpam_m', 'M_dhlam_m'
    added_supplies = []

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
    
    # catabolic_intermediates.append('M_o2_l')
    # catabolic_intermediates.append('M_h2o_l')
    # catabolic_intermediates.append('M_o2_c')
    # catabolic_intermediates.append('M_photon_e')
    # catabolic_intermediates.append('M_o2_e')
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
            
        if catabolic_reac is True and 'ADK1' not in r:
            print(r)
            model.setReactionBounds(r, 0, 0)
            kos.append(r)
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
        'M_atp_c', 'M_adp_c', 
        'M_nad_c', 'M_nadh_c', 
        'M_nadph_c', 'M_nadp_c', 
         'M_fdxr_42_c', 'M_fdxo_42_c', 
        'M_fad_c', 'M_fadh2_c', 'M_q8_c', 'M_q8h2_c',
        'M_atp_m', 'M_adp_m', 'M_nad_m', 'M_nadh_m',
        'M_nadph_m', 'M_nadp_m', 'M_f420_DASH_2h2_c',
        'M_f420_DASH_2_c', 'M_fdred_c', 'M_fdox_c',
        'M_mphenh2_c', 'M_mphen_c', 'M_omchr_e', 
        'M_omcho_e', 'M_focytC_c', 'M_ficytC_c',
        'M_flxso_c', 'M_flxr_c', 'M_mql8_c', 'M_mqn8_c',
        'M_focytc_m', 'M_ficytc_m', 'M_q6h2_m', 'M_q6_m',
        'M_fad_m', 'M_fadh2_m', 'M_lpam_m', 'M_dhlam_m', 
        'M_gdp_c', 'M_gtp_c', 'M_ttp_c',
        'M_tdp_c', 'M_utp_c', 'M_udp_c', 'M_fdxo_2_2_c', 'M_fdxr_c',
        'M_nmn_c', 'M_nmnh_c'
    ]

    # Step 1: Process inputs to identify species
    catabolic_intermediates, catabolic_external, energy_carriers = process_inputs(model, result_df, possible_energy_carriers)
    print(catabolic_external)
    print(catabolic_intermediates)
    print(energy_carriers)

    # Step 2: Create species and reactions
    ex_added = create_species_and_reactions(model, catabolic_intermediates, energy_carriers)

    # Step 3: Add supply reactions
    added_supplies = add_supply_reactions(model, energy_carriers)
    print(added_supplies)
    # Step 4: Adjust reaction bounds
    kos = adjust_reaction_bounds(model, catabolic_external, catabolic_intermediates, energy_carriers)
    
    # Set the objective for optimization
    model.setObjectiveFlux(biomass_reac, osense='minimize')
    # model.setObjectiveFlux('R_ATP_supply_fwd', osense='maximize')
    
    model_notsplit = model.clone()
    # model.setReactionBounds(substrate, 0, 0)
    # Analyze the model using CBGLPK
    model = cbmpy.CBTools.splitReversibleReactions(model)
    # model.setReactionBounds('R_EX_oaa_e_fwd', 0, 0)
    # model.setReactionBounds('R_ADK1_fwd', 0, 1000)
    # model_notsplit.setReactionBounds('R_ADK1_fwd', 0, 1000)
    # model.setObjectiveFlux('R_ATP_supply_fwd', osense='maximize')
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



# file_dict =  {r'readfile_yeast_NMN_GAPN.xlsx':'R_GROWTH', r'readfile_yeast_NMN_heterolactic.xlsx':'R_GROWTH', r'readfile_yeast_NMN_calvin.xlsx':'R_GROWTH',
#               r'readfile_yeast_aerobic_GAPN.xlsx':'R_GROWTH',r'readfile_yeast_aerobic_heterolactic.xlsx':'R_GROWTH',r'readfile_yeast_NMN_calvin.xlsx':'R_GROWTH',}

# file_dict = {r'readfile_yeast_aerobic_calvin.xlsx':'R_GROWTH',}

# file_dict = {r'readfile_ecoli_NMN_calvin_anaerobic.xlsx':'R_EX_BIOMASS',r'readfile_ecoli_NMN_ED_anaerobic.xlsx':'R_EX_BIOMASS',
#              r'readfile_ecoli_NMN_GAPN_anaerobic.xlsx':'R_EX_BIOMASS',}
# file_dict= {r'readfile_ecoli_anaerobic_ED.xlsx':'R_EX_BIOMASS',
#              r'readfile_ecoli_anaerobic_GAPN.xlsx':'R_EX_BIOMASS',r'readfile_ecoli_anaerobic_heterolactic.xlsx':'R_EX_BIOMASS',
#              r'readfile_ecoli_anaerobic_lacprod.xlsx':'R_EX_BIOMASS',r'readfile_ecoli_WT_aerobic.xlsx':'R_EX_BIOMASS',
#              r'readfile_ecoli_WT_anaerobic.xlsx':'R_EX_BIOMASS',r'readfile_ecoli_anaerobic_calvin.xlsx':'R_EX_BIOMASS',}

file_dict = {r'readfile_ecoli_NMN_GAPN_anaerobic.xlsx':'R_EX_BIOMASS', r'readfile_ecoli_NMN_heterolactic_calvin.xlsx':'R_EX_BIOMASS',}

# to do 2/4
# r'readfile_ecoli_NMN_heterolactic_anaerobic.xlsx':'R_EX_BIOMASS',
# r'readfile_ecoli_NMN_heterolactic_anaerobic_ACKko.xlsx':'R_EX_BIOMASS',

df_dict_anabolic_mceq = {}
result_path = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application\Result files'


# file_list = [file_list[0]]
for file, biomass_reac in file_dict.items():
    gluc =  sca.separate_cat_ana(file, all_cofactors=False)
    
    gluc.do_cat_ana_split()
    
    dfs = gluc.generate_dfs()
    
    result_df, atp_df, cofactor_df = dfs
    
    mceq_dict_biomass = gluc.generate_MCEQs(result_df, additional_active_metabolites=['M_atp_c', 'M_nadh_c', 'M_adp_c', 'M_nad_c', 'M_nadph_c', 'M_nadp_c'])
    df_mceqs = pd.DataFrame.from_records(mceq_dict_biomass, index=['overall']).T
    
    model = gluc.active_model.clone()
    # Example usage:
    fluxes, mceq, df_fva, kos = main(model, result_df, biomass_reac=biomass_reac)
    
    energy_obj = sca.energy_generation(model)
    if 'yeast' in file or 'synechocystis' in file:
        split_names = False
    else:
        split_names = True
    atp_consumption, atp_production, atp_special_consumption = energy_obj.calculate_energy(fluxes, split_names=split_names)
    nadh_prod, nadh_cons, transhydrogenase_prod, transhydrogenase_cons = energy_obj.calculate_NADH_prod(fluxes, split_names=split_names)
    nadph_prod, nadph_cons, transhydrogenase_prod, transhydrogenase_cons = energy_obj.calculate_NADPH_prod(fluxes, split_names=split_names)
    
    info_dict = {}
    info_dict['mceq'] = mceq
    info_dict['atp production'] = atp_production
    info_dict['atp consumption'] = atp_consumption
    info_dict['nadh production'] = nadh_prod
    info_dict['nadh consumption'] = nadh_cons
    info_dict['nadph production'] = nadph_prod
    info_dict['nadph consumption'] = nadph_cons
    
    with pd.ExcelWriter(result_path+"\\"+file[9:-5]+".xlsx") as writer:
        
        result_df.to_excel(writer, sheet_name='fluxes_separation1')
        atp_df.to_excel(writer, sheet_name="atp_separation1")
        cofactor_df.to_excel(writer, sheet_name="cofactor_separation1")
        df_mceqs.to_excel(writer, sheet_name="mceqs_separation1")
        
        
        # use to_excel function and specify the sheet_name and index 
        # to store the dataframe in specified sheet
        pd.Series(fluxes).to_excel(writer, sheet_name="fluxes_BP")
        pd.Series(info_dict).to_excel(writer, sheet_name="info_BP")
        df_fva.to_excel(writer, sheet_name="fva_BP")

        
        
    df_dict_anabolic_mceq[file[9:-5]] = str(mceq)