# -*- coding: utf-8 -*-
"""
Created on Fri Mar 21 10:46:43 2025

@author: mre283
"""

import cbmpy
import numpy as np
import scipy as sp
import os
import pandas as pd
import sys
import pickle
import cobra


#%%


modelPath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\89 Genome scale models\iML1515_macroBOF_core.xml'
model = cbmpy.loadModel(modelPath)
model.setReactionBounds('R_POR5', 0, 1000)

#%%

# Alternative redox cofactors
# model.createSpecies('M_nmn_c', compartment='c', name='oxidized nmn')
model.createSpecies('M_nmnh_c', compartment='c', name='reduced nmn')


#%%

model.createReaction('R_GAPDz', reversible=True, name='Glyceraldehyde-3-phosphate dehydrogenase NRC coupled')
model.createReactionReagent('R_GAPDz', metabolite='M_nmn_c', coefficient=-1)
model.createReactionReagent('R_GAPDz', metabolite='M_g3p_c', coefficient=-1)
model.createReactionReagent('R_GAPDz', metabolite='M_pi_c', coefficient=-1)
model.createReactionReagent('R_GAPDz', metabolite='M_13dpg_c', coefficient=1)
model.createReactionReagent('R_GAPDz', metabolite='M_nmnh_c', coefficient=1)
model.createReactionReagent('R_GAPDz', metabolite='M_h_c', coefficient=1)
model.setReactionBounds('R_GAPDz', lower=-1000, upper=1000)

# LDH with alternative redox factor
model.createReaction('R_LDH_Dz', reversible=True, name='Lactate dehydrogenase NRC coupled')
model.createReactionReagent('R_LDH_Dz', metabolite='M_lac__D_c', coefficient=-1)
model.createReactionReagent('R_LDH_Dz', metabolite='M_nmn_c', coefficient=-1)
model.createReactionReagent('R_LDH_Dz', metabolite='M_pyr_c', coefficient=1)
model.createReactionReagent('R_LDH_Dz', metabolite='M_nmnh_c', coefficient=1)
model.createReactionReagent('R_LDH_Dz', metabolite='M_h_c', coefficient=1)
model.setReactionBounds('R_LDH_Dz', lower=-1000, upper=1000)


#%%

model.setReactionBounds('R_EX_o2_e', 0, 0)
model.setReactionBounds('R_GAPD', 0, 0)
model.setReactionBounds('R_LDH_D', 0, 0)
model.setReactionBounds('R_LDH_D2', 0, 0)
model.getReaction('R_GLYCL').reversible = True
model.getReaction('R_EX_for_e').reversible = True

result =  cbmpy.CBGLPK.glpk_analyzeModel(model)

fluxes = model.getReactionValues()

#%%

savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application'

cbmpy.CBWrite.writeSBML3FBCV2(model, savepath+'/iML1515_NMN_ED_anaerobic.xml')

#%%

model.setReactionBounds('R_EDD', 0, 0)


model.createReaction('R_GAPN', name='glyceraldehyde-3-phosphate dehydroganes, non-phosphorylating', reversible=False)
model.createReactionReagent('R_GAPN', 'M_g3p_c', -1.0)
model.createReactionReagent('R_GAPN', 'M_h2o_c', -1.0)
model.createReactionReagent('R_GAPN', 'M_nadp_c', -1.0)
model.createReactionReagent('R_GAPN', 'M_3pg_c', 1.0)
model.createReactionReagent('R_GAPN', 'M_nadph_c', 1.0)
model.createReactionReagent('R_GAPN', 'M_h_c', 1.0)
model.setReactionBounds('R_GAPN', 0, 1000)

model.createReaction('R_GAPNx', name='glyceraldehyde-3-phosphate dehydroganes, non-phosphorylating', reversible=False)
model.createReactionReagent('R_GAPNx', 'M_g3p_c', -1.0)
model.createReactionReagent('R_GAPNx', 'M_h2o_c', -1.0)
model.createReactionReagent('R_GAPNx', 'M_nad_c', -1.0)
model.createReactionReagent('R_GAPNx', 'M_3pg_c', 1.0)
model.createReactionReagent('R_GAPNx', 'M_nadh_c', 1.0)
model.createReactionReagent('R_GAPNx', 'M_h_c', 1.0)
model.setReactionBounds('R_GAPNx', 0, 1000)


result =  cbmpy.CBGLPK.glpk_analyzeModel(model)

fluxes = model.getReactionValues()

#%%

savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application'

cbmpy.CBWrite.writeSBML3FBCV2(model, savepath+'/iML1515_NMN_GAPN_anaerobic.xml')

#%%

model.setReactionBounds('R_GAPN', 0, 0)

model.createReaction('R_XU5P_PKT', reversible=False, name='Phosphoketolase (xu5p utilizing)')
model.createReactionReagent('R_XU5P_PKT', metabolite='M_xu5p__D_c', coefficient=-1)
model.createReactionReagent('R_XU5P_PKT', metabolite='M_pi_c', coefficient=-1)
model.createReactionReagent('R_XU5P_PKT', metabolite='M_actp_c', coefficient=1)
model.createReactionReagent('R_XU5P_PKT', metabolite='M_g3p_c', coefficient=1)
model.createReactionReagent('R_XU5P_PKT', metabolite='M_h2o_c', coefficient=1)
model.setReactionBounds('R_XU5P_PKT', lower=0, upper=1000)


result =  cbmpy.CBGLPK.glpk_analyzeModel(model)

fluxes = model.getReactionValues()

#%%

savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application'

cbmpy.CBWrite.writeSBML3FBCV2(model, savepath+'/iML1515_NMN_heterolactic_anaerobic.xml')

#%%

model.setReactionBounds('R_XU5P_PKT', 0, 0)

model.createSpecies('M_rb15bp_c', compartment='c', name='D-Ribulose-1,5-bisphosphate', chemFormula='C5H8O11P2')

model.createReaction('R_PRUK', reversible=False, name='Phosphoribulokinase')
model.createReactionReagent('R_PRUK', metabolite='M_ru5p__D_c', coefficient=-1)
model.createReactionReagent('R_PRUK', metabolite='M_atp_c', coefficient=-1)
model.createReactionReagent('R_PRUK', metabolite='M_rb15bp_c', coefficient=1)
model.createReactionReagent('R_PRUK', metabolite='M_adp_c', coefficient=1)
model.createReactionReagent('R_PRUK', metabolite='M_h_c', coefficient=1)
model.setReactionBounds('R_PRUK', lower=0, upper=1000)

model.createReaction('R_RBPC', reversible=False, name='Ribulose-bisphosphate carboxylase')
model.createReactionReagent('R_RBPC', metabolite='M_rb15bp_c', coefficient=-1)
model.createReactionReagent('R_RBPC', metabolite='M_co2_c', coefficient=-1)
model.createReactionReagent('R_RBPC', metabolite='M_h2o_c', coefficient=-1)
model.createReactionReagent('R_RBPC', metabolite='M_3pg_c', coefficient=2)
model.createReactionReagent('R_RBPC', metabolite='M_h_c', coefficient=2)
model.setReactionBounds('R_RBPC', lower=0, upper=1000)


result = cbmpy.CBGLPK.glpk_analyzeModel(model)

fluxes = model.getReactionValues()

#%%

savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application'

cbmpy.CBWrite.writeSBML3FBCV2(model, savepath+'/iML1515_NMN_calvin_anaerobic.xml')


#%%

savepath = r'C:\Users\mre283\OneDrive - Vrije Universiteit Amsterdam\01 Separating catabolism and anabolism\python files\generalized_scripts\Metabolic engineering application'

model = cbmpy.loadModel(savepath+'/iML1515_NMN_calvin_anaerobic.xml')



# model.setReactionBounds('R_EX_o2_e', 0, 0)
model.setReactionBounds('R_GAPDz', 0, 0)
model.setReactionBounds('R_LDH_Dz', 0, 0)
model.setReactionBounds('R_GAPD', 0, 0)
model.setReactionBounds('R_PGI', 0, 0)
model.setReactionBounds('R_HEX7', 0, 0)
model.setReactionBounds('R_DHAPT', 0, 0)


result = cbmpy.CBGLPK.glpk_analyzeModel(model)

fluxes = model.getReactionValues()



