# -*- coding: utf-8 -*-
"""
Created on Tue Nov 18 15:27:12 2025

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

#%%

model.setReactionBounds('R_EDD', 0, 0)


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

cbmpy.CBWrite.writeSBML3FBCV2(model, savepath+'/iML1515_NMN_GAPNx_anaerobic.xml')