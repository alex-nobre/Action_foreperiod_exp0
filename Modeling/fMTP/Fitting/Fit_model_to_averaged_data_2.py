# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:02:49 2023

Changes: 
    - Included fitting for sequential effects

@author: alpno

Passos:
a.	Escolher uma métrica
    i.	Gerar valores de razão ativação e inibição (ver como Salet gera isso, incluindo a transformação linear) e gerar correlação com os RTS (não precisa transformar); minimizar 1-correlação (para maximizar correlação)
    ii.	Ou utilizar RMSE do z-score da razão ativação-inibição e do z-score do RT.
b.	Para cada conjunto de valores de parâmetro, gerar valores de k, r e c que minimizem o erro por sujeito
c.	Ou usar uma função de otimização pronta no código do Los


Salet et al. (2022) avaliam o ajuste do modelo utilizando R2.
"""

import numpy as np
import pandas as pd
import os
from scipy import optimize as op
from scipy import stats

os.chdir('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Fitting')

# Import for plotting
import matplotlib.pyplot as plt

# Import classes
from fmtp import fMTP, FPexp, FPgonogo
from hazard import fMTPhz
from fit import sort_fit, show_fit, get_fit 

#=================================== functions =============================

def temp2RT (prepValue, a, b):
    return (a*prepValue + b)

def sse (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    SSE = sum(pow(emp.RT - simRTs, 2))
    return SSE

def rmse (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    simzRTs = stats.zscore(simRTs)
    RMSE = np.sqrt(((emp.zRT - simzRTs)**2).mean())
    return RMSE

def anticorr (parm, emp, sim):
    simRTs = temp2RT(sim.prep, parm[0], parm[1])
    antiCorr = 1 - np.corrcoef(emp.RT, simRTs)[0][1]
    return antiCorr
    
    
#============================== Fit hazard =============================
# Dados empíricos: média por condição (externa vs ação) e por FP, apenas para valores 'go' (RTs em segundos)
empData = pd.read_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/dataActionFPAll.csv')

# Convert RT to ms
empData.RT=empData.RT*1000

# Get means of RTs and 
# empData = empData.groupby(['foreperiod','condition'],
#                             as_index=False)[['RT']].mean()
empData = empData.groupby(['foreperiod', 'oneBackFP', 'condition'],
                             as_index=False)[['RT']].mean()

empDataAv = empData.groupby(['foreperiod', 'condition'], as_index=False)[['RT']].mean()

# Create columns for z-score of RT in each data frame
empData['zRT'] = stats.zscore(empData.RT)
empDataAv['zRT'] = stats.zscore(empDataAv.RT)

# Split by condition
empDataAction = empData.loc[empData.condition=='action'].reset_index()
empDataExternal = empData.loc[empData.condition=='external'].reset_index()

empDataActionAv = empDataAv.loc[empDataAv.condition=='action'].reset_index()
empDataExternalAv = empDataAv.loc[empDataAv.condition=='external'].reset_index()


# Generate experiment from model
FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)
#exp = FPgonogo(FPs = FPs, distribution = distr, tr_per_block = 150, relax = 0)

# Parameters for simulations
kList=[i for i in range(1,11)]
rList=np.linspace(-4,-1,num=10)
cList=np.linspace(0,0.0003,num=10)

startVals=[4., 375.]

# Lists to store results
hazardResults=[None]*(len(kList)*len(rList)*len(cList))
seqEffResults=[None]*(len(kList)*len(rList)*len(cList))

# Simulate results by combinations of parameter values
fitEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        for ic in range(len(cList)):
            # Gerar modelo
            k=kList[ik]
            r=rList[ir]
            c=cList[ic]
            fmtp = fMTP(r, c, k)
            
            # Simulate experiment
            sim, prep_fmtp = exp.run_exp(fmtp)
            
            # Average preparation at discrete FPs
            sim = sim.iloc[1:, :] # first trial has no preparation
            sim = sim.groupby(['FP', 'FPn_1']).mean().reset_index()
            sim['distrib'] = distr
            simAv = sim.groupby('FP').mean().reset_index()
            
            # Fit model for hazard effect
               
            # Compute and minimize error measure
            # Action
            sseAction = op.minimize(sse, startVals, args = (empDataActionAv, simAv), method = 'L-BFGS-B')
            rmseAction = op.minimize(rmse, startVals, args = (empDataActionAv, simAv), method = 'L-BFGS-B')
            antiCorrAction = op.minimize(anticorr, startVals, args = (empDataActionAv, simAv), method = 'L-BFGS-B')
            RsseAction = np.corrcoef((sseAction.x[0] * simAv.prep + sseAction.x[1]), empDataActionAv.RT)[0][1]**2
            RrmseAction = np.corrcoef((rmseAction.x[0] * simAv.prep + rmseAction.x[1]), empDataActionAv.RT)[0][1]**2
            RantiCorrAction = np.corrcoef((antiCorrAction.x[0] * simAv.prep + antiCorrAction.x[1]), empDataActionAv.RT)[0][1]**2
            
            # External
            sseExternal = op.minimize(sse, startVals, args = (empDataExternalAv, simAv), method = 'L-BFGS-B')
            rmseExternal = op.minimize(rmse, startVals, args = (empDataExternalAv, simAv), method = 'L-BFGS-B')
            antiCorrExternal = op.minimize(anticorr, startVals, args = (empDataExternalAv, simAv), method = 'L-BFGS-B')
            RsseExternal = np.corrcoef((sseExternal.x[0] * simAv.prep + sseExternal.x[1]), empDataExternalAv.RT)[0][1]**2
            RrmseExternal = np.corrcoef((rmseExternal.x[0] * simAv.prep + rmseExternal.x[1]), empDataExternalAv.RT)[0][1]**2
            RantiCorrExternal = np.corrcoef((antiCorrExternal.x[0] * simAv.prep + antiCorrExternal.x[1]), empDataExternalAv.RT)[0][1]**2
            
            # Save to list
            hazardResults[fitEl]=(k, r, c,
                               sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                               rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                               antiCorrAction.x[0], antiCorrAction.x[1], antiCorrAction.fun, RantiCorrAction,
                               sseAction.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                               rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                               antiCorrExternal.x[0], antiCorrExternal.x[1], antiCorrExternal.fun, RantiCorrExternal)
                               
            # Fit model for sequential effect
               
            # Compute and minimize error measure
            # Action
            sseAction = op.minimize(sse, startVals, args = (empDataActionAv, simAv), method = 'L-BFGS-B')
            rmseAction = op.minimize(rmse, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
            antiCorrAction = op.minimize(anticorr, startVals, args = (empDataActionAv, simAv), method = 'L-BFGS-B')
            RsseAction = np.corrcoef((sseAction.x[0] * simAv.prep + sseAction.x[1]), empDataActionAv.RT)[0][1]**2
            RrmseAction = np.corrcoef((rmseAction.x[0] * sim.prep + rmseAction.x[1]), empDataAction.RT)[0][1]**2
            RantiCorrAction = np.corrcoef((antiCorrAction.x[0] * simAv.prep + antiCorrAction.x[1]), empDataActionAv.RT)[0][1]**2
            
            # External
            sseExternal = op.minimize(sse, startVals, args = (empDataExternalAv, simAv), method = 'L-BFGS-B')
            rmseExternal = op.minimize(rmse, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
            antiCorrExternal = op.minimize(anticorr, startVals, args = (empDataExternalAv, simAv), method = 'L-BFGS-B')
            RsseExternal = np.corrcoef((sseExternal.x[0] * simAv.prep + sseExternal.x[1]), empDataExternalAv.RT)[0][1]**2
            RrmseExternal = np.corrcoef((rmseExternal.x[0] * sim.prep + rmseExternal.x[1]), empDataExternal.RT)[0][1]**2
            RantiCorrExternal = np.corrcoef((antiCorrExternal.x[0] * simAv.prep + antiCorrExternal.x[1]), empDataExternalAv.RT)[0][1]**2
            
            # Save to list
            seqEffResults[fitEl]=(k, r, c,
                               sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                               rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                               antiCorrAction.x[0], antiCorrAction.x[1], antiCorrAction.fun, RantiCorrAction,
                               sseAction.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                               rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                               antiCorrExternal.x[0], antiCorrExternal.x[1], antiCorrExternal.fun, RantiCorrExternal)
            
            print(fitEl)
            fitEl+=1

# Join in a data frame
hazardResultsDF = pd.DataFrame(hazardResults,
                            columns = ['k', 'r', 'c', 
                                       'a SSE Action', 'b SSE Action', 'SSEAction', 'R2_SSE_Action', 
                                       'a RMSE Action', 'b RMSE Action', 'RMSE_Action', 'R2_RMSE_Action', 
                                       'a 1-corr Action', 'b 1-corr Action', 'One_minus_corr_Action', 'R2_1-corr_Action',
                                       'a SSE External', 'b SSE External', 'SSEExternal', 'R2_SSE_External', 
                                       'a RMSE External', 'b RMSE External', 'RMSE_External', 'R2_RMSE_External', 
                                       'a 1-corr External', 'b 1-corr External', 'One_minus_corr_External', 'R2_1-corr_External'])
                            #columns = ['k', 'r', 'c', 'a', 'b', 'SSE', 'R2'])
                            
seqEffResultsDF = pd.DataFrame(seqEffResults,
                            columns = ['k', 'r', 'c', 
                                       'a SSE Action', 'b SSE Action', 'SSEAction', 'R2_SSE_Action', 
                                       'a RMSE Action', 'b RMSE Action', 'RMSE_Action', 'R2_RMSE_Action', 
                                       'a 1-corr Action', 'b 1-corr Action', 'One_minus_corr_Action', 'R2_1-corr_Action',
                                       'a SSE External', 'b SSE External', 'SSEExternal', 'R2_SSE_External', 
                                       'a RMSE External', 'b RMSE External', 'RMSE_External', 'R2_RMSE_External', 
                                       'a 1-corr External', 'b 1-corr External', 'One_minus_corr_External', 'R2_1-corr_External'])
                            #columns = ['k', 'r', 'c', 'a', 'b', 'SSE', 'R2'])
                            

# Find row with best-fitting parameters
hazardResultsDF[hazardResultsDF['R2_RMSE_Action']==hazardResultsDF['R2_RMSE_Action'].max()]
seqEffResultsDF[seqEffResultsDF['R2_RMSE_Action']==seqEffResultsDF['R2_RMSE_Action'].max()]



    
    
    


    
    
    
    
    
    
    
    
    
    