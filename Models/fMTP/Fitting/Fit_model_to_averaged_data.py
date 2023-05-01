# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:02:49 2023

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

# Get means of RTs and # Create column for z-score of RT
empData = empData.groupby(['foreperiod','condition'],
                            as_index=False)[['RT']].mean()
empData['zRT'] = stats.zscore(empData.RT)

# Split by condition
empDataAction = empData.loc[empData.condition=='action'].reset_index()
empDataExternal = empData.loc[empData.condition=='external'].reset_index()


# Gerar RTs a partir dos valores
FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)
#exp = FPgonogo(FPs = FPs, distribution = distr, tr_per_block = 150, relax = 0)

# Criar listas de valores de k, r e c
k = 4 # temporal precision
r = -2.81 # rate of forgettting
c = 0.0001 # memory persistence

kList=[i for i in range(1,11)]  #np.linspace(2,8,num=10)
rList=np.linspace(-4,-1,num=10)
cList=np.linspace(0,0.0003,num=10)


startVals=[4., 375.]

hazardResults=[None]*(len(kList)*len(rList)*len(cList))

fitEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        for ic in range(len(cList)):
            # Gerar modelo
            k=kList[ik]
            r=rList[ir]
            c=cList[ic]
            fmtp = fMTP(r, c, k)
            
            # gerar RTs
            sim, prep_fmtp = exp.run_exp(fmtp)
            
            # Average preparation at discrete FPs
            sim = sim.iloc[1:, :] # first trial has no preparation
            sim = sim.groupby('FP').mean().reset_index()
            sim['distrib'] = distr
            
            # Average preparation curves
            prep = pd.DataFrame()
            prep['prep'] = np.mean(prep_fmtp.iloc[1:, :]) # again remove 1st trial
            prep['distrib'] = distr
            prep['FP'] = np.unique(sim.FP)[0]
            
            # Ajustar modelo
            
            # Generate RTs from preparation values
            medRT = temp2RT(sim.prep, 4, 375)
            
            # Compute and minimize sse
            # Action
            sseAction = op.minimize(sse, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
            rmseAction = op.minimize(rmse, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
            antiCorrAction = op.minimize(anticorr, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
            RsseAction = np.corrcoef((sseAction.x[0] * sim.prep + sseAction.x[1]), empDataAction.RT)[0][1]**2
            RrmseAction = np.corrcoef((rmseAction.x[0] * sim.prep + rmseAction.x[1]), empDataAction.RT)[0][1]**2
            RantiCorrAction = np.corrcoef((antiCorrAction.x[0] * sim.prep + antiCorrAction.x[1]), empDataAction.RT)[0][1]**2
            
            # External
            sseExternal = op.minimize(sse, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
            rmseExternal = op.minimize(rmse, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
            antiCorrExternal = op.minimize(anticorr, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
            RsseExternal = np.corrcoef((sseExternal.x[0] * sim.prep + sseExternal.x[1]), empDataExternal.RT)[0][1]**2
            RrmseExternal = np.corrcoef((rmseExternal.x[0] * sim.prep + rmseExternal.x[1]), empDataExternal.RT)[0][1]**2
            RantiCorrExternal = np.corrcoef((antiCorrExternal.x[0] * sim.prep + antiCorrExternal.x[1]), empDataExternal.RT)[0][1]**2
            
            # Save to list
            #fitResults[fitEl]=(k, r, c, fitAction.x[0], fitAction.x[1], fitAction.fun, RAction)
            hazardResults[fitEl]=(k, r, c,
                               sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                               rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                               antiCorrAction.x[0], antiCorrAction.x[1], antiCorrAction.fun, RantiCorrAction,
                               sseAction.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                               rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                               antiCorrExternal.x[0], antiCorrExternal.x[1], antiCorrExternal.fun, RantiCorrExternal)
            
            print(fitEl)
            fitEl+=1

# Join in a data frame
fitResultsDF = pd.DataFrame(fitResults,
                            columns = ['k', 'r', 'c', 'SSE', 'R2 SSE', 'RMSE', 'R2 SMSE', 'One minus corr', 'R2 1-corr'])
                            #columns = ['k', 'r', 'c', 'a', 'b', 'SSE', 'R2'])
                            

# Find row with best-fitting parameters
fitResultsDF[fitResultsDF['R2 SSE']==fitResultsDF['R2 SSE'].max()]
                            
#================================== fit sequential effects ============================


    
    
    


    
    
    
    
    
    
    
    
    
    