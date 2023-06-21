# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 11:21:35 2023

@author: alpno

Changes:
    - changed fmtp code to pass go/no-go trial sequence as argument to FPgonogo function


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
from fmtp_single_trial import fMTP, FPexp, FPgonogo
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
    

#============================ Fit single trial data =========================
empData = pd.read_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/dataActionFPAll.csv')

IDList=np.unique(empData['ID']).tolist()

#empData.RT=empData.RT*1000

# Create columns for z-score of RT in each data frame
empData['zRT'] = stats.zscore(empData.RT)

empDataAction = empData.loc[empData.condition=='action'].reset_index()
trialIndexAction = [[i for i in range(1,151)] for sub in range(int(len(empDataAction)/150))] #hard-coded n of trials, may change later
trialIndexAction = [index for ID in trialIndexAction for index in ID]
empDataAction['trial']=trialIndexAction

empDataExternal = empData.loc[empData.condition=='external'].reset_index()
trialIndexExternal = [[i for i in range(1,151)] for sub in range(int(len(empDataExternal)/150))] #hard-coded n of trials, may change later
trialIndexExternal = [index for ID in trialIndexExternal for index in ID]
empDataExternal['trial']=trialIndexExternal

# Generate experiment from model
FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
# exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 150)
#exp = FPgonogo(FPs = FPs, distribution = distr, tr_per_block = 150, relax = 0)

# k = 1 #4 # temporal precision
# r = -4 #-2.81 # rate of forgettting
# c = 0.0002 #0.0001 # memory persistence
# i = 0
# iID = 2

# Parameters for simulations
kList=[i for i in range(1,11)]
rList=np.linspace(-4,-1,num=10)
cList=np.linspace(0,0.0003,num=10)

startVals=[4., 375.]

# Lists to store results
hazardResults=[None]*(len(kList)*len(rList)*len(cList)*len(IDList))
seqEffResults=[None]*(len(kList)*len(rList)*len(cList)*len(IDList))

# Simulate results by combinations of parameter values
# testList=[None]*len(IDList)

# for i, iID in enumerate(IDList):
#     IDData=empData.loc[empData.ID==iID]
#     IDFPList=IDData['foreperiod'].values
#     print(IDFPList)

fitEl=0  
for i, iID in enumerate(IDList):
    print(i)
    
    # Action and external datasets for part
    IDDataAction=empDataAction.loc[empDataAction.ID==iID].reset_index()
    IDDataExternal=empDataExternal.loc[empDataExternal.ID==iID].reset_index()
    
    # FP lists for part
    IDFPAction=IDDataAction['foreperiod'].values.tolist()
    IDFPExternal=IDDataExternal['foreperiod'].values.tolist()
    
    # Trial type list for part
    IDTrialTypeAction=IDDataAction['trialType'].values.tolist()
    IDTrialTypeExternal=IDDataExternal['trialType'].values.tolist()
    
    # Remove first trial of each dataset for simulation
    IDDataAction = IDDataAction.iloc[1:,:]
    IDDataExternal = IDDataExternal.iloc[1:,:]
    
    # Remove no-go trials for estimation (replace dataframe name if goxno-go comparison is needed)
    IDDataAction = IDDataAction.loc[IDDataAction.trialType=='go']
    IDDataExternal = IDDataExternal.loc[IDDataExternal.trialType=='go']
    
    # Remove nan's
    actionNaIdx = IDDataAction[IDDataAction['RT'].isna()].index.tolist()
    IDDataAction = IDDataAction.dropna(subset=['RT'])
    
    externalNaIdx = IDDataExternal[IDDataExternal['RT'].isna()].index.tolist()
    IDDataExternal = IDDataExternal.dropna(subset=['RT'])
    
    # Exps for action and external conditions for part
    actionExp = FPgonogo(FPs = FPs, FPList = IDFPAction, distribution = distr, tr_per_block = 150, gngList = IDTrialTypeAction, relax = 0)
    externalExp = FPgonogo(FPs = FPs, FPList = IDFPExternal, distribution = distr, tr_per_block = 150, gngList = IDTrialTypeExternal, relax = 0)
     
    for ik in range(len(kList)):
        for ir in range(len(rList)):
            for ic in range(len(cList)):
                # Gerar modelo
                k=kList[ik]
                r=rList[ir]
                c=cList[ic]
                fmtp = fMTP(r, c, k)
                
                # Simulate experiments for action and external conditions
                simAction, prep_fmtpAction = actionExp.run_exp(fmtp)
                simExternal, prep_fmtpExternal = externalExp.run_exp(fmtp)
                
                # Average preparation at discrete FPs
                simAction = simAction.iloc[1:, :] # first trial has no preparation
                simExternal = simExternal.iloc[1:, :] # first trial has no preparation
                
                simAction['distrib'] = distr
                simExternal['distrib'] = distr
                
                
                # Remove no-go trials for estimation (replace dataframe name if goxno-go comparison is needed)
                # IDDataAction = IDDataAction.loc[IDDataAction.trialType=='go']
                # IDDataExternal = IDDataExternal.loc[IDDataExternal.trialType=='go']
                simAction = simAction.loc[simAction.gng=='go']
                simExternal = simExternal.loc[simExternal.gng=='go']
                
                # Remove nan's
                #actionNaIdx = IDDataAction[IDDataAction['RT'].isna()].index.tolist()
                simAction = simAction.drop(index=actionNaIdx)
                #IDDataActionTrim = IDDataAction.dropna(subset=['RT'])
                
                #externalNaIdx = IDDataExternal[IDDataExternal['RT'].isna()].index.tolist()
                simExternal = simExternal.drop(index=externalNaIdx)
                #IDDataExternal = IDDataExternal.dropna(subset=['RT'])
    
                # Fit model for hazard effect
                   
                # Compute and minimize error measure
                # Action
                sseAction = op.minimize(sse, startVals, args = (IDDataAction, simAction), method = 'L-BFGS-B')
                rmseAction = op.minimize(rmse, startVals, args = (IDDataAction, simAction), method = 'L-BFGS-B')
                antiCorrAction = op.minimize(anticorr, startVals, args = (IDDataAction, simAction), method = 'L-BFGS-B')
                RsseAction = np.corrcoef((sseAction.x[0] * simAction.prep + sseAction.x[1]), IDDataAction.RT)[0][1]**2
                RrmseAction = np.corrcoef((rmseAction.x[0] * simAction.prep + rmseAction.x[1]), IDDataAction.RT)[0][1]**2
                RantiCorrAction = np.corrcoef((antiCorrAction.x[0] * simAction.prep + antiCorrAction.x[1]), IDDataAction.RT)[0][1]**2
                
                # External
                sseExternal = op.minimize(sse, startVals, args = (IDDataExternal, simExternal), method = 'L-BFGS-B')
                rmseExternal = op.minimize(rmse, startVals, args = (IDDataExternal, simExternal), method = 'L-BFGS-B')
                antiCorrExternal = op.minimize(anticorr, startVals, args = (IDDataExternal, simExternal), method = 'L-BFGS-B')
                RsseExternal = np.corrcoef((sseExternal.x[0] * simExternal.prep + sseExternal.x[1]), IDDataExternal.RT)[0][1]**2
                RrmseExternal = np.corrcoef((rmseExternal.x[0] * simExternal.prep + rmseExternal.x[1]), IDDataExternal.RT)[0][1]**2
                RantiCorrExternal = np.corrcoef((antiCorrExternal.x[0] * simExternal.prep + antiCorrExternal.x[1]), IDDataExternal.RT)[0][1]**2
                
                # Save to list
                hazardResults[fitEl]=(iID, k, r, c,
                                   sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                                   rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                                   antiCorrAction.x[0], antiCorrAction.x[1], antiCorrAction.fun, RantiCorrAction,
                                   sseAction.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                                   rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                                   antiCorrExternal.x[0], antiCorrExternal.x[1], antiCorrExternal.fun, RantiCorrExternal)
                                   
                # Fit model for sequential effect
                   
                # # Compute and minimize error measure
                # # Action
                # sseAction = op.minimize(sse, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
                # rmseAction = op.minimize(rmse, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
                # antiCorrAction = op.minimize(anticorr, startVals, args = (empDataAction, sim), method = 'L-BFGS-B')
                # RsseAction = np.corrcoef((sseAction.x[0] * sim.prep + sseAction.x[1]), empDataAction.RT)[0][1]**2
                # RrmseAction = np.corrcoef((rmseAction.x[0] * sim.prep + rmseAction.x[1]), empDataAction.RT)[0][1]**2
                # RantiCorrAction = np.corrcoef((antiCorrAction.x[0] * sim.prep + antiCorrAction.x[1]), empDataAction.RT)[0][1]**2
                
                # # External
                # sseExternal = op.minimize(sse, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
                # rmseExternal = op.minimize(rmse, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
                # antiCorrExternal = op.minimize(anticorr, startVals, args = (empDataExternal, sim), method = 'L-BFGS-B')
                # RsseExternal = np.corrcoef((sseExternal.x[0] * sim.prep + sseExternal.x[1]), empDataExternal.RT)[0][1]**2
                # RrmseExternal = np.corrcoef((rmseExternal.x[0] * sim.prep + rmseExternal.x[1]), empDataExternal.RT)[0][1]**2
                # RantiCorrExternal = np.corrcoef((antiCorrExternal.x[0] * sim.prep + antiCorrExternal.x[1]), empDataExternal.RT)[0][1]**2
                
                # Save to list
                # seqEffResults[fitEl]=(k, r, c,
                #                    sseAction.x[0], sseAction.x[1], sseAction.fun, RsseAction, 
                #                    rmseAction.x[0], rmseAction.x[1], rmseAction.fun, RrmseAction,
                #                    antiCorrAction.x[0], antiCorrAction.x[1], antiCorrAction.fun, RantiCorrAction,
                #                    sseAction.x[0], sseExternal.x[1], sseExternal.fun, RsseExternal, 
                #                    rmseExternal.x[0], rmseExternal.x[1], rmseExternal.fun, RrmseExternal,
                #                    antiCorrExternal.x[0], antiCorrExternal.x[1], antiCorrExternal.fun, RantiCorrExternal)
                
                fitEl+=1


# Join in a data frame
hazardResultsDF = pd.DataFrame(hazardResults,
                            columns = ['ID', 'k', 'r', 'c', 
                                       'a SSE Action', 'b SSE Action', 'SSEAction', 'R2_SSE_Action', 
                                       'a RMSE Action', 'b RMSE Action', 'RMSE_Action', 'R2_RMSE_Action', 
                                       'a 1-corr Action', 'b 1-corr Action', 'One_minus_corr_Action', 'R2_1-corr_Action',
                                       'a SSE External', 'b SSE External', 'SSEExternal', 'R2_SSE_External', 
                                       'a RMSE External', 'b RMSE External', 'RMSE_External', 'R2_RMSE_External', 
                                       'a 1-corr External', 'b 1-corr External', 'One_minus_corr_External', 'R2_1-corr_External'])
                            #columns = ['k', 'r', 'c', 'a', 'b', 'SSE', 'R2'])
                            
# seqEffResultsDF = pd.DataFrame(seqEffResults,
#                             columns = ['k', 'r', 'c', 
#                                        'a SSE Action', 'b SSE Action', 'SSEAction', 'R2_SSE_Action', 
#                                        'a RMSE Action', 'b RMSE Action', 'RMSE_Action', 'R2_RMSE_Action', 
#                                        'a 1-corr Action', 'b 1-corr Action', 'One_minus_corr_Action', 'R2_1-corr_Action',
#                                        'a SSE External', 'b SSE External', 'SSEExternal', 'R2_SSE_External', 
#                                        'a RMSE External', 'b RMSE External', 'RMSE_External', 'R2_RMSE_External', 
#                                        'a 1-corr External', 'b 1-corr External', 'One_minus_corr_External', 'R2_1-corr_External'])
                            #columns = ['k', 'r', 'c', 'a', 'b', 'SSE', 'R2'])

hazardResultsDF.to_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Fitting/'+'fittingResults.csv')                            
                            