# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:46:34 2023

@author: alpno

Plots parameter values obtained from fitting process to check individual fitting.

Parameters for sequential effects are plotted in separate panel for each foreperiod for better visualization
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import itertools

import os

from fmtp import fMTP, FPexp, FPgonogo

def temp2RT (prepValue, a, b):
    return (a*prepValue + b)

hazardResultsDF = pd.read_csv('./Models/fMTP/Fitting/fittingResults.csv')

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEAction.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a SSE Action', 'b SSE Action', 'SSEAction']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').SSEExternal.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a SSE External', 'b SSE External', 'SSEExternal']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a SSE Action':'a SSE Action', 'b SSE Action':'b SSE Action', 
                                                                  'SSEAction':'SSE_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a SSE External':'a SSE External', 
                                                                    'b SSE External':'b SSE External', 
                                                                    'SSEExternal': 'SSE_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a SSE Action', 'b SSE Action', 'SSE_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a SSE External', 'b SSE External', 'SSE_External'])

# Plot SSE values
plt.scatter(stCoefs.index, stCoefs.SSE_Action)
plt.show()

# Check if they match (prob not)
fig = plt.figure()
ax = fig.subplots(3,2, sharey=True)
ax = ax.ravel()

titles = ['kAction', 'rAction', 'cAction', 'kExternal', 'rExternal', 'cExternal']

for i, panel in enumerate(ax):
    panel.hist(stCoefs.iloc[:,i+1])
    panel.set_title(titles[i])

plt.show()    

# Build plots with individual data and predictions
FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
exp = FPgonogo(FPs = FPs, distribution = distr, tr_per_block = 150, relax = 0)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure()
ax = fig.subplots(5,7,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, panel in enumerate(ax):
    sub = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = sub['kAction']
    rAction = sub['rAction']
    cAction = sub['cAction']
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, sub['a SSE Action'], sub['b SSE Action'])
    
    # External condition
    kExternal = sub['kExternal']
    rExternal = sub['rExternal']
    cExternal = sub['cExternal']
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, sub['a SSE External'], sub['b SSE External'])
    
    # Plot on single panel
    panel.plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    panel.plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    panel.plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    panel.plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    panel.set_title(titles[iSub])
    
fig.suptitle('ID fits')

figname='./Models/fMTP/Fitting/Plots/id_fit.png'
plt.savefig(figname,format='png')
plt.show()
    

# Plot by sequential effects
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True)
subfigs=fig.subfigures(5,7)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5,6]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    sub = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = sub['kAction']
    rAction = sub['rAction']
    cAction = sub['cAction']
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, sub['a SSE Action'], sub['b SSE Action'])
    
    # External condition
    kExternal = sub['kExternal']
    rExternal = sub['rExternal']
    cExternal = sub['cExternal']
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, sub['a SSE External'], sub['b SSE External'])
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,3,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='.Models/fMTP/Fitting/Plots/seq_effects_id_fit.png'
plt.savefig(figname,format='png')

#========================================================================================================================================#
# Based on RMSE

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').RMSE_Action.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a RMSE Action', 'b RMSE Action', 'RMSE_Action']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').RMSE_External.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a RMSE External', 'b RMSE External']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 'a RMSE Action':'a RMSE Action', 'b RMSE Action':'b RMSE Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 'a RMSE External':'a RMSE External', 'b RMSE External':'b RMSE External'})], axis = 1),
                       columns = ['ID', 'kAction', 'rAction', 'cAction', 'a RMSE Action', 'b RMSE Action',
                                  'kExternal', 'rExternal', 'cExternal', 'a RMSE External', 'b RMSE External'])

# Plot SSE values
plt.scatter(stCoefs.index, stCoefs.SSE_Action)
                                                   
#========================================================================================================================================#
# Based on 1-correlation

# Find min for SSE in each condition
stCoefsAction = hazardResultsDF.loc[hazardResultsDF.groupby('ID').One_minus_corr_Action.idxmin()].reset_index()
stCoefsAction = stCoefsAction[['ID', 'k', 'r', 'c', 'a 1-corr Action', 'b 1-corr Action', 'One_minus_corr_Action']]
stCoefsExternal = hazardResultsDF.loc[hazardResultsDF.groupby('ID').One_minus_corr_External.idxmin()].reset_index()
stCoefsExternal = stCoefsExternal[['k', 'r', 'c', 'a 1-corr External', 'b 1-corr External', 'One_minus_corr_External']]

stCoefs = pd.DataFrame(pd.concat([stCoefsAction.rename(columns = {'ID':'ID', 'k':'kAction', 'r':'rAction', 'c':'cAction', 
                                                                  'a 1-corr Action':'a 1-corr Action',
                                                                  'b 1-corr Action':'b 1-corr Action', 
                                                                  'One_minus_corr_Action':'One_minus_corr_Action'}), 
                                  stCoefsExternal.rename(columns = {'k':'kExternal', 'r':'rExternal', 'c':'cExternal', 
                                                                    'a 1-corr External':'a 1-corr External', 
                                                                    'b 1-corr External':'b 1-corr External', 
                                                                    'One_minus_corr_External': 'One_minus_corr_External'})], 
                                 axis = 1), columns = ['ID', 'kAction', 'rAction', 'cAction', 'a 1-corr Action', 'b 1-corr Action', 
                                                       'One_minus_corr_Action',
                                                       'kExternal', 'rExternal', 'cExternal', 'a 1-corr External', 'b 1-corr External',
                                                       'One_minus_corr_External'])

# Plot SSE values
plt.scatter(stCoefs.index, stCoefs.One_minus_corr_Action)
                                                
# Check if they match (prob not)
fig = plt.figure()
ax = fig.subplots(3,2, sharey=True)
ax = ax.ravel()

titles = ['kAction', 'rAction', 'cAction', 'kExternal', 'rExternal', 'cExternal']

for i, panel in enumerate(ax):
    panel.hist(stCoefs.iloc[:,i+1])
    panel.set_title(titles[i])
    
# Build plots with individual data and predictions
FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
exp = FPgonogo(FPs = FPs, distribution = distr, tr_per_block = 150, relax = 0)


# Plot by foreperiod only
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData = empData.groupby(['ID', 'foreperiod', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure()
ax = fig.subplots(5,7,sharey=False)
ax = ax.ravel()

titles = stCoefs.ID.tolist()

for iSub, panel in enumerate(ax):
    sub = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = sub['kAction']
    rAction = sub['rAction']
    cAction = sub['cAction']
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, sub['a 1-corr Action'], sub['b 1-corr Action'])
    
    # External condition
    kExternal = sub['kExternal']
    rExternal = sub['rExternal']
    cExternal = sub['cExternal']
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, sub['a 1-corr External'], sub['b 1-corr External'])
    
    # Plot on single panel
    panel.plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
    panel.plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
    panel.plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
    panel.plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')
    panel.set_title(titles[iSub])
    
fig.suptitle('ID fits')

figname='./ModelS/fMTP/Fitting/Plots/1-corr_id_fit.png'
plt.savefig(figname,format='png')
    

# Plot by sequential effects
empData = pd.read_csv('./Analysis/dataActionFPAll.csv')
empData = empData.groupby(['ID', 'foreperiod', 'oneBackFP', 'condition'], as_index=False)[['RT']].mean()

fig = plt.figure(constrained_layout=True)
subfigs=fig.subfigures(5,7)
subfigsCoords=list(itertools.product([0,1,2,3,4],[0,1,2,3,4,5,6]))

#subfigs = subfigs.ravel()
titles = stCoefs.ID.tolist()


for iSub in range(len(titles)):
    sub = stCoefs.iloc[iSub,:]
    empSub = empData[empData['ID']==titles[iSub]]
    empSubAction = empSub[empSub['condition']=='action']
    empSubExternal = empSub[empSub['condition']=='external']
    
    # Action condition
    kAction = sub['kAction']
    rAction = sub['rAction']
    cAction = sub['cAction']
    fmtpAction=fMTP(rAction,cAction,kAction)
    state_discr_action, state_con_action = exp.run_exp(fmtpAction)
    state_discr_action=state_discr_action[1:]
    mean_state_action=state_discr_action.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_action['RT'] = temp2RT(mean_state_action.prep, sub['a SSE Action'], sub['b SSE Action'])
    
    # External condition
    kExternal = sub['kExternal']
    rExternal = sub['rExternal']
    cExternal = sub['cExternal']
    fmtpExternal=fMTP(rExternal,cExternal,kExternal)
    state_discr_external, state_con_external = exp.run_exp(fmtpExternal)
    state_discr_external=state_discr_external[1:]
    mean_state_external=state_discr_external.groupby(['FP', 'FPn_1']).mean().reset_index()
    mean_state_external['RT'] = temp2RT(mean_state_external.prep, sub['a SSE External'], sub['b SSE External'])
    
    # Plot on single panel
    ax=subfigs[subfigsCoords[iSub]].subplots(1,3,sharey=True)
    
    for idx, iFP in enumerate(np.unique(mean_state_action.FPn_1)):
        FP = round(iFP, 2)
        state_n1_action = mean_state_action[mean_state_action.FPn_1 == iFP]
        emp_n1_action = empSubAction[empSubAction.oneBackFP == FP]
        ax[idx].plot(state_n1_action.FP, state_n1_action.RT, '.--', color = 'blue')
        ax[idx].plot(emp_n1_action.foreperiod, emp_n1_action.RT, '.-', color = 'blue')
        
        state_n1_external = mean_state_external[mean_state_external.FPn_1 == iFP]
        emp_n1_external = empSubExternal[empSubExternal.oneBackFP == FP]
        ax[idx].plot(state_n1_external.FP, state_n1_external.RT, '.--', color = 'orange')
        ax[idx].plot(emp_n1_external.foreperiod, emp_n1_external.RT, '.-', color = 'orange')
        
        ax[idx].set_title(FP)
    
    subfigs[subfigsCoords[iSub]].suptitle(titles[iSub])

fig.suptitle('Sequential effects')

figname='.Models/fMTP/Fitting/Plots/seq_effects_id_fit.png'
plt.savefig(figname,format='png')

#========================================================================================================================================#
#ax[iSub].plot(mean_state_action.FP, mean_state_action.RT,'.--', color = 'blue')
#ax[iSub].plot(empSubAction.foreperiod, empSubAction.RT, '.-', color = 'blue')
#ax[iSub].plot(mean_state_external.FP, mean_state_external.RT,'.--', color = 'orange')
#ax[iSub].plot(empSubExternal.foreperiod, empSubExternal.RT, '.-', color = 'orange')





