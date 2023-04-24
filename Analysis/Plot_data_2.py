# -*- coding: utf-8 -*-
"""
Changes:
    - Renamed 'prevFP' to 'n-1_FP'
    - Included plots for n-2 FP
    - Renamed 'task' to 'trialType'
    
Created on Mon Aug 22 10:40:00 2022
@author: alpno
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning


import os
os.chdir('G:/My Drive/Post-doc/Experimento 1')

sns.set_context('paper')

# If we want to save data from individual participants
SavePlots=True

# Open data 
data = pd.read_csv('./Analysis/dataActionFPAll.csv')

# Keep only go trials with correct responses to analyze RT
dataForModel = data[(data['trialType'] == 'go') & (data['Response.rt'].notnull()) & (data['Response.corr'] == 1)]

# Remove trials without n-1 FP values (i.e., first of each block)
dataForModel = dataForModel[dataForModel['n-1_FP'].notnull()]
dataForModel = dataForModel[dataForModel['n-2_FP'].notnull()]

# Remove outliers
print('notclean: ' + str(len(dataForModel)))
dataForModel = dataForModel[(dataForModel['Response.rt'] < 1) & (dataForModel['Response.rt'] > 0.1)]
print('clean: ' + str(len(dataForModel)))

# Individual plots
if SavePlots==True:    
    unique_IDs=np.unique(dataForModel['ID'])
    for iSub,Sub in enumerate(unique_IDs):
        dataSub=dataForModel.loc[dataForModel['ID']==Sub]
        
        # Boxplots           
        plt.figure()
        box = sns.boxplot(x='foreperiod', y='Response.rt',hue='condition',data=dataSub,
                    palette=sns.color_palette(['grey','red']))
        
        
        figname='./Plots/Individual_plots/'+'ActionFP_boxplot'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Point plots
        plt.figure()
        point = sns.pointplot(x='foreperiod', y='Response.rt', hue='condition',
                              data=dataSub, palette=sns.color_palette(['blue', 'orange','green']))
        
        figname='./Plots/Individual_plots/'+'ActionFP_pointplot'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Sequential effect
        plt.figure()
        seqeff = sns.factorplot(x='foreperiod', y='Response.rt', hue='n-1_FP', data=dataSub,
                              col = 'condition', palette=sns.color_palette(['blue', 'orange','green']))
        
        figname='./Plots/Individual_plots/'+'ActionFP_seqeffect'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
    
    
summaryData=dataForModel.groupby(['ID','foreperiod','condition'],as_index=False)[['Response.rt']].mean()

# Summary boxplot
g=sns.boxplot(x="foreperiod", y="Response.rt",hue="condition",data=summaryData.sort_values(by='condition'),
                palette=sns.color_palette(['grey','red']))

# Summary pointplot
plt.figure()
summaryPlot =sns.pointplot(x="foreperiod", y="Response.rt",hue="condition",data=summaryData.sort_values(by='condition'),
                palette=sns.color_palette(['blue','orange','green']))

# Sequential effects
summaryData=dataForModel.groupby(['ID','foreperiod','condition','n-1_FP'],as_index=False)[['Response.rt']].mean()
plt.figure()
sns.factorplot(x='foreperiod', y='Response.rt', hue='n-1_FP', data=summaryData, col='condition',
               palette=sns.color_palette(['blue','orange','green']))


# n-2 sequential effects
summaryData=dataForModel.groupby(['ID','foreperiod','condition','n-2_FP'],as_index=False)[['Response.rt']].mean()
plt.figure()
sns.factorplot(x='foreperiod', y='Response.rt', hue='n-2_FP', data=summaryData, col='condition',
               palette=sns.color_palette(['blue','orange','green']))


#g.set_ylim([0.2, 0.5])

dataForModel['CueTypeNum']=-0.5
dataForModel['CueTypeNum'].loc[dataForModel['CueType']=='Valid']=0.5


md = smf.mixedlm("RT ~ CueTypeNum*TgOnsetCond", dataForModel, groups=dataForModel["ID"],re_formula="~CueTypeNum+TgOnsetCond")
mdf = md.fit(method=["powell","lbfgs","cg"])
print(mdf.summary())

