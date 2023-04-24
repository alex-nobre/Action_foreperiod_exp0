 # -*- coding: utf-8 -*-
"""
Changes:
    - Created sections for data import, descriptive analysis
    - Compute trial-by-trial magnitude 
    - Included plots by trialType
    - Increased low cutoff limit for RTs to 150 ms, as per Los and Agter (2005), Los (2013) and Los et al. (2017)
    - Included plots for n-2 sequential effects and magnitude of sequential effect
    by trialType
    - Included FP x RT plots by participant
    - Include Acc x RT plots to check for speed-accuracy trade-off
    
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
SavePlots=False

#===================== 1. Prepare data ========================#
# Open data 
data = pd.read_csv('./Analysis/dataActionFPAll.csv')

# Keep only go trials with correct responses to analyze RT
goData = data[(data['trialType'] == 'go') & (data['RT'].notnull()) & (data['Response.corr'] == 1)]

#======================= 2. Go data ===========================#
# Remove trials without n-1 FP values (i.e., first of each block)
goData = goData[goData['n-1_FP'].notnull()]
#goData = goData[goData['n-2_FP'].notnull()]

# Remove outliers
goData = goData[(goData['RT'] < 1.0) & (goData['RT'] > 0.15)]


# Individual plots
if SavePlots==True:    
    unique_IDs=np.unique(goData['ID'])
    for iSub,Sub in enumerate(unique_IDs):
        dataSub=goData.loc[goData['ID']==Sub]
        
        # Boxplots           
        plt.figure()
        box = sns.boxplot(x='foreperiod', y='RT',hue='condition',data=dataSub,
                    palette=sns.color_palette(['grey','red']))    
        
        figname='./Plots/Individual_plots/'+'ActionFP_boxplot'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Histograms of RTs
        plt.figure()
        histRT = sns.distplot(dataSub['RT'],bins=20,norm_hist=True)
        
        figname='./Plots/Individual_plots/'+'ActionFP_RThistogram'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Point plots
        plt.figure()
        point = sns.pointplot(x='foreperiod', y='RT', hue='condition',
                              data=dataSub, palette=sns.color_palette(['blue', 'orange','green']))
        
        figname='./Plots/Individual_plots/'+'ActionFP_pointplot'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Sequential effect
        plt.figure()
        seqeff = sns.factorplot(x='foreperiod', y='RT', hue='n-1_FP', data=dataSub,
                              col = 'condition', palette=sns.color_palette(['blue', 'orange','green']))
        
        figname='./Plots/Individual_plots/'+'ActionFP_seqeffect'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
    

#========================== 2.2 Group plots ===================
summaryData=goData.groupby(['ID','foreperiod','condition','n-1_FP','n-1_trialType'],
                           as_index=False)[['RT']].mean()


# Summary boxplot
g=sns.boxplot(x="foreperiod", y="RT",hue="condition",data=summaryData.sort_values(by='condition'),
                palette=sns.color_palette(['grey','red']))

# Summary FP x RT plot
plt.figure()
summaryPlot=sns.pointplot(x="foreperiod", y="RT",hue="condition",
                          data=summaryData.sort_values(by='condition'),
                palette=sns.color_palette(['blue','orange','green']))

# FP x RT plot by participant
sns.pointplot(x='foreperiod',y='RT',hue='ID',
               data=goData)

sns.factorplot(x='foreperiod',y='RT',col='ID',col_wrap=4,
               data=goData)

# Sequential effects
summaryData=goData.groupby(['ID','foreperiod','condition','n-1_FP'],as_index=False)[['RT']].mean()
sns.factorplot(x='foreperiod',y='RT',hue='n-1_FP',
               data=summaryData,col='condition',
               palette=sns.color_palette(['blue','orange','green']))

# Effect of previous trial type
trialTypePlot =sns.pointplot(x="n-1_trialType", y="RT",
                             #hue="condition",
                             data=summaryData.sort_values(by='condition'),
                             palette=sns.color_palette(['blue','orange','green']))

trialTypePlot=sns.barplot(x="n-1_trialType", y="RT",
                             #hue="condition",
                             data=summaryData.sort_values(by='condition'),
                             palette=sns.color_palette(['blue','orange','green']))

# Sequential effects by trial type
sns.factorplot(x='foreperiod',y='RT',hue='n-1_FP',
               data=summaryData,col='n-1_trialType',
               #row='condition',
               palette=sns.color_palette(['blue','orange','green']))


# Magnitude of sequential effects by trial type
summaryDataWide=summaryData.pivot(['ID','foreperiod','condition','n-1_FP'],
                                  'n-1_trialType',
                                  'RT').reset_index()

summaryDataWide['gonogoDiff']=summaryDataWide['no-go']-summaryDataWide['go']

sns.barplot(x='foreperiod',y='gonogoDiff',hue='n-1_FP',
               data=summaryDataWide,
               palette=sns.color_palette(['blue','orange','green']))

# 2.3. n-2 sequential effects
summaryData=goData.groupby(['ID','foreperiod','condition','n-1_FP','n-2_FP','n-1_trialType'],
                           as_index=False)[['RT']].mean()

plt.figure()
sns.factorplot(x='foreperiod', y='RT', hue='n-2_FP', 
               data=summaryData, col='condition',
               palette=sns.color_palette(['blue','orange','green']))



#======================== 2.2. Both go and no-go data ===============================#
# # Remove trials without n-1 FP values (i.e., first of each block)
goNoGoData = data[data['n-1_FP'].notnull()]

summaryData=goNoGoData.groupby(['ID'],as_index=False)[['RT', 'Response.corr']].mean()
summaryFARate=goNoGoData[goNoGoData['trialType']=='no-go'].groupby(['ID'], as_index=False)[['Response.corr']].mean()
summaryData=pd.merge(summaryData, summaryFARate, on='ID')
summaryData=summaryData.rename(columns={'Response.corr_x':'Acc', 
                                        'Response.corr_y':'FARate'})


# Plot RT x ACC
RTAccPlot= sns.pointplot(x='Acc', y='RT', data=summaryData)
accValues=sorted(pd.Series.unique(summaryData['Acc']))
RTAccPlot.set_xticks(range(0,len(accValues),2))
RTAccPlot.set_xticklabels([round(item,2) for item in accValues[0:(len(accValues)-1):2]])

# Plot RT x FARate
RTFARatePlot=sns.pointplot(x='FARate', y='RT', data=summaryData)
sns.factorplot(x='FARate',y='RT',hue='n-1_FP',data=summaryData,col='n-1_trialType',row='condition',
               palette=sns.color_palette(['blue','orange','green']))

# Acc and false alarms per condition (action x sequential)
summaryData=goNoGoData.groupby(['ID','condition'],as_index=False)[['RT', 'Response.corr']].mean()
summaryFARate=goNoGoData[goNoGoData['trialType']=='no-go'].groupby(['ID','condition'], as_index=False)[['Response.corr']].mean()
summaryData=pd.merge(summaryData, summaryFARate, on=['ID','condition'])
summaryData=summaryData.rename(columns={'Response.corr_x':'Acc', 
                                        'Response.corr_y':'FARate'})


sns.pointplot(x='condition',y='Acc',data=summaryData)
sns.pointplot(x='condition',y='FARate',data=summaryData)


plt.figure()
sns.pointplot(x='condition',y='Acc',data=summaryData)
sns.pointplot(x='condition',y='FARate',data=summaryData)