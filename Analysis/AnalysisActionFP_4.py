 # -*- coding: utf-8 -*-
"""
Changes:
    - Renamed from 'Plot_data_4' to 'AnalysisActionFP_4'
    - Renamed n-1_FP to oneBackFP
    - Renamed n-2_FP to twoBackFP
    - Renamed n-1_trialType to oneBacktrialType
    - Renamed n-2_trialType to twoBacktrialType
    - Renamed n-1_effect to to oneBackEffect
    - Included ANOVAS
    
Created on Sep 30 18:12:00 2022
@author: Alexandre de Pontes Nobre
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
import seaborn as sns
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.tools.sm_exceptions import ConvergenceWarning
from statsmodels.stats.anova import AnovaRM


import os
os.chdir('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_1')

sns.set_context('paper')

# If we want to save data from individual participants
SavePlots=False

#===================== 1. Prepare data ========================#
# Open data 
data = pd.read_csv('./Analysis/dataActionFPAll.csv')

# Keep only go trials with correct responses to analyze RT
#goData = data[(data['trialType'] == 'go') & (data['RT'].notnull()) & (data['Response.corr'] == 1)]
goData = data[(data['trialType'] == 'go') & (data['RT'].notnull())]

#======================= 2. Go data ===========================#
# Remove trials without n-1 FP values (i.e., first of each block)
goData = goData[goData['oneBackFP'].notnull()]
goData = goData[goData['twoBackFP'].notnull()]

# Remove outliers
goData = goData[(goData['RT'] < 1.0) & (goData['RT'] > 0.15)]


# Individual plots
if SavePlots==True:    
    unique_IDs=np.unique(goData['ID'])
    for iSub,Sub in enumerate(unique_IDs):
        dataSub=goData.loc[goData['ID']==Sub]
        
        # Boxplots           
        # plt.figure()
        # box = sns.boxplot(x='foreperiod', y='RT',hue='condition',data=dataSub,
        #             palette=sns.color_palette(['grey','red']))    
        
        # figname='./Plots/Individual_plots/'+'ActionFP_boxplot'+str(dataSub['ID'].iloc[0])+'.png'
        # plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Histograms of RTs
        # plt.figure()
        # histRT = sns.distplot(dataSub['RT'],bins=20,norm_hist=True)
        # histRT.set(xlim=(0,1.0))
        
        # figname='./Plots/Individual_plots/'+'ActionFP_RThistogram'+str(dataSub['ID'].iloc[0])+'.png'
        # plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Point plots
        # plt.figure()
        # point = sns.pointplot(x='foreperiod', y='RT', hue='condition',
        #                       data=dataSub, palette=sns.color_palette(['blue', 'orange','green']))
        
        # figname='./Plots/Individual_plots/'+'ActionFP_pointplot'+str(dataSub['ID'].iloc[0])+'.png'
        # plt.savefig(figname,format='png',bbox_inches='tight')
        
        # Sequential effect
        plt.figure()
        seqeff = sns.catplot(x='foreperiod', y='RT', hue='oneBackFP',kind='point',
                              data=dataSub,
                              palette=sns.color_palette(['blue', 'orange','green']))
        
        figname='./Plots/Individual_plots/'+'ActionFP_seqeffect'+str(dataSub['ID'].iloc[0])+'.png'
        plt.savefig(figname,format='png',bbox_inches='tight')
    

#========================== 2.2 Group plots ===================
summaryData=goData.groupby(['ID','foreperiod','condition','oneBackFP','oneBackAcc','oneBacktrialType','block','counterbalance'],
                           as_index=False)[['RT','Acc']].mean()


# 2.2.1. Check for learning of fatigue effects
sns.pointplot(x='block',y='RT',
              data=goData)

sns.pointplot(x='block',y='RT',hue='condition',
              data=goData)

sns.pointplot(x='block',y='RT',hue='counterbalance',
              data=goData)

sns.catplot(x='block',y='RT',hue='counterbalance',
            kind='point',col='foreperiod',
            data=goData)

# sns.pointplot(x='block',y='RT',hue='ID',
#               data=goData)

sns.pointplot(x='foreperiod',y='RT',hue='block',
              data=goData)

# 2.2.2. Check for order effects
sns.catplot(x='condition',y='RT',col='counterbalance',kind='point',
            data=goData)

sns.catplot(x='condition',y='RT',hue='ID',col='counterbalance',kind='point',
            data=goData)

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

sns.catplot(x='foreperiod',y='RT',
            kind='point',
            col='ID',col_wrap=4,
               data=goData)

sns.factorplot(x='foreperiod',y='RT',hue='ID',col='condition',
               data=goData)

sns.catplot(x='foreperiod',y='RT',hue='condition',kind='point',
            col='ID',col_wrap=4,
            data=goData)

# Summary histograms
sns.histplot(data=summaryData,x='RT',hue='foreperiod',
                         element='poly',
                         palette=sns.color_palette(['blue','orange','green']))
sns.histplot(data=summaryData[summaryData['condition']=='external'],
                         x='RT',hue='foreperiod',element='poly',
                         palette=sns.color_palette(['blue','orange','green']))
sns.histplot(data=summaryData[summaryData['condition']=='action'],
                         x='RT',hue='foreperiod',element='poly',
                         palette=sns.color_palette(['blue','orange','green']))

sns.catplot(x='RT',kind='violin',
             data=summaryData,col='condition',
             palette=sns.color_palette(['blue','orange','green']))


# Sequential effects
sns.catplot(x='foreperiod',y='RT',hue='oneBackFP',kind='point',
               data=summaryData,col='condition',
               palette=sns.color_palette(['blue','orange','green']))

sns.catplot(x='foreperiod',y='RT',hue='oneBackFP',kind='point',
               data=summaryData,
               palette=sns.color_palette(['blue','orange','green']))

# Effect of previous trial type on RT
trialTypePlot=sns.pointplot(x="oneBacktrialType", y="RT",
                             hue="condition",
                             data=summaryData.sort_values(by='condition'),
                             palette=sns.color_palette(['blue','orange','green']))

sns.pointplot(x="foreperiod", y="RT",
              hue="oneBacktrialType",
              data=summaryData.sort_values(by='condition'),
              palette=sns.color_palette(['lightblue','magenta']))

sns.barplot(x='oneBacktrialType',y='Acc',hue='condition',
            data=summaryData)

sns.catplot(x='oneBacktrialType', y='RT', hue='condition',
                             col='foreperiod',
                             kind='point',
                             data=summaryData.sort_values(by='condition'),
                             palette=sns.color_palette(['blue','orange','green']))

# By hit or miss on trial n-1
sns.catplot(x='oneBacktrialType',y='RT',hue='condition',
            col='oneBackAcc',
            kind='point',
            data=summaryData,
            palette=sns.color_palette(['blue','orange','green']))

sns.pointplot(x='oneBackAcc',y='RT',hue='condition',
            data=summaryData[summaryData['oneBacktrialType']=='no-go'],
            palette=sns.color_palette(['blue','orange','green']))

# Sequential effects by n-1 trial type
sns.catplot(x='foreperiod',y='RT',hue='oneBackFP',kind='point',
               data=summaryData[summaryData['condition']=='external'],
               col='oneBacktrialType',
               #row='condition',
               palette=sns.color_palette(['blue','orange','green']))

sns.catplot(x='foreperiod',y='RT',hue='oneBackFP',kind='point',
               data=summaryData,
               col='condition',
               palette=sns.color_palette(['blue','orange','green']))


# Magnitude of sequential effects by trial type
summaryData=goData.groupby(['ID','foreperiod','condition','oneBackFP','oneBackAcc','oneBacktrialType','block','counterbalance'],
                           as_index=False)[['RT']].mean()

summaryDataWide=summaryData.pivot(['ID','foreperiod','condition','oneBackFP'],
                                  'oneBacktrialType',
                                  'RT').reset_index()

summaryDataWide['gonogoDiff']=summaryDataWide['no-go']-summaryDataWide['go']

sns.barplot(x='foreperiod',y='conditionDiff',hue='oneBackFP',
               data=summaryDataWide,
               palette=sns.color_palette(['blue','orange','green']))


summaryDataWide['conditionDiff']=summaryDataWide['action']-summaryDataWide['external']
sns.pointplot(x='foreperiod',y='conditionDiff',
               data=summaryDataWide)


# Magnitude of difference between conditions
summaryData=goData.groupby(['ID','condition'],
                           as_index=False)[['RT']].mean()
summaryDataWide=summaryData.pivot('ID',
                                  'condition',
                                  'RT').reset_index()
summaryDataWide['conditionDiff']=summaryDataWide['action']-summaryDataWide['external']

histConditionDiff=sns.distplot(summaryDataWide['conditionDiff'],
                               norm_hist=True)


# 2.3. n-2 sequential effects
summaryData=goData.groupby(['ID','foreperiod','condition','oneBackFP','twoBackFP','oneBacktrialType'],
                           as_index=False)[['RT']].mean()

plt.figure()
sns.factorplot(x='foreperiod', y='RT', hue='twoBackFP', 
               data=summaryData, col='condition',
               palette=sns.color_palette(['blue','orange','green']))



#======================== 2.2. Both go and no-go data ===============================#
# # Remove trials without n-1 FP values (i.e., first of each block)
goNoGoData = data[data['oneBackFP'].notnull()]


# Plot Acc by condition
summaryData=goNoGoData.groupby(['ID','foreperiod','condition','trialType',
                                'oneBackFP','oneBacktrialType','block', 'counterbalance'],
                               as_index=False)[['RT', 'Acc']].mean()
sns.barplot(x='trialType',y='Acc',hue='condition',
            data=summaryData)

sns.barplot(x='oneBacktrialType',y='Acc',hue='condition',
            data=summaryData)

# Plot RT x Acc
RTAccPlot= sns.pointplot(x='Acc', y='RT', data=summaryData)
accValues=sorted(pd.Series.unique(summaryData['Acc']))
RTAccPlot.set_xticks(range(0,len(accValues),2))
RTAccPlot.set_xticklabels([round(item,2) for item in accValues[0:(len(accValues)-1):2]])


# FA rate
summaryData=goNoGoData.groupby(['ID'],as_index=False)[['RT', 'Response.corr']].mean()
summaryFARate=goNoGoData[goNoGoData['trialType']=='no-go'].groupby(['ID'], as_index=False)[['Response.corr']].mean()
summaryData=pd.merge(summaryData, summaryFARate, on='ID')
summaryData=summaryData.rename(columns={'Response.corr_x':'Acc', 
                                        'Response.corr_y':'FARate'})


# Plot RT x FARate
RTFARatePlot=sns.pointplot(x='FARate', y='RT', data=summaryData)
sns.factorplot(x='FARate',y='RT',hue='oneBackFP',data=summaryData,col='oneBacktrialType',row='condition',
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

#=============== 3. Inferential analysis ======================#
# ANOVAs for RT by condition and by FP
FPxRTAnova=AnovaRM(summaryData,depvar='RT', subject='ID', 
                     within=['foreperiod','condition'],
                     aggregate_func='mean')
AnovaOut=FPxRTAnova.fit()
AnovaOut.summary()

FPxRTAnova=AnovaRM(summaryData,depvar='RT', subject='ID', 
                     within=['foreperiod','condition','order'],
                     aggregate_func='mean')
AnovaOut=FPxRTAnova.fit()
AnovaOut.summary()


# ANOVAs for sequential effects (aggregated by condition)
seqEffectsAnova=AnovaRM(summaryData,depvar='RT', subject='ID', 
                     within=['foreperiod','oneBackFP'],
                     aggregate_func='mean')
AnovaOut=seqEffectsAnova.fit()
AnovaOut.summary()


# ANOVAs for RT by condition, FP and FPn-1
seqEffectsAnova=AnovaRM(summaryData,depvar='RT', subject='ID', 
                     within=['foreperiod','oneBackFP','condition'],
                     aggregate_func='mean')
AnovaOut=seqEffectsAnova.fit()
AnovaOut.summary()


# Effect of trialType on trial n-1 (go x no-go)
trialTypexconditionAnova=AnovaRM(summaryData,depvar='RT', subject='ID', 
                     within=['oneBacktrialType','condition'],
                     aggregate_func='mean')
AnovaOut=trialTypexconditionAnova.fit()
AnovaOut.summary()

# Effect of FP, FPn-1 and FP n-2
summaryData=goData.groupby(['ID','foreperiod','condition','oneBackFP','twoBackFP'],
                           as_index=False)[['RT']].mean()

seqModel=smf.ols('RT ~ C(foreperiod) + C(oneBackFP) + C(twoBackFP)', summaryData)
seqModelOut=seqModel.fit()
seqModelOut.summary()

md = smf.mixedlm("RT ~ CueTypeNum*TgOnsetCond", goData, groups=goData["ID"],re_formula="~CueTypeNum+TgOnsetCond")
mdf = md.fit(method=["powell","lbfgs","cg"])
print(mdf.summary())