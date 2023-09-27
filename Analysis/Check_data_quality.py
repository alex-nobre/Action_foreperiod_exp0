# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 13:33:24 2022

@author: alpno
"""

import pandas as pd
import glob
import os

filesPath = '.\Data'

fileList=glob.glob(filesPath + './*.csv')
fileList.sort()

nFiles=int(len(fileList))

qualityTable = pd.DataFrame(columns=['participant','errorRateAll','errorRateGo','errorRateNoGo',
'errorRateAction', 'errorRateExternal', 'trimmedRTRate'])

for file in fileList:
     subData=pd.read_csv(file)
     
     # Remove unnecessary columns
     subData = subData[['participant', 'date', 'Response.corr', 'blockCondition', 'block', 'condition', 'foreperiod', 'corrAns', 'Response.rt', 'action_trigger.rt', 'Response.keys', 'counterbalance', 'extFixationDuration']]

     # Rename condition columns for clarity
     subData=subData.rename(columns={'condition':'trialType'})
     subData=subData.rename(columns={'blockCondition':'condition'})
     subData=subData.rename(columns={'Response.rt':'RT'})
     
     # Remove practice trials
     subData = subData[(subData['condition'] != 'practice') & (subData['condition'].notnull())]
    
    # ID
     participantID=set(subData['participant'].tolist()).pop()
     
     # Error rate
     errorRateAll=(len(subData[subData['Response.corr']==0])/len(subData)) * 100 # overall
     errorRateGo=(len(subData[(subData['Response.corr']==0)&(subData['trialType']=='go')])/len(subData[subData['trialType']=='go'])) * 100 # go trials
     errorRateNoGo=(len(subData[(subData['Response.corr']==0)&(subData['trialType']=='no-go')])/len(subData[subData['trialType']=='no-go'])) * 100 # no-go trials
     errorRateAction=(len(subData[(subData['Response.corr']==0)&(subData['condition']=='action')])/len(subData[subData['condition']=='action'])) * 100 # action condition
     errorRateExternal=(len(subData[(subData['Response.corr']==0)&(subData['condition']=='external')])/len(subData[subData['condition']=='external'])) * 100 # external condition
     
     # Percentage of trials excluded due to RT trimming
     subData=subData[(subData['trialType']=='go')&(subData['RT'].notnull())&(subData['Response.corr']==1)]
     trimmedRTs=subData[(subData['RT']>1)&(subData['RT']<0.1)]     
     trimmedRTRate=(len(trimmedRTs)/len(subData)) * 100
     
     subQualityData=[participantID,errorRateAll,errorRateGo,
                     errorRateNoGo,errorRateAction, errorRateExternal, trimmedRTRate]
     
     qualityTable.loc[len(qualityTable.index)]=subQualityData
 
qualityTable.to_csv('./Analysis/data_quality.csv')
