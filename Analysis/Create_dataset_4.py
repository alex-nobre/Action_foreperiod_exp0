# -*- coding: utf-8 -*-
"""
Changes:
    - Renamed n-1_FP to oneBackFP
    - Renamed n-2_FP to twoBckFP
    - Renamed n-1_trialType to oneBacktrialType
    - Renamed n-2_trialType to twoBacktrialType
    - Renamed n-1_effect to to oneBackEffect
    - Removed import of unused libraries
    
    
Created on Sep 30 18:24:00 2022
@author: Alexandre de Pontes Nobre
"""

import pandas as pd
import glob
import os

os.chdir('G:/My Drive/Post-doc/Experimento 1/Data/Pavlovia')

FileList=glob.glob('./*.csv')
FileList.sort()

nFiles=int(len(FileList))

dataActionFPAll=pd.DataFrame()

for iFile,FileName in enumerate(FileList):
    
    dataActionFP = pd.read_csv(FileName)  
    
    #Get info for this file
    ID=FileName[2:5]
         
    # Remove unnecessary columns
    dataActionFP = dataActionFP[['participant', 'date', 'Response.corr', 'blockCondition', 'block', 'condition', 'foreperiod', 'corrAns', 'Response.rt', 'action_trigger.rt', 'Response.keys', 'counterbalance', 'extFixationDuration']]
    
    # Rename columns for clarity
    dataActionFP = dataActionFP.rename(columns={'condition':'trialType'})
    dataActionFP = dataActionFP.rename(columns={'blockCondition':'condition'})
    dataActionFP = dataActionFP.rename(columns={'Response.rt':'RT'})
    
    # Remove practice trials
    dataActionFP = dataActionFP[(dataActionFP['condition'] != 'practice') & (dataActionFP['condition'].notnull())]
    
    # Create columns for n-1 and n-2 foreperiods by block
    dataActionFP['oneBackFP'] = dataActionFP.groupby(['block'])['foreperiod'].shift(1)
    dataActionFP['twoBackFP'] = dataActionFP.groupby(['block'])['foreperiod'].shift(2)
    
    # Create columns for n-1 and n-2 trial types
    dataActionFP['oneBacktrialType'] = dataActionFP.groupby(['block'])['trialType'].shift(1)
    dataActionFP['twoBacktrialType'] = dataActionFP.groupby(['block'])['trialType'].shift(2)
    
    # Compute trial-by-trial magnitude of sequential effects
    dataActionFP['oneBackEffect']=dataActionFP.groupby(['block'])['RT'].diff()

    # Replace participant's ID by three digits ID from file name
    dataActionFP['ID']=ID
    cols = dataActionFP.columns.tolist()
    cols = cols[-1:] + cols[1:-1]
    dataActionFP = dataActionFP[cols]
       
    dataActionFPAll=pd.concat([dataActionFPAll, dataActionFP], axis=0)    

dataActionFPAll.to_csv('G:/My Drive/Post-doc/Experimento 1/Analysis/'+'dataActionFPAll.csv')
