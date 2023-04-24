# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 16:29:45 2022

@author: Alexandre Nobre
"""

import pandas as pd
import numpy as np
import numpy.matlib
import matplotlib.pyplot as plt
from os import listdir
import glob
from scipy import stats

import seaborn as sns


from statsmodels.formula.api import ols


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
    
    # Rename condition columns for clarity
    dataActionFP = dataActionFP.rename(columns={'condition':'task'})
    dataActionFP = dataActionFP.rename(columns={'blockCondition':'condition'})
    
    # Create column for n-1 foreperiod by block
    dataActionFP['prevFP'] = dataActionFP.groupby(['block'])['foreperiod'].shift(1)
    
    # Remove practice trials
    dataActionFP = dataActionFP[dataActionFP['condition'] != 'practice']
    
    # Replace participant's ID by three digits ID from file name
    dataActionFP['ID']=ID
    cols = dataActionFP.columns.tolist()
    cols = cols[-1:] + cols[1:-1]
    dataActionFP = dataActionFP[cols]
       
    dataActionFPAll=pd.concat([dataActionFPAll, dataActionFP], axis=0)    

dataActionFPAll.to_csv('G:/My Drive/Post-doc/Experimento 1/Analysis/'+'dataActionFPAll.csv')
