# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:23:42 2022

@author: alpno
"""

import pandas as pd
import os

os.chdir('D:/Sync/Experimento_0_v2/Pavlovia/')


subFileName = input('Enter File Name:\n')
# if len(subnumber) < 2:
#     subFileName = '00' + str(subnumber) + '_action_fp_gabor_disc_2022-04-29_11h39.23.149.csv'
# else:
#     subFileName = '00' + str(subnumber) + '_action_fp_gabor_disc_2022-04-29_11h39.23.149.csv'

# Read 
subData = pd.read_csv(subFileName)

# Remove unnecessary columns
subData = subData[['participant', 'date', 'Response.corr', 'blockCondition', 'block', 'condition', 'foreperiod', 'corrAns', 'Response.rt', 'action_trigger.rt', 'Response.keys', 'counterbalance', 'extFixationDuration']]

# Rename condition columns for clarity
subData=subData.rename(columns={'condition':'trialType'})
subData=subData.rename(columns={'blockCondition':'condition'})
subData=subData.rename(columns={'Response.rt':'RT'})

# Create columns for n-1 and n-2 foreperiods by block
subData['oneBackFP'] = subData.groupby(['block'])['foreperiod'].shift(1)
subData['twoBackFP'] = subData.groupby(['block'])['foreperiod'].shift(2)

# Create columns for n-1 and n-2 trial types
subData['oneBacktrialType'] = subData.groupby(['block'])['trialType'].shift(1)
subData['twoBacktrialType'] = subData.groupby(['block'])['trialType'].shift(2)

# Remove practice trials and lines with nan values
subData = subData[(subData['condition'] != 'practice') & (subData['condition'].notnull())]

# Check n and % of errors
print(len(subData[subData['Response.corr'] == 0]))
print((len(subData[subData['Response.corr'] == 0])/len(subData)) * 100)

# Errors by condition
print(len(subData[(subData['trialType']=='no-go') & (subData['Response.corr']==0)])/
      len(subData[subData['trialType']=='no-go'])*100)

print(len(subData[(subData['trialType']=='go') & (subData['Response.corr']==0)]))
print(len(subData[(subData['trialType']=='go') & (subData['Response.corr']==0)])/
      len(subData[subData['trialType']=='go'])*100)


# Keep only go trials with correct responses to analyze RT
subData = subData[(subData['trialType'] == 'go') & (subData['RT'].notnull()) & (subData['Response.corr'] == 1)]

# Remove trials without n-1 FP values (i.e., first of each block)
subData = subData[subData['oneBackFP'].notnull()]

# Remove outliers
print('notclean: ' + str(len(subData)))
subData = subData[(subData['RT'] < 1) & (subData['RT'] > 0.1)]
print('clean: ' + str(len(subData)))



