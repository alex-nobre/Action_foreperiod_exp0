# -*- coding: utf-8 -*-
"""
Created on Fri Dec  9 11:58:22 2022

@author: alpno
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import seaborn as sns

import os
os.chdir("G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling")

from fmtp import fMTP, FPexp, FPgonogo

#====================================================================================    
#============================= Using slope function =================================
#====================================================================================    
# Function to compute slope
def slope(x1, x2, y1, y2):
    m=(y2-y1)/(x2-x1)
    return m

# Initialize values for simulation
FP = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
# Set-up experiment using "FPexp"
exp = FPexp(FPs = FP, distribution = distr, tr_per_block = 150)
#exp = FPgonogo(FPs = FP, distribution = distr, tr_per_block = 150, relax = 0)      
                                              

# create list of parameter values
kList=[i for i in range(1,11)]  #np.linspace(2,8,num=10)
rList=np.linspace(1,4,num=10)
cList=np.linspace(0,0.0003,num=10)


# Create list to store outputs:
simResults=[None]*(len(kList)*len(rList)*len(cList))   
    
simEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        for ic in range(len(cList)):
            k=kList[ik]
            r=rList[ir]
            c=cList[ic]
            fmtp=fMTP(r,c,k)
            state_discr, state_con = exp.run_exp(fmtp)
            state_discr=state_discr[1:]
            
            mean_state=state_discr.groupby(['FP']).mean().reset_index()

            lowFP=min(np.unique(mean_state.FP))
            highFP=max(np.unique(mean_state.FP))
            
            # Hazard slope
            hazardSlope=slope(lowFP,
                              highFP,
                              mean_state.loc[mean_state['FP']==lowFP,'prep'].values[0],
                              mean_state.loc[mean_state['FP']==highFP,'prep'].values[0])
            
            #### Sequential effects difference of slopes
            mean_state = state_discr.groupby(['FP', 'FPn_1']).mean().reset_index()
            lowFPn_1=min(np.unique(mean_state.FPn_1))
            highFPn_1=max(np.unique(mean_state.FPn_1))
            
            # low FP n-1
            state_lowFPn_1=mean_state[mean_state.FPn_1==lowFPn_1]
            lowSeqSlope=slope(lowFP,
                              highFP,
                              state_lowFPn_1.loc[state_lowFPn_1['FP']==lowFP,'prep'].values[0],
                              state_lowFPn_1.loc[state_lowFPn_1['FP']==highFP,'prep'].values[0])
            
            # High FP n-1
            state_highFPn_1=mean_state[mean_state.FPn_1==highFPn_1]
            highSeqSlope=slope(lowFP,
                               highFP,
                               state_highFPn_1.loc[state_highFPn_1['FP']==lowFP,'prep'].values[0],
                               state_highFPn_1.loc[state_highFPn_1['FP']==highFP,'prep'].values[0])
            
            seqSlopeDiff=highSeqSlope-lowSeqSlope
            simResults[simEl]=(hazardSlope,seqSlopeDiff)
            print(simEl)
            simEl+=1


plt.scatter(*zip(*simResults))
plt.show()

# By parameters value
chunkedList=list()
chunkSize=int(len(simResults)/len(kList))

for i in range(0,len(simResults),chunkSize):
    chunkedList.append(simResults[i:i+chunkSize])
    
 
fig,ax=plt.subplots(nrows=2,ncols=5)

kSublist=0
for row in ax:
    for col in row:
        plotSublist=chunkedList[kSublist]
        col.scatter(*zip(*plotSublist))
        kSublist+=1
        
        
fig,ax=plt.subplots(nrows=2,ncols=5)
        
chunkedChunkedList=list()
chunkSize=int(len(cList))
kChunk=chunkedList[0]

for i in range(0,len(kChunk),chunkSize):
    chunkedChunkedList.append(kChunk[i:i+chunkSize])

rSublist=0
for row in ax:
    for col in row:
        plotSublist=chunkedChunkedList[rSublist]
        col.scatter(*zip(*plotSublist))
        rSublist+=1

#===================================================================================#
#===================================================================================#
#===================================================================================#

################## Sketches ###########################


k = 4
r = -2.81 # rate of forgetting
c = 1e-4 # memory persistence

# loop through the parameters and fit model to obtain curves
# Initialize fmtp class
fmtp=fMTP(r,c,k)

# Run experiment using object "exp" and "fmtp"
state_discr, state_con = exp.run_exp(fmtp)

state_discr=state_discr[1:]

#=============== extract indices for hazard and sequential effects =================#
mean_state=state_discr.groupby(['FP']).mean().reset_index()

lowFP=min(np.unique(mean_state.FP))
highFP=max(np.unique(mean_state.FP))

### compute slope of hazard curve

#hazardLm=smf.ols(formula='prep~FP',data=mean_state).fit()
hazardLm=smf.ols(formula='prep~FP',data=state_discr).fit()
hazardSlope=round(hazardLm.params[1],4)

### compute difference of slopes for FPn-1
mean_state = state_discr.groupby(['FP', 'FPn_1']).mean().reset_index()


lowFPn_1=min(np.unique(mean_state.FPn_1))
highFPn_1=max(np.unique(mean_state.FPn_1))


state_lowFPn_1=state_discr[state_discr.FPn_1==lowFPn_1]
lowSeqLm=smf.ols(formula='prep~FP',data=state_lowFPn_1).fit()
lowSeqSlope=lowSeqLm.params[1]

state_highFPn_1=state_discr[state_discr.FPn_1==highFPn_1]
highSeqLm=smf.ols(formula='prep~FP',data=state_highFPn_1).fit()
highSeqSlope=highSeqLm.params[1]

seqSlopeDiff=highSeqSlope-lowSeqSlope

# Use slope formula (would need to adapt to select cells automatically)
#hazardSlope=slope(0.6,1.8,7.71,2.47) # hazard

state_lowFPn_1=mean_state[mean_state.FPn_1==lowFPn_1] # sequential
lowSeqSlope=slope(lowFP,
                  highFP,
                  state_lowFPn_1.loc[state_lowFPn_1['FP']==lowFP,'prep'].values[0],
                  state_lowFPn_1.loc[state_lowFPn_1['FP']==highFP,'prep'].values[0]) #slope(0.6,1.8,3.18,1.46)

state_highFPn_1=mean_state[mean_state.FPn_1==highFPn_1]
highSeqSlope=slope(lowFP,
                   highFP,
                   state_highFPn_1.loc[state_highFPn_1['FP']==lowFP,'prep'].values[0],
                   state_highFPn_1.loc[state_highFPn_1['FP']==highFP,'prep'].values[0]) # slope(0.6,1.8,12.90,3.55)



