# -*- coding: utf-8 -*-
"""
Created on Wed Jan 25 14:24:41 2023

@author: alpno

Simulate several ways in which internal references may be created which can
generate sequential effects:
    - Weighted average:
        - geometric average
        - Others
    - Individual trials
"""


import matplotlib.pyplot as plt
import numpy as np

#=============================================================================
#=================== 1. Point values of means (no noise) =====================
#=============================================================================
# This is the typical experimental context, where experimental parameters are
# randomly draw in each trial
FPs = np.arange(1.0, 3.0, 0.6)

nTrials = np.random.randint(low=1, high=6, size=999)
nTrials = np.random.choice(FPs, size=999)

#======================= Equal weights moving average =========================
# Simple averaging of average in previous trial and current trial value; this
# is a particular case of the geometric moving average (see below) with g = 0.5
ewma = [None] * len(nTrials)


for nv in range(len(nTrials)):
    if nv == 0:
        ewma[nv] = nTrials[nv]
        prevValue = nTrials[nv]
    else:
        updateValue = (prevValue+nTrials[nv])/2
        ewma[nv] = updateValue
        prevValue = updateValue


plt.plot([i for i in range(len(nTrials))],
         ewma)

#===================== Geometric moving average ================================
# This is the model employed by Dyjas et al. (2012)

# First compute a single gma (for later comparison with other models)
gma = [None] * len(nTrials)

g=0.9
for nv in range(len(nTrials)):
    if nv == 0:
        gma[nv] = nTrials[nv]
        prevValue = nTrials[nv]
    else:
        updateValue = prevValue * g + nTrials[nv] * (1-g)
        gma[nv] = updateValue
        prevValue = updateValue

plt.plot([i for i in range(len(nTrials))],
          gma)

# Function to compute gmas
def buildGma(gValue, listTrials):
    gma = [None] * len(listTrials)
    for nv in range(len(listTrials)):
        if nv == 0:
            gma[nv] = listTrials[nv]
            prevValue = listTrials[nv]
        else:
            updateValue = prevValue * gValue + listTrials[nv] * (1-gValue)
            gma[nv] = updateValue
            prevValue = updateValue
    
    return(gma)


# Plot multiple values of g
gList = np.arange(0.8, 0.99, 0.02)
fig = plt.figure()
ax = fig.subplots(2, 5)
ax = ax.ravel()

for iG, panel in enumerate(ax):
    gValue = gList[iG]
    gma = buildGma(gValue, nTrials)
    panel.plot([i for i in range(len(nTrials))],
             gma)
    

#=================== Equal weights with multiple traces (values) =============
ewmt = [None] * (len(nTrials))

for nv in range(len(nTrials)):
    if nv == 0:
        ewmt[nv] = nTrials[nv]
        prevValue = nTrials[nv]
    else:
        #updateValue = (sum(nTrials[0:nv])+nTrials[nv])/(1+nv)
        updateValue = (sum(ewmt[0:nv])+nTrials[nv])/(1+nv)
        ewmt[nv] = updateValue
        prevValue = updateValue


plt.plot([i for i in range(len(nTrials))],
         ewmt)

# Plot one against the other
plt.plot([i for i in range(len(nTrials))], ewma, label='ewma')
plt.plot([i for i in range(len(nTrials))], gma, label='gma')
plt.plot([i for i in range(len(nTrials))], ewmt, label='ewmt')
plt.legend()

# Plot differences
diffs = [v1 - v2 for (v1, v2) in zip(ewma, ewmt)]

plt.plot([i for i in range(len(nValues))],
         diffs)



#=============== Decaying weights with multiple traces (values) ================
    




#=============================================================================
#=============== 2. Values draw from distributions (no noise) ================
#=============================================================================
# These models draw duration values from distributions. It does not estimate
# the variance; it is only used to produce noisy estimates

#======================= Equal weights moving average =========================
ewma = [None] * len(nTrials)

sd = 0.1
for nv in range(len(nTrials)):
    tValue = np.random.normal(nTrials[nv], sd)
    if nv == 0:
        ewma[nv] = tValue
        prevValue = tValue
    else:
        updateValue = (prevValue+tValue)/2
        ewma[nv] = updateValue
        prevValue = updateValue


plt.plot([i for i in range(len(nTrials))],
         ewma)


#===================== Geometric moving average ================================
# First compute a single gma (for later comparison with other models)
gma = [None] * len(nTrials)

sd = 0.1
g=0.9
for nv in range(len(nTrials)):
    tValue = np.random.normal(nTrials[nv], sd)
    if nv == 0:    
        gma[nv] = tValue
        prevValue = tValue
    else:
        updateValue = prevValue * g + tValue * (1-g)
        gma[nv] = updateValue
        prevValue = updateValue

plt.plot([i for i in range(len(nTrials))],
          gma)

# Function to compute gmas
def buildGma(gValue, sdValue, listTrials):
    gma = [None] * len(listTrials)
    for nv in range(len(listTrials)):
        tValue = np.random.normal(listTrials[nv], sdValue)
        if nv == 0:
            gma[nv] = tValue
            prevValue = tValue
        else:
            updateValue = prevValue * gValue + tValue * (1-gValue)
            gma[nv] = updateValue
            prevValue = updateValue
    
    return(gma)


# Plot multiple values of g
gList = np.arange(0.8, 0.99, 0.02)
fig = plt.figure()
ax = fig.subplots(2, 5)
ax = ax.ravel()

for iG, panel in enumerate(ax):
    gValue = gList[iG]
    gma = buildGma(gValue, 0.1, nTrials)
    panel.plot([i for i in range(len(nTrials))],
             gma)
    ax[iG].set_title('g =' + str(round(gValue, 2)))
    
# Plot multiple values of sd
sdList = np.arange(0.05, 0.55, 0.05)
fig = plt.figure()
ax = fig.subplots(2, 5)
ax = ax.ravel()

for iSd, panel in enumerate(ax):
    sdValue = sdList[iSd]
    gma = buildGma(0.98, sdValue, nTrials)
    panel.plot([i for i in range(len(nTrials))],
             gma)
    ax[iSd].set_title('sd =' + str(round(sdValue, 2)))


#=================== Equal weights with multiple traces (values) =============
ewmt = [None] * (len(nTrials))

sd = 0.5

for nv in range(len(nTrials)):
    tValue = np.random.normal(nTrials[nv], sd)
    if nv == 0:
        ewmt[nv] = tValue
        prevValue = tValue
    else:
        #updateValue = (sum(nTrials[0:nv])+nTrials[nv])/(1+nv)
        updateValue = (sum(ewmt[0:nv])+tValue)/(1+nv)
        ewmt[nv] = updateValue
        prevValue = updateValue


plt.plot([i for i in range(len(nTrials))],
         ewmt)


# Plot one against the other
plt.plot([i for i in range(len(nTrials))], ewma, label='ewma')
plt.plot([i for i in range(len(nTrials))], gma, label='gma')
plt.plot([i for i in range(len(nTrials))], ewmt, label='ewmt')
plt.legend()

# Plot differences
diffs = [v1 - v2 for (v1, v2) in zip(ewma, ewmt)]

plt.plot([i for i in range(len(nValues))],
         diffs)











# Random Walk




#=======================================================================================
#============================== Jones et al. (2013) ====================================
#=======================================================================================

def base (trial, epsilon):
    


