# -*- coding: utf-8 -*-
"""
Created on Tue Dec  6 14:10:32 2022
- Por enquanto, utilizar apenas parâmetros investigados no artigo (k, c e r)
- 10 valores para cada parâmetro: linearmente espaçados entre o mínimo e o máximo
daqueles utilizados no artigo
- para cada combinação, gerar curvas para efeito hazard e efeitos sequenciais
- quantificar o efeito hazard (calcular slope fp x rt)
- quantificar o efeito sequencial (subtrair slopes fp n-1 longo e curto)
- plotar scatterplots de (índices de) efeito hazard x efeito sequencial;
- observar se há:
    - uma linha (o efeito sequencial e o efeito hazard estão atrelados para todos
                 os valores de cada parâmetro);
    - pontos completamente dispersos (cada combinação gera um efeito sem relação
                                      com o outro - improvável)
    - clusters (correlação para alguns valores, ausência para outros)
- Hipótese: como se trata de um modelo que descende do condicionamento de traço,
é improvável que ele possa gerar efeitos dissociados entre hazard e sequencial:
no condicionamento de traço, o hazard é uma consequência do sequencial; na 
fMTP, o hazard e o sequencial são produtos de um mesmo mecanismo de memória,
então provavelmente sempre estarão correlacionados
- Neste caso, outro modelo, como o de duplo-processo - p. ex., na implementação 
de Grabenhorst et al., 2019/2021 - ou um modelo bayesiano talvez expliquem melhor
os resultados.
- Por enquanto, utilizar os dados que existem, sabendo que estão por ser
confirmados pelo segundo experimento

# Changes
- Rounded slope values for better visualization
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
import statsmodels.formula.api as smf
import seaborn as sns
import itertools

import os
os.chdir("G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Grid_search")

from fmtp import fMTP, FPexp, FPgonogo


    
#================================== Constants and parameters ================================#
# Initialize values for simulation
FP = np.arange(0.6, 1.8, 0.6)
distr = 'uni'
# Set-up experiment using "FPexp"
#exp = FPexp(FPs = FP, distribution = distr, tr_per_block = 150)
exp = FPgonogo(FPs = FP, distribution = distr, tr_per_block = 150, relax = 0)
#exp = FPgonogo(exp)      
                                              

# create list of parameter values
kList=[i for i in range(1,11)]  #np.linspace(2,8,num=10)
#rList=np.linspace(1,4,num=10)
rList=np.linspace(-4,-1,num=10)
cList=np.linspace(0,0.0003,num=10)

# Simulate experiments and store in lists
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
            
            # Hazard slope
            hazardLm=smf.ols(formula='prep~FP', data=state_discr).fit()
            hazardSlope=round(hazardLm.params[1],3)
            
            # Sequential effects difference of slopes
            lowFPn_1=min(np.unique(state_discr.FPn_1))
            highFPn_1=max(np.unique(state_discr.FPn_1))
            # low FP n-1
            state_lowFPn_1=state_discr[state_discr.FPn_1==lowFPn_1]
            lowSeqLm=smf.ols(formula='prep~FP',data=state_lowFPn_1).fit()
            lowSeqSlope=lowSeqLm.params[1]
            # High FP n-1
            state_highFPn_1=state_discr[state_discr.FPn_1==highFPn_1]
            highSeqLm=smf.ols(formula='prep~FP',data=state_highFPn_1).fit()
            highSeqSlope=highSeqLm.params[1]
            
            seqSlopeDiff=round(highSeqSlope-lowSeqSlope,3)
            simResults[simEl]=(hazardSlope,seqSlopeDiff)
            # simResults[simEl]=(hazardSlope,
            #                    round(highSeqSlope,3),
            #                    round(lowSeqSlope,3),
            #                    seqSlopeDiff)
            print(simEl)
            simEl+=1

simResults2=simResults
            

#=================== plot indices on a scatterplot ======================#
plt.scatter(*zip(*simResults2))
plt.show()

# Histograms of sequential effects
sepSimResults=list(set(zip(*simResults)))

sns.histplot(sepSimResults[1],
             binwidth=0.1)

positiveSeqEff=[ef for ef in sepSimResults[0] if ef >0] # None are positive

# Plot on same panel, but with colors and annotations to indicate value of k        
chunkedList=list()
chunkSize=int(len(simResults)/len(kList))

for i in range(0,len(simResults),chunkSize):
    chunkedList.append(simResults[i:i+chunkSize])

fig=plt.figure()

kSublist=0
for ik in range(len(chunkedList)):
    k=kList[kSublist]
    plotSublist=chunkedList[kSublist]
    plt.scatter(*zip(*plotSublist))
    plt.annotate(('k='+str(k)),xy=(plotSublist[-1][0],plotSublist[-1][1]+0.002))
    kSublist+=1

figname='G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Plots/Modelling/hazardSeqSlopeScatterplot.svg'
plt.savefig(figname,format='svg')


# Plot by value of k (each plot with all combinations of r and c)
fig,ax=plt.subplots(nrows=2,ncols=5)

kSublist=0
for row in ax:
    for col in row:
        k=kList[kSublist]
        plotSublist=chunkedList[kSublist]
        col.scatter(*zip(*plotSublist))
        col.set_title('k='+str(k))
        kSublist+=1

plt.savefig('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Plots/Modelling/hazardSeqSlopeScatterplot_facet.png',
            format='png')


# Separate plots for each r value (within a given value of k)
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

plt.savefig('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Plots/Modelling/rxc_k1_facet.png',
            format='png',bbox_inches='tight')
        
        
# Compute correlation
# hazardSlopeList=[value[0] for value in simResults]
# seqSlopeDiffList=[value[1] for value in simResults]

# np.corrcoef(hazardSlopeList, seqSlopeDiffList)
np.corrcoef(*zip(*simResults2))

corrList=list()
for ik in range(len(chunkedList)):
    corrSublist=chunkedList[ik]
    corrList.append(np.corrcoef(*zip(*corrSublist)))
    
    
        
# Without c
c = cList[0]
# List to store values
simResultsNoC=[None]*(len(kList)*len(rList))

simEl=0   
for ik in range(len(kList)):
    for ir in range(len(rList)):
        k=kList[ik]
        r=rList[ir]
        fmtp=fMTP(r,c,k)
        state_discr, state_con = exp.run_exp(fmtp)
        state_discr=state_discr[1:]
        
        # Hazard slope
        hazardLm=smf.ols(formula='prep~FP', data=state_discr).fit()
        hazardSlope=round(hazardLm.params[1],3)
        
        # Sequential effects difference of slopes
        lowFPn_1=min(np.unique(state_discr.FPn_1))
        highFPn_1=max(np.unique(state_discr.FPn_1))
        # low FP n-1
        state_lowFPn_1=state_discr[state_discr.FPn_1==lowFPn_1]
        lowSeqLm=smf.ols(formula='prep~FP',data=state_lowFPn_1).fit()
        lowSeqSlope=lowSeqLm.params[1]
        # High FP n-1
        state_highFPn_1=state_discr[state_discr.FPn_1==highFPn_1]
        highSeqLm=smf.ols(formula='prep~FP',data=state_highFPn_1).fit()
        highSeqSlope=highSeqLm.params[1]
        
        seqSlopeDiff=round(highSeqSlope-lowSeqSlope,3)
        simResultsNoC[simEl]=(hazardSlope,seqSlopeDiff)
        
        print(simEl)
        simEl+=1

sepSimResultsNoC=list(set(zip(*simResultsNoC)))
hazardList=list(sepSimResultsNoC[0])
reshapedHazardList=np.reshape(hazardList, (10,10))

hazardData=pd.DataFrame(reshapedHazardList.transpose(),
                        columns=kList,
                        index=[round(item,3) for item in rList])

hazardData.to_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Grid_search/hazardData_gonogo.csv')

seqEffList=list(sepSimResultsNoC[1])
reshapedSeqEffList=np.reshape(seqEffList, (10,10))
seqEffData=pd.DataFrame(reshapedSeqEffList.transpose(),
                        columns=kList,
                        index=[round(item,3) for item in rList])

seqEffData.to_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Grid_search/SeqEffData_gonogo.csv')

#=========================================================================================#

paramPairs=list(itertools.product(kList, [round(item,3) for item in rList]))

paramCombinations=list(itertools.combinations(paramPairs, 2))

abs(seqEffData.loc[paramCombinations[4789][0][1],paramCombinations[4789][0][0]]-seqEffData.loc[paramCombinations[4789][1][1],paramCombinations[4789][1][0]])<0.03

relParams={pair[0]:pair[1] for pair in paramCombinations if 
              abs(seqEffData.loc[pair[0][1],pair[0][0]] - seqEffData.loc[pair[1][1],pair[1][0]]) < 0.03}

len(relParams)

# Simulate experiment
for params in relParams:
    k=params[0]
    r=params[1]
    fmtp=fMTP(r,c,k)
    state_discr, state_con = exp.run_exp(fmtp)
    state_discr=state_discr[1:]
    mean_state=state_discr.groupby(['FP']).mean().reset_index()
    f, ax = plt.subplots(1,2)
    ax[0].plot(mean_state.FP, mean_state.prep, '.-')
    ax[0].set_title('k='+str(k)+', r='+str(r))
    
    mean_state = state_discr.groupby(['FP', 'FPn_1']).mean().reset_index()        
    for iFP in np.unique(mean_state.FPn_1):
        state_n1 = mean_state[mean_state.FPn_1 == iFP]
        ax[1].plot(state_n1.FP, state_n1.prep, '.-')
    ax[1].legend([0.6, 1.2, 1.8])




# Build plots
fig=plt.figure(constrained_layout=True)
subfigs=fig.subfigures(4,5)

subfigsCoords=list(itertools.product([0,1,2,3],[0,1,2,3,4]))

iSubfig=0
for params in relParams:
    # Simulate experiment
    k=params[0]
    r=params[1]
    fmtp=fMTP(r,c,k)
    state_discr, state_con = exp.run_exp(fmtp)
    state_discr=state_discr[1:]
    
    # Create subplots within subfigures
    ax=subfigs[subfigsCoords[iSubfig]].subplots(1,2,sharey=True)
    
    mean_state=state_discr.groupby(['FP']).mean().reset_index()
    ax[0].plot(mean_state.FP, mean_state.prep, '.-')
    #ax[0].set_title('k='+str(k)+', r='+str(r))
    
    mean_state = state_discr.groupby(['FP', 'FPn_1']).mean().reset_index()        
    for iFP in np.unique(mean_state.FPn_1):
        state_n1 = mean_state[mean_state.FPn_1 == iFP]
        ax[1].plot(state_n1.FP, state_n1.prep, '.-')
    ax[1].legend([0.6, 1.2, 1.8])
    subfigs[subfigsCoords[iSubfig]].suptitle('k='+str(k)+', r='+str(r))
    iSubfig+=1

    
    
    


