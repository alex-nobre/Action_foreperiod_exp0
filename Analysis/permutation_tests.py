####### For the permutation test, the python packages mne, pandas and numpy are used:
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mne

datAction=pd.read_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento 1/Analysis/lm_coeff_action.csv')

actionPermtestSign=pd.DataFrame(columns=['varname','pval','cval'])

for varname in ['intercept','fpMain','oneBackFPMain','interaction']:
    X = datAction[[varname]].values
    X = X.reshape((34, X.size//34))
    X-= X[:,0][:, np.newaxis]
    tobs, c,p,H0 = mne.stats.permutation_cluster_1samp_test(X,n_permutations=10000)
    for cval,pval in zip(c,p):
        if pval<0.1:
            testels=[varname,pval,cval]
            actionPermtestSign=actionPermtestSign.append(pd.DataFrame(testels,
                                                                      columns=['varname','pval','cval']),
                                                         ignore_index=True)
            print('teste significativo')
    print('===')
    

actionPermtestSign.to_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento 1/Analysis/'+'action_permtest_sign.csv')
 