# -*- coding: utf-8 -*-
"""
Created on Wed Jan 11 16:55:54 2023

@author: alpno
"""
import numpy as np
import pandas as pd
import os

os.chdir('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Fitting')

# Import for plotting
import matplotlib.pyplot as plt

# Import classes
from fmtp import fMTP, FPexp, FPgonogo
from hazard import fMTPhz
from fit import sort_fit, show_fit, get_fit 


FPs = np.arange(0.6, 1.8, 0.6)
distr = 'uni'

# Parameterization fMTP
k = 4 # temporal precision
r = -2.81 # rate of forgettting
c = 0.0001 # memory persistence
fmtp = fMTP(r, c, k)
    
    
# Define fMTPhz
fMTPhz_lin = fMTPhz(-3.0, 3e-4,  k)
fMTPhz_inv = fMTPhz(-2.8, 1e-4, k)


# Initialize the fitting procedure
sorter = sort_fit(fmtp, fMTPhz_lin, fMTPhz_inv)
    
# Empirical data
emp_dat = pd.read_csv('G:/My Drive/Post-doc/Projetos/Action_foreperiod/Experimento_0/Analysis/Modelling/Code_Salet_et_al/model_code/empirical_data/niemi1979.csv')
# Discrete FPs 
FP_con = np.array([0.5, 5.5]) # constant
FP_uni = np.arange(0.5, 3.5, 0.5) # uniform
# Sort empirical data
emp_con = emp_dat.loc[emp_dat.condition == 'bright_pure'].RT.values        
emp_uni = emp_dat.loc[emp_dat.condition == 'bright'].RT.values  
emp = pd.DataFrame({'RT' : np.r_[emp_con, emp_uni], 'distrib': np.r_[
['constant'] * len(emp_con), ['uni'] * len(emp_uni)]})
emp['FP'] = np.r_[FP_con, FP_uni]

# Simulate experiment
sim_con, prep_con = sorter.run_dist(FP_con, None, np.repeat(['constant'], 
                                                          len(FP_con)))
sim_uni, prep_uni = sorter.run_dist([None], FP_uni, ['uni'])
sim = pd.concat([sim_con, sim_uni]).reset_index()

# Plot settings
sim.loc[sim.FP == 5.5, 'FP'] = 4.25 # change x-axis
emp.loc[emp.FP == 5.5, 'FP'] = 4.25
xlim = [0, 4250 + 500]
xticks = np.arange(500, 5001, 1250)
ylim = [200, 350]
yticks = np.arange(225, 351, 50)

# Fit models and display preparation curves together with empirical data
models = [['fMTP'], ('fMTPhz', 'fMTPhz_inv'), 
      ('sub_hz', 'sub_hz_inv'), ('c_hz', 'c_hz_inv')]
# Initialise plot
show = show_fit(xlim, ylim, xticks, yticks)
# Fit constant and uniform paradigm separately
emp_con = emp[emp.distrib == 'constant']
emp_uni = emp[emp.distrib == 'uni']
# Start fitting procedure looping through all the models
show.print_results(title = 'Niemi (1979)', model = '', coef = '', R = '')
for mod in models: 
    show.init_plot(title = mod)
    for sub_mod in mod:
        # Fit models
        fit = get_fit(emp_con, sim_con, prep_con, 'distrib') # constant
        coef_con, R_con = fit.run_fit('Niemi (1979)', sub_mod)
    
        fit = get_fit(emp_uni, sim_uni, prep_uni, 'distrib')  # uniform
        coef_uni, R_uni = fit.run_fit('Niemi (1979)', sub_mod)
        # Show plots 
        print('Constant FP paradigm: ')
        show.show_dist(sub_mod, coef_con, emp_con, sim_con, prep_con, 
                       'distrib', R_con, clr_it = 0)
        print('Uniform FP paradigm: ')
        show.show_dist(sub_mod, coef_uni, emp_uni, sim_uni, prep_uni, 
                       'distrib', R_uni, clr_it = 1)
        print('')
        

# Run fMTP 
exp = FPexp(FPs = FPs, distribution = distr, tr_per_block = 500)                          
sim, prep_fmtp = exp.run_exp(fmtp) # fMTP 
sim.rename(columns = {'prep' : 'fMTP'}, inplace = True)
sim_hz, prep_hz = exp.run_exp(hz_lin) # fMTPhz
sim['fMTPhz'] = sim_hz.prep
sim_hz_inv, prep_hz_inv = exp.run_exp(self.hz_inv, inv_map = True) 
sim['fMTPhz_inv'] = sim_hz_inv.prep # inverse fMTPhz

# Average preparation at discrete FPs
sim = sim.iloc[1:, :] # first trial has no preparation
sim = sim.groupby('FP').mean().reset_index()
sim['distrib'] = distr

# Average preparation curves
prep = pd.DataFrame()
prep['fMTP'] = np.mean(prep_fmtp.iloc[1:, :]) # again remove 1st trial
prep['fMTPhz'] = np.mean(prep_hz.iloc[1:, :])  
prep['fMTPhz_inv'] = np.mean(prep_hz_inv.iloc[1:, :])
prep['distrib'] = distr
prep['FP'] = np.unique(sim.FP)[0]

return(sim, prep)

sim_uni, prep_uni =sorter.sim_dist(FPs, 'uni')
