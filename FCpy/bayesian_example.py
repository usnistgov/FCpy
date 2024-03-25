# -*- coding: utf-8 -*-
"""
@author: Evan Groopman
@affiliation: National Institute of Standards and Technology (NIST)
@email: evan.groopman@nist.gov

Bayesian example for calculating confidence intervals on a Poisson process 
with a known background count rate.
"""

import arviz as az
az.style.use("arviz-darkgrid")
import matplotlib.pyplot as plt
import numpy as np
import pymc as pm

#%%

bkgd_rate = (1000.*10)/(3600.*10)
bkgd_sig = np.sqrt(1000*10)/(3600*10)
time = 3600 # seconds
n = 1200 #number of observed counts in [time]

with pm.Model() as model:
    # since we have observed counts and known background rate, convert everything to counts to start
    # otherwise the sampling for potentially very small count rates can cause lots of divergences in the MCMC
    # convert back to count rate at the end
    
    ### Priors ###
    bkgd_counts = pm.TruncatedNormal('bkgd_counts', mu=bkgd_rate*time, sigma= bkgd_sig*time, lower=0.0) #lower limit at 0
    #use TruncatedNormal prior for counts that is relatively wide (flat)
    counts = pm.TruncatedNormal('counts', mu=0.5, sigma=500, lower=0.0) #lower limit at 0. select appropriate mu and sigma as needed.
    
    # Total counts from both processes
    counts_tot = counts + bkgd_counts

    pois = pm.Poisson('counts_total', mu=counts_tot, observed=n)
    
    # convert counts of interest to count rate.
    # use pm.Deterministic() to ensure posterior sampling
    count_rate = pm.Deterministic('count_rate', counts/time)
    
    # perform sampling
    # need cores=1 for Spyder use, otherwise crashes
    # If running standalone, more cores can be used
    data = pm.sample(4000, init='advi+adapt_diag', chains=8, tune=2000, target_accept=0.95, progressbar=True, cores=1)
    

#%% Calculate CIs and plot trace results
#highest density interval = CI
hdi95 = az.hdi(data.posterior, hdi_prob=0.95) 
print(hdi95)
az.plot_trace(data, filter_vars="regex", var_names=["count_rate"], show=True, combined=True)

#%% Plot posterior
# az.plot_posterior(data.posterior.count_rate, hdi_prob=0.95, show=True)
az.plot_posterior(data.posterior.counts, hdi_prob=0.95, show=True)

