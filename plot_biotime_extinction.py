import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from scipy.stats import spearmanr, pearsonr


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

def corr_invasion_extinction(pres, inv_ref = "t1", ext_ref = "t0"
                             , scaled = "_scaled"):
    
    # present next time point, but not now
    invasions = np.sum(pres[1:] & ~pres[:-1], axis = 1)
    
    if inv_ref == "t0": # of all present, which did invade?
        invasions = np.append(np.nan, invasions)
    else:
        invasions = np.append(invasions, np.nan)
        
    # not present last year, but present this year
    extinctions = np.sum(pres[:-1] & ~pres[1:], axis = 1)
    if ext_ref == "t0": # of all present, which did invade?
        extinctions = np.append(np.nan, extinctions)
    else:
        extinctions = np.append(extinctions, np.nan)  

    # divide by total richness
    if scaled:
        richness = np.sum(pres, axis = 1)
        invasions = invasions/richness                    
        extinctions = extinctions/richness
    ind = np.isfinite(invasions*extinctions)
    with warnings.catch_warnings(record = True):
        return np.corrcoef(invasions[ind], extinctions[ind])[0,1]


itera = 200

study_ids = presences["study_ids"]
study_ids = study_ids[[len(presences[study_id]) > 6 for study_id in study_ids]]



combs = [[x,y,z] for y in ["t0", "t1"] for x in ["t0", "t1"] for z in ["", "_scaled"]]

corrs_real = np.full((len(combs), len(study_ids)), np.nan)
corrs_virtual = np.full((len(combs), len(study_ids), itera), np.nan)

for i, comb in enumerate(combs):
    ps = np.full(len(study_ids), np.nan)
    for j, study_id in enumerate(study_ids):
        # first compute real correlation bewteen invasion and extinction
        pres = presences[study_id] == 1 # conversion to boolean
        corrs_real[i,j] = corr_invasion_extinction(pres, *comb)
        if np.isnan(corrs_real[i,j]):
            continue
        for k in range(itera):
            # shuffle years
            pres_virtual = pres[np.argsort(np.random.rand(len(pres)))]
            corrs_virtual[i,j,k] = corr_invasion_extinction(pres_virtual, *comb)


ps = np.sum(corrs_real[...,np.newaxis] > corrs_virtual, axis = -1)/itera
#"""

"""
fig, ax = plt.subplots(4,2, sharex = True, sharey = True, figsize = (9,9))
ax = ax.flatten()
for i in range(len(combs)):
    
    ax[i].hist(np.median(corrs_virtual[i], axis = -1),
               bins = np.linspace(-1,1,31), label = "simulated")
    ax[i].hist(corrs_real[i],
               bins = np.linspace(-1,1,31), alpha = .5, label = "real")
    ax[i].set_title(combs[i][0] + "_" + combs[i][1] + combs[i][2])
    ax[i].axvline(0, color = "red")

#ax[0].set_ylim([-1,1])
ax[-1].set_xlabel("correlation")
ax[-2].set_xlabel("correlation")
#ax[0].set_xlim([0,1])
ax[0].legend()

fig.savefig("Figure_extinction_biotime_correlation.pdf")

# p value of real correlations
fig, ax = plt.subplots(4,2, sharex = True, sharey = True, figsize = (9,9))
ax = ax.flatten()
for i in range(len(combs)):
    
    ax[i].hist(ps[i], bins = np.linspace(0,1,31))
    ax[i].set_title(combs[i][0] + "_" + combs[i][1] + combs[i][2])


#ax[0].set_ylim([-1,1])
ax[-1].set_xlabel("percentile")
ax[-2].set_xlabel("percentile")
#ax[0].set_xlim([0,1])
fig.savefig("Figure_extinction_biotime_percentile.pdf")
"""
study_len = np.array([len(presences[study]) for study in study_ids])

fig, ax = plt.subplots(4,2, sharex = True, sharey = True, figsize = (9,9))
ax = ax.flatten()

for i in range(len(combs)):
    
    ind = study_len>20
    ax[i].scatter(study_len, ps[i])
    ax[i].set_title(combs[i][0] + "_" + combs[i][1] + combs[i][2])
    ax[i].set_xlabel(np.round(pearsonr(study_len[ind], ps[i, ind])[0], 3))
    
fig.tight_layout()


    
    
    
   

    
