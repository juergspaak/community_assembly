

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from functions_for_plotting import biotime_colors, keys

# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

try:
    d_prey
    d_pred
except NameError:
    d_prey = np.load(path + "data_LV_assembly.npz")
    d_pred = np.load(path + "data_LV_assembly_predators.npz")
    d_prey = {key: d_prey[key] for key in d_prey.files}
    d_pred = {key: d_pred[key] for key in d_pred.files}


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

study_ids = presences["study_ids"]
# only work with long datasets
study_ids = study_ids[[presences[study_id].shape[0] >= 30
                      for study_id in study_ids]]

meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape', index_col="STUDY_ID")


meta_data = meta_data.loc[[int(std_id) for std_id in study_ids]]
meta_data["identity"] = [biotime_colors[key][0]
                         if key in biotime_colors.keys() else np.nan
                         for key in meta_data.TAXA]
meta_data["corr"] = np.nan # store correlation between invasion and extinction
meta_data["p_value"] = np.nan # p_value for correlation
meta_data["rand_corr"] = np.nan # correlation of perturbed data


def inv_ext_correlation(pres):
    # compute correlation between invasion and extinction
    
    # total species richness
    richness = np.sum(pres, axis = 1)
    
    # not present last year, but present this year
    inv = np.sum(pres[1:] & ~pres[:-1], axis = 1)
    # how many of the present species are invaders?
    inv = inv/richness[1:]
    
    # present last year, but not present this year
    ext = np.sum(~pres[1:] & pres[:-1], axis = 1)
    # how many of the present species are invaders?
    ext = ext/richness[:-1]
    
    #d = 1
    # return correlation, extinctions cause invasions
    #return np.corrcoef(inv[d:], ext[:-d])[0,1]
    # invasions cause extinctions
    #return np.corrcoef(inv[:-1], ext[1:])[0,1]
    # 
    return np.corrcoef(inv, ext)[0,1]


# compute correlation between invasion and extinction
itera = 1000
p_distribution = np.empty((len(meta_data), itera))
for j, study_id in enumerate(meta_data.index):
    pres = presences[str(study_id)] == 1 # conversion to boolean
    
    meta_data.loc[study_id, "corr"] = inv_ext_correlation(pres)
    
    for i in range(itera):
        order = np.argsort(np.random.rand(len(pres)))
        p_distribution[j, i] = inv_ext_correlation(pres[order])
        
    meta_data.loc[study_id, "p_value"] = np.sum(meta_data.loc[study_id, "corr"]>p_distribution[j])/itera
    meta_data.loc[study_id, "rand_corr"] = np.mean(p_distribution[j])
    
#"""
identities = [biotime_colors[key][0] for key in keys]

fig = plt.figure()
plt.hist([1-meta_data.loc[meta_data.identity == ident, "p_value"]
          for ident in identities],
         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
         label = keys)
plt.xlabel("p value of correlation")
plt.ylabel("Frequency")

plt.xticks([0,0.05, 0.5, 1],  ["0",0.05, 0.5, 1])
plt.xlim([0,1])
fig.savefig("Figure_invasion_extinctions_correlation.pdf")

fig, ax = plt.subplots(1,2,figsize = (9,7), sharey = True)

bins = np.linspace(-1,1,41)
ax[0].hist(meta_data["corr"], bins = bins, label = "Real correlation")
ax[0].hist(meta_data.rand_corr, bins = bins, label = "Reference", alpha = 0.5)

ax[1].hist([1-meta_data.loc[meta_data.identity == ident, "p_value"]
          for ident in identities],
         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
         label = keys)

ax[0].legend()
ax[0].axvline(0, color = "k")

ax[0].set_ylabel("Frequency")
ax[0].set_xlabel("Correlation of\ninvasion and extinction")
ax[1].set_xlabel("p-value of correlation")



fig.savefig("Figure_ap_correlation.pdf")
#"""
##############################################################################
# compare actual correlation with correlation based on randomness

fig, ax = plt.subplots(5,5, sharex = True, sharey = True, figsize = (9,9))
for i in range(len(ax)):
    ax[i,0].set_ylabel("Frequency")
    ax[-1,i].set_xlabel("Correlation")

ax = ax.flatten()
ind = np.argsort(1-meta_data.p_value.values)

for i, a in enumerate(ax):
    a.hist(p_distribution[ind[i]], bins = bins, color = biotime_colors[meta_data.TAXA.values[ind[i]]],
           density = True)
    a.axvline(meta_data["corr"].values[ind[i]], color = "r")
    a.axvline(0, color = "grey")
    a.set_title(meta_data.TAXA.values[ind[i]])
    
fig.tight_layout()
fig.savefig("Figure_ap_comparison_correlatoin.pdf")



