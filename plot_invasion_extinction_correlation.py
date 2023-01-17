import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from functions_for_plotting import biotime_colors, keys



path = "C:/Users/Jurg/Documents/science_backup/P14_community_assembly/"
presences = np.load("biotime_converted.npz")

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
meta_data["p_value"] = np.nan # p_value for 
meta_data["p_value_high_rich"] = np.nan # p_value for correlation for high richness
meta_data["rand_corr"] = np.nan # correlation of perturbed data

meta_data["include"] = True


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
itera = 100
p_distribution = np.empty((2, len(meta_data), itera))
for j, study_id in enumerate(meta_data.index):
    pres = presences[str(study_id)] == 1 # conversion to boolean
    richness = np.sum(pres, axis = 1)
    
    # for exclusion of communities with to strong fluctuations
    if min(richness)/max(richness)<0.25:
        meta_data.loc[study_id, "include"] = False
    
    # for exclusion of timepoints with very low species richness
    ind_high_rich = richness>=5
    
    meta_data.loc[study_id, "corr"] = inv_ext_correlation(pres)
    meta_data.loc[study_id, "corr_high_rich"] = inv_ext_correlation(pres[ind_high_rich])
    
    
    for i in range(itera):
        order = np.argsort(np.random.rand(len(pres)))
        p_distribution[0, j, i] = inv_ext_correlation(pres[order])
    
    pres_high_rich = pres[ind_high_rich]
    for i in range(itera):
        order = np.argsort(np.random.rand(len(pres_high_rich)))
        p_distribution[1, j, i] = inv_ext_correlation(pres_high_rich[order])
        
    meta_data.loc[study_id, "p_value"] = np.sum(meta_data.loc[study_id, "corr"]>p_distribution[0, j])/itera
    meta_data.loc[study_id, "rand_corr"] = np.mean(p_distribution[0, j])
    meta_data.loc[study_id, "p_value_high_rich"] = np.sum(meta_data.loc[study_id, "corr_high_rich"]>p_distribution[1, j])/itera
    if np.isnan(meta_data.loc[study_id, "p_value"]):
        raise
    
#"""
keys = ["Birds","Inverte-\nbrates",
                        "Terrestrial\nplants", "Fish", "Other"]
keys = list(set(meta_data["REALM"]))

identities = [biotime_colors[key][0] for key in keys]

fig = plt.figure()
#plt.hist([1-meta_data.loc[meta_data.identity == ident, "p_value"]
#          for ident in identities],
#         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
#         label = keys)
plt.hist([1-meta_data.loc[meta_data.REALM == ident, "p_value"]
          for ident in set(meta_data.REALM)],
         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
         label = keys)
fs = 16
plt.xlabel("p value of correlation", fontsize = fs)
plt.ylabel("Frequency", fontsize = fs)

plt.legend(fontsize = fs)
plt.xticks([0,0.05, 0.5, 1],  ["0",0.05, 0.5, 1])
plt.xlim([0,1])
fig.savefig("Figure_invasion_extinctions_correlation.pdf")

##############################################################################
# exclude datasets with to large fluctuations in richness

fig, ax = plt.subplots(1,2,sharex = True, figsize = (9,7), sharey = True)

ax[0].hist([1-meta_data.loc[meta_data.identity == ident, "p_value_high_rich"]
          for ident in identities],
         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
         label = keys)

m_data = meta_data[meta_data.include]
ax[1].hist([1-m_data.loc[m_data.identity == ident, "p_value"]
          for ident in identities],
         bins = np.linspace(0,1,41), stacked = True, color = [biotime_colors[key] for key in keys],
         label = keys)
fs = 16
ax[0].set_xlabel("p value of correlation", fontsize = fs)
ax[1].set_xlabel("p value of correlation", fontsize = fs)
ax[0].set_ylabel("Frequency", fontsize = fs)

ax[0].legend(fontsize = fs)
plt.xticks([0, 0.5, 1],  [0, 0.5, 1])
plt.xlim([0,1])

ax[0].set_title("A", loc = "left", fontsize = fs)
ax[1].set_title("B", loc = "left", fontsize = fs)

fig.savefig("Figure_ap_correlation_exclude_datasets.pdf")

#"""
##############################################################################
# compare actual correlation with correlation based on randomness

fig, ax = plt.subplots(5,5, sharex = True, sharey = True, figsize = (9,9))
bins = np.linspace(-1,1,51)
for i in range(len(ax)):
    ax[i,0].set_ylabel("Frequency")
    ax[-1,i].set_xlabel("Correlation")

ax = ax.flatten()
ind = np.argsort(1-meta_data.p_value.values)

for i, a in enumerate(ax):
    a.hist(p_distribution[0, ind[i]], bins = bins, color = biotime_colors[meta_data.TAXA.values[ind[i]]],
           density = True)
    a.text(0.08, 0.85, "p = {}".format(np.round(1-meta_data.p_value.values[ind[i]],2)), 
           transform = a.transAxes)
    a.axvline(meta_data["corr"].values[ind[i]], color = "r")
    a.axvline(0, color = "grey")
    a.set_title(meta_data.TAXA.values[ind[i]])
    
fig.tight_layout()
fig.savefig("Figure_ap_comparison_correlatoin.pdf")
#"""
