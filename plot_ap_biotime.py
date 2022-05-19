import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from scipy.stats import pearsonr

from functions_for_plotting import biotime_colors as colors, keys


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz", allow_pickle=True)

study_ids = presences["study_ids"]
# only work with long datasets
study_ids = study_ids[[presences[study_id].shape[0]>=30
                      for study_id in study_ids]]
study_ids_int = [int(s) for s in study_ids]

meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape', index_col="STUDY_ID")
meta_data = meta_data.loc[study_ids_int]
meta_data["identity"] = [colors[key][0]
                         if key in colors.keys() else np.nan
                         for key in meta_data.TAXA]

meta_data["p_inv"] = np.nan
meta_data["p_ext"] = np.nan
meta_data["R2_inv"] = np.nan
meta_data["R2_ext"] = np.nan

##############################################################################
# compute correlations

for study_id in study_ids:
    # load dataset, shape (years, species)
    pres = presences[study_id] == 1
    
    years = presences["year_{}".format(study_id)][:-1]
    
    # proportion of invaders
    invasion = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[1:], axis = 1)
    
    meta_data.loc[int(study_id), ["R2_inv", "p_inv"]] = pearsonr(years, invasion)

    # proportion of extinctions
    extinction = np.sum(pres[:-1] & ~pres[1:], axis = 1)/np.sum(pres[:-1], axis = 1)
    meta_data.loc[int(study_id), ["R2_ext", "p_ext"]] = pearsonr(years, extinction)
    

meta_data["R2_ext"] = meta_data["R2_ext"]**2
meta_data["R2_inv"] = meta_data["R2_inv"]**2
    
identities = [colors[key][0] for key in keys]

fig, ax = plt.subplots(2,2, figsize = (9,9), sharex = True, sharey = True)

ax[0,0].set_ylabel("Frequency")
ax[1,0].set_ylabel("Frequency")
ax[1,0].set_xlabel("P-value")
ax[1,1].set_xlabel("$R^2$")

ax[0,0].set_xlim([0,1])

ax[0,0].set_title("Invasions")
ax[0,1].set_title("Extinctions")

ax = ax.flatten()

for i, par in enumerate(["p_inv", "R2_inv", "p_ext", "R2_ext"]):
    ax[i].hist([meta_data.loc[meta_data.identity == ident, par]
          for ident in identities],
         bins = np.linspace(0,1,41), stacked = True, color = [colors[key] for key in keys],
         label = keys)
    ax[i].set_title("ABCD"[i], loc = "left")
    

fig.savefig("Figure_ap_biotime.pdf")