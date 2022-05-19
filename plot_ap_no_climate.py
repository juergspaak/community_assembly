import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from scipy.stats import pearsonr, linregress

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

##############################################################################
# compute correlations
ps = []
rs = []
occurences = []
for study_id in study_ids:
    # load dataset, shape (years, species)
    pres = presences[study_id] == 1
    
    years = presences["year_{}".format(study_id)]
    
    for i in range(pres.shape[1]):
        s, i, r, p, std = linregress(years, pres[:,i])
        ps.append(p)
        rs.append(r)
    
    # count how often a species invades and goes extinct
    occur = np.sum(pres[1:] & ~pres[:-1], axis = 0)
    occur += pres[0]
    occurences.extend(occur)
   

ps = np.array(ps)
rs = np.array(rs)
occurences = np.array(occurences)
rs = rs**2

fig, ax = plt.subplots(1,2, sharex = True, sharey = True)

ax[0].hist(ps, bins = np.linspace(0,1,41), density = True,
           label = "One occurence")
ax[0].hist(ps[occurences >= 2], bins = np.linspace(0,1,41), density = True,
           label = "Multiple occurences", alpha = 0.5)
ax[1].hist(rs, bins = np.linspace(0,1,41), density = True)
ax[1].hist(rs[occurences >= 2], bins = np.linspace(0,1,41), density = True,
           alpha = 0.5)

ax[0].set_xlabel("P-value")
ax[1].set_xlabel("$R^2$")

ax[0].set_ylabel("Frequency")

ax[0].set_title("A", loc = "left")
ax[1].set_title("B", loc = "left")

ax[0].set_xlim([0,1])
ax[0].legend()

fig.savefig("Figure_ap_no_climate.pdf")
