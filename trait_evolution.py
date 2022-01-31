import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from function_biotime_trends import plotting
import warnings
"""
path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')
meta_data["duration"] = meta_data["END_YEAR"] - meta_data["START_YEAR"]

ind_meta = ((meta_data.TAXA == "Marine plants")
            | (meta_data.TAXA == "Freshwater plants"))
#ind_meta = meta_data.TAXA == "Freshwater plants"
ind_meta = ind_meta & (meta_data["duration"] >= 10)
study_ids = [str(std_id) for std_id in meta_data.STUDY_ID[ind_meta].values]

presences = np.load(path + "biotime_converted.npz", allow_pickle=True)

tr_base = pd.read_csv("Phytoplankton_trait_data_rimet.csv",
                     encoding = 'unicode_escape')
litch = pd.read_csv("Phytoplankton_size.csv", encoding = "unicode_escape")
traits_all = {}
years = {}

richness = {}

for study_id in study_ids:
    years[study_id] = presences["year_" + study_id]
    traits = np.full(presences[study_id].shape[1], np.nan)
    traits2 = traits.copy()
    trait = "Cell biovolume mum3"
    for i, sp in enumerate(presences["species {}".format(study_id)]):
        traits[i] = np.nanmean(np.log(np.append(
            tr_base.loc[tr_base["Genus + species name"] == sp,
                                trait].values,
            litch.loc[litch["species"] == sp, "volume"].values)))
        #traits[i] = np.nanmean(np.log(litch.loc[litch["species"] == sp, "volume"].values))

    print(study_id, np.sum(np.isfinite(traits)), len(traits), np.sum(np.isfinite(traits2)))
    #traits[traits == 0] = np.nan
    traits_all[study_id] = traits.copy()
    
    # add richness counts
    pres = presences[study_id] == 1 # conversion to boolean
    # total species richness
    richness[study_id] = np.sum(pres, axis = 1)
  
fig, ax = plt.subplots(4,2, sharex = "row", figsize = (10,10))
trait_means = {}
trait_numbers = {}
trait_fraction = {}
for key in study_ids:
    with warnings.catch_warnings(record = True):
        dummy = traits_all[key]*presences[key]
        dummy[dummy == 0] = np.nan
        trait_means[key] = np.nanmean(dummy, axis = 1)
        trait_numbers[key] = np.sum(np.isfinite(dummy), axis = 1)
        trait_fraction[key] = np.mean(np.isfinite(dummy), axis = 1)
    trait_means[key][trait_means[key] == 0] = np.nan
met_traits = plotting(years, trait_means, ax[:,0], presences)
#met_numbs = plotting(years, trait_numbers, ax[:,1], presences)
#plotting(years, trait_fraction, ax[:,2], presences)
plotting(years, richness, ax[:,1], presences)
fig.tight_layout()
ax[0,0].legend()
fig.savefig("Figure_phytoplankton_trait_evolutin.pdf")


fig, ax = plt.subplots(2,3, figsize = (9,9), sharex = True)
ax = ax.flatten()
for i, study_id in enumerate(study_ids):
    
    if not np.any(np.isfinite(traits_all[study_id])):
        continue
    ax[i].hist(traits_all[study_id], bins = np.sum(np.isfinite(traits_all[study_id]))//5)
    ax[i].text(0.9,0.9,np.sum(np.isfinite(traits_all[study_id])),
               transform=ax[i].transAxes, va = "bottom", ha = "left")
    ax[i].set_title(study_id)
    ax[i].set_xlabel("Phytoplankton size")
    ax[i].set_ylabel("Frequency")
fig.tight_layout()

fig.savefig("Figure_phytoplankton_trait_distribution.pdf")
"""
key = "254"
dummy = traits_all[key]*presences[key]
dummy[dummy == 0] = np.nan
bins = np.linspace(np.nanmin(dummy), np.nanmax(dummy), 31)

fig, ax = plt.subplots(4,5, sharex = True, sharey = True)
ax = ax.flatten()
density = False
for i in range(len(dummy)):
    ax[i].hist(traits_all[key], density = density, bins = bins)
    ax[i].hist(dummy[i], density = density, bins = bins, alpha = 0.5)