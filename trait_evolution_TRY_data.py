import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress
from function_biotime_trends import plotting
import warnings

"""
path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')
meta_data["duration"] = meta_data["END_YEAR"] - meta_data["START_YEAR"]

ind_meta = ((meta_data.TAXA == "Terrestrial plants")
             & (meta_data["duration"] >= 10)).values

study_ids = [str(std_id) for std_id in meta_data.STUDY_ID[ind_meta].values]

presences = np.load(path + "biotime_converted.npz", allow_pickle=True)
years = {study_id: presences["year_" + study_id] for study_id in study_ids}
"""
trait_ids = range(3106,3118)
#trait_ids = [3106] # only plant height, most data for it
for trait_id in trait_ids:
    
    ps = [[],[],[]]
    tr_base = pd.read_csv(path + "TRY_subdata_{}.csv".format(trait_id))
    traits_all = {}
    
    for study_id in study_ids:
        traits = np.full(presences[study_id].shape[1], np.nan)
        
        for i, sp in enumerate(presences["species {}".format(study_id)]):
            traits[i] = np.nanmedian(tr_base.loc[tr_base["SpeciesName"] == sp, "StdValue"])
            #traits[i] = np.nanmean(np.log(litch.loc[litch["species"] == sp, "volume"].values))
    
        print(study_id, np.sum(np.isfinite(traits)), len(traits))
        traits_all[study_id] = traits.copy()
    #''' 
    fig, ax = plt.subplots(4,3, sharex = "row", figsize = (10,10))
    
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
    met_numbs = plotting(years, trait_numbers, ax[:,1], presences)
    plotting(years, trait_fraction, ax[:,2], presences)
    fig.tight_layout()        
    
    ax[0,0].set_title("Trait evolution")
    ax[0,1].set_title("Number of species\nwith trait measurmetns")
    ax[0,2].set_title("Percentages of species\nwith trait measurments")
    
    ax[0,0].set_ylabel("Average trait")
    ax[0,1].set_ylabel("Number of species with traits")
    ax[0,2].set_ylabel("Percentages of species with traits")
    
    for i in range(3):
        ax[0,i].set_xlabel("Year")
    
    fig.suptitle(tr_base["TraitName"].iloc[0])
    fig.tight_layout()
    fig.savefig("Figure_trait_evolution_{}".format(trait_id))

    