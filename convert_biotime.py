import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
file = "BioTIMEQuery_24_06_2021.csv"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')
meta_data["duration"] = meta_data["END_YEAR"] - meta_data["START_YEAR"]


ind_meta = ((meta_data.TAXA == "Marine plants")
            | (meta_data.TAXA == "Freshwater plants"))
ind_meta = ind_meta & (meta_data.duration >=10)
ind_meta = meta_data.duration >= 10

studies = meta_data.STUDY_ID[ind_meta].values
try:
    biotime_data
except NameError:
    biotime_data = pd.read_csv(path + file)
    
study = biotime_data[biotime_data.STUDY_ID == studies[0]]



presences = {}
years_dict = {}
species_dict = {}

for i, study_id in enumerate(studies):
    print(i, len(studies))
    study = biotime_data[biotime_data.STUDY_ID == study_id]
    species = np.unique(study.GENUS_SPECIES)
    years = np.unique(study.YEAR)
    pres =  pd.DataFrame(0, index = years, columns=species)
    
    for j, year in enumerate(years):
        pres.loc[year, set(study.loc[study.YEAR == year, "GENUS_SPECIES"])] = 1

    if np.any(np.isnan(pres)):
        raise
    presences[str(study_id)] = pres
    years_dict["year_"+str(study_id)] = pres.index.values
    species_dict["species {}".format(study_id)] = pres.columns.values
    
#"""

np.savez(path + "biotime_converted.npz", **presences, **years_dict,
         study_ids = list(presences.keys()),
         **species_dict)