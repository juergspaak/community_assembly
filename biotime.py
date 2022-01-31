import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""
path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
file = "BioTIMEQuery_24_06_2021.csv"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')
meta_data["duration"] = meta_data["END_YEAR"] - meta_data["START_YEAR"]
ind_meta = ((meta_data.TAXA == "Marine plants")
            | (meta_data.TAXA == "Freshwater plants"))
ind_meta = ind_meta & (meta_data.duration >=10)
#ind_meta = meta_data.duration >= 10

studies = meta_data.STUDY_ID[ind_meta].values
try:
    biotime_data
except NameError:
    biotime_data = pd.read_csv(path + file)
    
study = biotime_data[biotime_data.STUDY_ID == studies[0]]



presences = {}
changes = {}

for study_id in studies:
    study = biotime_data[biotime_data.STUDY_ID == study_id]
    species = set(study.SPECIES)
    presences[study_id] = pd.DataFrame(0, index = set(study.YEAR),
                                       columns=set(study.SPECIES))
    
    changes[study_id] = np.empty((2, len(presences[study_id])-1))
    for sp in set(study.SPECIES):
        presences[study_id].loc[set(study.loc[study.SPECIES == sp, "YEAR"]), sp] = 1
    
    pres = presences[study_id].values
    changes[study_id][0] = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[:-1], axis = 1)
    changes[study_id][1] = np.sum(~pres[1:] & pres[:-1], axis = 1)/np.sum(pres[:-1], axis = 1)

    np.round(changes[study_id], 10)
    
#"""
studies = studies
changes = changes    

fig, ax = plt.subplots(2,2, figsize = (9,9))
for study_id in studies:
    ax[0,0].plot(changes[study_id][0])
    ax[0,1].plot(changes[study_id][1])

ax[0,0].set_title("Invasion")
ax[0,1].set_title("Extincions")

corrs = np.empty(len(studies))
for i, study_id in enumerate(studies):
    ax[1,0].plot(changes[study_id][0][1:], changes[study_id][1][:-1])
    corrs[i] = np.corrcoef(changes[study_id][0][1:], changes[study_id][1][:-1])[0,1]
    ax[1,1].plot(changes[study_id].shape[1], corrs[i], 'o')
    #ax[1,1].plot(presences[study_id].shape[1], corrs[i], 'o')
    
invasion_all = []
extinction_all = []

for study_id in studies:
    invasion_all.extend(changes[study_id][0])
    extinction_all.extend(changes[study_id][1])
    
ax[1,1].axhline(np.corrcoef(invasion_all, extinction_all)[0,1], color = "r")

ax[0,0].set_xlabel("Time")
ax[0,1].set_xlabel("TIme")
ax[1,1].set_xlabel("Dataset length")
ax[1,0].set_xlabel("Invasions next year")
ax[1,0].set_ylabel("Extincions this year")
ax[1,1].set_ylabel("Correlation")
ax[0,0].set_ylabel("Probability")

fig.savefig("Figure_extinction_biotime.pdf")
