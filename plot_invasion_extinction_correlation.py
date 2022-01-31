import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from function_biotime_trends import plotting


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

invasions_t = {}
invasions_t1 = {}
extinctions_t = {}
extinctions_t1 = {}
richness = {}
years = {}



study_ids = presences["study_ids"]
for study_id in study_ids:
    pres = presences[study_id] == 1 # conversion to boolean
    # total species richness
    richness[study_id] = np.sum(pres, axis = 1)
    
    # not present last year, but present this year
    inv = np.sum(pres[1:] & ~pres[:-1], axis = 1)

    invasions_t1[study_id] = np.append(np.nan, inv)
    invasions_t[study_id] = np.append(inv, np.nan)
    
    
    
    # not present last year, but present this year
    ext = np.sum(pres[:-1] & ~pres[1:], axis = 1)
    extinctions_t1[study_id] = np.append(np.nan, ext)
    extinctions_t[study_id] = np.append(ext, np.nan)  

    # divide by total richness
    years[study_id] = presences["year_" + study_id]
    
#"""    
    
    
fig, ax = plt.subplots(4,6, figsize = (14,12), sharey = "row")

met_inv = plotting(invasions_t, invasions_t1, ax[:,0], presences, color = "k")
ax[0,0].set_xlabel("Invasions t")
ax[0,0].set_ylabel("Invasions t1")
met_inv = plotting(extinctions_t, extinctions_t1, ax[:,1], presences, color = "k")
ax[0,1].set_xlabel("Extinctions t")
ax[0,1].set_ylabel("Extinctions t1")
met_inv = plotting(invasions_t, extinctions_t, ax[:,2], presences, color = "k")
ax[0,2].set_xlabel("Invasions t")
ax[0,2].set_ylabel("Extinctions t")
met_inv = plotting(invasions_t1, extinctions_t1, ax[:,3], presences, color = "k")
ax[0,3].set_xlabel("Invasions t1")
ax[0,3].set_ylabel("Extinctions t1")
met_inv = plotting(invasions_t, extinctions_t1, ax[:,4], presences, color = "k")
ax[0,4].set_xlabel("Invasions t")
ax[0,4].set_ylabel("Extinctions t1")
met_inv = plotting(invasions_t1, extinctions_t, ax[:,5], presences, color = "k")
ax[0,5].set_xlabel("Invasions t1")
ax[0,5].set_ylabel("Extinctions t")


for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKLMNOPQRSTUVWXYZ"[i], loc = "left")
fig.tight_layout()

fig.savefig("Figure_invasion_extinction_correlation.pdf")
