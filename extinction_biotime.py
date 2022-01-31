import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from function_biotime_trends import plotting


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

invasions = {}
extinctions = {}
richness = {}
years = {}

inv_all = []
ext_all = []
years_all = []

inv_ref = "current"
ext_ref = "future"
scaled = "_scaled"

study_ids = presences["study_ids"]
for study_id in study_ids:
    pres = presences[study_id] == 1 # conversion to boolean
    # total species richness
    richness[study_id] = np.sum(pres, axis = 1)
    
    # not present last year, but present this year
    invasions[study_id] = np.sum(pres[1:] & ~pres[:-1], axis = 1)
    if inv_ref == "current": # of all present, which did invade?
        invasions[study_id] = np.append(np.nan, invasions[study_id])
    else:
        invasions[study_id] = np.append(invasions[study_id], np.nan)
    
    
    
    # not present last year, but present this year
    extinctions[study_id] = np.sum(pres[:-1] & ~pres[1:], axis = 1)
    if ext_ref == "current": # of all present, which did invade?
        extinctions[study_id] = np.append(np.nan, extinctions[study_id])
    else:
        extinctions[study_id] = np.append(extinctions[study_id], np.nan)  

    # divide by total richness
    if scaled:
        invasions[study_id] = invasions[study_id]/richness[study_id]                    
        extinctions[study_id] = extinctions[study_id]/richness[study_id]
    years[study_id] = presences["year_" + study_id]
    
    years_all.extend(years[study_id])
    inv_all.extend(invasions[study_id])
    ext_all.extend(extinctions[study_id])

years_all = np.array(years_all)
inv_all = np.array(inv_all)
ext_all = np.array(ext_all)
#"""    
    
    
fig, ax = plt.subplots(4,3, figsize = (10,10), sharey = "row")

met_inv = plotting(years, invasions, ax[:,0], presences, color = "k")
met_ext = plotting(years, extinctions, ax[:,1], presences, color = "k")
met_inv_ext = plotting(invasions, extinctions, ax[:,2], presences, color = "k")

##############################################################################
# layout

ax[0,0].set_title(inv_ref)
ax[0,1].set_title(ext_ref)

ax[0,0].set_xlabel("Year")
ax[0,1].set_xlabel("Year")
ax[0,2].set_xlabel("Invasion probability")

ax[0,0].set_ylabel("Invasion probability")
ax[0,1].set_ylabel("Extinction probability")
ax[0,2].set_ylabel("Extinction probability")

for i in range(3):
    ax[2,i].semilogx()


for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKL"[i], loc = "left")
    
fig.tight_layout()

fig.savefig("Figure_biotime_changes_{}_{}{}.pdf".format(inv_ref, ext_ref, scaled))