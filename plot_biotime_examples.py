import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings
from matplotlib.cm import Set1

colors = Set1(np.linspace(0,1,9))

# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

try:
    
    datas
except NameError:
    d_prey = np.load(path + "data_LV_assembly.npz")
    d_pred = np.load(path + "data_LV_assembly_predators.npz")
    d_prey = {key: d_prey[key] for key in d_prey.files}
    d_pred = {key: d_pred[key] for key in d_pred.files}


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

study_ids = presences["study_ids"]
# only work with long datasets
study_ids = study_ids[[presences[study_id].shape[0]>=40
                      for study_id in study_ids]]

path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')
organisms = [meta_data[meta_data.STUDY_ID == int(sid)].TAXA.values[0] for sid in study_ids]

keys = ["Birds","Marine\ninvertebrates",
                        "Terrestrial\nplants", "Fish", "Other"]
colors = {key: colors[i] for i, key in enumerate(keys)}
colors["Mammals"] = colors["Other"]
colors["Benthos"] = colors["Other"]
colors["All"] = colors["Other"]
colors["Marine invertebrates"] = colors["Marine\ninvertebrates"]
colors["Terrestrial plants"] = colors["Terrestrial\nplants"]


fig, ax = plt.subplots(4,3,figsize = (12,12), sharey = "row")
persist_max = 100
bins = np.arange(1, 100, 5)

##############################################################################
# case for prey only
surv = d_prey["equi_all"]>0
time = np.arange(surv.shape[1])
com = 0
ax[0,0].plot(time, np.sum(surv[com], axis = -1), 'b--', alpha = 0.2,
             label = "Example")
ax[0,0].plot(time, np.mean(np.sum(surv, axis = -1), axis = 0), 'b',
             label = "Average case")
ax[0,0].legend()

# invasion probability
invasion = np.sum(surv[:,1:] & (~surv[:,:-1]), axis = -1)/np.sum(surv[:,:-1], axis = -1)
ax[1,0].plot(time[:-1], invasion[com], 'b--', alpha = 0.2)
ax[1,0].plot(time[:-1], np.mean(invasion, axis = 0), 'b')

# extinction porbability
extinction = np.sum(~surv[:,1:] & surv[:,:-1], axis = -1)/np.sum(surv[:,1:], axis = -1)
ax[2,0].plot(time[:-1], extinction[com], 'b--', alpha = 0.2)
ax[2,0].plot(time[:-1], np.mean(extinction, axis = 0), 'b')

# persistence time
persistence_time = np.sum(surv[:,100:], axis = 1)
hist, bins = np.histogram(persistence_time.flatten(), bins = bins,
                          density = True)
ax[3,0].plot(bins[:-1], hist, 'b')

##############################################################################
# case for prey and pred
surv = d_pred["equi_all"]>0
time = np.arange(surv.shape[1])
surv = np.array([surv & (d_pred["level"][:,np.newaxis] == 0),
                 surv & (d_pred["level"][:,np.newaxis]== 1)])

com = 0
label = ["Prey", "Predator"]
for i in range(2):
    ax[0,1].plot(time, np.sum(surv[i, com], axis = -1),
                 "br"[i] + '--', alpha = 0.2)
    ax[0,1].plot(time, np.mean(np.sum(surv[i], axis = -1), axis = 0), 'br'[i]
                 , label = label[i])
    
    # invasion probability
    with warnings.catch_warnings(record = True):
        invasion = np.sum(surv[i,:,1:] & (~surv[i,:,:-1]), axis = -1)/np.sum(surv[i,:,1:], axis = -1)
        ax[1,1].plot(time[:-1], invasion[com], 'br'[i]+'--', alpha = 0.2)
        ax[1,1].plot(time[:-1], np.nanmean(invasion, axis = 0), 'br'[i])
    
    # extinction porbability
    with warnings.catch_warnings(record = True):
        extinction = np.sum(~surv[i,:,1:] & surv[i,:,:-1], axis = -1)/np.sum(surv[i,:,:-1], axis = -1)
        #extinction[np.isfinite(extinction)] = np.nan
        ax[2,1].plot(time[:-1], extinction[com], 'br'[i] +'--', alpha = 0.2)
        ax[2,1].plot(time[:-1], np.nanmean(extinction, axis = 0), 'br'[i])

ax[0,1].legend()
# persistence time
persistence_time = np.sum(surv[:,:,100:], axis = 2)
hist, bins = np.histogram(persistence_time[0], bins = bins,
                          density = True)
ax[3,1].plot(bins[:-1], hist, 'b')
hist, bins = np.histogram(persistence_time[1], bins = bins,
                          density = True)
ax[3,1].plot(bins[:-1], hist, 'r')

##############################################################################
# biotime data

alpha = 0.5
for study_id in study_ids:
    pres = presences[study_id]
    col = colors[meta_data[meta_data.STUDY_ID == int(study_id)].TAXA.values[0]]
    ax[0,2].plot(presences["year_{}".format(study_id)],
                 np.sum(pres, axis = 1), '-', alpha = alpha, color = col)
    invasion = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[1:], axis = 1)
    ax[1,2].plot(presences["year_{}".format(study_id)][1:],
                 invasion, color = col, alpha = alpha)
    extinction = np.sum(pres[:-1] & ~pres[1:], axis = 1)/np.sum(pres[:-1], axis = 1)
    ax[2,2].plot(presences["year_{}".format(study_id)][:-1],
                 extinction, color = col, alpha = alpha)
    persistence_time = np.nansum(presences[study_id], axis = 0)
    hist, bins = np.histogram(persistence_time, bins = bins,
                              density = True)
    #hist[bins[:-1]>=len(pres)] = 0
    ax[3,2].plot(bins[:-1], hist, '.', color = col, alpha = alpha)

for key in keys:
    ax[3,2].plot(np.nan, np.nan, color = colors[key], label = key)
ax[3,2].legend(ncol = 2, fontsize = 10)
ax[3,2].set_ylim([1e-3, 8])

##############################################################################
# add layout

ax[3,0].semilogy()
ax[0,0].semilogy()

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Proportion of invaders")
ax[2,0].set_ylabel("Proportion of extinctions")
ax[3,0].set_ylabel("Frequency")

# xlabel
for a in ax[:-1].flatten():
    a.set_xlabel("Time")

ax[3,0].set_xlabel("Persistence time")
ax[3,1].set_xlabel("Persistence time")
ax[3,2].set_xlabel("Persistence time")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKLMNOP"[i], loc = "left")
    
for i, a in enumerate(ax[:3,:2].flatten()):
    a.set_xlim(time[[0,-1]])
for i, a in enumerate(ax[-1]):
    a.set_xlim(bins[[0,-1]])

ax[1,0].set_ylim([0,1])
ax[2,0].set_ylim([0,1])
    
ax[0,0].set_title("Prey species\nonly")
ax[0,1].set_title("Predator and\nprey speceis")
ax[0,2].set_title("BioTime data")

fig.tight_layout()
fig.savefig("Figure_biotime_examples.pdf")
#"""

