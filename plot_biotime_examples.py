import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from functions_for_plotting import biotime_colors as colors, keys

# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

try:
    d_prey
    d_pred
except NameError:
    d_prey = np.load(path + "data_LV_assembly.npz")
    d_pred = np.load(path + "data_LV_assembly_predators.npz")
    d_prey = {key: d_prey[key] for key in d_prey.files}
    d_pred = {key: d_pred[key] for key in d_pred.files}


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz", allow_pickle=True)

study_ids = presences["study_ids"]
# only work with long datasets
study_ids = study_ids[[presences[study_id].shape[0]>=30
                      for study_id in study_ids]]

path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape')

fig, ax = plt.subplots(4,3,figsize = (12,12))
persist_max = 100
bins = np.arange(1, 100, 2)

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
ax[3,0].set_frame_on(False)
ax[3,0].set_xticks([])

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
    # load dataset, shape (years, species)
    pres = presences[study_id]
    # get color of taxa
    col = colors[meta_data[meta_data.STUDY_ID == int(study_id)].TAXA.values[0]]
    
    # plot species richness
    ax[0,2].plot(presences["year_{}".format(study_id)],
                 np.sum(pres, axis = 1), '-', alpha = alpha, color = col)
    # plot portion of invaders
    invasion = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[1:], axis = 1)
    ax[1,2].plot(presences["year_{}".format(study_id)][1:],
                 invasion, color = col, alpha = alpha)
    # plot portion of extinctions
    extinction = np.sum(pres[:-1] & ~pres[1:], axis = 1)/np.sum(pres[:-1], axis = 1)
    ax[2,2].plot(presences["year_{}".format(study_id)][:-1],
                 extinction, color = col, alpha = alpha)
    # plot persistence time
    persistence_time = np.nansum(presences[study_id], axis = 0)
    ext = pres[:-1] & ~pres[1:]
    ext[-1] = 1 # include last time point for all communities
    persist_total = np.cumsum(pres, axis = 0)
    persistence_t = []
    for i in range(pres.shape[1]):
        temp = persist_total[1:, i][ext[:,i] == 1]
        temp[1:] -= temp[:-1]
        persistence_t.extend(temp)
    hist, bins = np.histogram(persistence_t, bins = bins,
                              density = True)
    #hist[bins[:-1]>=len(pres)] = 0
    ax[3,2].plot(bins[:-1], hist, '.', color = col, alpha = alpha)

for key in keys:
    ax[3,2].plot(np.nan, np.nan, color = colors[key], label = key)
ax[3,2].legend(ncol = 2, fontsize = 10)

##############################################################################
# add layout

# y-ticks
for i in range(3):
    ax[0,i].set_ylim([1,200])
    ax[0,i].semilogy()
    ax[1,i].set_ylim([0,1])
    ax[2,i].set_ylim([0,1])
    ax[3,i].set_ylim([1e-3, 8])
    

for i in range(4):
    for j in range(3):
        if j == 0:
            continue
        ax[i,j].set_yticklabels([])

ax[3,0].set_yticks([])
ax[3,1].set_yticks([1e-3, 1e-2, 1e-1, 1])
ax[3,1].semilogy()
ax[3,2].semilogy()
    
ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Proportion of invaders")
ax[2,0].set_ylabel("Proportion of extinctions")
ax[3,1].set_ylabel("Frequency")



# xlabel
for a in ax[:-1].flatten():
    a.set_xlabel("Time")

ax[3,1].set_xlabel("Persistence time")
ax[3,2].set_xlabel("Persistence time")

for i, a in enumerate(ax[:3].flatten()):
    a.set_title("ABCDEFGHIJKLMNOP"[i], loc = "left")

ax[3,1].set_title("J", loc = "left")
ax[3,2].set_title("K", loc = "left")
    
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

