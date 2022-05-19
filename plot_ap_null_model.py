import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings


# load simulated data
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
try:
    d_prey
    d_pred
except NameError:
    d_pred = np.load(path + "data_LV_assembly_predators.npz")
    d_pred = {key: d_pred[key] for key in d_pred.files}

n_coms, years, n_spec = d_pred["equi_all"].shape


# persistence time, ignore initial communities
persist_time = 1.0*np.sum(d_pred["equi_all"][:,years//2:, years//2:]>0, axis = 1)
persist_time[persist_time == 0] = np.nan

# average persistence time
av_persist = np.nanmean(persist_time)

# get invasion success, ignore initial 100 communities
prob_invasion = np.mean(np.isfinite(persist_time))
prob_invasion = 0.8

###############################################################################
# generate communities with non-interacting species
# community assembly based on pure randomness

# persistence time is exponential distribution
persist = np.round(np.random.exponential(av_persist, size = (n_coms, 1, n_spec)))
# certain amount of species can't invade at all
persist[np.random.uniform(size = persist.shape)>prob_invasion] = 0
# when do species invade
inv_time = np.random.randint(2*years, size = (n_coms, 1, n_spec))

# to store species presences
present = np.full((n_coms, years, n_spec), False, dtype = bool)

# each species invades at its timepoint
present[:, np.arange(n_spec)<=np.arange(years)[:,np.newaxis]] = True
#species go extinct after their persistence time
present[(np.arange(n_spec) + persist)<=np.arange(years)[:,np.newaxis]] = False

# each species invades at its timepoint
present[inv_time<=np.arange(years)[:,np.newaxis]] = True
#species go extinct after their persistence time
present[(inv_time + persist)<=np.arange(years)[:,np.newaxis]] = False

# compute richness
richness = np.sum(present, axis = -1)

###############################################################################
# plot results
fig, ax = plt.subplots(4,2,figsize = (12,12), sharey = "row", sharex = "row")
persist_max = 100
bins = np.arange(1, 100, 2)

# percentiles for plotting
q = 1.98

##############################################################################
# Null model
surv = present
time = np.arange(surv.shape[1])

# compute richness
rich = np.sum(surv, axis = 1)
mean_rich = np.nanmean(rich, axis = 0)
std_rich = np.nanstd(rich, axis = 0)
ax[0,1].plot(time[:-1], mean_rich, 'b')
ax[0,1].fill_between(time[:-1], mean_rich - std_rich*q, mean_rich + std_rich*q) 

# invasion probability
with warnings.catch_warnings(record = True):
    invasion_null = np.sum(surv[:,1:] & (~surv[:,:-1]), axis = -1)/np.sum(surv[:,1:], axis = -1)
    mean_inv = np.nanmean(invasion_null, axis = 0)
    std_inv = np.nanstd(invasion_null, axis = 0)
    ax[1,1].plot(time[:-1], mean_inv, 'b')
    ax[1,1].fill_between(time[:-1], mean_inv - std_inv*q, mean_inv + std_inv*q)

# extinction porbability
with warnings.catch_warnings(record = True):
    extinction_null = np.sum(~surv[:,1:] & surv[:,:-1], axis = -1)/np.sum(surv[:,:-1], axis = -1)
    #extinction[np.isfinite(extinction_null)] = np.nan
    mean_ext = np.nanmean(extinction_null, axis = 0)
    std_ext = np.nanstd(extinction_null, axis = 0)
    ax[2,1].plot(time[:-1], mean_ext, 'b')
    ax[2,1].fill_between(time[:-1], mean_ext - std_ext*q, mean_ext + std_ext*q)

persistence_time = np.sum(surv[:,100:], axis = 1)
hist, bins = np.histogram(persistence_time, bins = bins,
                          density = True)
ax[3,1].plot(bins[:-1], hist, 'b')

##############################################################################
# actual simulation
surv = d_pred["equi_all"]>0
time = np.arange(surv.shape[1])

# compute richness
rich = np.sum(surv, axis = 1)
mean_rich = np.nanmean(rich, axis = 0)
std_rich = np.nanstd(rich, axis = 0)
ax[0,0].plot(time[:-1], mean_rich, 'b')
ax[0,0].fill_between(time[:-1], mean_rich - std_rich*q, mean_rich + std_rich*q) 

# invasion probability
with warnings.catch_warnings(record = True):
    invasion = np.sum(surv[:,1:] & (~surv[:,:-1]), axis = -1)/np.sum(surv[:,1:], axis = -1)
    mean_inv = np.nanmean(invasion, axis = 0)
    std_inv = np.nanstd(invasion, axis = 0)
    ax[1,0].plot(time[:-1], mean_inv, 'b')
    ax[1,0].fill_between(time[:-1], mean_inv - std_inv*q, mean_inv + std_inv*q)

# extinction porbability
with warnings.catch_warnings(record = True):
    extinction = np.sum(~surv[:,1:] & surv[:,:-1], axis = -1)/np.sum(surv[:,:-1], axis = -1)
    #extinction[np.isfinite(extinction)] = np.nan
    mean_ext = np.nanmean(extinction, axis = 0)
    std_ext = np.nanstd(extinction, axis = 0)
    ax[2,0].plot(time[:-1], mean_ext, 'b')
    ax[2,0].fill_between(time[:-1], mean_ext - std_ext*q, mean_ext + std_ext*q)

persistence_time = np.sum(surv[:,100:], axis = 1)
hist, bins = np.histogram(persistence_time, bins = bins,
                          density = True)
ax[3,0].plot(bins[:-1], hist, 'b')

##############################################################################
# add layout

ax[3,1].set_yticks([1e-3, 1e-2, 1e-1, 1])
ax[3,0].semilogy()
    
ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Proportion of invaders")
ax[2,0].set_ylabel("Proportion of extinctions")
ax[3,1].set_ylabel("Frequency")

# xlabel
for a in ax[:-1].flatten():
    a.set_xlabel("Time")

ax[3,1].set_xlabel("Persistence time")
ax[3,0].set_xlabel("Persistence time")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKLMNOP"[i], loc = "left")
    
for i, a in enumerate(ax[:3,:2].flatten()):
    a.set_xlim(time[[100,-100]])
for i, a in enumerate(ax[-1]):
    a.set_xlim(bins[[0,-1]])

ax[0,0].set_ylim([0, 40])
ax[1,0].set_ylim([0, 0.5])
ax[2,0].set_ylim([0, 0.5])

    
ax[0,0].set_title("Null model")
ax[0,1].set_title("Simulation model")

fig.tight_layout()
fig.savefig("Figure_ap_null_model.pdf")
#"""


##############################################################################
# do we get a correlation between invasion and extinction
plt.figure()
corr_null = np.empty(len(invasion_null))
for i in range(len(invasion_null)):
    ind = np.isfinite(invasion_null[i])*np.isfinite(extinction_null[i])
    ind[:100] = False
    ind[-100:] = False
    corr_null[i] = np.corrcoef(invasion_null[i, ind], extinction_null[i, ind])[0,1]
    
plt.hist(corr_null,
         label = "Null model")
plt.hist([np.corrcoef(invasion[i,100:-100], extinction[i,100:-100])[0,1] for i in range(len(invasion_null))],
         label = "Community model")
plt.legend()
