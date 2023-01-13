import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from scipy.stats import linregress

from functions_for_plotting import biotime_colors as colors, keys

# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data


hist_plot = False

try:
    d_prey
    d_pred
except NameError:
    d_prey = np.load(path + "data_LV_assembly2.npz")
    d_pred = np.load(path + "data_LV_assembly_predators2.npz")
    d_prey = {key: d_prey[key] for key in d_prey.files}
    d_pred = {key: d_pred[key] for key in d_pred.files}

if False:
    # command to regenerate the stored files
    import assembly_functions as ap
    years = 800
    n_coms = 100
    range_omega = [2,5]
    omega = np.random.uniform(*range_omega,n_coms)
    omega[:5] = np.linspace(*range_omega, 5)

    # prey only     
    species_id = ap.generate_species(n_coms, years, level = [1,0], omega = omega)
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    np.savez(path + "data_LV_assembly2.npz", equi_all = equi_all, present = present,
             **species_id)
    
    species_id = ap.generate_species(n_coms, years, level = [1,1], omega = omega)
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    np.savez(path + "data_LV_assembly_predators2.npz", equi_all = equi_all, present = present,
             **species_id)
    
    d_prey = np.load(path + "data_LV_assembly2.npz")
    d_pred = np.load(path + "data_LV_assembly_predators2.npz")
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

fig, ax = plt.subplots(4,3,figsize = (14,12), sharey = "row")
fs = 16
persist_max = 100
bins = np.arange(1, 200, 2)
#study_ids = study_ids[:2]

##############################################################################
# case for prey only
surv = d_prey["equi_all"]>0
time = np.arange(surv.shape[1])
com = 0
ax[0,0].plot(time, np.sum(surv[com], axis = -1), 'b--', alpha = 0.2,
             label = "Example")
ax[0,0].plot(time, np.mean(np.sum(surv, axis = -1), axis = 0), 'b',
             label = "Average case")
ax[0,0].legend(fontsize = fs)

# invasion probability
invasion = np.sum(surv[:,1:] & (~surv[:,:-1]), axis = -1)/np.sum(surv[:,:-1], axis = -1)
ax[1,0].plot(time[:-1], invasion[com], 'b--', alpha = 0.2)
ax[1,0].plot(time[:-1], np.mean(invasion, axis = 0), 'b')

# extinction porbability
extinction = np.sum(~surv[:,1:] & surv[:,:-1], axis = -1)/np.sum(surv[:,1:], axis = -1)
ax[2,0].plot(time[:-1], extinction[com], 'b--', alpha = 0.2)
ax[2,0].plot(time[:-1], np.mean(extinction, axis = 0), 'b')

# persistence time
persistence_time = np.sum(surv[:,600:], axis = -2).reshape(-1)
persist = np.arange(np.amax(persistence_time))+1
percent_persist = np.sum(persistence_time[...,np.newaxis] >= persist,
                         axis = 0)

if hist_plot:
    # version with histograms
    hist, bins = np.histogram(persistence_time.flatten(), bins = bins,
                              density = True)
    ax[3,0].plot(bins[:-1], hist, 'b')
else:
    ax[3,0].plot(persist, percent_persist/percent_persist[0], 'b')

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

ax[0,1].legend(fontsize = fs)
# persistence time
persistence_time = np.sum(surv[:,:,400:], axis = -2)
persist = np.arange(np.amax(persistence_time))+1
percent_persist = np.sum(persistence_time[...,np.newaxis] >= persist,
                         axis = -2)

# persistence times for different communities
for i in range(5):
    ax[3,1].plot(persist, percent_persist[0,i]/percent_persist[0,i,0], 'b--'
                 , alpha = 0.5)
    ax[3,1].plot(persist, percent_persist[1,i]/percent_persist[1,i,0], 'r--'
                 , alpha = 0.5)

# average over multiple communities
percent_persist = np.sum(percent_persist, axis = 1)
ax[3,1].plot(persist, percent_persist[0]/percent_persist[0,0], 'b')
ax[3,1].plot(persist, percent_persist[1]/percent_persist[1,0], 'r')

##############################################################################
# biotime data
R2 = []
R2_log = []
r_max = 1
alpha = 0.5
for study_id in study_ids:
    # load dataset, shape (years, species)
    pres = presences[study_id]
    # get color of taxa
    col = colors[meta_data[meta_data.STUDY_ID == int(study_id)].TAXA.values[0]]
    # for plotting by realm
    #col = colors[meta_data[meta_data.STUDY_ID == int(study_id)].REALM.values[0]]
    
    # plot species richness
    ax[0,2].plot(presences["year_{}".format(study_id)],
                 np.sum(pres, axis = 1), '-', alpha = alpha, color = col)
    r_max = max(np.append(np.sum(pres, axis = 1), r_max))
    # plot portion of invaders
    invasion = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[1:], axis = 1)
    ax[1,2].plot(presences["year_{}".format(study_id)][1:],
                 invasion, color = col, alpha = alpha)
    # plot portion of extinctions
    extinction = np.sum(pres[:-1] & ~pres[1:], axis = 1)/np.sum(pres[:-1], axis = 1)
    ax[2,2].plot(presences["year_{}".format(study_id)][:-1],
                 extinction, color = col, alpha = alpha)
    ### plot persistence time
    pres = pres == 1
    # when did species go extinct?
    ext = pres[:-1] & ~pres[1:]
    # when did species invade?
    inv = pres[1:] & ~pres[:-1]
    
    persist_total = np.cumsum(pres, axis = 0)
    deaths = [] # when death occured
    all_data = [] # includes censored data for species with unknown total survival
    for i in range(pres.shape[1]):
        if not ext[:,i].any():
            # species never goes extinct, only note maximum persistence time
            all_data.append(np.sum(pres[:,i]))
            continue
        temp = persist_total[1:, i][ext[:,i]]
        # add potential observation where we did not observe the extinction
        if pres[-1,i]:
            all_data.append(persist_total[-1,i] - np.max(temp))
        
        # conmpute survival time of each individual occurence
        temp[1:] -= temp[:-1]
        if pres[0,i]:
            # censor first deat observation, could have been there earlier
            deaths.extend(temp[1:])
            # shorten first occurence by 1, as death point is unknown
            temp[0] -= 1
            all_data.extend(temp)
        else:
            deaths.extend(temp)
            all_data.extend(temp)
        
    
    if True:
        # when did species go extinct?
        ext = pres[:-1] & ~pres[1:]
        ext[-1] = 1 # include last time point for all communities
        persist_total = np.cumsum(pres, axis = 0)
        persistence_t = []
        for i in range(pres.shape[1]):
            temp = persist_total[1:, i][ext[:,i] == 1]
            temp[1:] -= temp[:-1]
            persistence_t.extend(temp)
    
    
    
    if hist_plot:
        # old histogram plots
        hist, bins = np.histogram(persistence_t, bins = bins,
                                  density = True)
        s,i,r2,p,std = linregress(bins[:-1], hist)
        R2.append(r2**2)
        s,i,r2,p,std = linregress(np.log(bins[:-1]), hist)
        R2_log.append(r2**2)
        
        #hist[bins[:-1]>=len(pres)] = 0
        ax[3,2].plot(bins[:-1], hist, '.', color = col, alpha = alpha)
    else:
        persist = np.arange(max(persistence_t))+1
        d_bins = np.arange(max(all_data)+1)
        d_hist,empty = np.histogram(deaths, bins = d_bins)
        a_hist,empty = np.histogram(all_data, bins = d_bins)
        #a_hist,empty = np.histogram(deaths, bins = d_bins)
        d_cum = np.cumsum(d_hist)
        a_cum = np.cumsum(a_hist[::-1])[::-1]
        prob_death = d_hist/a_cum
        survival = np.cumprod(1-prob_death)
        ax[3,2].plot(d_bins[1:], survival, color = col, linestyle = "--") 
    

for key in keys:
    ax[3,2].plot(np.nan, np.nan, color = colors[key], label = key)
ax[3,2].legend(fontsize = 10, loc = "upper right")
ax[3,2].set_xlim([1,200])

##############################################################################
# add layout

# y-ticks
ax[0,0].set_ylim([1,500])
ax[1,0].set_ylim([1e-3,1])
ax[2,0].set_ylim([1e-3,1])
ax[3,0].set_ylim([1e-2,1.1])


if hist_plot:
    ax[3,0].semilogx()
    ax[3,1].semilogx()
    ax[3,2].semilogx()

for a in ax[[0,3],0]:
    a.semilogy()

ax[1,0].set_yticks([0,1])
ax[2,0].set_yticks([0,1])

#ax[3,0].set_yticks([])
ax[3,1].set_yticks([1e-3, 1e-2, 1e-1, 1])
#ax[3,2].semilogx()

ax[0,0].set_ylabel("Species\nrichness", fontsize = fs)
ax[1,0].set_ylabel("Proportion of\ninvaders", fontsize = fs)
ax[2,0].set_ylabel("Proportion of\nextinctions", fontsize = fs)
ax[3,0].set_ylabel("Frequency", fontsize = fs)

# xlabel
for a in ax[:-1].flatten():
    a.set_xlabel("Time", fontsize = fs)

for i in range(3):
    ax[3,i].set_xlabel("Persistence Time", fontsize = fs)
    ax[3,i].semilogx()

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKLMNOP"[i], loc = "left", fontsize = fs)



for i, a in enumerate(ax[:3,:2].flatten()):
    a.set_xlim([0,500])
    

    


    
ax[0,0].set_title("Prey species\nonly", fontsize = fs)
ax[0,1].set_title("Predator and\nprey species", fontsize = fs)
ax[0,2].set_title("BioTime data", fontsize = fs)

fig.tight_layout()
fig.savefig("Figure_biotime_examples_hist.pdf")
#"""

