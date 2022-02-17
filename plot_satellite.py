import numpy as np
import matplotlib.pyplot as plt
import assembly_functions as af

import warnings


# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

species_id = af.generate_species(20,500, sigs = [0.5, 0.2, .8],
                                  max_utils = [0.5, 0.8, 0.08], locs = [0,0,0],
                                  sig_locs = [1,1,0], ms = [0.6, 0.3, 0])
    
present, equi_all, surv = af.community_assembly(species_id)

try:
    d_free
except NameError:
    d_free = np.load(path + "Data_LV_satelitte.npz")
    d_free = {key: d_free[key] for key in d_free.files}

d_free = species_id
d_free["equi_all"] = equi_all

    
fig, ax = plt.subplots(2,2)
com = 4
surv = d_free["equi_all"][com]>0
time = np.arange(d_free["equi_all"].shape[1])
n_spec = d_free["loc"].shape[1]
for i in range(n_spec):
    pres = surv[:,i]
    if not pres.any():
        continue # species is never present
    ax[0,0].plot(time[pres], np.full(np.sum(pres), d_free["loc"][com, i]),
                 "br"[d_free["level"][com, i]])
    
ax[0,0].set_ylabel("Species mean trait")
ax[0,0].set_xlabel("Time")

n_bins = 20
trait_bins, db = np.linspace(-2.2,2.2,n_bins-1, retstep = True)
trait_id = np.zeros(d_free["loc"].shape)
for b in trait_bins:
    trait_id[d_free["loc"]>b] += 1

invasions, extinctions = np.empty((2, n_bins))
persistence_time = np.full(n_bins, np.nan)
for i in range(n_bins):
    surv_bin = (d_free["equi_all"]>0) & ((d_free["level"] == 0) & (trait_id == i))[:,np.newaxis]
    persistence_temp = 1.0*np.sum(surv_bin, axis = 1)
    persistence_temp[persistence_temp == 0] = np.nan
    # only take species that arrived in the middle time, i.e. extinciont actually observed
    persistence_time[i] = np.nanmean(persistence_temp[:,200:-200])
    #inv_prob = np.mean(persi)
    """invasions[i] = np.nanmean((np.sum(surv_bin[:,1:] & ~surv_bin[:,:-1], axis = -1)
                    /np.sum(surv_bin[:,1:], axis = -1))[0,1000:])
    extinctions[i] = np.nanmean((np.sum(surv_bin[:,:-1] & ~surv_bin[:,1:], axis = -1) \
                    /np.sum(surv_bin[:,:-1], axis = -1))[:,1000:])"""

bin_locations = np.append(trait_bins[0]-db/2, trait_bins + db/2)

ax[0,1].plot(bin_locations, persistence_time, 'b.')
ax[0,1].set_xlabel("Trait location")
ax[0,1].set_ylabel("Persistence time")
ax[0,1].set_ylim([0,None])

fig.tight_layout()
fig.savefig("Figure_satelitte.pdf")