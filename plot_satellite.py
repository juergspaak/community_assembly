import numpy as np
import matplotlib.pyplot as plt
import assembly_functions as af

import warnings


# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

omega = 3
species_id = af.generate_species(2,10000, omega = omega, ms = [0.01,0.5],
                                 sigs = [0.5,0.2], level = [1,2],
                                 mu_max = 1.2)
#species_id = af.generate_species(50,1000)
    
present, equi_all, surv = af.community_assembly(species_id)

try:
    d_free
except NameError:
    d_free = np.load(path + "Data_LV_satelitte.npz")
    d_free = {key: d_free[key] for key in d_free.files}

d_free = species_id
d_free["equi_all"] = equi_all

    
fig, ax = plt.subplots(2,2)
com = 1
surv = d_free["equi_all"][com]>0
time = np.arange(d_free["equi_all"].shape[1])
n_spec = d_free["loc"].shape[1]
tresh = 9000
time = time[tresh:]
for i in range(n_spec):
    pres = surv[tresh:,i]
    if not pres.any():
        continue # species is never present
    ax[0,0].plot(time[pres], np.full(np.sum(pres), d_free["loc"][com, i]),
                 "br"[d_free["level"][com, i]])
    
ax[0,0].set_ylabel("Species mean trait")
ax[0,0].set_xlabel("Time")

n_bins = 40
trait_bins, db = np.linspace(-2*omega, 2*omega,n_bins-1, retstep = True)
trait_id = np.zeros(d_free["loc"].shape)
for b in trait_bins:
    trait_id[d_free["loc"]>b] += 1

invasions, extinctions = np.empty((2, n_bins))
persistence_time = np.full(n_bins, np.nan)
persistence_per_species = 1.0*np.sum(d_free["equi_all"]>0, axis = 1)
persistence_per_species[persistence_per_species == 0] = np.nan
persistence_per_species[d_free["level"] == 1] = np.nan

plt.hist(np.log(persistence_per_species[:,500:-200]).flatten(), bins = 20)
a, b = 100, 200
for i in range(n_bins):
    # only take species that arrived in the middle time, i.e. extinciont actually observed
    persistence_time[i] = np.nanmean(persistence_per_species[:,a:-b][trait_id[:,a:-b] == i])

bin_locations = np.append(trait_bins[0]-db/2, trait_bins + db/2)

ax[0,1].plot(bin_locations, persistence_time, 'b.')
ax[0,1].set_xlabel("Trait location")
ax[0,1].set_ylabel("Persistence time")
ax[0,1].set_ylim([0,None])

fig.tight_layout()
fig.savefig("Figure_satelitte.pdf")
#'''