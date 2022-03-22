import numpy as np
import matplotlib.pyplot as plt
import warnings

import assembly_functions as af

try:
    species_id = np.load("Data_average_trait_distance.npz")
    surv = species_id["surv"]
except FileNotFoundError:
    species_id = af.generate_species(n_coms = 600, years = 200)
    
    present, equi_all, surv = af.community_assembly(species_id)
    max_richness = np.amax(np.sum(surv, axis = -1))
    
    np.savez("Data_average_trait_distance.npz", **species_id,
             surv = surv)
    



time = np.arange(surv.shape[1])

fig, ax = plt.subplots(2,2, sharex = True, figsize = (9,9), sharey = "row")

fac = [np.sqrt(2), np.sqrt(2)]

for i in range(2):
    traits = species_id["loc"].copy()
    
    # remove other trophic level from data
    traits[species_id["level"] != i] = np.nan
    traits = np.repeat(traits[:,np.newaxis], surv.shape[1], axis = 1 )
    
    # remove speceis who are not present
    traits[~surv] = np.nan

    # compare to randomly drawn distribution from same meta-community
    ref = np.random.normal(0, 1, traits.shape)
    ref[np.isnan(traits)] = np.nan
    percents = [25,50,75]
    with warnings.catch_warnings(record = True):
        # plot trait mean over time
        percentiles = np.nanpercentile(np.nanmean(traits, axis = -1),
                                       percents, axis = 0)
        ax[0,i].plot(time, percentiles[1], color = "r", label = "Simulation")
        ax[0,i].fill_between(time, percentiles[0], percentiles[2], color = "r", alpha = 0.5)
        
        # compare to random draws
        percentiles = np.nanpercentile(np.nanmean(ref, axis = -1),
                                       percents, axis = 0)
        ax[0,i].plot(time, percentiles[1], color = "b", label = "Random")
        ax[0,i].fill_between(time, percentiles[0], percentiles[2], color = "b", alpha = 0.5)
        
        # plot trait variance over time
        percentiles = np.nanpercentile(np.nanvar(traits, axis = -1),
                                       percents, axis = 0)
        ax[1,i].plot(time, percentiles[1], color = "r")
        ax[1,i].fill_between(time, percentiles[0], percentiles[2], color = "r", alpha = 0.5)
        
        # compare to random draws
        percentiles = np.nanpercentile(np.nanvar(ref, axis = -1),
                                       percents, axis = 0)
        ax[1,i].plot(time, percentiles[1], color = "b")
        ax[1,i].fill_between(time, percentiles[0], percentiles[2], color = "b", alpha = 0.5)
    
ax[0,0].set_title("Basal species")
ax[0,1].set_title("Predator species")

ax[0,0].set_ylabel("Mean of tratis")
ax[1,0].set_ylabel("Variance of traits")
ax[1,0].set_xlabel("Time")
ax[1,1].set_xlabel("Time")
ax[1,0].set_xlim([100,time[-1]])
ax[1,0].set_ylim([0, None])

ax[0,0].legend()

fig.tight_layout()
fig.savefig("Figure_trait_differences.pdf")