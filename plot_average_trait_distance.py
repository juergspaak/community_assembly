import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
"""
species_id = af.generate_species(n_coms = 100, years = 1000, locs = np.zeros(3),
                     sigs = [1, 0.5, .8],
                     level = [0.5, 0.5, 0],
                     max_utils = [1.5, 3, 0.08])

present, equi_all, surv = af.community_assembly(species_id)
max_richness = np.amax(np.sum(surv, axis = -1))"""

time = np.arange(surv.shape[1])

fig, ax = plt.subplots(3,2, sharex = True, figsize = (9,9), sharey = "row")

fac = [np.sqrt(2), np.sqrt(2)]

for i in range(2):
    traits = species_id["loc"].copy()
    # remove other trophic level from data
    traits[species_id["level"] != i] = np.nan
    traits = np.repeat(traits[:,np.newaxis], surv.shape[1], axis = 1 )
    
    # remove speceis who are not present
    traits[~surv] = np.nan
    traits = np.sort(traits, axis = -1)[...,:max_richness]
    
    # compute trait differences
    trait_diffs = traits[...,1:]-traits[...,:-1]
    
    # plot percentiles of trait differences
    ax[0,i].plot(time, np.nanpercentile(trait_diffs, [1,25, 50, 75, 99], axis = (0,-1)).T,
             color = "r")
    ax[0,i].axhline(2*fac[i]*np.mean(species_id["sig"][i]), color = "g",
                    label = "limiting similarity")
    
    # compare to randomly drawn distribution from same meta-community
    ref = np.random.normal(0, 1, traits.shape)
    ref[np.isnan(traits)] = np.nan
    ref = np.sort(ref, axis = -1)
    diff_ref = ref[..., 1:] - ref[...,:-1]
    
    ax[0,i].plot(time, np.nanpercentile(diff_ref, [1,25, 50, 75, 99], axis = (0,-1)).T,
             color = "b", alpha = 0.5)
    
    # plot trait mean over time
    ax[1,i].plot(time, np.nanpercentile(np.nanmean(traits, axis = -1), [25,50,75], axis = 0).T,
                 color = "r")
    ax[1,i].plot(time, np.nanpercentile(np.nanmean(ref, axis = -1), [25,50,75], axis = 0).T,
                 color = "b", alpha = 0.5)
    
    # plot trait variance over time
    ax[2,i].plot(time, np.nanpercentile(np.nanstd(traits, axis = -1), [25,50,75], axis = 0).T,
                 color = "r", label = "Simulation")
    ax[2,i].plot(time, np.nanpercentile(np.nanstd(ref, axis = -1), [25,50,75], axis = 0).T,
                 color = "b", alpha = 0.5, label = "Random")
    
ax[0,0].set_ylabel("Trait difference")
ax[0,0].set_title("Basal species")
ax[0,1].set_title("Predator species")

ax[1,0].set_ylabel("Mean tratis")
ax[2,0].set_ylabel("Std(traits)")
ax[2,0].set_xlabel("Time")
ax[2,1].set_xlabel("Time")
ax[0,1].legend()
ax[0,0].legend()

ax[2,0].legend(["Sim 75 perc.", "Sim 50 perc", "Sim 25 perc", "Rand 75", "Rand 50", "Rand 25"])

fig.tight_layout()
fig.savefig("Figure_trait_differences.pdf")