import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
"""
species_id = af.generate_species(n_coms = 100, years = 1000, locs = np.zeros(3),
                     sigs = [0.5, 0.25, .8],
                     level = [0.5, 0.5, 0],
                     max_utils = [1.5, 3, 0.08])

present, equi_all, surv = af.community_assembly(species_id)
"""

species_id = species_id
present = present
surv = surv
max_richness = np.amax(np.sum(surv, axis = -1))

# number of families = richness/4
n_fam = int(np.round(np.mean(np.sum(surv[:,-1], axis = -1), axis = 0)/4))
range_fam = np.append(np.linspace(-2,2, n_fam), np.inf)
loc_fam = np.linspace(-2.5, 2.5, n_fam +1)

n_spec = 10*n_fam
range_spec = np.append(np.linspace(-2,2, n_spec), np.inf)
loc_spec = np.linspace(-2.5, 2.5, n_spec +1)

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
    traits = np.sort(traits, axis = -1)[...,:max_richness]
    traits = traits[:,100:] # cut off burning phase
    
    # richness per family
    rich_fam = np.empty((n_fam+1, ) + traits.shape[:2])
    for j in range(n_fam+1):
        rich_fam[j] = np.sum(traits<range_fam[j], axis = -1)
    rich_fam[1:] -= rich_fam[:-1]
    
    ax[0, i].plot(loc_fam, np.std(rich_fam, axis = (1,2))/np.mean(rich_fam, axis = (1,2)), "ro",
                  label = "Simulation")
    
    # compare to average
    ref = np.random.normal(0, 1, traits.shape)
    ref[np.isnan(traits)] = np.nan
    
    rich_fam = np.empty((n_fam+1, ) + traits.shape[:2])
    for j in range(n_fam+1):
        rich_fam[j] = np.sum(ref<range_fam[j], axis = -1)
    rich_fam[1:] -= rich_fam[:-1]
    ax[0, i].plot(loc_fam, np.std(rich_fam, axis = (1,2))/np.mean(rich_fam, axis = (1,2)), "bo",
                  label = "Randomness")
    
    
    # richness per species
    rich_spec = np.empty((n_spec+1, ) + traits.shape[:2])
    for j in range(n_spec+1):
        rich_spec[j] = np.sum(traits<range_spec[j], axis = -1)
    rich_spec[1:] -= rich_spec[:-1]
    
    ax[1, i].plot(loc_spec, np.std(rich_spec, axis = (1,2))/np.mean(rich_spec, axis = (1,2)), "ro")
    
    rich_spec = np.empty((n_spec+1, ) + traits.shape[:2])
    for j in range(n_spec+1):
        rich_spec[j] = np.sum(ref<range_spec[j], axis = -1)
    rich_spec[1:] -= rich_spec[:-1]
    
    ax[1, i].plot(loc_spec, np.std(rich_spec, axis = (1,2))/np.mean(rich_spec, axis = (1,2)), "bo")


ax[0,0].legend()

ax[0,0].set_title("Basal species")
ax[0,1].set_title("Predator species")

ax[0,0].set_ylabel("CV(Family richness)")
ax[1,0].set_ylabel("CV(Species richness)")

ax[1,0].set_xlabel("Trait location")
ax[1,1].set_xlabel("Trait location")

fig.savefig("Figure_family_vs_species_richness.pdf")
    
    
    