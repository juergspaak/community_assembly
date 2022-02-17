import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
import functions_for_plotting as fp
from scipy.stats import gaussian_kde

species_id = af.generate_species(n_coms = 2, sigs = [0.5, 0.2, 0.2],
                                  years = 1000)
present, equi_all, surv = af.community_assembly(species_id)
"""
n_reg_species = 50
pred_species = np.random.choice(species_id["loc"][species_id["level"] == 1],
                                size = n_reg_species)
basal_species = np.random.choice(species_id["loc"][species_id["level"] == 0],
                                size = n_reg_species)

species_id["loc"][species_id["level"] == 0] = np.random.choice(basal_species,
                                            size = species_id["loc"].shape)[species_id["level"] == 0]
species_id["loc"][species_id["level"] == 1] = np.random.choice(pred_species,
                                            size = species_id["loc"].shape)[species_id["level"] == 1]"""
"""
# create regional species pool
x = np.linspace(-2,2, 1000)
n_fam = 3
p = (2+np.cos(x*np.pi/x[-1]*n_fam))
p = p/np.sum(p)
reg_species = np.random.choice(np.linspace(-2,2, 1000), replace = False, p = p,
                               size = 50)
plt.hist(reg_species)
species_id["loc"] = np.random.choice(reg_species, size = species_id["loc"].shape)
"""

# convert species presence to families
n_fam = 6 # number of families
# boundaries of the families
fam_range, dfam = np.linspace(-2,2,n_fam-1, retstep=True)
com = 0

rich_fam = np.zeros((2, surv.shape[1], n_fam), dtype = int)

# associate to each species a family
species_id["fam"] = np.zeros(species_id["loc"].shape)
for i in range(n_fam-1):
    species_id["fam"][species_id["loc"]>fam_range[i]] = i+1
    
# compute richness per family    
for level in range(2):
    for fam in range(n_fam):
        rich_fam[level, :, fam] = np.sum(surv[com] & # present
                                       (species_id["fam"][com] == fam)
                                       & (species_id["level"][com] == level),
                                       axis = 1)
        
# convert richness into "trait location" of families
trait_fam = np.full(rich_fam.shape + (np.amax(rich_fam), ), False, dtype = bool)
for i in range(np.amax(rich_fam)):
    trait_fam[...,i][rich_fam>i] = True

# remove correspondence to families
trait_fam = np.reshape(trait_fam, trait_fam.shape[:2] + (-1,))
loc_fam = fam_range[0] + np.arange(n_fam)*dfam
loc_fam = loc_fam[:,np.newaxis] + np.linspace(-0.8,0.2, np.amax(rich_fam))*dfam
loc_fam = loc_fam.reshape(-1)
        
##############################################################################
# plot results
time = np.arange(equi_all.shape[1])
fig, ax = plt.subplots(2,2, figsize = (14,10), sharex = True, sharey = "row")
fp.plot_traits([ax[0,0]], surv, species_id)
ax[0,0].set_ylabel("Species traits")
ax[0,0].set_xlim(time[[0,-1]])

ax[0,0].set_title("Species focus")
ax[0,1].set_title("Family focus")


for f in fam_range:
    ax[0,0].axhline(f, color = "k", zorder = 0)
    ax[0,1].axhline(f, color = "k", zorder = 0)
    
# plot family richness
for i in range(len(loc_fam)):
    ax[0,1].plot(time[trait_fam[0,:,i]],
                 np.full(np.sum(trait_fam[0, :, i]), -0.05 + loc_fam[i]), 'b,')
    ax[0,1].plot(time[trait_fam[1,:,i]],
                 np.full(np.sum(trait_fam[1, :, i]), 0.05 + loc_fam[i]), 'r,')
        
###############################################################################
# plot jaccard similarity

# reference time
t = 500
trait_spec = np.array([surv[com] & (species_id["level"][com] == i)
                       for i in range(2)])
jaccard_spec = np.sum(trait_spec & trait_spec[:,[t]], axis = -1)  \
                / np.sum(trait_spec | trait_spec[:,[t]], axis = -1)
ax[1,0].plot(time, jaccard_spec[0], 'b')
ax[1,0].plot(time, jaccard_spec[1], 'r')

jaccard_fam = np.sum(trait_fam[:,[t]] & trait_fam, axis = -1) \
                / np.sum(trait_fam | trait_fam[:,[t]], axis = -1)
ax[1,1].plot(time, jaccard_fam[0], 'b')
ax[1,1].plot(time, jaccard_fam[1], 'r')

ax[1,0].set_ylabel("Jaccard similarity")
ax[1,1].set_xlabel("Time")
ax[1,0].set_xlabel("Time")

fig.tight_layout()
fig.savefig("Figure_jaccard_similarity.pdf")


