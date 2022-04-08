import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp
import assembly_functions as ap
"""
itera = 30
sigs = np.linspace(0.4, 3, itera)

richness = np.empty((2,itera,itera))
similarity = np.empty((itera, itera))

years = 2000
omega = 4
for i, sigi in enumerate(sigs):
    for j, sigj in enumerate(sigs):
        species_id = ap.generate_species(1, years, level = [1,1],
                                 omega = omega, sigs = [sigi, sigj],
                                 alpha_max = [1,1.2])
        present, equi_all, surv = ap.community_assembly(species_id, pr = False)
        
        lv = species_id["level"][:,np.newaxis]
        richness[0,i,j] = np.mean(np.sum(surv & (lv == 0), axis = -1)[:,years//2:])
        richness[1,i,j] = np.mean(np.sum(surv & (lv == 1), axis = -1)[:,years//2:])
        
        similarity[i,j] = np.sum(surv[:,-200] & surv[:,-1])/np.sum(surv[:,-200] | surv[:,-1])
        
        print(i,j, richness[:,i,j])
"""

fig, ax = plt.subplots(3,1, sharex = True, sharey = True, figsize = (9,9))

cmap = ax[0].imshow(similarity, extent = [sigs[0], sigs[-1], sigs[0], sigs[-1]], origin = "lower")
fig.colorbar(cmap, ax = ax[0])
ax[0].axvline(omega/1.87, color = "k")
ax[0].plot(sigs, sigs, 'r')

cmap = ax[1].imshow(richness[0], extent = [sigs[0], sigs[-1], sigs[0], sigs[-1]]
                    , origin = "lower", vmin = 0,
                    vmax = np.nanpercentile(richness[0], 95))
fig.colorbar(cmap, ax = ax[1])

cmap = ax[2].imshow(richness[1], extent = [sigs[0], sigs[-1], sigs[0], sigs[-1]]
                    , origin = "lower", vmin = 0,
                    vmax = np.nanpercentile(richness[1], 95))
fig.colorbar(cmap, ax = ax[2])

# add layout
ax[-1].set_xlabel("Niche width Predator")
for i, a in enumerate(ax):
    a.set_title("ABCDEFGHI"[i], loc = "left")
    a.set_ylabel("Niche width basal species")
    
ax[0].set_title("Jaccard similarity")
ax[1].set_title("Prey richness")
ax[2].set_title("Predator richness")

fig.savefig("Figure_ap_stability_vs_niche_width.pdf")