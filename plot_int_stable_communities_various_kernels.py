import numpy as np
import matplotlib.pyplot as plt

from timeit import default_timer as timer

import assembly_functions as ap
import functions_for_plotting as fp
import various_competition_kernels as vck

"""
itera = 10
keys = list(vck.kernels.keys())
years = 1000

richness = np.full((len(keys), 2, itera, itera), np.nan)
jaccard = np.full((len(keys), itera, itera), np.nan)

sigs = np.linspace(0.5, 2, itera)
species_id = ap.generate_species(1, years, omega = 2, level = [1,1])
lv = species_id["level"][:,np.newaxis]
start = timer()
for i, sigi in enumerate(sigs):
    for j, sigj in enumerate(sigs):
        
        species_id["sig"][0] = sigi
        species_id["sig"][1] = sigj
        
        for id_key, key in enumerate(keys):
            
            species_id["kernel"] = key
            try:
                present, equi_all, surv = ap.community_assembly(species_id, pr = False)
            except ValueError:
                continue
            
            
            richness[id_key, 0,i,j] = np.mean(np.sum(surv & (lv == 0),
                                                     axis = -1)[:,years//2:])
            richness[id_key, 1,i,j] = np.mean(np.sum(surv & (lv == 1),
                                                     axis = -1)[:,years//2:])
            
            jaccard[id_key, i,j] = (np.sum(surv[:,-200] & surv[:,-1])
                                    /np.sum(surv[:,-200] | surv[:,-1]))
        print(i,j, np.round([sigi, sigj],2), np.round(jaccard[:,i,j], 2),
              timer()-start)"""
        
fig, ax = plt.subplots(2,3,sharex = True, sharey = True, figsize = (9,9))

ax = ax.flatten()
extent = [sigs[0], sigs[-1], sigs[0], sigs[-1]]

for i, key in enumerate(keys):
    ax[i].imshow(jaccard[i], vmin = 0,vmax = 10, origin = "lower",
                 extent = extent)
    ax[i].set_title(key)