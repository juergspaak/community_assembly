import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp
import assembly_functions as ap
from timeit import default_timer as timer

start = timer()
itera = 20
sigs, dsig = np.linspace(1, 3, itera, retstep = True)

richness = np.empty((2,itera,itera))
similarity = np.empty((itera, itera))

years = 1000
omega = 4

n = 8
n_coms = 5
variances = np.empty((n_coms,itera,itera))
A = np.empty(n_coms, dtype = "object")

for i, sigi in enumerate(sigs):
    for j, sigj in enumerate(sigs):
        species_id = ap.generate_species(n_coms, years, level = [1,1],
                                 omega = sigj*2.5, sigs = [sigi, sigj])
        
        for k in range(n_coms):
            mu, A[k] = ap.compute_LV_param(species_id, i = 0,
                                    pres = np.full(years,True, dtype = bool))

            np.fill_diagonal(A[k], np.nan)
            #variances[k,i,j] = np.nanvar(np.nansum(A[k,i,j], axis = -1))
            variances[k,i,j] = np.nanvar(np.nansum(A[k], axis = -1))

    
        print(i,j, timer()-start)

extent = [sigs[0]-dsig/2,
          sigs[-1]+dsig/2,
          sigs[0]-dsig/2,
          sigs[-1]+dsig/2]        

fig = plt.figure()
plt.imshow(np.nanmean(variances, axis = 0), origin = "lower", extent = extent)
plt.colorbar()

plt.ylabel("Prey Niche")
plt.xlabel("Predator niche")
        
fig.savefig("Figure_stability_condition.pdf")