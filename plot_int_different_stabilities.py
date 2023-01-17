import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp
import assembly_functions as ap
from timeit import default_timer as timer

start = timer()
itera = 20
sigs, dsig = np.linspace(0.2, 3, itera, retstep = True)

richness = np.empty((2,itera,itera))
similarity = np.empty((itera, itera))

years = 500
omega = 4

n = 8
prey = np.arange(-n, n+1, 1)
pred = np.arange(-n,n, 1) +0.5
level = len(prey)*[0] + len(pred)*[1]
mu, A = np.empty((2,2,itera, itera), dtype = object)
variances = np.empty((2,2,itera,itera))

"""
for i, sigi in enumerate(sigs):
    for j, sigj in enumerate(sigs):
        species_id = ap.generate_species(1, years, level = [1,1],
                                 omega = sigj*2.5, sigs = [sigi, sigj])
        #species_id = ap.generate_species(1, years, level = [1,1],
        #                         omega = 5, sigs = [sigi, sigj])
        # position first species at potentiall optimal locations
        #species_id["loc"][:,:len(level)] = np.append(2*sigi*prey, 2*sigj*pred)
        #species_id["level"][:,:len(level)] = level
        
        present, equi_all, surv = ap.community_assembly(species_id, pr = False)
        
        mu[0,i,j], A[0,i,j] = ap.compute_LV_param(species_id, i = 0, pres = present[0,-1])
        mu[1,i,j], A[1,i,j] = ap.compute_LV_param(species_id, i = 0, pres = np.full(years, True, dtype = bool))
        
        #np.fill_diagonal(A[0,i,j], np.nan)
        #np.fill_diagonal(A[1,i,j], np.nan)
        
        for mat in range(2):
            A_copy = A[mat, i,j].copy()
            np.fill_diagonal(A_copy, np.nan)
            # variance of the sum
            variances[0, mat, i, j] = np.nanvar(np.nansum(A_copy, axis = 1))
            # overall variance
            variances[1, mat, i, j] = np.nanvar(A_copy)
    
        print(i,j, richness[:,i,j], timer()-start)

extent = [sigs[0]-dsig/2,
          sigs[-1]+dsig/2,
          sigs[0]-dsig/2,
          sigs[-1]+dsig/2]        

fig, ax = plt.subplots(2,2, sharex = True, sharey = True)
for i in range(2):
    for j in range(2):
        cmap = ax[i,j].imshow(variances[i,j], origin = "lower", extent = extent)
        fig.colorbar(cmap, ax = ax[i,j])

title = ["Local\nspecies", "Regional\nspecies"]
ylab = ["var(sum_j a_ij)", "var(a_ij)"]
for i in range(2):
    ax[i,0].set_ylabel(ylab[i] + "\nPrey Niche")
    ax[-1,i].set_xlabel("Predator niche")
    ax[0,i].set_title(title[i])
    
for i, a in enumerate(ax.flatten()):
    a.set_title("ABCD"[i], loc = "left")
        
fig.savefig("Figure_ap_stability_condition.pdf")
#"""
##############################################################################
# show spectrum of different matrices
fig, ax = plt.subplots(3,3, figsize = (9,9))

extent = [sigs[0]-dsig/2,
          sigs[-1]+dsig/2,
          sigs[0]-dsig/2,
          sigs[-1]+dsig/2]  

ij_comb = [[3, 17], [10,17], [16, 18],
           [3, 10], [0,0], [16, 10],
           [3, 2],  [10, 3], [16, 3]]
fig_loc = [[1/3, 2/3], [1/2, 2/3], [2/3, 2/3],
           [1/3, 1/2], [1/2, 1/2], [2/3, 1/2],
           [1/3, 1/3], [1/2, 1/3], [2/3, 1/3]]

#cmap = ax[1,1].imshow(variances[0,0], extent = extent, origin = "lower")
ax[1,1].plot(sigs[[0,-1]], sigs[[0,-1]], 'orange')
#fig.colorbar(cmap, ax = ax[1,1])


for k, a in enumerate(ax.flatten()):
    if k == 4: continue
    i,j = ij_comb[k]
    sigi, sigj = sigs[[i,j]]
    species_id = ap.generate_species(1, 2000, level = [1,1],
                                 omega = sigj*3, sigs = [sigi, sigj])

    # position first species at potentiall optimal locations
    species_id["loc"][:,:len(level)] = np.append(2*sigi*prey, 2*sigj*pred)
    species_id["level"][:,:len(level)] = level
    
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    mu[0,i,j], A[0,i,j] = ap.compute_LV_param(species_id, i = 0, pres = present[0,-1])
    
    # compute eigenvector decomposition
    eigvals = np.linalg.eigvals(-A[0, i,j])
    a.scatter(np.real(eigvals), np.imag(eigvals))
    a.set_title(str(len(A[0,i,j])) + " species")
    
    # annotate locations
    ax[1,1].annotate("", np.array([sigs[i], sigs[j]])+dsig/2, fig_loc[k]
                 , xycoords = "data", textcoords = "figure fraction",
                 arrowprops=dict(color = "r", shrink = 0))
    
    
fig.tight_layout()
"""
##############################################################################
# show spectrum of different matrices
fig, ax = plt.subplots(3,3, figsize = (9,9))

ij_comb = [[3, 17], [10,17], [16, 18],
           [3, 10], [0,0], [16, 10],
           [3, 2],  [10, 3], [16, 3]]
fig_loc = [[1/3, 2/3], [1/2, 2/3], [2/3, 2/3],
           [1/3, 1/2], [1/2, 1/2], [2/3, 1/2],
           [1/3, 1/3], [1/2, 1/3], [2/3, 1/3]]

cmap = ax[1,1].imshow(variances[0,0], extent = extent, origin = "lower")
ax[1,1].plot(sigs[[0,-1]], sigs[[0,-1]], 'orange')
fig.colorbar(cmap, ax = ax[1,1])


for k, a in enumerate(ax.flatten()):
    if k == 4: continue
    i,j = ij_comb[k]
    
    # compute eigenvector decomposition
    eigvals = np.linalg.eigvals(-A[1, i,j])
    a.scatter(np.real(eigvals), np.imag(eigvals))
    a.set_title(str(len(A[1,i,j])) + " species")
    
    # annotate locations
    ax[1,1].annotate("", np.array([sigs[i], sigs[j]])+dsig/2, fig_loc[k]
                 , xycoords = "data", textcoords = "figure fraction",
                 arrowprops=dict(color = "r", shrink = 0))
    
    
fig.tight_layout()"""