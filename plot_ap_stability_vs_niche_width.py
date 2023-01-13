import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp
import assembly_functions as ap
from timeit import default_timer as timer
"""
start = timer()
itera = 20
sigs, dsig = np.linspace(1, 3, itera, retstep = True)

richness = np.empty((2,itera,itera))
similarity = np.empty((itera, itera))

years = 2000
omega = 2.5

n = 8
prey = np.arange(-n, n+1, 1)
pred = np.arange(-n,n, 1) +0.5
level = len(prey)*[0] + len(pred)*[1]

for i, sigi in enumerate(sigs):
    for j, sigj in enumerate(sigs):
        species_id = ap.generate_species(1, years, level = [1,1],
                                 omega = omega, sigs = [sigi, sigj])
        # position first species at potentiall optimal locations
        species_id["loc"][:,:len(level)] = np.append(2*sigi*prey, 2*sigj*pred)
        species_id["level"][:,:len(level)] = level
        
        present, equi_all, surv = ap.community_assembly(species_id, pr = False)
        
        lv = species_id["level"][:,np.newaxis]
        richness[0,i,j] = np.mean(np.sum(surv & (lv == 0), axis = -1)[:,years//2:])
        richness[1,i,j] = np.mean(np.sum(surv & (lv == 1), axis = -1)[:,years//2:])
        
        similarity[i,j] = np.sum(surv[:,-200] & surv[:,-1])/np.sum(surv[:,-200] | surv[:,-1])
        
        print(i,j, richness[:,i,j], timer()-start)
#"""

years = 2000
fig = plt.figure(figsize = (10,11))
fs = 16
cmap = "bwr"
invasions = np.empty((len(ap.species_id_invader_basal["level"]), years))
extent =  [0, years,ap.species_id_invader_basal["loc"][0], 
                                   ap.species_id_invader_basal["loc"][-1]]
sigs_examples = [[1.6,1.5], [2.7, 2.5], [1.4,1.5], [1.3,2.5]]
sigs_examples = np.array([sigs[[7,6]], sigs[[-3,5]],
                          sigs[[5,6]], sigs[[3,-3]]]) + dsig/2
ax_id = [1,2,5,6]
# first show example cases
for j in range(4):
    cax = fig.add_subplot(3,2,ax_id[j])
    species_id = ap.generate_species(2, years, level = [1,1],
                                 omega = omega, sigs = sigs_examples[j])
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    
    fp.plot_traits([cax], surv, species_id)
    cax.set_title("ABCD"[j], loc = "left", fontsize = fs)

    cax.set_xlim([0,years])
    cax.set_xlabel("Year", fontsize = fs)
    cax.set_ylabel("Trait", fontsize = fs)#"""
cax = fig.add_axes([0.3,0.4,0.4,0.25])
cmap_ax = fig.add_axes([0.7,0.4,0.05,0.25])
extent = [sigs[0]-dsig/2,
          sigs[-1]+dsig/2,
          sigs[0]-dsig/2,
          sigs[-1]+dsig/2]
cmap = cax.imshow(similarity, extent = extent,
                  origin = "lower")
fig.colorbar(cmap, cax = cmap_ax)
cax.plot(sigs, sigs, 'r')
cmap_ax.set_ylabel("Jaccard similarity", fontsize = fs)

# annotate locations
cax.annotate("", sigs_examples[0][::-1]+dsig/2, [1/3,2/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "k", shrink = 0))
cax.annotate("", sigs_examples[1][::-1]+dsig/2, [2/3,2/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "k", shrink = 0))
cax.annotate("", sigs_examples[2][::-1]+dsig/2, [1/3,1/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "k", shrink = 0))
cax.annotate("", sigs_examples[3][::-1]+dsig/2, [2/3,1/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "k", shrink = 0))


cax.set_xlabel("Predator niche width\n$\sigma_P$", fontsize = fs)
cax.set_ylabel("Prey niche width\n$\sigma_B$", fontsize = fs)
cax.set_title("E", loc = "left", fontsize = fs)
fig.tight_layout()

fig.savefig("Figure_ap_stability_vs_niche_width.pdf")