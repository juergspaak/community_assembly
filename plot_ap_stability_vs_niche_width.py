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
        
        print(i,j, richness[:,i,j])
"""

fig = plt.figure(figsize = (10,11))
fs = 16
time = 2000
cmap = "bwr"
invasions = np.empty((len(ap.species_id_invader_basal["level"]), time))
extent =  [0, time,ap.species_id_invader_basal["loc"][0], 
                                   ap.species_id_invader_basal["loc"][-1]]
sigs_examples = [[2.5,1.5], [2.8,2.7], [0.6,1.5], [0.6,2.5]]
ax_id = [1,2,5,6]
# first show example cases
for j in range(4):
    cax = fig.add_subplot(3,2,ax_id[j])
    species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega, sigs = sigs_examples[j])
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    
    fp.plot_traits([cax], surv, species_id)
    cax.set_title("ABCD"[j], loc = "left", fontsize = fs)

    cax.set_xlim([0,time])
    cax.set_xlabel("Year", fontsize = fs)
    cax.set_ylabel("Trait", fontsize = fs)
cax = fig.add_axes([0.3,0.4,0.4,0.25])
cmap_ax = fig.add_axes([0.7,0.4,0.05,0.25])
cmap = cax.imshow(similarity, extent = [sigs[0], sigs[-1], sigs[0], sigs[-1]],
                  origin = "lower")
fig.colorbar(cmap, cax = cmap_ax)
cax.plot(sigs, sigs, 'r')
cax.axvline(omega/2, color = "k")
cmap_ax.set_ylabel("Jaccard similarity", fontsize = fs)

# annotate locations
cax.annotate("", sigs_examples[0][::-1], [1/3,2/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "red"))
cax.annotate("", sigs_examples[1][::-1], [2/3,2/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "red"))
cax.annotate("", sigs_examples[2][::-1], [1/3,1/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "red"))
cax.annotate("", sigs_examples[3][::-1], [2/3,1/3]
             , xycoords = "data", textcoords = "figure fraction",
             arrowprops=dict(color = "red"))

cax.set_xlabel("Predator niche width\n$\sigma_P$", fontsize = fs)
cax.set_ylabel("Prey niche width\n$\sigma_B$", fontsize = fs)
cax.set_title("E", loc = "left", fontsize = fs)
fig.tight_layout()

fig.savefig("Figure_ap_stability_vs_niche_width.pdf")