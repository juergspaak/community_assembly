import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp
import assembly_functions as ap

fig = plt.figure(figsize = (12, 15))
np.random.seed(0)
fs = 24

n_specs = np.arange(10,100, 10)
n_coms = 1
years = 1000
x_lim = [0, years]

ax = []
for i, n_spec in enumerate(n_specs):
    
    # genrate species
    species_id = ap.generate_species(n_coms, years)
    # select a random species identity for each species
    sp_id = np.random.randint(n_spec, size = (n_coms, years))
    # choose location and level
    locs = species_id["loc"][np.arange(n_coms),:n_spec]
    species_id["loc"] = locs[np.arange(n_coms)[:,np.newaxis], sp_id]
    species_id["level"][:] = sp_id%2
    
    # run assembly
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    
    # plot assembly
    ax.append(fig.add_subplot(6,3,i+1))
    fp.plot_traits([ax[-1]], surv, species_id)
    

n_coms = 100
jaccard = np.empty((len(n_specs), n_coms))
# how many lead to permanent community?
for i, n_spec in enumerate(n_specs):
    print(i)
    years = n_spec*10 # each species can invade 10 times on average
    # genrate species
    species_id = ap.generate_species(n_coms, years)
    # select a random species identity for each species
    sp_id = np.random.randint(n_spec, size = (n_coms, years))
    # choose location and level
    locs = species_id["loc"][np.arange(n_coms),:n_spec]
    species_id["loc"] = locs[np.arange(n_coms)[:,np.newaxis], sp_id]
    species_id["level"][:] = sp_id%2
    
    # run assembly
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    
    # are communities stable
    jaccard[i] = (np.sum(surv[:,-50] & surv[:,-1], axis = -1)
                  /np.sum(surv[:,-50] | surv[:,-1], axis = -1))
    
ax_jaccard = fig.add_subplot(2,1,2)

ax_jaccard.plot(n_specs, np.mean(jaccard, axis = -1), color = "purple"
                , label = "Jaccard similarity")
ax_jaccard.plot(n_specs, np.percentile(jaccard, [5,95], axis = -1).T, '--',
                color = "purple")

#ax_stable = ax_jaccard.twinx()
ax_jaccard.plot(n_specs, np.mean(jaccard >0.95, axis = -1),'o',
               color = "green", label = "Stable communities")

ax_jaccard.set_ylabel("Jaccard similarity", fontsize = fs)
ax_jaccard.set_title("K", loc = "left", fontsize = fs)
ax_jaccard.set_ylim([0,1])
ax_jaccard.legend(fontsize = fs)
#ax_stable.set_ylabel("Proportion of stable communities")
#ax_stable.set_ylim([0,1])

ax_jaccard.set_xlabel("Species richness\nProportion of stable communities"
                      , fontsize = fs)

y_lim = [-3,3]
for i, a in enumerate(ax):
    a.set_title("ABCDEFGHIJ"[i], loc = "left", fontsize = fs)
    a.set_title("Regional\nspecies richness {}".format(n_specs[i]), fontsize = 20)
    a.set_ylim(y_lim)
    a.set_xlim(x_lim)
    if i%3 != 0:
        a.set_yticks([])
    else:
        a.set_ylabel("Traits", fontsize = fs)
    if i>=6:
        a.set_xlabel("Year", fontsize = fs)
    else:
        a.set_xticks([])
        
fig.tight_layout()

fig.savefig("Figure_ap_finite_species_richness.pdf")