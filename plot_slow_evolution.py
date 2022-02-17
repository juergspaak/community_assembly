import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
import functions_for_plotting as fp

species_id = af.generate_species(n_coms = 2, years = 10000, locs = np.zeros(3),
                     sigs = [1, 2, .8], level = [0.5, 0.5, 0], max_utils = [1.5, 3, 0.08])

n_coms, years = species_id["loc"].shape

present = np.full((n_coms, years+1, years), False, dtype = bool)
present[:,np.arange(years), np.arange(years)] = True
n_specs = np.arange(years)

equi_all = np.zeros(present.shape)


for i in range(n_coms):
    ind_spec = np.arange(1)
    #equi_all[i,0,ind_spec] = mu[i,0]/A[i,0,0]
    equi = 0
    print(i)
    for j in range(years):
        # get species interactions
        mu, A = af.compute_LV_param(species_id, i, present[i,j],
                                 af.sig_res, af.tot_res)
        
        # can the new species invade?
        if ((mu - np.sum(A[:,:-1]*equi, axis = -1))<1e-10).all():
            present[i,j+1, ind_spec] = True
            equi_all[i,j] = equi_all[i,j-1]
            continue
        
        try:
            equi = af.community_equilibrium(mu, A)
        except:
            break
        ind = equi>0        
        ind_spec = n_specs[present[i, j]][ind]
        equi = equi[ind]
        equi_all[i,j, ind_spec] = equi
        present[i, j+1, ind_spec] = True
        
        # change identity of of next species
        # if all present species are basal species, introduce a random predator
        if j == species_id["level"].shape[-1]-1:
            continue
        if (species_id["level"][i, ind_spec] == 0).all():
            species_id["level"][i, j+1] = 1
        else:
            #continue
            # otherwise, a random species creates a copy
            id_copy = np.random.choice(ind_spec)
            species_id["level"][i, j+1] = species_id["level"][i, id_copy]
            species_id["loc"][i, j+1] = species_id["level"][i, id_copy] + np.random.normal(0, 0.1)
        
    equi_all[i,-1, ind_spec] = equi
    print(i, np.sum((equi_all[i, -1]>0) & (species_id["level"][i] == 0)),
          np.sum((equi_all[i,-1]>0) & (species_id["level"][i] == 1)))
    
surv = equi_all>0
fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)
#fp.plot_traits(ax[2,:2], surv, species_id)
#plot_invasion_prob(ax[-1,-1], surv, species_id)

#fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)
#fp.plot_mean_traits(ax[2,-1], surv, species_id, 1)