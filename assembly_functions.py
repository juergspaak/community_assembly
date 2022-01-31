import numpy as np
from itertools import combinations

sig_res = 2
tot_res = 10

def LV_model(t, N, A, mu):
    return N*(mu-A.dot(N))

def generate_species(n_coms = 200, years = 800, locs = np.zeros(3), sig_locs = np.ones(3),
                     sigs = [1, .3, .8], utils = [1.5, 3, 0.08],
                     ms = [.5, .2, .25], level = [0.5, 0.5, 0]):
    
    level = np.array(level)/np.sum(level)
    com = (n_coms, years)
    species_id = {"loc": np.random.normal(0,1, com),
              "sig": np.array([np.full(com, sig) for sig in sigs]),
              "util": np.array([np.full(com, util) for util in utils]),
              "m": np.array([np.full(com, m) for m in ms]),
              "level": np.random.choice(np.arange(3), p = level, size = com)}
    
    # potentially change location and variation of trophic levels
    for i in range(len(level)):
        species_id["loc"][species_id["level"] == i] *= sig_locs[i]
        species_id["loc"][species_id["level"] == i] += locs[i]
    
    species_id["level"][:,:2] = 0
    
    return species_id

###############################################################################
# one-dimensional resource axis

def compute_LV_param(species_id, i = 0, pres = [0,1],
                     sig_res = sig_res, tot_res = tot_res):
    # interaction matrix
    A = np.zeros((np.sum(pres), np.sum(pres)))
    
    # get data
    level = species_id["level"][i,pres]
    loc = species_id["loc"][i, pres]
    sig = species_id["sig"][level, i, pres]
    util = species_id["util"][level, i, pres]
    id_prey = np.arange(len(loc))[species_id["level"][i,pres] == 0]
    id_pred = np.arange(len(loc))[species_id["level"][i,pres] == 1]
    id_mut = np.arange(len(loc))[species_id["level"][i,pres] == 2]
    
    sig2 = sig[:,np.newaxis]**2 + sig**2
    # species interaction if all species were prey
    A = (util[:,np.newaxis]*util
            /np.sqrt(2*np.pi*sig2)
            *np.exp(-(loc[:,np.newaxis] - loc)**2/2/sig2))
    
    
    # effect of predator on prey
    A[id_prey[:,np.newaxis], id_pred] = (util/np.sqrt(2*np.pi*sig**2)
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2))[:,id_pred]
    # assuming uniform competition kernel
    """A = (util[:,np.newaxis]*util
            *np.amax([sig[:,np.newaxis] + sig - np.abs(loc[:,np.newaxis] - loc)
                      , np.zeros(A.shape)], axis = 0))
    A[id_prey[:,np.newaxis], id_pred] = np.where((sig<np.abs(loc[id_prey, np.newaxis]-loc)).T, util[id_prey], 0).T[:,id_pred]
    """             
    
    # effect of prey on predator
    A[id_pred[:,np.newaxis], id_prey] = -.2*A[id_prey[:,np.newaxis], id_pred].T 
    
    # predators don't interact with each other
    A[id_pred[:,np.newaxis], id_pred] = 0
    
    # effect of mutualists on plants
    A[id_prey[:,np.newaxis], id_mut] = -(util/np.sqrt(2*np.pi*sig**2)
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2))[:,id_mut]
    
    # effect of mutualists on plants
    A[id_mut[:,np.newaxis], id_prey] = 2*A[id_prey[:,np.newaxis], id_mut].T
    
    # mutualists don't interact with each other
    A[id_mut[:, np.newaxis], id_mut] = 0.2
    A[id_mut, id_mut] = 2
    
        
    # compute the intrinsic growth rates
    sig2_res = (sig**2 + sig_res**2)
    mu = tot_res*util/np.sqrt(2*np.pi*sig2_res)*np.exp(-loc**2/2/sig2_res) - species_id["m"][0, i, pres]
    mu[id_pred] = -species_id["m"][1, i, pres][id_pred]
    mu[id_mut] = -species_id["m"][2, i, pres][id_mut]
    
    return mu, np.round(A,8)

def community_equilibrium(mu, A):
    # compute maximum subcommunity that coexists
    
    n_spec = len(mu) # number of species
    
    # new species can invade compute new equilibrium
    try:
            equi = np.linalg.solve(A, mu)
    except np.linalg.LinAlgError:
        equi = -np.ones(n_spec)
    if (equi>0).all():
        return equi # all species coexist
    inv_gr = np.zeros(n_spec)  
    
    ind = np.full(n_spec, True)
    sp_loc = np.arange(n_spec)
    
    # try to find a coexisting community, assuming most species will coexist
    for i in range(n_spec):
        
        # remove species with lowest density
        ind[np.argmin(equi)] = False
        equi[~ind] = 0
        sp = sp_loc[ind]
        
        # compute new equilibrium
        try:
            equi[ind] = np.linalg.solve(A[sp[:,np.newaxis], sp], mu[sp])
        except np.linalg.LinAlgError: # too many predators, remove one by chance
            ind[np.random.choice(np.where(np.diag(A) == 0)[0])] = False

        if (equi<0).any():
            continue # still not all coexist, compute next
        else: # all present species coexist, copute invasion
            inv_gr = mu - A.dot(equi)
            
            if (inv_gr<1e-10).all():
                return equi # no other species can invade
            else: # add species with highest invasion growth rate
                ind[np.argmax(inv_gr)] = True
        
    # community that coexists and stable not found in fast version
    # start systematic search
    
    
    for i in range(1, n_spec):
        n_pres = n_spec - i # euilibrium of the resident community
        
        # create all possible subcommunities with n_pres communities
        sub_coms = np.array(list(combinations(np.arange(n_spec), n_pres)))
        # randomize community
        sub_coms = sub_coms[np.argsort(np.random.rand(len(sub_coms)))]
        for sub_com in sub_coms:
            equi = np.zeros(n_spec)
            try:
                equi[sub_com] = np.linalg.solve(A[sub_com[:,np.newaxis], sub_com],
                                   mu[sub_com])
            except np.linalg.LinAlgError:
                continue # not stable equilibrium possible
            
            inv_gr = mu - A.dot(equi)
            
            if (equi>=0).all() and (inv_gr<1e-10).all():
                return equi
    print(mu)
    print(A)
    raise

def community_assembly(species_id, sig_res = sig_res, tot_res = tot_res):
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
            mu, A = compute_LV_param(species_id, i, present[i,j],
                                     sig_res, tot_res)
            
            # can the new species invade?
            if ((mu - np.sum(A[:,:-1]*equi, axis = -1))<1e-10).all():
                present[i,j+1, ind_spec] = True
                equi_all[i,j] = equi_all[i,j-1]
                continue
            
            try:
                equi = community_equilibrium(mu, A)
            except:
                break
            ind = equi>0        
            ind_spec = n_specs[present[i, j]][ind]
            equi = equi[ind]
            equi_all[i,j, ind_spec] = equi
            present[i, j+1, ind_spec] = True
        equi_all[i,-1, ind_spec] = equi
        print(i, np.sum((equi_all[i, -1]>0) & (species_id["level"][i] == 0)),
              np.sum((equi_all[i,-1]>0) & (species_id["level"][i] == 1)))
    
    return present, equi_all, equi_all>0

if __name__ == "__main__":
    species_id = generate_species(3,5000, sigs = [0.5, 1, .8],
                                  utils = [1.5, 1.5, 0.08])
    
    present, equi_all, surv = community_assembly(species_id)
    
    """np.savez("Data_LV_reference.npz", **species_id, equi_all = equi_all,
             surv = surv, present = present)"""
    import functions_for_plotting as fp
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
    fp.plot_richness(ax[0], surv, species_id)
    fp.plot_traits(ax[1,:2], surv, species_id)
    #fp.plot_traits(ax[2,:2], surv, species_id)
    #plot_invasion_prob(ax[-1,-1], surv, species_id)
    
    fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)
    fp.plot_mean_traits(ax[2,-1], surv, species_id, 1)
