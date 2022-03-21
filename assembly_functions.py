import numpy as np
from itertools import combinations

def LV_model(t, N, A, mu):
    return N*(mu-A.dot(N))

# default traits
sigs = np.array([1, 0.8])
alpha_max = np.array([1.0,1.0])
omega = 2
mu_max = 1
ms = np.array([.1,.1])

def_specs = {
    "alpha_max_red": 1,
    "sig_red": 1,
    "m_red": 1,
    "pred_red": 1,
    "mu_red": 1,
    "pred_niche": 1}

def generate_species(n_coms = 200, years = 800, locs = np.zeros(2),
                     sig_locs = None, sigs = sigs,
                     alpha_max = alpha_max, omega = omega, mu_max = mu_max,
                     ms = ms, level = [0.5, 0.5],
                     p_defended = 0, def_specs = def_specs):
    
    # choose random trait values such that species cover all trait space
    if sig_locs is None:
        # about 2% of basal species have negative intrinsic growth rate
        sig = omega*np.sqrt(2*np.log(mu_max/ms[0]))/2.5
        sig_locs = [sig, sig]
    
    level = np.array(level)/np.sum(level)
    com = (n_coms, years)
    species_id = {"loc": np.random.normal(0,1, com), # trait location
                  # consumption variation
              "sig": np.array([np.full(com, sig) for sig in sigs]),
              # maximum consumption
              "alpha_max": np.array([np.full(com, a) for a in alpha_max]),
              # morality rate
              "m": np.array([np.full(com, m) for m in ms]),
              # trophic level
              "level": np.random.choice(np.arange(2), p = level, size = com),
              # whether basal species are defended
              "defended": np.random.choice([False, True],
                                    p = [1 - p_defended, p_defended],
                                    size = com),
              "def_specs": def_specs}
    
    # potentially change location and variation of trophic levels
    for i in range(len(level)):
        species_id["loc"][species_id["level"] == i] *= sig_locs[i]
        species_id["loc"][species_id["level"] == i] += locs[i]
    
    # first species to arrive must be basal to avoid issues
    species_id["level"][:,:2] = 0
    
    # add environmental data
    species_id["omega"] = omega
    species_id["mu_max"] = mu_max
    
    # alter certain traits of defended species
    def_basal = species_id["defended"] & (species_id["level"] == 0)
    for trait in ["alpha_max", "sig", "m"]:
        species_id[trait][0,def_basal] *= def_specs[trait + "_red"]
    
    return species_id

def compute_LV_param(species_id, i = 0, pres = [0,1]):
    
    # interaction matrix
    A = np.zeros((np.sum(pres), np.sum(pres)))
    
    # get data
    level = species_id["level"][i,pres]
    loc = species_id["loc"][i, pres]
    sig = species_id["sig"][level, i, pres]
    alpha_max = species_id["alpha_max"][level, i, pres]
    id_prey = np.where(species_id["level"][i,pres] == 0)[0]
    id_pred = np.where(species_id["level"][i,pres] == 1)[0]
    
    defended = np.where(species_id["defended"][i, pres])[0]
    
    # species interaction if all species were prey
    sig_basal = np.sqrt(sig[:,np.newaxis]**2 + sig**2)
    A = (alpha_max[:,np.newaxis]*alpha_max
            *np.exp(-(loc[:,np.newaxis] - loc)**2/2/sig_basal**2))
    
    # effect of predator on prey
    # defended species have reduced visibility for predators
    red_niche_def = np.where(species_id["defended"][i, pres],
                             species_id["def_specs"]["pred_niche"], 1)
    A[id_prey[:,np.newaxis], id_pred] = (alpha_max
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2/red_niche_def))[:,id_pred]

    # some species might be protected against predation
    A[defended[:,np.newaxis], id_pred] *= species_id["def_specs"]["pred_red"]        
    
    # effect of prey on predator
    A[id_pred[:,np.newaxis], id_prey] = -A[id_prey[:,np.newaxis], id_pred].T 
    
    # predators don't interact with each other
    A[id_pred[:,np.newaxis], id_pred] = 0
    
    # compute the intrinsic growth rates
    mu = species_id["mu_max"]*np.exp(-loc**2/2/species_id["omega"]**2) - species_id["m"][0, i, pres]
    
    # defended species have reduced intrinsic growth rates
    mu[defended] *= species_id["def_specs"]["mu_red"]
    mu[id_pred] = -species_id["m"][1, i, pres][id_pred] 
    
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
    print("printed from community_equilibrium function")
    print("mu", mu)
    print("A", A)
    raise ValueError

def community_assembly(species_id, pr = True):
    n_coms, years = species_id["loc"].shape

    present = np.full((n_coms, years+1, years), False, dtype = bool)
    present[:,np.arange(years), np.arange(years)] = True
    n_specs = np.arange(years)
    
    equi_all = np.zeros(present.shape)


    for i in range(n_coms):
        ind_spec = np.arange(1)
        #equi_all[i,0,ind_spec] = mu[i,0]/A[i,0,0]
        equi = 0
        if pr:
            print(i)
        for j in range(years):
            # get species interactions
            mu, A = compute_LV_param(species_id, i, present[i,j])
            
            # can the new species invade?
            if ((mu - np.sum(A[:,:-1]*equi, axis = -1))<1e-10).all():
                present[i,j+1, ind_spec] = True
                equi_all[i,j] = equi_all[i,j-1]
                continue
            
            try:
                equi = community_equilibrium(mu, A)
            except ValueError:
                raise
            ind = equi>0        
            ind_spec = n_specs[present[i, j]][ind]
            equi = equi[ind]
            equi_all[i,j, ind_spec] = equi
            present[i, j+1, ind_spec] = True
        equi_all[i,-1, ind_spec] = equi
        if pr: 
            print(i, np.sum((equi_all[i, -1]>0) & (species_id["level"][i] == 0)),
              np.sum((equi_all[i,-1]>0) & (species_id["level"][i] == 1)))
    
    return present, equi_all, equi_all>0
    

if __name__ == "__main__":
    species_id = generate_species(2, 5000)
    present, equi_all, surv = community_assembly(species_id, pr = False)

    import functions_for_plotting as fp
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
    fp.plot_richness(ax[0], surv, species_id)
    fp.plot_traits(ax[1,:2], surv, species_id)
