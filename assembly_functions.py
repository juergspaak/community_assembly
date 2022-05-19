import numpy as np
from itertools import combinations

from various_competition_kernels import compute_LV_from_traits_kernels

def LV_model(t, N, A, mu):
    return N*(mu-A.dot(N))

# default traits
sigs = np.array([1, 0.8])
alpha_max = np.array([1.0,10])
omega = 2
mu_max = 1
ms = np.array([.1,.1])

n_invs = 201
species_id_invader_basal = {"level": np.zeros(n_invs),
                  "loc": np.linspace(-4,4,n_invs),
                  "sig": np.full(n_invs, sigs[0]),
                  "alpha_max": np.full(n_invs, alpha_max[0]),
                  "m": np.full(n_invs, ms[0])}

species_id_invader_pred = {"level": np.ones(n_invs),
                  "loc": np.linspace(-3,3,n_invs),
                  "sig": np.full(n_invs, sigs[1]),
                  "alpha_max": np.full(n_invs, alpha_max[1]),
                  "m": np.full(n_invs, ms[1])}

def generate_species(n_coms = 200, years = 800, locs = np.zeros(2),
                     sig_locs = None, sigs = sigs,
                     alpha_max = alpha_max, omega = omega, mu_max = mu_max,
                     ms = ms, level = [0.5, 0.5], kernel = False):
    
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
              # competition kernel
              "kernel": kernel}
    
    # potentially change location and variation of trophic levels
    for i in range(len(level)):
        species_id["loc"][species_id["level"] == i] *= sig_locs[i]
        species_id["loc"][species_id["level"] == i] += locs[i]
    
    # first species to arrive must be basal to avoid issues
    species_id["level"][:,:2] = 0
    
    # add environmental data
    species_id["omega"] = omega
    species_id["mu_max"] = mu_max
    
    return species_id

def compute_LV_from_traits(level, loc, sig, alpha_max, m, omega = omega,
                           mu_max = mu_max):
    
    id_prey = np.where(level == 0)[0]
    id_pred = np.where(level == 1)[0]
    
    # species interaction if all species were prey
    sig_basal = np.sqrt(sig[:,np.newaxis]**2 + sig**2)
    A = (alpha_max[:,np.newaxis]*alpha_max
            *np.exp(-(loc[:,np.newaxis] - loc)**2/2/sig_basal**2))
    
    # effect of predator on prey
    A[id_prey[:,np.newaxis], id_pred] = (alpha_max
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2))[:,id_pred]      
    
    # effect of prey on predator
    A[id_pred[:,np.newaxis], id_prey] = -0.1*A[id_prey[:,np.newaxis], id_pred].T 
    
    # predators don't interact with each other
    A[id_pred[:,np.newaxis], id_pred] = 0
    
    # compute the intrinsic growth rates
    mu = mu_max*np.exp(-loc**2/2/omega**2) - m
    
    # mortality rate of predators
    mu[id_pred] = -m[id_pred] 
    
    return mu, np.round(A,8)

def compute_LV_param(species_id, i = 0, pres = [0,1]):
    
    # get data
    level = species_id["level"][i,pres]
    loc = species_id["loc"][i, pres]
    sig = species_id["sig"][level, i, pres]
    alpha_max = species_id["alpha_max"][level, i, pres]
    m = species_id["m"][0, i, pres]
    
    if not species_id["kernel"]:
        return compute_LV_from_traits(level, loc, sig, alpha_max, m,
                                  omega = species_id["omega"],
                                  mu_max = species_id["mu_max"])
    
    else:
        return compute_LV_from_traits_kernels(level, loc, sig, alpha_max, m,
                                  omega = species_id["omega"],
                                  mu_max = species_id["mu_max"],
                                  kernel = species_id["kernel"])

def invasion_success(species_id, i, equi, species_id_invader):
    
    pres = np.where(equi != 0)[0]
    
    level = species_id["level"][i,pres]
    loc = np.append(species_id["loc"][i, pres],
                    species_id_invader["loc"])
    sig = np.append(species_id["sig"][level, i, pres],
                    species_id_invader["sig"])
    alpha_max = np.append(species_id["alpha_max"][level, i, pres],
                          species_id_invader["alpha_max"])
    m = np.append(species_id["m"][0, i, pres], species_id_invader["m"])
    level = np.append(species_id["level"][i,pres],
                      species_id_invader["level"])
    
    mu, A = compute_LV_from_traits(level, loc, sig, alpha_max, m,
                                   mu_max = species_id["mu_max"],
                                  omega = species_id["omega"])
    equi = np.append(equi[pres], np.zeros(len(species_id_invader["level"])))
    inv = mu - A.dot(equi)
    return inv, loc
    return inv[len(species_id["level"][i,pres]):], loc[len(species_id["level"][i,pres]):]

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

def invasion_scheme(mu, A):
    n = len(mu)
    counter = 0
    equi_all = np.zeros(( 2**n, n))
    r_i = np.full(equi_all.shape, np.nan)
    # compute all possible equilibria
    for n_spec in range(1, n+1):
        for comb in np.array(list(combinations(range(n), n_spec))):
            counter += 1
            try:
                equi = np.linalg.solve(A[comb[:, np.newaxis], comb], mu[comb])
            except np.linalg.LinAlgError:
                equi_all[counter] = np.nan
                continue # not a stable equilibrium
            if (equi<0).any():
                equi_all[counter] = np.nan
                continue # not a stable equilibrium
            equi_all[counter, comb] = equi
    equi_all[0] = 0
    equi_all = equi_all[np.all(np.isfinite(equi_all), axis = 1)]
    r_i = mu - np.nansum(A*equi_all[:,np.newaxis], axis = -1)
    
    r_i = np.round(r_i, 10)
    r_i = r_i[np.isfinite(r_i[:,0])] # remove communities that are not present
    return r_i
    

if False and __name__ == "__main__":
    species_id = generate_species(2, 200, omega = 2)
    present, equi_all, surv = community_assembly(species_id, pr = False)

    import functions_for_plotting as fp
    import matplotlib.pyplot as plt
    
    fig, ax = plt.subplots(2,2, sharex = True, sharey = "row", figsize = (10,10))
    fp.plot_richness(ax[0], surv, species_id)
    fp.plot_traits(ax[1], surv, species_id)