import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.stats import linregress
from matplotlib.cm import viridis

def LV_model(t, N, A, mu):
    return N*(mu-A.dot(N))

n_coms = 500
years = 200
com = (n_coms, years)

col = viridis(np.linspace(0,1, years))
ls = ["-", ":", "--"]
###############################################################################
# one-dimensional resource axis
# create species

sig_res = 4
tot_res = 10
H_mut = 2 # halfsaturation constant for mutualists
predator_presence = "_predators"

pred = .5 if predator_presence else 0
p = np.array([3,1,3])
p = p/np.sum(p)
species_id = {"loc": np.random.normal(0,1, com),
              "sig": np.array([np.full(com, 1),
                               np.full(com, .2),
                               np.full(com, .2)]),
              "util": np.array([np.full(com, 1.5),
                                np.full(com, 2),
                                np.full(com, 2)]),
              "d": np.array([np.full(com, .5),
                             np.full(com, .2),
                             np.full(com, .25)]),
              "level": np.random.choice([0,1,2], replace = True,
                                        p = p, size = com)}

species_id["level"][:,0] = 0 # first 2 species are plants species
#species_id["level"][:,2] = 2

#species_id["loc"][0,:3] = [-1,1,1]


def compute_LV_param(species_id, i = 0, pres = [True, True]):
    # returns the community parameters for the model
    
    # get data
    level = species_id["level"][i,pres]
    loc = species_id["loc"][i, pres]
    sig = species_id["sig"][level, i, pres]
    util = species_id["util"][level, i, pres]
    d = species_id["d"][level, i, pres]
    id_prey = np.arange(len(loc))[species_id["level"][i,pres] == 0]
    id_pred = np.arange(len(loc))[species_id["level"][i,pres] == 1]
    id_mut = np.arange(len(loc))[species_id["level"][i,pres] == 2]
    
    sig2 = sig[:,np.newaxis]**2 + sig**2
    # linear species interactions between of same guild
    A = (util[:,np.newaxis]*util
            /np.sqrt(2*np.pi*sig2)
            *np.exp(-(loc[:,np.newaxis] - loc)**2/2/sig2))
    
    # predators don't interact with each other
    A[id_pred[:,np.newaxis], id_pred] = 0
    

    # generally, species from different trophic levels don't interact
    A[level != level[:,np.newaxis]] = 0
    
    # effect of predator on prey
    A[id_prey[:,np.newaxis], id_pred] = (util/np.sqrt(2*np.pi*sig**2)
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2))[:,id_pred]
    
    # effect of prey on predator
    A[id_pred[:,np.newaxis], id_prey] = -.2*A[id_prey[:,np.newaxis], id_pred].T 
    
    
    # saturating effect of mutualism
    B = np.zeros(A.shape)
    
    # effect of mutualists on plants
    B[id_prey[:,np.newaxis], id_mut] = (util/np.sqrt(2*np.pi*sig**2)
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2))[:,id_mut]
    
    # effect of mutualists on plants
    B[id_mut[:,np.newaxis], id_prey] = B[id_prey[:,np.newaxis], id_mut].T
        
    # compute the intrinsic growth rates
    sig2_res = (sig**2 + sig_res**2)
    mu = tot_res*util/np.sqrt(2*np.pi*sig2_res)*np.exp(-loc**2/2/sig2_res)
    # 
    mu[id_pred] = 0
    mu[id_mut] = np.sum(util/np.sqrt(2*np.pi*sig**2)
                    *np.exp(-(loc[id_prey, np.newaxis]-loc)**2/2/sig**2)
                    , axis = 0)[id_mut]
    #mu[id_mut] = 0
    mu -= d
    
    return mu, np.round(A,8), np.round(B, 8) 

def model(t, N, mu, A, B):
    pollen = np.sum(B*N/(.01+np.sum(B*N, axis = 1)), axis = 0)
    return N*(mu - A.dot(N) + pollen/(1+pollen))
    return N*(mu - A.dot(N) + B.dot(N)/(H_mut + B.dot(N)))

present = np.full((n_coms, years+1, years), False, dtype = bool)
present[:,np.arange(years), np.arange(years)] = True
n_specs = np.arange(years)

equi_all = np.zeros(present.shape)

mu, A, B = compute_LV_param(species_id, 0, np.arange(3))
print(species_id["level"][0,:3])
print(species_id["loc"][0,:3])
print(mu)
print(A)
print(B)

N_start = .1*np.ones(3)

sol = solve_ivp(model, [0,20], N_start, args = (mu, A, B))

#plt.plot(sol.t, sol.y.T)

time = np.array([0,20])
plot_until = 5
for i in range(n_coms):
    if i<plot_until:
        fig, ax = plt.subplots(2,1)
    ind_spec = np.arange(1)
    #equi_all[i,0,ind_spec] = mu[i,0]/A[i,0,0]
    equi = []
    print(i)
    for j in range(years):
        # get species interactions
        mu, A, B = compute_LV_param(species_id, i, present[i,j])
        
        equi = np.append(equi,1e-3)
        # can the new species invade?
        """if (model(0,equi, mu, A, B)<1e-10).all():
            present[i,j+1, ind_spec] = True
            equi_all[i,j] = equi_all[i,j-1]
            equi = equi[:-1]
            continue"""
        
        sol = solve_ivp(model, time + j*time[-1], equi, args = (mu, A, B))
        
        
        ind = sol.y[:,-1]>1e-3        
        ind_spec = n_specs[present[i, j]][ind]
        if i<plot_until:
            ax[0].axvline(j*time[-1])
            for id_k, k in enumerate(ind_spec):
                ax[0].plot(sol.t, sol.y[id_k], color = col[k],
                           ls = ls[species_id["level"][i,k]])
                ax[1].plot(sol.t[[0,-1]], 2*[species_id["loc"][i,k]]
                           , color = col[k], ls = ls[species_id["level"][i,k]])
        equi = sol.y[:,-1][ind]
        equi_all[i,j, ind_spec] = equi
        present[i, j+1, ind_spec] = True
    equi_all[i,-1, ind_spec] = equi
    print(i, [np.sum((equi_all[i, -1]>0) & (species_id["level"][i] == l)) for l in [0,1,2]])
    #raise
    if i < plot_until:
        plt.show()

richness = np.mean(np.sum(present, axis = -1), axis = 0)
richness[1:-1] -= 1#"""

np.savez("data_LV_assembly_mut_pred.npz", **species_id, present = present,
         equi_all = equi_all, richness = richness)

import test_plot_all
plt.show()
