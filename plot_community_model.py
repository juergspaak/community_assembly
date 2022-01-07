import numpy as np
import matplotlib.pyplot as plt

sig_res = 4
tot_res = 10
predator_presence = "_predators"

n_coms = 50
years = 500
com = (n_coms, years)

pred = .5 if predator_presence else 0

species_id = {"loc": np.random.normal(0,1, com),
              "sig": np.array([np.full(com, .5),
                               np.full(com, .25),
                               np.full(com, .05)]),
              "util": np.array([np.full(com, 1.5),
                                np.full(com, 3),
                                np.full(com, .08)]),
              "m": np.array([np.full(com, .5),
                             np.full(com, .2),
                             np.full(com, .25)]),
              "level": np.random.binomial(1,pred,com)}

species_id["level"][:,:6] = [0,0,0,1,1,1]

def compute_LV_param(species_id, i = 0, pres = [True, True]):
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

mu, A = compute_LV_param(species_id, pres = np.arange(3))
equi = np.linalg.solve(A, mu)

res = np.linspace(-5,5,1001)

fig, ax = plt.subplots(2,1,sharex = True, figsize = (7,7))
ax[0].plot(res, tot_res*np.exp(-res**2/sig_res), "k", label = "Resource level")

for i in range(3):
    ax[0].plot(res, species_id["util"][0,0,i]*np.exp(-(res - species_id["loc"][0,i])**2/species_id["sig"][0,0,i]),
             label = "Prey {}".format(i))
    ax[1].bar(species_id["loc"][0,i], 5*equi[i], width = .1)
    



color = ["red", "blue", "purple"]
for i in range(4,7):
    ax[1].plot(res, species_id["util"][1,0,i]*np.exp(-(res - species_id["loc"][0,i])**2/species_id["sig"][1,0,i]),
             label = "Predator {}".format(i-4),
             color = color[i-4])

ax[0].legend()
ax[1].legend()
ax[0].set_xlabel("Resource/\nTrait location")
ax[0].set_ylabel("Resource availability\nresource consumption")
ax[1].set_ylabel("Prey densities\n(in absence of predators)\nPredator attack rates")
ax[1].set_xlabel("Prey trait location\nPredator trait location")

fig.savefig("Figure_community_model.pdf")


