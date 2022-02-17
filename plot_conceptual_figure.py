import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
import functions_for_plotting as fp


np.random.seed(0)

species_id = af.generate_species(2, level = [1,0,0],
                                 sigs = [1.9,0.8,0], years=1000)
species_id["loc"][0,:3] = [2, -0.5, -2.5]
present, equi_all, surv = af.community_assembly(species_id)

fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)

def plot_species(time, ax, com = 0, species_id = species_id, equi_all = equi_all):
    for i in range(species_id["loc"].shape[1]):
        if equi_all[com, time, i] == 0:
            continue
        lv = species_id["level"][com, i]
        ax.plot(res, equi_all[com, time, i]*species_id["max_util"][lv, com, i]*
                np.exp(-(res-species_id["loc"][com, i])**2
                            /2/species_id["sig"][lv, com, i]**2), color = "br"[lv])

res = np.linspace(-7, 7, 1001)

omega = af.sig_res
K = 2
sig_basal = 1
sig_pred = 0.5

fig, ax = plt.subplots(3,5, sharex = True, sharey = True, figsize = (14,9))


ax[0,0].plot(res, af.max_res*np.exp(-res**2/2/omega**2), color = "k")
ax[0,1].plot(res, af.max_res*np.exp(-res**2/2/omega**2), color = "k")
ax[0,2].plot(res, af.max_res*np.exp(-res**2/2/omega**2), color = "k")
ax[0,-1].plot(res, af.max_res*np.exp(-res**2/2/omega**2), color = "k")


plot_species(0, ax[0,0])
plot_species(1, ax[0,1])
plot_species(2, ax[0,2])
plot_species(800, ax[0,-1])


species_pred = af.generate_species(2, level = [0.5,0.5,0],
                                 sigs = [1.5,0.8,0])
species_pred["loc"][0,:5] = species_id["loc"][0][np.argsort(equi_all[0,-1])][-5:]
species_pred["level"][0,:5] = 0
#species_pred["level"][0,5] = 1
present, equi_pred, surv = af.community_assembly(species_pred)

for i, t in enumerate(np.arange(4, 14)):
    plot_species(t, ax[1:].flatten()[i], species_id = species_pred, equi_all = equi_pred)
    
for a in ax[-1]:
    a.set_xlabel("Trait axis")

for a in ax[:,0]:
    a.set_ylabel("Consumption rate")
    
fig.tight_layout()    
fig.savefig("Figure_conceptual.pdf")
