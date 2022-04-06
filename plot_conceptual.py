import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap
import functions_for_plotting as fp

np.random.seed(0)
omega = 2
time = 250
ylim = [-3,3]
s = 1

# basal species only
species_id = ap.generate_species(2, time, level = [1,0],
                                 omega = omega)
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

fig, ax = plt.subplots(2,1, sharex = True, sharey = True, figsize = (9,7))
fp.plot_traits(ax[[0]], surv, species_id)

invasions = np.empty((len(ap.species_id_invader_basal["level"]), time))
for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[:,i] = inv[-ap.n_invs:]
    
# change to binary
invasions = np.where(invasions>0, 1.0, np.nan)

extent =  [0, time,ap.species_id_invader_basal["loc"][0], 
                                   ap.species_id_invader_basal["loc"][-1]]
cmap = "bwr"
ax[0].plot(np.nan, np.nan, 'bs', label = "Possible invaison", alpha = 0.5)
ax[0].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 1, vmax = 1.5, extent = extent,
             aspect = "auto", cmap = cmap)
ax[0].legend()

species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega)
species_id["level"][:,:50] = 0
present, equi_all, surv = ap.community_assembly(species_id, pr = False)
fp.plot_traits(ax[[1]], surv, species_id)

for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[:, i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan)
ax[1].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 1, vmax = 1.5, extent = extent,
             aspect = "auto", cmap = cmap)


for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_pred)
    invasions[:,i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan) 
ax[1].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 0.5, vmax = 1, extent = extent,
             aspect = "auto", cmap = cmap)

"""ax[1].scatter(np.repeat(np.arange(time)[:,np.newaxis], ap.n_invs, axis = 1)[invasions>0],
        np.repeat(ap.species_id_invader_pred["loc"][np.newaxis], time, axis = 0)[invasions>0]
            , color = 'r', s = s, alpha = 0.2, marker = "s")"""

# layout

ax[1].plot(np.nan, np.nan, "r", label = "Predator")
ax[1].plot(np.nan, np.nan, "b", label = "Prey")
ax[1].legend(loc = "upper left")

ax[1].set_xlim([0, time])
ax[1].set_ylim(ylim)

ax[1].set_xlabel("Time")
ax[0].set_ylabel("Trait")
ax[1].set_ylabel("Trait")

ax[0].set_title("A", loc = "left")
ax[0].set_title("Prey only")

ax[1].set_title("B", loc = "left")
ax[1].set_title("Predator and Prey")

fig.savefig("Figure_conceptual.pdf")