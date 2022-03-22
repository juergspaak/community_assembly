import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap
import functions_for_plotting as fp

np.random.seed(0)
omega = 2
time = 200
ylim = [-3.5,3.5]

# basal species only
species_id = ap.generate_species(2, time, level = [1,0],
                                 omega = omega)
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

fig, ax = plt.subplots(2,1, sharex = True, sharey = True, figsize = (9,7))
fp.plot_traits(ax[[0]], surv, species_id)

invasions = np.empty((time, len(ap.species_id_invader_basal["level"])))
for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[i] = inv[-ap.n_invs:]

ax[0].scatter(np.repeat(np.arange(time)[:,np.newaxis], ap.n_invs, axis = 1)[invasions>0],
        np.repeat(ap.species_id_invader_basal["loc"][np.newaxis], time, axis = 0)[invasions>0]
            , color = 'b', s = 4, alpha = 0.2, label = "Possible invasion")
ax[0].legend()

species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega)
species_id["level"][:,:50] = 0
present, equi_all, surv = ap.community_assembly(species_id, pr = False)
fp.plot_traits(ax[[1]], surv, species_id)

invasions = np.empty((time, len(ap.species_id_invader_basal["level"])))
for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[i] = inv[-ap.n_invs:]

ax[1].scatter(np.repeat(np.arange(time)[:,np.newaxis], ap.n_invs, axis = 1)[invasions>0],
        np.repeat(ap.species_id_invader_basal["loc"][np.newaxis], time, axis = 0)[invasions>0]
            , color = 'b', s = 4, alpha = 0.2)


for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_pred)
    invasions[i] = inv[-ap.n_invs:]

ax[1].scatter(np.repeat(np.arange(time)[:,np.newaxis], ap.n_invs, axis = 1)[invasions>0],
        np.repeat(ap.species_id_invader_pred["loc"][np.newaxis], time, axis = 0)[invasions>0]
            , color = 'r', s = 4, alpha = 0.2)


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
ax[1].set_title("Preadtor and Prey")

fig.savefig("Figure_conceptual.pdf")