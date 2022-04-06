import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap
import functions_for_plotting as fp

np.random.seed(0)
omega = 2
time = 2000
ylim = [-3,3]

fig, ax = plt.subplots(3,1, sharex = True, sharey = True, figsize = (9,7))
extent =  [0, time,ap.species_id_invader_basal["loc"][0], 
                                   ap.species_id_invader_basal["loc"][-1]]
invasions = np.empty((len(ap.species_id_invader_basal["level"]), time))
cmap = "bwr"

###############################################################################
# continuous invasion
species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega, sigs = [1,0.8])
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

fp.plot_traits(ax[[0]], surv, species_id)

for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[:, i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan)
ax[0].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 1, vmax = 1.5, extent = extent,
             aspect = "auto", cmap = cmap)

for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_pred)
    invasions[:,i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan) 
ax[0].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 0.5, vmax = 1, extent = extent,
             aspect = "auto", cmap = cmap)

###############################################################################
# intermediate case
species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega, sigs = [1,1.2])
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

###############################################################################
# stable case
species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega, sigs = [1,1.4])
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

fp.plot_traits(ax[[2]], surv, species_id)

for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_basal)
    invasions[:, i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan)
ax[2].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 1, vmax = 1.5, extent = extent,
             aspect = "auto", cmap = cmap)

for i in range(time):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, i],
                                   ap.species_id_invader_pred)
    invasions[:,i] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan) 
ax[2].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 0.5, vmax = 1, extent = extent,
             aspect = "auto", cmap = cmap)

# layout

ax[0].plot(np.nan, np.nan, "r", label = "Predator")
ax[0].plot(np.nan, np.nan, "b", label = "Prey")
ax[0].plot(np.nan, np.nan, 'bs', label = "Possible invaison", alpha = 0.5)
ax[0].legend(loc = "upper left")

ax[1].set_xlim([0, time])
ax[1].set_ylim(ylim)

ax[-1].set_xlabel("Time")
ax[0].set_ylabel("Trait")
ax[1].set_ylabel("Trait")
ax[2].set_ylabel("Trait")

ax[0].set_title("A", loc = "left")
ax[0].set_title("Continuous invasions, small predator niche")

ax[1].set_title("B", loc = "left")
ax[1].set_title("Eventual stability, intermediate predator niche")

ax[2].set_title("C", loc = "left")
ax[2].set_title("Stable community, large predator niche")

fig.savefig("Figure_ap_stable_community.pdf")