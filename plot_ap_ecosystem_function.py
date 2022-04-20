import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap

np.random.seed(0)

omega = 5
time = 500
years = np.arange(time+1)
ylim = [-3,3]
s = 1
n_coms = 200

# basal species only
species_id = ap.generate_species(50, time, level = [1,1],
                                 omega = omega)
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

# prey and predator identities
id_pred = np.where(species_id["level"] == 1)
id_prey = np.where(species_id["level"] == 0)

# create predator and prey equilibira
equi_pred = equi_all.copy()
equi_pred[id_prey[0], :, id_prey[1]] = np.nan

equi_prey = equi_all.copy()
equi_prey[id_pred[0], :, id_pred[1]] = np.nan
equi = np.array([equi_prey, equi_pred])

EF = np.nansum(equi, axis = -1)



##############################################################################
# add a random walk with similar properties

change_EF = (EF[...,1:] - EF[...,:-1]).reshape(2,-1)

change_EF_random_walk = np.array([np.random.choice(change_EF[0], size = (len(years), n_coms)),
                                  np.random.choice(change_EF[1], size = (len(years), n_coms))])

EF_random_walk = np.nanmean(EF[...,100:], axis = (-1,-2)) + np.cumsum(change_EF_random_walk[:,100:], axis = 1).T

##############################################################################
# plot results
fig, ax = plt.subplots(1,2,sharex = True, sharey = False, figsize = (9,9))
titles = ["Basal species", "Predator species"]
fs = 16
q = [5, 95]
for i in range(2):
    ax[i].plot(years, EF[i, 0], label = "Simulation", color = "blue")
    ax[i].fill_between(years, np.nanpercentile(EF[i], q[0], axis = 0),
                     np.nanpercentile(EF[i], q[1], axis = 0),
                     color = "blue", alpha = 0.5)
    ax[i].plot(years[100:], EF_random_walk[0,:,i], label = "Random walk", color = "orange")
    ax[i].fill_between(years[100:], np.nanpercentile(EF_random_walk[...,i], q[0], axis = 0),
                     np.nanpercentile(EF_random_walk[...,i], q[1], axis = 0),
                     color = "orange", alpha = 0.5)

    ax[i].set_xlabel("Year", fontsize = fs)
    ax[i].set_title("AB"[i], loc = "left", fontsize = fs)
    ax[i].set_title(titles[i], fontsize = fs + 2)
    

ax[0].set_ylabel("Total biomass", fontsize = fs)
ax[0].legend(fontsize = fs)

ax[0].set_xlim([100, time])

fig.savefig("Figure_ap_EF_over_time.pdf")


