import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap
import functions_for_plotting as fp

np.random.seed(0)
omega = 5
time = 500
ylim = [-3,3]
s = 1

# basal species only
species_id = ap.generate_species(2, time, level = [1,1],
                                 omega = omega)
present, equi_all, surv = ap.community_assembly(species_id, pr = False)



equi_pred = equi_all.copy()
id_pred = np.where(species_id["level"] == 1)
id_prey = np.where(species_id["level"] == 0)
equi_pred[id_prey[0], :, id_prey[1]] = np.nan
equi_prey = equi_all.copy()

equi_prey[id_pred[0], :, id_pred[1]] = np.nan

fig = plt.figure()

plt.plot(np.nansum(equi_prey, axis = -1).T)
plt.plot(np.nansum(equi_pred, axis = -1).T, '--')
plt.xlabel("Time")
plt.ylabel("Total biomass")

fig.savefig("Figure_ap_EF_over_time.pdf")

fig = plt.figure()
dens = equi_pred[np.isfinite(np.log(equi_pred))]
plt.hist(np.log(dens), bins = 30)
