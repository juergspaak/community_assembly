import numpy as np
import matplotlib.pyplot as plt

import functions_for_plotting as fp

path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
file = path + "Data_LV_{}.npz"

names = ["only_prey", "pred_and_prey", "uniform_pred_prey", "pred_free"]
fig, ax = plt.subplots(2,2, sharex = True, sharey = True)
ax[0,0].set_ylabel("Frequency")
ax[1,0].set_ylabel("Frequency")
ax[1,0].set_xlabel("Life time")
ax[1,1].set_xlabel("Life time")

ax = ax.flatten()
cutoff = 100

for i, name in enumerate(names):
    data = np.load(file.format(name), allow_pickle = True)
    trait_condition = (data["loc"] > -1) & (data["loc"] <1)
    species_id = {"loc": data["loc"], "level": data["level"]}
    surv = data["surv"]
    
    fp.plot_persistence_time(ax[i], surv, species_id,
                             cutoff)
    ax[i].set_title(name)
    



