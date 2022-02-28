import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as af
import functions_for_plotting as fp
from scipy.stats import gaussian_kde

n_year = 1000
species_id = af.generate_species(n_coms = 2, sigs = [0.5, 0.25],
                                  years = n_year)
present, equi_all, surv = af.community_assembly(species_id)
# cut of initial 500 years
surv = surv[:,500:]

time = np.arange(surv.shape[1])
com = 0
fig, ax = plt.subplots(2,3, figsize = (14,10), sharex = "col",
                       sharey = "col")

for i in range(n_year):
    if not np.any(surv[com,:,i]):
        continue
    lv = species_id["level"][com,i]
    ax[lv,0].plot(time[surv[com, :, i]],
                 np.full(np.sum(surv[com, :, i]), species_id["loc"][com, i]),
                         "br"[lv])

ax[0,0].set_ylabel("Species traits")
ax[1,0].set_ylabel("Species traits")

##############################################################################
# trait distribution
x = 1.2*np.linspace(*ax[0,0].get_ylim(), 1000)
for j, bw in enumerate([0.5, 0.2]):
    for t in np.arange(200,500,10):
    
        lv = species_id["level"][com, surv[com, t]]
        traits = species_id["loc"][com, surv[com, t]]
        kernel = gaussian_kde(traits[lv == 0], bw_method = bw)
        ax[0,j+1].plot(x, kernel.evaluate(x), 'b', alpha = 0.1)
        ax[0,j+1].plot([-bw, bw], [0.01, 0.01], 'b')
        ax[0,j+1].text(0, 0.01, "$\sigma$", ha = "center", va = "bottom")
        kernel_pred = gaussian_kde(traits[lv == 1], bw_method = bw)
        ax[1,j+1].plot(x, kernel_pred.evaluate(x), 'r', alpha = 0.1)
        ax[1,j+1].plot([-bw, bw], [0.01, 0.01], 'r')
        ax[1,j+1].text(0, 0.01, "$\sigma$", ha = "center", va = "bottom")
    
ax[0,1].set_title("Large kernel (Family focus)")
ax[0,2].set_title("Medium kernel")

ax[0,1].set_ylim([0,None])
ax[0,2].set_ylim([0,None])

for a in ax[:,1:].flatten():
    a.set_ylabel("Frequency")
for a in ax[1,1:].flatten():
    a.set_xlabel("Species traits")

fig.tight_layout()
fig.savefig("Figure_trait_distribution.pdf")


