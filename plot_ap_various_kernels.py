import numpy as np
import matplotlib.pyplot as plt

import assembly_functions as ap
import functions_for_plotting as fp
import various_competition_kernels as vck

np.random.seed(0)

fig, ax = plt.subplots(2,3, sharex = True, sharey = True,
                       figsize = (9,9))
fs = 16
ax[0,0].set_ylabel("Trait", fontsize = fs)
ax[1,0].set_ylabel("Trait", fontsize = fs)
ax[1,0].set_xlabel("Year", fontsize = fs)
ax[1,1].set_xlabel("Year", fontsize = fs)
ax[1,2].set_xlabel("Year", fontsize = fs)

ax = ax.flatten()

itera = 101
x = np.linspace(-2,2,itera)
keys = list(vck.kernels.keys())

species_ref = ap.generate_species(1, itera, omega = 2, level = [1,0])
species_ref["loc"][0] = x

years  = 200
species_id = ap.generate_species(1, years, level = [1,1])

for i, key in enumerate(keys):
    ax_inset = ax[i].inset_axes([0.05,0.65,0.4,0.3])
    u_il, res_loc = vck.kernel_functions(np.zeros(1), 2*np.ones(1), x, key)
    ax_inset.plot(x, u_il[:,0])
    ax_inset.set_xlim(x[[0,-1]])
    ax_inset.set_ylim([0, 1.2])
    
    species_ref["kernel"] = key
    mu, A = ap.compute_LV_param(species_ref, 0, np.full(itera, True))
    ax_inset.plot(x, A[:,itera//2])
    ax_inset.set_yticks([])
    ax_inset.set_xticks([])
    
    species_id["kernel"] = key
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    fp.plot_traits(ax[[i]], surv, species_id)
    
    ax[i].set_title("ABCDEF"[i], loc = "left", fontsize = fs)
    ax[i].set_title(key, fontsize = fs)
    
ax[0].set_xlim([0,years])

fig.savefig("Figure_ap_various_kernels.pdf")
    