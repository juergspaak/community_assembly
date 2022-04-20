import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

import assembly_functions as ap
import functions_for_plotting as fp

np.random.seed(2)

fig, ax = plt.subplots(2,1, sharex = False, sharey = False, figsize = (9,7))
fs = 16

n_year_shown = 20 # number of years shown with densities over time

omega = 2
years = 250
colors = plt.cm.viridis(np.random.rand(years))
ylim = [-3,3]
s = 1
n_coms = 2

# basal species only
species_id = ap.generate_species(n_coms, years, level = [1,1],
                                 omega = omega)
present, equi_all, surv = ap.community_assembly(species_id, pr = False)

time_sim = np.linspace(0,1,30)

n_coms, years = species_id["loc"].shape

present = np.full((n_coms, years+1, years), False, dtype = bool)
present[:,np.arange(years), np.arange(years)] = True
n_specs = np.arange(years)

equi_all = np.zeros(present.shape)
equi = []
n_spec = np.arange(years)

# only do it for one community
i = 0
ind_spec = np.arange(1)
for j in range(years):
    #print(j)
    mu, A = ap.compute_LV_param(species_id, i, present[i,j])
    
    # can the new species invade?
    if ((mu - np.sum(A[:,:-1]*equi, axis = -1))<1e-10).all():
        present[i,j+1, ind_spec] = True
        equi_all[i,j] = equi_all[i,j-1]
        continue
    
    equi = np.append(equi, 1e-2)
    
    # solve densities over time
    sol = solve_ivp(lambda t,N: 100*N*(mu-A.dot(N)), time_sim[[0,-1]],
                    equi, t_eval = time_sim)
    if j>years-n_year_shown:
        print(j)
        for n_pres in range(len(sol.y)):
            lv = species_id["level"][i,present[i,j]][n_pres]
            ax[0].plot(sol.t + j*time_sim[-1], sol.y[n_pres]
                      , color = colors[present[i,j]][n_pres],
                      linestyle = '-' if lv else '--')
    #plt.plot(sol.t + j*time_sim[-1], sol.y.T)
    ind = sol.y[:,-1]/np.sum(sol.y[:,-1])>1e-2
    ind_spec = n_specs[present[i, j]][ind]
    equi = sol.y[ind,-1]
    equi_all[i,j, ind_spec] = equi
    present[i, j+1, ind_spec] = True
    if sum(ind) == 0:
        raise
    

ax[0].semilogy()
ax[0].set_ylim([1e-3,2])


fp.plot_traits(ax[[1]], surv, species_id)

invasions = np.empty((len(ap.species_id_invader_basal["level"]), years))
for j in range(years):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, j],
                                   ap.species_id_invader_basal)
    invasions[:, j] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan)
cmap = "bwr"
extent =  [0, years,ap.species_id_invader_basal["loc"][0], 
                                   ap.species_id_invader_basal["loc"][-1]]
ax[1].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 1, vmax = 1.5, extent = extent,
             aspect = "auto", cmap = cmap)


for j in range(years):
    inv, loc = ap.invasion_success(species_id, 0, equi_all[0, j],
                                   ap.species_id_invader_pred)
    invasions[:,j] = inv[-ap.n_invs:]

invasions = np.where(invasions>0, 1.0, np.nan) 
ax[1].imshow(invasions[::-1], interpolation = "none", alpha = 0.5,
             vmin = 0.5, vmax = 1, extent = extent,
             aspect = "auto", cmap = cmap)

#### add layout

ax[0].set_xlim([years-n_year_shown, years])
ax[0].set_xlabel("Year", fontsize = fs)
ax[1].set_xlabel("Year", fontsize = fs)
ax[1].set_ylabel("Species trait", fontsize = fs)
ax[0].set_ylabel("Species densities", fontsize = fs)

ax[0].set_title("A", loc = "left", fontsize = fs)
ax[1].set_title("B", loc = "left", fontsize = fs)
fig.tight_layout()

fig.savefig("Figure_ap_continuous_times.pdf")