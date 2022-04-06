import numpy as np
import matplotlib.pyplot as plt

from matplotlib.cm import Set1

import assembly_functions as ap

colors = ["b", "r", "g"]
labels = ["Phytoplankton", "Zooplankton", "Mutualists"]

def plot_traits(axc, equi, species_id, plot_invasion = False):
    
    surv = equi>0
    years = np.arange(surv.shape[-2])
    for i, a in enumerate(axc):
        for y in years[:-1]:
            # does this species survive for a certain time?
            value = species_id["level"][i,y] + 2*(1-species_id["level"][i,y])*species_id["defended"][i,y]
            if np.sum(surv[i,:,y])==1:
                a.scatter(years[surv[i,:,y]],
                       np.repeat(species_id["loc"][i,y], np.sum(surv[i,:,y]))
                       , color = colors[value], s = 2*1.5**2)
            elif np.sum(surv[i,:,y])>1:
                a.plot(years[surv[i,:,y]],
                       np.repeat(species_id["loc"][i,y], np.sum(surv[i,:,y])),
                       color = colors[value], linewidth = 3)
    if plot_invasion:
        plot_invasion_prob(axc, equi, species_id)
                
def plot_mean_traits(axc, surv, species_id, level):
    years = np.arange(surv.shape[-2])
    time = np.ones(surv.shape)*years[:,np.newaxis]
    loc = np.ones(surv.shape)*species_id["loc"][:,np.newaxis]
    sp_pres = surv & (species_id["level"][:,np.newaxis] == level)
    
    traits, binx, biny = np.histogram2d(loc[sp_pres],
                                    time[sp_pres],
                                    bins = [50, len(years)-1])
    #traits = traits/np.sum(traits, axis = 0, keepdims = True)
    #traits[traits == 0] = np.nan
    axc.imshow(traits, aspect = "auto", origin = "lower",
                    extent = [biny[0], biny[-1], binx[0], binx[-1]], vmin = 1)
                
def plot_richness(axc, surv, species_id):
    years = np.arange(surv.shape[-2])
    for i, a in enumerate(axc[:-1]):
        for level in range(len(colors)):
            if not (species_id["level"] == level).any():
                continue
            
            a.plot(years, np.sum((species_id["level"][i] == level) & surv[i], axis = -1),
                                 colors[level])
    for level in range(len(colors)):
        if not (species_id["level"] == level).any():
                continue
        axc[-1].plot(years, np.mean(np.sum((species_id["level"][:,np.newaxis] == level)
                                           & surv, axis = -1 ), axis = 0), colors[level],
                     label = labels[level])
    axc[-1].legend()
    
def plot_invasion_prob(axc, equi, species_id):
    itera = 101
    time = np.arange(equi.shape[1])
    inv_basal = {"loc": np.linspace(*np.nanpercentile(species_id["loc"], [1, 99]), itera),
                 "level": np.zeros(itera),
                 "sig": np.full(itera, species_id["sig"][0,0,0]),
                 "m": np.full(itera, species_id["m"][0,0,0]),
                 "alpha_max": np.full(itera, species_id["alpha_max"][0,0,0]),
                 "defended": np.full(itera, False, dtype = bool)}
    inv_pred = {"loc": np.linspace(*np.nanpercentile(species_id["loc"], [1, 99]), itera),
                 "level": np.ones(itera),
                 "sig": np.full(itera, species_id["sig"][1,0,0]),
                 "m": np.full(itera, species_id["m"][1,0,0]),
                 "alpha_max": np.full(itera, species_id["alpha_max"][1,0,0]),
                 "defended": np.full(itera, False, dtype = bool)}
    invasions = np.empty((2, itera, len(time)))
    for run, ax in enumerate(axc):
        for i in time:
            inv, loc = ap.invasion_success(species_id, 0, equi[0, i],
                                           inv_basal)
            invasions[0, :,i] = inv[-itera:] # initial entries are from species
            inv, loc = ap.invasion_success(species_id, 0, equi[0, i],
                                           inv_pred)
            invasions[1, :,i] = inv[-itera:] # initial entries are from species
        
        # change to binary
        invasions = np.where(invasions>0, 1.0, np.nan)
        
        extent =  [0, time[-1],ap.species_id_invader_basal["loc"][0], 
                                           ap.species_id_invader_basal["loc"][-1]]
        cmap = "bwr"
        ax.imshow(invasions[0, ::-1], interpolation = "none", alpha = 0.5,
                     vmin = 1, vmax = 1.5, extent = extent,
                     aspect = "auto", cmap = cmap)
        ax.imshow(invasions[1, ::-1], interpolation = "none", alpha = 0.5,
                     vmin = 0.5, vmax = 1, extent = extent,
                     aspect = "auto", cmap = cmap) 
        
        
            
def plot_persistence_time(axc, surv, species_id, cutoff = 50
                          , trait_condition = True):
    persistence_length = np.sum(surv, axis = 1)[:,cutoff:]
    years = np.arange(surv.shape[-2]+1)
    for level in range(len(colors)):
        if not (species_id["level"] == level).any():
                continue
        hist, bins = np.histogram(persistence_length[trait_condition & 
                                    species_id["level"][:,cutoff:] == level],
                                  bins = years[::5])
        hist[hist == 0] = -1
        axc.plot(bins[1:], hist, '.', color = colors[level])
    axc.semilogy()
    
biotime_colors = Set1(np.linspace(0,1,9))
keys = ["Birds","Invertebrates",
                        "Terrestrial\nplants", "Fish", "Other"]
biotime_colors = {key: biotime_colors[i] for i, key in enumerate(keys)}
biotime_colors["Mammals"] = biotime_colors["Other"]
biotime_colors["Benthos"] = biotime_colors["Other"]
biotime_colors["All"] = biotime_colors["Other"]
biotime_colors["Amphibians"] = biotime_colors["Other"]
biotime_colors["Marine invertebrates"] = biotime_colors["Invertebrates"]
biotime_colors["Terrestrial invertebrates"] = biotime_colors["Invertebrates"]
biotime_colors["Freshwater invertebrates"] = biotime_colors["Invertebrates"]
biotime_colors["Terrestrial plants"] = biotime_colors["Terrestrial\nplants"]

if __name__ == "__main__":
    species_id = ap.generate_species(2, 5000, omega = 4, sigs = [3, 2])
    species_id["level"][:,:20] = 0
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    plot_traits([plt.gca()], equi_all, species_id, plot_invasion=True)
