import numpy as np
import matplotlib.pyplot as plt

colors = ["b", "r", "g"]
labels = ["Phytoplankton", "Zooplankton", "Mutualists"]

def plot_traits(axc, surv, species_id):
    years = np.arange(surv.shape[-2])
    for i, a in enumerate(axc):
        for y in years[:-1]:
            # does this species survive for a certain time?
            if surv[i,:,y].any():
                a.plot(years[surv[i,:,y]],
                       np.repeat(species_id["loc"][i,y], np.sum(surv[i,:,y])),
                       colors[species_id["level"][i,y]])
                
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
    
def plot_invasion_prob(axc, surv, species_id):
    years = np.arange(surv.shape[-1])
    persistence_length = np.sum(surv, axis = 1)
    for level in range(len(colors)):
        if not (species_id["level"] == level).any():
                continue
        axc.plot(years, np.mean((persistence_length > 0) & (species_id["level"] == level), axis = 0),
                 colors[level])
        
def plot_persistence_time(axc, surv, species_id, cutoff = 50):
    persistence_length = np.sum(surv, axis = 1)[:,cutoff:]
    years = np.arange(surv.shape[-2]+1)
    for level in range(len(colors)):
        if not (species_id["level"] == level).any():
                continue
        hist, bins = np.histogram(persistence_length[species_id["level"][:,cutoff:] == level],
                                  bins = years[::5])
        hist[hist == 0] = -1
        axc.plot(bins[1:], hist, '.', color = colors[level])
    axc.semilogy()
    
    
if __name__ == "__main__":
    data = np.load("Data_LV_reference.npz", allow_pickle = True)
    species_id = {**data}
    surv = data["surv"]
    
    
    plot_persistence_time(plt, surv, species_id, 50)
    
    fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
    plot_richness(ax[0], surv, species_id)
    plot_traits(ax[1,:2], surv, species_id)
    plot_traits(ax[2,:2], surv, species_id)
    plot_invasion_prob(ax[-1,-1], surv, species_id)
    
    plot_mean_traits(ax[1,-1], surv, species_id, 0)
    plot_mean_traits(ax[2,-1], surv, species_id, 1)