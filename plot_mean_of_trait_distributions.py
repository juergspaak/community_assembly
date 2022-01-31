import numpy as np
import matplotlib.pyplot as plt

path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
data = np.load(path  + "data_LV_assembly_predators.npz")
#data = {**data}



cutoff = 99 # data before is burn in time
years = np.arange(data["present"].shape[1])[cutoff:]

fig, ax = plt.subplots(3,4, figsize = (9,9))

# remove traits of species that are not present
traits_base = ((data["level"] == 1)*data["loc"])[:,np.newaxis]*data["present"]
traits_base[traits_base == 0] = np.nan
traits_base = traits_base[:,cutoff:]

combs = [[fun, level] for fun in ["mean", "std"] for level in range(2)]

for i, comb in enumerate(combs):
    fun = np.nanmean if comb[0] == "mean" else np.nanstd
    
    traits_base = ((data["level"] == comb[1])*data["loc"])[:,np.newaxis]*data["present"]
    traits_base[traits_base == 0] = np.nan
    traits_base = traits_base[:,cutoff:]

    a = fun(traits_base[0], axis = -1)
    a[a == 0] = np.nan # in case of variance over only one species
    site_mean = np.nanmean(a)
    ax[0,i].plot(years, a)
    reg_mean = fun(data["loc"])
    ax[0,i].axhline(reg_mean, color = "r", label = "regional") # regional 
    ax[0,i].axhline(site_mean, color = "b", label = "site",
                    linestyle = "--") # regional
    
    ax[1,i].hist(a)
    ax[1,i].axvline(reg_mean, color = "r", label = "regional") # regional 
    ax[1,i].axvline(site_mean, color = "b", label = "site",
                    linestyle = "--") # regional
    
    ax[2,i].hist(fun(traits_base[:,-1], axis = -1))
    ax[2,i].axvline(reg_mean, color = "r", label = "regional") # regional 
    ax[2,i].axvline(np.nanmean(fun(traits_base[:,-1], axis = -1))
                    , color = "b", label = "site", linestyle = "--") # regional
    
    ax[0,i].set_title(str(comb))
    
fig.tight_layout()