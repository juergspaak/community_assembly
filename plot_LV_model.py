import matplotlib.pyplot as plt
import numpy as np

path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
data = np.load(path + "data_LV_assembly.npz")
data = {**data}

years = np.arange(data["present"].shape[1])

present = data["present"]
surv = data["equi_all"]>0
persistence_length = np.sum(present, axis = 1)


res = np.nanpercentile(data["loc"], [0,100])
sig_ref = np.nanmean(data["sig"][0])
util_ref = np.nanmean(data["util"][0])
m_ref = np.nanmean(data["m"][0])


def invasion_prob(equi, loc, util, sig):
    
    # generate interaction matrix
    loc_inv = np.linspace(*res, 1000)
    
    sig_comb = (sig_ref**2 + sig**2)
    A_invader = (util_ref*util/np.sqrt(2*np.pi*sig_comb)*
          np.exp(-(loc_inv[:,np.newaxis]-loc)**2/2/sig_comb))
    
    sig_res = (data["sig_res"]**2 + sig_ref**2)
    mu_inv = data["tot_res"]*util_ref/np.sqrt(2*np.pi*sig_res)*np.exp(-loc_inv**2/2/sig_res) - m_ref
    invasion = mu_inv - np.sum(A_invader*equi, axis = 1)
    
    return loc_inv, invasion



fig, ax = plt.subplots(3,3, sharex = "col", sharey = "row", figsize = (10,10))
invasion_all = np.empty((len(years), 1000))
for i in range(2):
    ax[0,i].plot(years, np.sum(present[i], axis = -1), 'b.-')
    
    
    for j in years[:-1]:
        if persistence_length[i,j]>1:
            ax[1,i].plot(years[present[i,:,j]],
                         np.repeat(data["loc"][i,j], persistence_length[i,j]),
                         color = "blue")
            
    # add invasion probability
    for j in years:
        loc_inv, invasion = invasion_prob(data["equi_all"][i,j][data["present"][i,j]],
                                 data["loc"][i][data["present"][i,j]],
                                 data["util"][0,i][data["present"][i,j]],
                                 data["sig"][0,i][data["present"][i,j]])
        invasion_all[j] = invasion
    locs, times = np.meshgrid(loc_inv, years)
    ax[1,i].scatter(times[invasion_all>0], locs[invasion_all>0], s = 1, color = "grey")
        
    ax[2,i].plot(years, np.mean((invasion_all>0), axis = 1), 'b')
            
ax[0,-1].plot(years, np.mean(np.sum(present, axis = -1), axis = 0), 'b')


time = np.ones(surv.shape)*years[:,np.newaxis]
loc = np.ones(surv.shape)*data["loc"][:,np.newaxis]

traits, binx, biny = np.histogram2d(loc[surv], time[surv], bins = [50, len(years)-1])
#traits = traits/np.sum(traits, axis = 0, keepdims = True)
#traits[traits == 0] = np.nan
ax[1,-1].imshow(traits, aspect = "auto", origin = "lower",
                extent = [biny[0], biny[-1], binx[0], binx[-1]], vmin = 1)

ax[2,-1].plot(years[1:], np.mean(persistence_length>1, axis = 0), 'b')
ax[0,0].set_xlim([0, 100])
ax[0,1].set_xlim([0, 100])

ax[-1,0].set_ylim([0,1])

# add layout
ax[-1,0].set_xlabel("Time")
ax[-1,0].set_xlabel("Time")
ax[-1,0].set_xlabel("Time")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Trait distribution")
ax[2,0].set_ylabel("Invasion probability")

ax[0,0].set_title("Example 1")
ax[0,1].set_title("Example 2")
ax[0,2].set_title("Average case")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHI"[i], loc = "left")
    
fig.savefig("Figure_LV_no_predators.pdf")
fig.savefig("Figure_LV_no_predators.png")

#"""
