import matplotlib.pyplot as plt
import numpy as np


"""
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
data = np.load(path  + "data_LV_assembly_predators.npz")
data = {**data}

years = np.arange(data["present"].shape[1])


surv = data["equi_all"]>0
persistence_length = np.sum(data["present"], axis = 1)"""
"""
path = ""
data = np.load(path  + "data_LV_assembly_predators.npz")
data = {**data}"""

years = np.arange(data["present"].shape[1])


surv = data["equi_all"]>0
persistence_length = np.sum(data["present"], axis = 1)

res = np.nanpercentile(data["loc"], [0,100])
sig_ref = np.nanmean(data["sig"], axis = (1,2))
util_ref = np.nanmean(data["util"], axis = (1,2))
m_ref = np.nanmean(data["m"], axis = (1,2))


def invasion_prob_prey(equi, loc, util, sig, level):
    
    # generate interaction matrix
    loc_inv = np.linspace(*res, 1000)
    #loc_inv = np.append(loc_inv, loc[level == 0])
    
    sig_comb = (sig_ref[0]**2 + sig[0]**2)
    A_invader = (util_ref[0]*util[0]/np.sqrt(2*np.pi*sig_comb)*
          np.exp(-(loc_inv[:,np.newaxis]-loc)**2/2/sig_comb))
    
    # change interaction for predator species
    A_invader[:,level == 1] = (util[1]/np.sqrt(2*np.pi*sig[1]**2)
                    *np.exp(-(loc_inv[:,np.newaxis]-loc)**2/2/sig[1]**2))[:,level == 1]
    
    sig_res = (data["sig_res"]**2 + sig_ref[0]**2)
    mu_inv = data["tot_res"]*util_ref[0]/np.sqrt(2*np.pi*sig_res)*np.exp(-loc_inv**2/2/sig_res) - m_ref[0]
    invasion = mu_inv - np.sum(A_invader*equi, axis = 1)
    
    return loc_inv, invasion

def invasion_prob_pred(equi, loc, util, sig, level):
    
    # generate interaction matrix
    loc_inv = np.linspace(*res, 1000)
    
    # change interaction for predator species
    A_invader = -.2*(util_ref[1]/np.sqrt(2*np.pi*sig_ref[1]**2)
                    *np.exp(-(loc_inv[:,np.newaxis]-loc)**2/2/sig_ref[1]**2))
    A_invader[:,level == 1] = 0
    
    mu_inv = -m_ref[1]
    invasion = mu_inv - np.sum(A_invader*equi, axis = 1)
    
    return loc_inv, invasion

i = 0
j = 2
j = 80
level = data["level"][i][surv[i,j]]
equi = data["equi_all"][i,j][surv[i,j]]
loc = data["loc"][i,surv[i,j]]
util = data["util"][:, i, surv[i,j]]
sig = data["sig"][:, i, surv[i,j]]

print(level)

loc_inv, invasion = invasion_prob_prey(equi, loc, util, sig, level)
plt.plot(loc_inv, invasion)

loc_inv, invasion = invasion_prob_pred(equi, loc, util, sig, level)
plt.plot(loc_inv, invasion)

plt.plot(loc[level == 0], np.zeros(sum(level == 0)), 'o')
plt.plot(loc[level == 1], np.zeros(sum(level == 1)), 'o')

fig, ax = plt.subplots(4,3, sharex = "col", sharey = "row", figsize = (10,10))
invasion_all = np.empty((len(years), 1000))
for i in range(2):
    ax[0,i].plot(years, np.sum((data["level"][i] == 0)
                               & (surv[i]), axis = -1), 'b.-')
    ax[0,i].plot(years, np.sum((data["level"][i] == 1)
                               & (surv[i]), axis = -1), 'r.-')
    
    
    for j in years[:-1]:
        if surv[i,:,j].any():
            ax[1,i].plot(years[surv[i,:,j]],
                         np.repeat(data["loc"][i,j], np.sum(surv[i,:,j])),
                         color = "r" if data["level"][i,j] else "b")
            ax[2,i].plot(years[surv[i,:,j]],
                         np.repeat(data["loc"][i,j], np.sum(surv[i,:,j])),
                         color = "r" if data["level"][i,j] else "b")
            
    # add invasion probability
    for j in years:
        loc_inv, invasion = invasion_prob_prey(data["equi_all"][i,j][surv[i,j]],
                                 data["loc"][i][surv[i,j]],
                                 data["util"][:,i,surv[i,j]],
                                 data["sig"][:,i,surv[i,j]],
                                 data["level"][i,surv[i,j]])
        
        invasion_all[j] = invasion
    dl = loc_inv[1]-loc_inv[0]
    locs, times = np.meshgrid(loc_inv, years)
    #ax[1,i].scatter(times[invasion_all>0], locs[invasion_all>0], s = 1, color = "lightblue", alpha = .1)
    ax[3,i].plot(years, dl*np.sum(1/np.sqrt(2*np.pi)*np.exp(-loc_inv**2)
                                *(invasion_all>0), axis = 1), 'b')
    
    # add invasion probability for predators
    for j in years:
        loc_inv, invasion = invasion_prob_pred(data["equi_all"][i,j][surv[i,j]],
                                 data["loc"][i][surv[i,j]],
                                 data["util"][:,i,surv[i,j]],
                                 data["sig"][:,i,surv[i,j]],
                                 data["level"][i,surv[i,j]])
        
        invasion_all[j] = invasion
    locs, times = np.meshgrid(loc_inv, years)
    #ax[2,i].scatter(times[invasion_all>0], locs[invasion_all>0], s = 1, color = "pink", alpha = .1)
        
    ax[3,i].plot(years, dl*np.sum(1/np.sqrt(2*np.pi)*np.exp(-loc_inv**2/2)
                                *(invasion_all>0), axis = 1), 'r')

richness = np.array([np.mean(np.sum(surv
                        & (data["level"][:,np.newaxis] == 0), axis = -1),
                             axis = 0),
                     np.mean(np.sum(surv
                        & (data["level"][:,np.newaxis] == 1), axis = -1),
                             axis = 0)])

ax[0,-1].plot(years, richness[0], 'b', label = "phytoplankton")
ax[0,-1].plot(years, richness[1], 'r', label = "zooplankton")

time = np.ones(surv.shape)*years[:,np.newaxis]
loc = np.ones(surv.shape)*data["loc"][:,np.newaxis]

for i in range(2):
    traits, binx, biny = np.histogram2d(loc[surv & (data["level"][:,np.newaxis] == i)],
                                        time[surv & (data["level"][:,np.newaxis] == i)],
                                        bins = [50, len(years)-1])
    #traits = traits/np.sum(traits, axis = 0, keepdims = True)
    #traits[traits == 0] = np.nan
    ax[1+i,-1].imshow(traits, aspect = "auto", origin = "lower",
                    extent = [biny[0], biny[-1], binx[0], binx[-1]], vmin = 1)

ax[3,-1].plot(years[1:], np.mean((persistence_length>1)
                                 & (data["level"] == 0), axis = 0), 'b')
ax[3,-1].plot(years[1:], np.mean((persistence_length>1)
                                 & (data["level"] == 1), axis = 0), 'r')
ax[0,0].set_xlim([0, 100])
ax[0,1].set_xlim([0, 100])

ax[-1,0].set_ylim([0,1])
ax[1,0].set_ylim([-2.5, 2.5])
ax[2,0].set_ylim([-2.5, 2.5])

# add layout
ax[-1,0].set_xlabel("Time")
ax[-1,0].set_xlabel("Time")
ax[-1,0].set_xlabel("Time")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Trait distribution\nPhytoplankton")
ax[2,0].set_ylabel("Trait distribution\nZooplankton")
ax[3,0].set_ylabel("Invasion probability")

ax[0,0].set_title("Example 1")
ax[0,1].set_title("Example 2")
ax[0,2].set_title("Average case")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEFGHIJKLMNOP"[i], loc = "left")
    
ax[0,-1].legend()
fig.savefig("Figure_LV_predators.pdf")
fig.savefig("Figure_LV_predators.png")
#"""


fig = plt.figure()

a = persistence_length[:,50:] # cut off initial conditions
bins = plt.hist(np.log(a[a>1]), bins = 15)
plt.semilogy()
#plt.semilogx()