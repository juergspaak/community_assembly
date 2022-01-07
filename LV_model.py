import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import viridis
from scipy.integrate import solve_ivp
from scipy.stats import linregress

def LV_model(t, N, A, mu):
    return N*(mu-A.dot(N))

n_coms = 200
years = 500
colors = viridis(np.linspace(0,1, years))
"""
##############################################################################
# random LV matrix

A = np.random.uniform(0,1,(n_coms, years, years))
A[:,np.arange(years), np.arange(years)] = 1
mu = np.random.uniform(.2, 2, (n_coms, years))

equi = np.full((n_coms, years), np.nan)
present = np.full((n_coms, years+1, years), False, dtype = bool)
present[:,np.arange(years), np.arange(years)] = True
n_specs = np.arange(years)

time = [0,20]

for i in range(n_coms):
    #plt.figure()
    N_start = [mu[i, 0]/A[i,0,0]]
    for j in range(years):
        ind_spec = n_specs[present[i, j]]
        A_c = A[i, ind_spec[:,np.newaxis], ind_spec]
        mu_c = mu[i, ind_spec]
        
        equi = solve_ivp(LV_model, time, N_start, args = (A_c, mu_c))
        
        present[i, j+1, ind_spec[equi.y[:,-1]>1e-4]] = True
        N_start = np.append(equi.y[:,-1], 1e-2)
        N_start = N_start[N_start>1e-4]
    print(i)
        #for ind_equi, k in enumerate(ind_spec):
        #    plt.plot(time[-1] * j + equi.t, equi.y[ind_equi], color = colors[k])
        
richness = np.mean(np.sum(present, axis = -1), axis = 0)
richness[1:-1] -= 1
#"""


###############################################################################
# one-dimensional resource axis

# species traits
res, dr = np.linspace(-5,5,101, retstep = True)
sig_k = 4
max_K = 4
sig_ref = .8
tot_r_ref = 2
mort_ref = .5

def generate_species(com, plot = False):
    mu_r = np.random.normal(0, 1, com)
    #mu_r = np.random.uniform(-4,4, com)
    sigma_r = np.random.normal(.8, .1, com)
    #sigma_r = np.random.uniform(1,1.2, com)
    sigma_r[:] = sig_ref
    tot_r = np.random.normal(2, .2, com)
    tot_r[:] = tot_r_ref
    mort = np.random.normal(.5, .1, com)
    mort[:] = mort_ref
    
    if plot:
        consum = tot_r/np.sqrt(2*np.pi)/sigma_r*np.exp(-(mu_r-res.reshape(-1,1,1))**2/2/sigma_r**2)
        for i in range(com[1]):
            plt.plot(consum[:,0,i])
        plt.show()
    
    sig = (sigma_r[...,np.newaxis]**2 + sigma_r[:,np.newaxis]**2)
    A = (tot_r[...,np.newaxis]*tot_r[:,np.newaxis]/np.sqrt(2*np.pi*sig)*
          np.exp(-(mu_r[...,np.newaxis]-mu_r[:,np.newaxis])**2/2/sig))
    
    # resource availability is given by gaussian
    
    sig_res = (sig_k**2 + sigma_r**2)
    mu = max_K*tot_r/np.sqrt(2*np.pi*sig_res)*np.exp(-mu_r**2/2/sig_res) - mort
    
    return mu, A, mu_r, sigma_r, mort

mu, A,  mu_r, sigma_r, mort = generate_species((n_coms, years), plot = False)
mu[mu<0] = 0
"""A = np.random.uniform(0,1,(n_coms, years, years))
A[:,np.arange(years), np.arange(years)] = 1
mu = np.random.uniform(.2, 2, (n_coms, years))"""

present = np.full((n_coms, years+1, years), False, dtype = bool)
present[:,np.arange(years), np.arange(years)] = True
n_specs = np.arange(years)

equi_all = np.zeros(present.shape)


for i in range(n_coms):
    ind_spec = np.arange(1)
    #equi_all[i,0,ind_spec] = mu[i,0]/A[i,0,0]
    equi = 0
    print(i)
    for j in range(years):
        # check whether all present species are at equlibrium
        """if (equi<=0).any() and j!= 0:
            raise
        if j>=1:
            r_i = mu[i,ind_spec] - A[i, ind_spec[:,np.newaxis], ind_spec].dot(equi_all[i,j-1, ind_spec])
            if np.any(np.abs(r_i)>1e-10):
                raise"""        
        if (mu[i,j] - np.sum(A[i,j, ind_spec]*equi))<1e-10:
            present[i,j+1, ind_spec] = True
            equi_all[i,j] = equi_all[i,j-1]
            continue
        
        ind_spec = n_specs[present[i, j]]
        
        A_c = A[i, ind_spec[:,np.newaxis], ind_spec]
        mu_c = mu[i, ind_spec]
        
        
        equi = -np.ones(len(ind_spec))
        ind = np.full(len(ind_spec), True)
        sp_loc = np.arange(len(ind_spec))
        counter = 0
        while (equi<0).any() or (inv_gr>1e-10).any():
            sp = sp_loc[ind]
            equi[~ind] = 0
            equi[ind] = np.linalg.solve(A_c[sp[:,np.newaxis], sp],
                                            mu_c[sp])
            
            
            if (equi<0).any():
                # remove species with lowest density
                ind[np.argmin(equi)] = False
            else: # to avoid infinite loops:
                # add back species that have positive growth rates
                inv_gr = mu_c - A_c.dot(equi)
                ind[np.argmax(inv_gr)] = True
            counter += 1
            if counter > 2**len(ind):
                # going towards infinite loop, simple community assembly
                # failed, break loop
                present[i,j:] = False
                print("break", i)
                break
        
        ind_spec = ind_spec[ind]
        equi = equi[ind]
        equi_all[i,j, ind_spec] = equi
        present[i, j+1, ind_spec] = True
    equi_all[i,-1, ind_spec] = equi

# check whether community assembly was done correctly
if (equi_all<0).any():
    raise RuntimeError("Computed some equilibrium wrong")
    
# compute whether invasion was done correctly
i_check = np.random.randint(n_coms)
r_i = mu[i_check,np.newaxis] - np.einsum("ij,tj->ti", A[i_check], equi_all[i_check])
surv = present[i_check].copy()
surv[np.arange(years), np.arange(years)] = False
if (r_i[surv]>1e-10).any():
    raise RuntimeError("Computed invasion growth rates wrong")


richness = np.mean(np.sum(present, axis = -1), axis = 0)
richness[1:-1] -= 1# 

##############################################################################
# plot results

fig, ax = plt.subplots(2,3, sharex = False,  figsize = (9,9))
ax[0,0].plot(richness, '.-')

persistence_length = np.sum(present, axis = 1)
ax[0,1].plot(np.mean(persistence_length, axis = 0))

ax[0,2].plot(np.mean(persistence_length>1, axis = 0))
ax[0,2].plot(richness[-1]/(1+np.arange(years)), label = "naive")



x = 1 + np.arange(years)
y = np.mean(persistence_length>1, axis = 0)
t_cut = 400
s, inter, r, p, std = linregress(x[:t_cut], (x*y)[:t_cut])
print(s, std, s-2*std)
t_cut = years
s, inter, r, p, std = linregress(x[:t_cut], (x*y)[:t_cut])
print(s, std, s-2*std)
ax[0,2].plot(x, inter/x + s)
ax[0,2].axhline(s, color = "green", linestyle = "--")

ax[0,2].set_ylim([0,.5])



ax[1,0].hist(mu_r.flatten(), bins=  25, density = True)
ax[1,0].hist(mu_r[present[:,-1]], bins = 25, alpha = .5, density = True)

ax[1,1].hist(sigma_r.flatten(), bins=  25, density = True)
ax[1,1].hist(sigma_r[present[:,-1]], bins = 25, alpha = .5, density = True)


#ax[1,2].hist(mu.flatten(), bins=  25, density = True)
#ax[1,2].hist(mu[present[:,-1]], bins = 25, alpha = .5, density = True)
ax[1,2].hist(mort.flatten(), bins=  25, density = True)
ax[1,2].hist(mort[present[:,-1]], bins = 25, alpha = .5, density = True)

ax[0,0].set_xlabel("Year")
ax[0,0].set_ylabel("Species richness")

ax[0,1].set_title("persistence time")
ax[0,2].set_title("Invasion probability")

ax[1,0].set_title("Location")
ax[1,1].set_title("Variance")
ax[1,2].set_title("Total consumption")

#fig.savefig("Figure_example_LV_resources.pdf")
def plot(i, t = -1):
    ind = present[i,t]
    plt.figure()
    consum = (mu[i, ind]/np.sqrt(2*np.pi)/sigma_r[i, ind]*
                np.exp(-(mu_r[i, ind]-res.reshape(-1,1,1))**2/2/sigma_r[i, ind]**2))
        
    for i in range(consum.shape[-1]):
        plt.plot(res, consum[:,0,i])

def plot2(i):
    fig = plt.figure()
    for j in np.arange(years)[np.sum(present[i], axis = 0)>1]:
        plt.plot(np.arange(years+1)[present[i,:,j]],
                 np.full(np.sum(present[i,:,j]), mu_r[i,j]), zorder = 3,
                 linewidth=3)
    plt.axis([None, None, -4,4])
    return fig
        
def invasion_fun(i, t = years, plot = False):
    ind = np.arange(t)[present[i,t,:t]]
    equi = np.linalg.solve(A[i, ind[:,np.newaxis], ind], mu[i, ind])
    if (equi<0).any():
        raise
    
    mu_invader = np.linspace(*res[[0,-1]], 1000)
    sig = 2*sig_ref**2
    tr = 2 # total consumptoin
    A_invader = (tr*tr/np.sqrt(2*np.pi*sig)*
          np.exp(-(mu_invader[:,np.newaxis]-mu_r[i, ind])**2/2/sig))
    
    sig_res = (sig_k**2 + sig_ref**2)
    mu_inv = max_K*tot_r_ref/np.sqrt(2*np.pi*sig_res)*np.exp(-mu_invader**2/2/sig_res) - mort_ref
    invasion = mu_inv - np.sum(A_invader*equi, axis = 1)
    
    if plot:
        plt.figure()
        plt.plot(mu_invader, invasion)
        plt.axhline(0, color = "grey")
        
        plt.plot(mu_r[i, ind], np.zeros(len(ind)), 'r.')
        plt.ylim([-.1, .1])
        plt.title("{} {}".format(i, t))
    return mu_invader, invasion

def plot3(i):
    fig = plot2(i)
    for year in range(years):
        loc, invasion  = invasion_fun(i, year)
        plt.scatter(np.full(np.sum(invasion>0), year), loc[invasion>0],
                    color = "grey", s = .5, alpha = .1,
                    zorder = 0)

#"""