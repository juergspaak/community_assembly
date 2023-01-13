import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import solve_ivp

# each species has a mass
itera = 100

m_i = np.array([1, 100]) # species pass
s_i = np.array([1e-10, 1]) # selectivity center
f_i = np.array([1e-10, 1]) # feeding center

specs = {"m_i": np.log(m_i),
         "s_i": s_i, "f_i": np.log(f_i),
         "n_spec": 2,
         "K": 100,
         "R": 1}
epsilon = 2e-4 # extinciton threshhold
B_start = np.array([100, epsilon*specs["m_i"][-1]])

"""
m_i = np.random.normal(15,3, itera)
s_i = np.random.uniform(0.5, 1.5, itera)
f_i = np.random.uniform(-np.log(1000), -np.log(3), itera)
f_i = m_i + f_i
m_i[0] = 0
f_i[0] = np.log(1e-10)
s_i[0] = 1e-10
specs = {"m_i": m_i,
         "s_i": s_i, "f_i": f_i,
         "n_spec": itera,
         "K": 100,
         "R": 1}
B_start = np.full(itera, 1e-3)
B_start[0] = specs["K"]"""



speciation = 1e-3 # proability of speciation



# speciation
def new_species(specs):
    i = np.random.randint(1, specs["n_spec"], 1)
    m_i_mutant = specs["m_i"][i]
    m_i_mutant = np.random.uniform(m_i_mutant - np.log(2),
                                   m_i_mutant + np.log(2))
    specs["m_i"] = np.append(specs["m_i"],m_i_mutant)
    specs["s_i"] = np.append(specs["s_i"], np.random.uniform(0.5,1.5))
    specs["s_i"][-1] = 0.5 + m_i_mutant/(3+m_i_mutant)
    f_i_mutant = np.random.uniform(m_i_mutant - np.log(1000),
                                   m_i_mutant - np.log(3))
    specs["f_i"] = np.append(specs["f_i"], f_i_mutant)
    specs["n_spec"] += 1
    return specs
    
def extinction(specs, extinct):
    survive = ~extinct
    for key in ["m_i", "f_i", "s_i"]:
        specs[key] = specs[key][survive]
    
    specs["n_spec"] = sum(survive)
    return specs
    
    
def LV_params2(specs):
    
    # base consumption speed
    a_i = np.exp(specs["m_i"]*0.75)   
    
    # depending on optimal trait size
    a_ij = np.exp(-(specs["m_i"][:,np.newaxis] - specs["f_i"])**2/
                  2/specs["s_i"]**2)
    a_ij = (a_i*a_ij/(specs["s_i"]*np.sqrt(2*np.pi))).T
    
    e_j = np.full(specs["n_spec"], 0.15)
    e_j[0] = 0.1 # plant species can be assimilated less efficiently
    
    
    # actual linear interaction depends on consumption and predation
    a_ij = (e_j/np.exp(specs["m_i"])[:,np.newaxis]*a_ij # predation
            - 1/(np.exp(specs["m_i"])*a_ij).T) # consumption
    
    mu = np.zeros(specs["n_spec"])
    mu[:] = -0.1
    mu[0] = specs["R"]
    
    return mu, a_ij

def LV_params(specs):
    m_i = specs["m_i"]
    f_i = specs["f_i"]
    
    a_ij = np.exp(-(m_i[:,np.newaxis] - f_i)**2/(2*specs["s_i"]**2))
    m_i = np.exp(m_i)
    a_ij *= m_i**0.75/(specs["s_i"]*np.sqrt(2*np.pi))
    a_ij = a_ij.T
    
    g_ij = a_ij/m_i[:,np.newaxis]
    
    e_j = np.full(specs["n_spec"], 0.85)
    e_j[0] = 0.45 # plant species can be assimilated less efficiently
    
    a_ij = (e_j*g_ij - g_ij.T)
    a_ij[np.arange(specs["n_spec"]), np.arange(specs["n_spec"])] -= 0.0
    a_ij[0,0] = -specs["R"]/specs["K"]
    
    mu = -0.314*m_i**(-0.25) # mass respiration
    mu[0] = specs["R"]
    return mu, a_ij
    
    

n_time = int(1e6)
speciation = np.random.binomial(1, speciation, n_time) == 1
event_speciation = np.append(0,np.arange(n_time)[speciation])

dens = np.full((len(event_speciation)-1, len(event_speciation)), np.nan)
m_is = [specs["m_i"][-1]]
f_is = [specs["f_i"][-1]]
s_is = [specs["s_i"][-1]]
present = np.array([0,1])


for i in range(len(event_speciation)-1):
    print(event_speciation[i], specs["n_spec"])
    if i != 0:
        specs = new_species(specs)
        B_start = np.append(B_start, specs["m_i"][-1]*epsilon)
        m_is.append(specs["m_i"][-1])
        present = np.append(present,i+1)
    time = [event_speciation[i], event_speciation[i+1]]
    #time = [0,10000]
    mu, a_ij = LV_params(specs)
    
    sol = solve_ivp(lambda t, N: N*(mu + a_ij.dot(N)),
                    time, B_start)
    
    dens[i,present] = sol.y[:,-1]
    
    extinct = sol.y[:,-1]<epsilon    
    extinct[0] = False # plant species can't go extinct
    present = present[~extinct]
    specs = extinction(specs, extinct)
    B_start = sol.y[~extinct, -1]
#"""
    
m_is = np.array(m_is)

fig, ax = plt.subplots(2,2, sharex = True, figsize = (9,9))
ax[0,0].plot(dens[:,0])
ax[0,0].set_ylabel("Resource density")

ax[1,0].plot(np.sum(dens>0, axis = 1))
ax[1,0].set_ylabel("Species richness")

pres = dens>epsilon
ax[0,1].plot(np.sum(pres[1:] & ~pres[:-1], axis = 1), label = "invasion")
ax[0,1].plot(np.sum(pres[:-1] & ~pres[1:], axis = 1), label = "extinciton")
ax[0,1].legend()

step = max(1, len(pres)//1000)
for i, p in enumerate(pres[::step]):
    ax[1,1].scatter(np.full(np.sum(p[1:]), step*i), m_is[p[1:]], color = "k", s = 1)
