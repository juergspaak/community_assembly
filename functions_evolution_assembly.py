import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

pars = dict(sig_niche = 5,
            sig_basal = 2,
            sig_pred = 0.8,
            sig_eps = 1, mu_eps = 0.1)


m_basal = -0.1
m_pred = -0.2
pars["niche_width"] = pars["sig_niche"]*np.sqrt(np.log(1/m_basal**2))


def LV_parameters(x, level, pars = pars):
    
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # intrinsic growth rates follow a gaussian
    mu = np.exp(-x**2/(2*pars["sig_niche"]**2)) + m_basal
    mu += pars["mu_eps"]*np.exp(-x**2/(2*pars["sig_eps"]**2))
    mu[pred] = m_pred # predators have negative intrinsic growth rate

    # interaction matrix
    # basal basal interaction
    A = np.exp(-(x-x[:,np.newaxis])**2/(2*pars["sig_basal"]**2))
    A[pred, basal[:,np.newaxis]] = -np.exp(-(x-x[:,np.newaxis])**2/(2*pars["sig_pred"]**2))[pred, basal[:,np.newaxis]]
    A[basal[:,np.newaxis], pred] = np.exp(-(x-x[:,np.newaxis])**2/(2*pars["sig_pred"]**2))[basal[:,np.newaxis], pred]
    A[pred, pred[:, np.newaxis]] = 0
    
    return mu, A

def growth_rates(N, x, level, pars = pars):
    
    mu, A = LV_parameters(x, level, pars)
    return N*(mu - A.dot(N))

def evolution(N, x, level, pars = pars):
    
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # current parameters
    mu, A  = LV_parameters(x, level, pars)
    
    # change in intrinsic growth rate
    dmu = -x/pars["sig_niche"]**2*np.exp(-x**2/(2*pars["sig_niche"]**2)) \
         -x/pars["sig_eps"]**2*pars["mu_eps"]*np.exp(-x**2/(2*pars["sig_eps"]**2))
    dmu[level == 1] = 0
    
    # how basal species interaction change
    dA = (x - x[:,np.newaxis])*A
    dA[basal[:,np.newaxis], basal] /= pars["sig_basal"]**2
    dA[pred, basal[:,np.newaxis]] /= pars["sig_pred"]**2
    dA[basal[:,np.newaxis], pred] /= pars["sig_pred"]**2
    

    evol_speed = np.full(len(level), 5.0)
    #evol_speed[pred] = 10
    #evol_speed[:] = 10
    evol_speed[:] = .1
    return evol_speed*(dmu - dA.dot(N))

def ode(t, Nx, level, pars = pars):
    N = Nx[:len(level)]
    x = Nx[len(level):]
    
    return np.append(growth_rates(N, x, level, pars),
                     evolution(N, x, level, pars))

def fitness(N, x, level, pars = pars, y_range = None):
    
    if y_range is None:
        y_range = pars["sig_niche"]*np.sqrt(np.log(1/m_basal**2))
        y_range = [-y_range, y_range]
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # all possible invader locations
    y = np.linspace(*y_range, 1001)[:,np.newaxis]
    
    fit_ref = (np.exp(-y**2/(2*pars["sig_niche"]**2)) + m_basal
               + pars["mu_eps"]*np.exp(-y**2/(2*pars["sig_eps"]**2)))[:,0]
    
    # intrinsic growth rate
    fit_basal = fit_ref.copy()
    # minus competition
    fit_basal -= np.sum((N*np.exp(-(y - x)**2/(2*pars["sig_basal"]**2)))[:,basal], axis = 1)
    # minus predation
    fit_basal -= np.sum((N*np.exp(-(y - x)**2/(2*pars["sig_pred"]**2)))[:,pred], axis = 1)
    
    # fitness of predators
    fit_pred = np.sum((N*np.exp(-(y - x)**2/(2*pars["sig_pred"]**2)))[:,basal], axis = 1)+m_pred
    
    return y, fit_basal, fit_pred, fit_ref

def plot_fitness(N, x, level, ax, pars = pars, y_range = None):
        
    traits, fit_basal, fit_pred, fit_ref = fitness(N,x, level, pars, y_range)
    ax.plot(traits, fit_basal, '--b')
    ax.plot(traits, fit_pred, "--r")
    ax.plot(traits, fit_ref, "--g")
    ax.axhline(0, color = "grey", zorder = 0)
    
    growth = growth_rates(N, x, level)/N
    for i in range(len(level)):
        ax.plot(x[i], growth[i], ".", color = "r" if level[i] else "b")
        
    ax.set_xlabel("Species traits")
    ax.set_ylabel("Growth rate")
    ax.set_xlim(traits[[0,-1]])
    

if __name__ == "__main__":
    iteras = 1
    label = ["Basal", "Predator"]
    dt = 3
    time = np.arange(0,100, 1/dt)
    
    x_start = np.array([0.1,-0.1])
    #x_start = np.linspace(-10,10, 5)
    level = np.zeros(len(x_start))
    level[-1] = 1
    N_start = np.full(len(level), 0.1)
    N_start[level == 1] /= 10
    c = ["r" if l else "b" for l in level]
    
    fig, ax = plt.subplots(2,4, figsize = (12,9), sharey = False,
                           sharex = False)
    

    c = ["r" if l else "b" for l in level]
    sol = solve_ivp(ode, time[[0,-1]], np.append(N_start, x_start),
                t_eval = time, args = (level, pars))
    
    for i in range(len(level)):
        ax[0,0].plot(time, sol.y[i], color = c[i])
        ax[1,0].plot(time, sol.y[i+ len(level)], color = c[i], label = label[i])
        ax[1,0].axhline(0, color = "grey", zorder = 0, linewidth =3)
        
    ax[1,0].legend() 
    ax[0,0].set_ylabel("Species densities")
    ax[1,0].set_ylabel("Species trait location")
    ax[1,0].set_xlabel("Time")
    
    ax_f = ax[:,1:].flatten()
    t = [13, 14, 30, 35, 45, 99]
    y_ranges = [[-pars["niche_width"], pars["niche_width"]],
                [-pars["niche_width"], pars["niche_width"]],
                [5,15],[5,15], [5, 15],
                [-pars["niche_width"], pars["niche_width"]]]
    for j in range(len(t)):
        plot_fitness(sol.y[:len(level), dt*t[j]],
                 sol.y[len(level):, dt*t[j]], level, ax_f[j], y_range = y_ranges[j])
        ax_f[j].text(0.8,0.9, "t={}".format(t[j]), transform=ax_f[j].transAxes)
        
    fig.tight_layout()
    
    #fig.savefig("Figure_example_evolution1.pdf")
    
    ###########################################################################
    # investigate effect of sig_pred
    
    fig, ax = plt.subplots(2,3, figsize = (12,9), sharey = "row",
                           sharex = True)
    time = np.arange(0,1000, 1/dt)
    x_start = np.array([0.1,-0.1])
    #x_start = np.linspace(-10,10, 5)
    level = np.zeros(len(x_start))
    level[-1] = 1
    N_start = np.full(len(level), 0.1)
    N_start[level == 1] /= 10
    c = ["r" if l else "b" for l in level]
    
    sigs = [0.8, 1.5, 2.5]
    
    
    for j, sig in enumerate(sigs):
        pars["sig_pred"] = sig
        sol = solve_ivp(ode, time[[0,-1]], np.append(N_start, x_start),
                t_eval = time, args = (level, pars))
        for i in range(len(level)):
            ax[0,j].plot(time, sol.y[i], color = c[i])
            ax[1,j].plot(time, sol.y[i+ len(level)], color = c[i],
                         label = label[i])
            ax[1,j].axhline(0, color = "grey", zorder = 0, linewidth =3)
            
    ax[0,0].set_ylabel("Densities")
    ax[1,0].set_ylabel("Traits")
    ax[1,0].set_xlabel("Time")
    ax[1,1].set_xlabel("Time")
    ax[1,2].set_xlabel("Time")
    
    ax[0,0].set_title("Low predation width")
    ax[0,1].set_title("Middle predation width")
    ax[0,2].set_title("wide predation width")
    ax[1,0].legend()
            
    fig.tight_layout()
            
    fig.savefig("Figure_evolution_effect_of_sig_pred.pdf")
    
    
