import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

sig_niche = 5
sig_basal = 2
m_basal = -0.1
m_pred = -0.2
sig_pred = 0.8
niche_width = sig_niche*np.sqrt(np.log(1/m_basal**2))

def LV_parameters(x, level):
    
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # intrinsic growth rates follow a gaussian
    mu = np.exp(-x**2/(2*sig_niche**2)) + m_basal
    mu[pred] = m_pred # predators have negative intrinsic growth rate

    # interaction matrix
    # basal basal interaction
    A = np.exp(-(x-x[:,np.newaxis])**2/(2*sig_basal**2))
    A[pred, basal[:,np.newaxis]] = -np.exp(-(x-x[:,np.newaxis])**2/(2*sig_pred**2))[pred, basal[:,np.newaxis]]
    A[basal[:,np.newaxis], pred] = np.exp(-(x-x[:,np.newaxis])**2/(2*sig_pred**2))[basal[:,np.newaxis], pred]
    A[pred, pred[:, np.newaxis]] = 0
    
    return mu, A

def growth_rates(N, x, level):
    
    mu, A = LV_parameters(x, level)
    return N*(mu - A.dot(N))

def evolution(N, x, level):
    
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # current parameters
    mu, A  = LV_parameters(x, level)
    
    # change in intrinsic growth rate
    dmu = -x/sig_niche**2*(mu - m_basal)
    dmu[level == 1] = 0
    
    # how basal species interaction change
    dA = (x - x[:,np.newaxis])*A
    dA[basal[:,np.newaxis], basal] /= sig_basal**2
    dA[pred, basal[:,np.newaxis]] /= sig_pred**2
    dA[basal[:,np.newaxis], pred] /= sig_pred**2
    
    """
    # compare to numerical differentiation, works good
    a = dmu - dA.dot(N)
    b = np.empty(len(a))
    dx_org = 1e-6
    for i in range(len(level)):
        dx = x.copy()
        dx[i] += dx_org
        dmu2, dA2 = LV_parameters(dx, level)
        b[i] = ((dmu2-dA2.dot(N) - (mu - A.dot(N)))/dx_org)[i]
    
    if not np.allclose(a,b, rtol=1e-3, atol = 1e-8):
        print(N, x, level)
        print(a, b)
        print((a-b)/a)
        raise"""
    evol_speed = np.full(len(level), 5)
    evol_speed[pred] = 10
    return evol_speed*(dmu - dA.dot(N))

def ode(t, Nx, level):
    N = Nx[:len(level)]
    x = Nx[len(level):]
    
    return np.append(growth_rates(N, x, level), evolution(N, x, level))

def fitness(N, x, level, y_range = [-niche_width, niche_width]):
    
    # predator and basal species
    pred = np.where(level == 1)[0]
    basal = np.where(level == 0)[0]
    
    # all possible invader locations
    y = np.linspace(*y_range, 1001)[:,np.newaxis]
    
    fit_ref = (np.exp(-y**2/(2*sig_niche**2)) + m_basal)[:,0]
    
    # intrinsic growth rate
    fit_basal = (np.exp(-y**2/(2*sig_niche**2)) + m_basal)[:,0]
    # minus competition
    fit_basal -= np.sum((N*np.exp(-(y - x)**2/(2*sig_basal)))[:,basal], axis = 1)
    # minus predation
    fit_basal -= np.sum((N*np.exp(-(y - x)**2/(2*sig_pred)))[:,pred], axis = 1)
    
    # fitness of predators
    fit_pred = np.sum((N*np.exp(-(y - x)**2/(2*sig_pred)))[:,basal], axis = 1)+m_pred
    
    return y, fit_basal, fit_pred, fit_ref

def plot_fitness(N, x, level, ax, y_range = [-niche_width, niche_width]):
    
    traits, fit_basal, fit_pred, fit_ref = fitness(N,x, level, y_range = y_range)
    ax.plot(traits, fit_basal, '--b')
    ax.plot(traits, fit_pred, "--r")
    ax.plot(traits, fit_ref, "--g")
    ax.axhline(0, color = "grey", zorder = 0)
    
    growth = growth_rates(N, x, level)/N
    for i in range(len(level)):
        ax.plot(x[i], growth[i], ".", color = "r" if level[i] else "b")
        
    ax.set_xlabel("Species traits")
    ax.set_ylabel("Growth rate")
    ax.set_xlim(y_range)
    

if __name__ == "__main__":
    iteras = 1
    
    dt = 3
    time = np.arange(0,100, 1/dt)
    t_org = time[-1]
    
    x_start = np.array([0.1,-0.1])
    level = np.zeros(len(x_start))
    level[-1] = 1
    N_start = np.full(len(level), 0.1)
    N_start[level == 1] /= 10
    c = ["r" if l else "b" for l in level]
    
    fig, ax = plt.subplots(2,4, figsize = (12,9), sharey = False,
                           sharex = False)
    
    for i in range(iteras):
        if i != 0:
            N_start = np.append(N_start, 0.01)
            x_start = np.append(x_start, np.random.normal(0,sig_niche/2))
        
            level = np.append(level, np.random.binomial(1, 0.5))
        c = ["r" if l else "b" for l in level]
        sol = solve_ivp(ode, time[[0,-1]], np.append(N_start, x_start),
                    t_eval = time, args = (level,))
        
        for i in range(len(level)):
            ax[0,0].plot(time, sol.y[i], color = c[i])
            ax[1,0].plot(time, sol.y[i+ len(level)], color = c[i])
            ax[1,0].axhline(0, color = "grey")
        """N_start = sol.y[:len(level), -1]
        x_start = sol.y[len(level):, -1]
        surv = N_start>1e-3
        N_start = N_start[surv]
        x_start = x_start[surv]
        level = level[surv]
        time += t_org"""
        
        
    ax[0,0].set_ylabel("Species densities")
    ax[1,0].set_ylabel("Species trait location")
    ax[1,0].set_xlabel("Time")
    
    ax_f = ax[:,1:].flatten()
    t = [11, 12, 30, 35, 50, 99]
    y_ranges = [[-niche_width, niche_width],
                [-niche_width, niche_width],
                [-15,-5],[-15,-5], [-15, -5], [-niche_width, niche_width]]
    for j in range(len(t)):
        plot_fitness(sol.y[:len(level), dt*t[j]],
                 sol.y[len(level):, dt*t[j]], level, ax_f[j], y_range = y_ranges[j])
        ax_f[j].text(0.8,0.9, "t={}".format(t[j]), transform=ax_f[j].transAxes)
        
    fig.tight_layout()
    
