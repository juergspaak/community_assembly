import numpy as np
import matplotlib.pyplot as plt

import functions_evolution_assembly as fea
from scipy.integrate import solve_ivp

pars = fea.pars
pars["sig_niche"] = 4

iteras = 5
    
dt = 3
time = np.arange(0,100, 1/dt)
t_org = time[-1]

x_start = np.array([0.1,-0.1])
level = np.zeros(len(x_start))
N_start = np.full(len(level), 0.1)
N_start[level == 1] /= 10
c = ["r" if l else "b" for l in level]

fig, ax = plt.subplots(3,iteras, figsize = (12,9), sharey = "row",
                       sharex = False)

for i in range(iteras):
    if i != 0:
        N_start = np.append(N_start, 0.01)
        x_start = np.append(x_start, np.random.normal(0,fea.pars["sig_niche"]))
    
        level = np.append(level, np.random.binomial(1, 0.0))
    c = ["r" if l else "b" for l in level]
    sol = solve_ivp(fea.ode, time[[0,-1]], np.append(N_start, x_start),
                t_eval = time, args = (level,pars))
    
    for j in range(len(level)):
        ax[0,i].plot(time, sol.y[j], color = c[j])
        ax[1,i].plot(time, sol.y[j+ len(level)], color = c[j])
        ax[1,i].axhline(0, color = "grey", zorder = 0)
        fea.plot_fitness(sol.y[:len(level), -1],
                         sol.y[len(level):, -1], level, ax[2, i],
                         pars = pars)
    N_start = sol.y[:len(level), -1]
    x_start = sol.y[len(level):, -1]
    surv = N_start>1e-3
    N_start = N_start[surv]
    x_start = x_start[surv]
    level = level[surv]
    time += t_org

ax[0,0].set_ylabel("Densities")
ax[1,0].set_ylabel("Traits")
ax[2,0].set_ylabel("Growth rate")

for i in range(iteras):
    ax[0,i].set_xlabel("Time")
    ax[1,i].set_xlabel("Time")
    ax[2,i].set_xlabel("Traits")
    
fig.tight_layout()
fig.savefig("Figure_evolution_basal_only.pdf")