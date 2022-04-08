import numpy as np
import matplotlib.pyplot as plt

from scipy.integrate import simps

# default traits
sigs = np.array([1, 0.8])
alpha_max = np.array([1.0,1.0])
omega = 2
mu_max = 1
ms = np.array([.1,.1])

def res_axis(loc, sig):
    # resource axis for integration
    return np.linspace(np.amin(loc)-3*np.amax(sig),
                        np.amax(loc)+3*np.amax(sig), 1001)

keys = ["Gaussian", "Flat Gaussian", "Flat", "Triangular", "Quadratic",
        "Asymetric"]
scaling = {}

kernel_fun = {
    "Gaussian": lambda theta: np.exp(-theta**2),
    "Flat Gaussian": lambda theta: np.exp(-theta**4),
    "Flat": lambda theta: np.where(np.abs(theta)<1, 1, 0),
    "Triangular": lambda theta: np.where(np.abs(theta)<1, 1-np.abs(theta),0),
    "Quadratic": lambda theta: np.where(np.abs(theta)<1, 1-np.abs(theta)**2, 0)
    }

def asymetric_fun(theta):
    u_il = np.where(theta>0, 1-theta, 1 + theta/3)
    u_il[u_il<0] = 0
    return u_il

kernel_fun["Asymetric"] = asymetric_fun

for key in keys:
    theta = np.linspace(-3,3,10001)
    u_il = kernel_fun[key](theta)
    u_2 = simps(u_il**2, theta, axis = 0)
    u_1 = simps(u_il, theta, axis = 0)
    scaling[key] = [u_1/u_2, u_2/u_1**2]

def kernel_functions(loc, sig, x, kernel = "Gaussian", scaling = scaling):
    x = x.reshape(-1,1)
    theta = (loc-x)/sig/scaling[kernel][1]
    u_il = 1/np.sqrt(sig)*scaling[kernel][0]*kernel_fun[kernel](theta)
    return u_il, x[:,0]

def gaussian_kernel(loc, sig, x, scaling = scaling["Gaussian"]):
    # gaussian competition kernel    
    return np.exp(-1/2*((loc-x[:,np.newaxis])/sig)**2), x

def flat_gaussian_kernel(loc, sig, x, scaling = scaling["Flat Gaussian"]):
    # gaussian with exponent of 4
    return np.exp(-1/2*((loc-x[:,np.newaxis])/sig)**4), x

def flat_kernel(loc, sig, x, scaling = scaling["Flat"]):
    # indicator funciton
    return np.where(np.abs(loc-x[:,np.newaxis])<2*sig, 1,0), x

def triangular_kernel(loc, sig, x, scaling = scaling["Triangular"]):
    # triangular kernel shape
    x = x.reshape(-1,1)
    u_il = 1-np.abs(loc-x)/(2*sig)
    u_il[u_il<0] = 0
    return u_il, x[:,0]

def quadratic_kernel(loc, sig, x, scaling = scaling["Quadratic"]):
    # quadratic polynomial
    x = x.reshape(-1,1)
    u_il = 1-(loc-x)**2/(2*sig)**2
    u_il[u_il<0] = 0
    return u_il, x[:,0]

def asymetric_kernel(loc, sig, x, scaling = scaling["Asymetric"]):
    # assymetric triangular kernel
    x = x.reshape(-1,1)
    u_il = np.where(loc>x, 1-np.abs(loc-x)/(1.5*sig), 1-np.abs(loc-x)/(3*sig))
    u_il[u_il<0] = 0
    return u_il, x[:,0]



kernels = {"Gaussian": gaussian_kernel,
           "Flat Gaussian": flat_gaussian_kernel,
           "Flat": flat_kernel,
           "Triangular": triangular_kernel,
           "Quadratic": quadratic_kernel,
           "Asymetric": asymetric_kernel}

def compute_LV_from_traits_kernels(level, loc, sig, alpha_max, m, omega = omega,
                           mu_max = mu_max, kernel = "Gaussian"):
    
    id_prey = np.where(level == 0)[0]
    id_pred = np.where(level == 1)[0]
    
    # get resource consumption vectors
    #u_il, res_loc = kernels[kernel](loc[id_prey], sig[id_prey],
    #                                res_axis(loc[id_prey], sig[id_prey]))
    u_il, res_loc = kernel_functions(loc[id_prey], sig[id_prey],
                                    res_axis(loc[id_prey], sig[id_prey]),
                                    kernel)
    
    # entire interaction matrix
    A = np.full((len(level), len(level)), np.nan)
    # species interaction if all species were prey
    A[id_prey, id_prey[:,np.newaxis]] = simps(
        u_il[:,np.newaxis]*u_il[...,np.newaxis], axis = 0, x = res_loc)
    # scale
    #A /= np.mean(np.diag(A)[id_prey])
    A *= alpha_max[:,np.newaxis]*alpha_max
    
    # effect of predator on prey
    u_il, res_loc = kernels[kernel](loc, sig, loc)
    u_il, res_loc = kernel_functions(loc, sig, loc, kernel)
    A[id_prey[:,np.newaxis], id_pred] = alpha_max[id_pred]*u_il[id_prey[:,np.newaxis], id_pred]
    
    # effect of prey on predator
    A[id_pred[:,np.newaxis], id_prey] = -A[id_prey[:,np.newaxis], id_pred].T 
    
    # predators don't interact with each other
    A[id_pred[:,np.newaxis], id_pred] = 0
    
    # compute the intrinsic growth rates
    mu = mu_max*np.exp(-loc**2/2/omega**2) - m
    
    # mortality rate of predators
    mu[id_pred] = -m[id_pred] 
    
    return mu, np.round(A,8)

if __name__ == "__main__":
    
    for key in keys:
        u_il, x = kernel_functions(0, 2, theta, key)
        plt.plot(x, u_il[:,0], label = key)
        print(key, simps(u_il, x = theta, axis = 0), simps(u_il**2, x = theta, axis = 0))
    plt.legend()