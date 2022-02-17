import numpy as np
import matplotlib.pyplot as plt
import assembly_functions as af

fig, ax = plt.subplots(3,2, sharex = False,
                       sharey = "row", figsize = (8,8))

# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

try:
    
    datas
except NameError:
    d_prey = np.load(path + "data_LV_assembly.npz")
    d_pred = np.load(path + "data_LV_assembly_predators.npz")
    datas = [{key: d_prey[key] for key in d_prey.files},
             {key: d_pred[key] for key in d_pred.files}]


##############################################################################
# plot presences over time
com = 0 # take the first community

for j in range(2):
    data = datas[j]
    surv = data["equi_all"][com]>0
    time = np.arange(data["equi_all"].shape[1])
    n_spec = data["loc"].shape[1]
    for i in range(n_spec):
        pres = surv[:,i]
        if not pres.any():
            continue # species is never present
        ax[0,j].plot(time[pres], np.full(np.sum(pres), data["loc"][com, i]),
                     "br"[data["level"][com, i]])
        
ax[0,0].set_ylabel("Species mean trait")
ax[0,0].set_xlabel("Time")
ax[0,1].set_xlabel("Time")

##############################################################################
# plot species richness
ax[1,0].plot(time, np.sum(datas[0]["equi_all"][com]>0, axis = 1), 'b--',
             label = "Example")
ax[1,0].plot(time, np.mean(np.sum(datas[0]["equi_all"]>0, axis = -1),axis = 0),
             'b-', label = "Average")
ax[1,0].legend()

# plot for case with predators
id_pred = datas[1]["level"] == 1


ax[1,1].plot(time, np.sum((datas[1]["equi_all"]>0) &
                           (datas[1]["level"][:,np.newaxis] == 0)
                             ,axis = -1)[0], 'b--')
ax[1,1].plot(time, np.sum((datas[1]["equi_all"]>0) &
                           (datas[1]["level"][:,np.newaxis] == 1)
                             ,axis = -1)[0], 'r--')

ax[1,1].plot(time, np.mean(np.sum((datas[1]["equi_all"]>0) &
                           (datas[1]["level"][:,np.newaxis] == 0), axis = -1)
                             ,axis = 0), 'b-')
ax[1,1].plot(time, np.mean(np.sum((datas[1]["equi_all"]>0) &
                           (datas[1]["level"][:,np.newaxis] == 1), axis = -1)
                             ,axis = 0), 'r-')

# axis labeling
ax[1,0].set_ylabel("Species richness")
ax[1,0].set_xlabel("Time")
ax[1,1].set_xlabel("Time")

##############################################################################
# species persistence time
treshhold = 150 # initial phase not counted towards persistence time
bins = np.arange(1, 400, 5)
persistence_length = np.sum(datas[0]["equi_all"][:,treshhold:]>0, axis = 1)
hist, bins = np.histogram(persistence_length.flatten(), bins = bins, density = True)
x = (bins[1:] + bins[:-1])/2
ax[2,0].plot(x, hist, 'b')

persistence_length = np.sum(datas[1]["equi_all"][:,treshhold:]>0, axis = 1)
hist, bins = np.histogram(persistence_length[datas[1]["level"] == 0].flatten(),
                          bins = bins, density = True)
ax[2,1].plot(x, hist, 'b', label = "Prey")
hist, bins = np.histogram(persistence_length[datas[1]["level"] == 1].flatten(),
                          bins = bins, density = True)
ax[2,1].plot(x, hist, 'r', label = "Predator")
ax[2,1].legend()

ax[2,1].set_xlim([0,100])

ax[2,0].semilogy()

ax[2,0].set_xlabel("Persistence time")
ax[2,1].set_xlabel("Persistence time")
ax[2,0].set_ylabel("Frequency")

for i, a in enumerate(ax.flatten()):
    a.set_title("ABCDEF"[i], loc = "left")
    
ax[0,0].set_title("Prey only")
ax[0,1].set_title("Predator\nand Prey")

fig.tight_layout()
fig.savefig("Figure_example_dynamics.pdf")

