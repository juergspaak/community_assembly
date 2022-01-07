import numpy as np
import matplotlib.pyplot as plt
from string import ascii_uppercase as ABC
import assembly_functions as af
import functions_for_plotting as fp

path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
file = path + "Data_LV_{}.npz"
fig_file = "Figure_LV_{}.pdf"
##############################################################################
# lower trophic level only

name = "only_prey"
try:
    data = np.load(file.format(name), allow_pickle = True)
    species_id = {**data}
    surv = data["surv"]
except FileNotFoundError:
    
    # generate species only for the first trophic level
    species_id = af.generate_species(level = [1,0,0])
    
    # run community assembly
    present, equi_all, surv = af.community_assembly(species_id)
    
    # save data 
    np.savez(file.format(name), **species_id, equi_all = equi_all,
            surv = surv, present = present)
    
fig, ax = plt.subplots(3,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)
fp.plot_invasion_prob(ax[-1,-1], surv, species_id)

fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)

for i in range(3):
    ax[-1,i].set_xlabel("Year")
    
for i, a in enumerate(ax.flatten()):
    a.set_title(ABC[i], loc = "left")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Species traits")
ax[2,0].set_ylabel("Invasion probability")
ax[0,0].set_title("Example community 1")
ax[0,1].set_title("Exmple community 2")
ax[0,-1].set_title("Average case")

fig.savefig(fig_file.format(name))

##############################################################################
# both trophic levels

name = "pred_and_prey"
try:
    data = np.load(file.format(name), allow_pickle = True)
    species_id = {**data}
    surv = data["surv"]
except FileNotFoundError:
    
    # generate species only for the first trophic level
    species_id = af.generate_species()
    
    # run community assembly
    present, equi_all, surv = af.community_assembly(species_id)
    
    # save data 
    np.savez(file.format(name), **species_id, equi_all = equi_all,
            surv = surv, present = present)
    
fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)
fp.plot_traits(ax[2,:2], surv, species_id)
fp.plot_invasion_prob(ax[-1,-1], surv, species_id)

fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)
fp.plot_mean_traits(ax[2,-1], surv, species_id, 1)

for i in range(3):
    ax[-1,i].set_xlabel("Year")
    
for i, a in enumerate(ax.flatten()):
    a.set_title(ABC[i], loc = "left")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Species traits")
ax[2,0].set_ylabel("Invasion probability")
ax[0,0].set_title("Example community 1")
ax[0,1].set_title("Exmple community 2")
ax[0,-1].set_title("Average case")

fig.savefig(fig_file.format(name))

##############################################################################
# uniform trait distributions

name = "uniform_pred_prey"
try:
    data = np.load(file.format(name), allow_pickle = True)
    species_id = {**data}
    surv = data["surv"]
except FileNotFoundError:
    
    # generate species only for the first trophic level
    species_id = af.generate_species()
    species_id["loc"] = np.random.uniform(-2,2, species_id["loc"].shape)
    
    # run community assembly
    present, equi_all, surv = af.community_assembly(species_id)
    
    # save data 
    np.savez(file.format(name), **species_id, equi_all = equi_all,
            surv = surv, present = present)
    
fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)
fp.plot_traits(ax[2,:2], surv, species_id)
fp.plot_invasion_prob(ax[-1,-1], surv, species_id)

fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)
fp.plot_mean_traits(ax[2,-1], surv, species_id, 1)

for i in range(3):
    ax[-1,i].set_xlabel("Year")
    
for i, a in enumerate(ax.flatten()):
    a.set_title(ABC[i], loc = "left")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Species traits")
ax[2,0].set_ylabel("Invasion probability")
ax[0,0].set_title("Example community 1")
ax[0,1].set_title("Exmple community 2")
ax[0,-1].set_title("Average case")

fig.savefig(fig_file.format(name))

##############################################################################
# uniform trait distributions

name = "pred_free"
try:
    data = np.load(file.format(name), allow_pickle = True)
    species_id = {**data}
    surv = data["surv"]
except FileNotFoundError:
    
    # generate species only for the first trophic level
    species_id = af.generate_species(locs = [0,-1.5,0])
    
    
    # run community assembly
    present, equi_all, surv = af.community_assembly(species_id)
    
    # save data 
    np.savez(file.format(name), **species_id, equi_all = equi_all,
            surv = surv, present = present)
    
fig, ax = plt.subplots(4,3, sharex = True, sharey = "row", figsize = (10,10))
fp.plot_richness(ax[0], surv, species_id)
fp.plot_traits(ax[1,:2], surv, species_id)
fp.plot_traits(ax[2,:2], surv, species_id)
fp.plot_invasion_prob(ax[-1,-1], surv, species_id)

fp.plot_mean_traits(ax[1,-1], surv, species_id, 0)
fp.plot_mean_traits(ax[2,-1], surv, species_id, 1)

for i in range(3):
    ax[-1,i].set_xlabel("Year")
    
for i, a in enumerate(ax.flatten()):
    a.set_title(ABC[i], loc = "left")

ax[0,0].set_ylabel("Species richness")
ax[1,0].set_ylabel("Species traits")
ax[2,0].set_ylabel("Invasion probability")
ax[0,0].set_title("Example community 1")
ax[0,1].set_title("Exmple community 2")
ax[0,-1].set_title("Average case")

fig.savefig(fig_file.format(name))