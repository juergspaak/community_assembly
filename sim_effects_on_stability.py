import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

import assembly_functions as ap
import functions_for_plotting as fp
"""
inv_basal = ap.species_id_invader_basal
inv_pred = ap.species_id_invader_pred

itera = 200
cols = [key + level for key in ["n_mean", "inv_t0", "inv_t1", "inv_t2"]
        for level in ["_pred", "_prey"]] + ["omega", "change_comp"]

data = pd.DataFrame(np.nan, index = range(itera), columns = cols)
# length of run
years = 2000
data["years"] = years
# niche width of prey and predators
data["sig_pred"] = np.random.uniform(0.5, 2, itera)
data["sig_prey"] = np.random.uniform(0.5, 2, itera)
# niche breadth
data["omega"] = np.random.uniform(3, 5, itera)
# maximum interaction strength
data["alpha_pred"] = np.random.uniform(0.5, 2, itera)
data["alpha_prey"] = np.random.uniform(0.5, 2, itera)
data["alpha_pred"] = data["alpha_prey"]
# mortality rate
data["m_pred"] = np.random.uniform(0.01, 0.2, itera)
data["m_prey"] = np.random.uniform(0.01, 0.2, itera)
data["m_pred"] = data["m_prey"]
# maximum growth rate
data["mu_prey"] = np.random.uniform(0.5,2, itera)


for i, row in data.iterrows():
    # generate species according to generation
    species_id = ap.generate_species(1, int(row["years"]), level = [1,1],
                omega = row["omega"], sigs = row[["sig_prey", "sig_pred"]],
                alpha_max = row[["alpha_prey", "alpha_pred"]],
                mu_max = row["mu_prey"], ms = row[["m_prey", "m_pred"]])

    # compute community dynamics
    present, equi_all, surv = ap.community_assembly(species_id, pr = False)
    lv = species_id["level"][:,np.newaxis]
    # save average species richness
    data.loc[i, "n_mean_prey"] = np.mean(np.sum(surv & (lv == 0),
                                                axis = -1)[:, years//2:])
    data.loc[i, "n_mean_pred"] = np.mean(np.sum(surv & (lv == 1),
                                                axis = -1)[:, years//2:])
    # adjust invading species
    inv_basal["loc"] = np.linspace(*np.nanpercentile(species_id["loc"], [5, 95]),
                                   len(inv_basal["loc"]))
    inv_basal["sig"][:] = row["sig_prey"]
    inv_basal["alpha_max"][:] = row["alpha_prey"]
    inv_basal["m"][:] = row["m_prey"]
    # invading predators
    inv_pred["loc"] = inv_basal["loc"]
    inv_pred["sig"][:] = row["sig_pred"]
    inv_pred["alpha_max"][:] = row["alpha_pred"]
    inv_pred["m"][:] = row["m_pred"]

    # save invasion probabilities
    for j, timepoint in enumerate([years//2, 3*years//4, years]):
        inv, loc = ap.invasion_success(species_id, 0, equi_all[0, timepoint],
                                   inv_basal)
        data.loc[i, "inv_t{}_basal".format(j)] = np.mean(inv>0)
        
        inv, loc = ap.invasion_success(species_id, 0, equi_all[0, timepoint],
                                   inv_pred)
        data.loc[i, "inv_t{}_pred".format(j)] = np.mean(inv>0)
        
    # compute changes in composition
    data.loc[i, "change_comp"] = np.sum(present[:,-200] & present[:,-1])/np.sum(present[:,-200] | present[:,-1])


    if data.loc[i, "change_comp"]>0.3:
        plt.figure()
        fp.plot_traits([plt.gca()], surv, species_id)
        plt.xlim([1800,None])
        plt.xlabel(data.loc[i, "change_comp"])
        plt.title(str(np.round(data.loc[i, ["sig_prey", "sig_pred"]].values, 3)))
        plt.show()
    print(i, data.loc[i, ["n_mean_prey", "n_mean_pred"]].values)
data = data.loc[:i]

#"""
data_cop = data.copy()
data["frac_sig"] = data["sig_pred"]/data["sig_prey"]
data["diff_sig"] = data["sig_pred"]-data["sig_prey"]

plt.scatter(data["sig_pred"]-data["sig_prey"], data["change_comp"],
            c = data["n_mean_prey"].values)
plt.colorbar()
#plt.scatter(data["sig_pred"]-data["sig_prey"], data.inv_t2_pred)

responses = ["change_comp", "n_mean_pred", "n_mean_prey", "inv_t2_basal"]
predictors = ["omega", "sig_prey", "sig_pred", "mu_prey", "alpha_pred",
              "alpha_prey", "change_comp", "frac_sig", "diff_sig"]

fig, ax = plt.subplots(len(responses), len(predictors), figsize = (13,13),
                       sharex = "col", sharey = "row")

for i, predictor in enumerate(predictors):
    ax[-1,i].set_xlabel(predictor)
    for j, response in enumerate(responses):
        ax[j,0].set_ylabel(response)
        cmap = ax[j,i].scatter(data[predictor], data[response], c = data["n_mean_prey"], vmax = 10)
        
fig.colorbar(cmap)

