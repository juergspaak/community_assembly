import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from scipy.stats import pearsonr, linregress

from functions_for_plotting import biotime_colors, keys

"""
# load simulations
path = "C:/Users/Juerg Spaak/Documents/Science backup/TND/"
# load prey only data

path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

study_ids = presences["study_ids"]
# only work with long datasets
study_ids = study_ids[[presences[study_id].shape[0] >= 30
                      for study_id in study_ids]]

meta_data = "BioTIMEMetadata_24_06_2021.csv"
meta_data = pd.read_csv(path + meta_data, usecols=np.arange(27),
                            encoding = 'unicode_escape', index_col="STUDY_ID")


meta_data = meta_data.loc[[int(std_id) for std_id in study_ids]]
meta_data["identity"] = [biotime_colors[key][0]
                         if key in biotime_colors.keys() else np.nan
                         for key in meta_data.TAXA]
meta_data["corr"] = np.nan # store correlation between invasion and extinction
meta_data["p_value"] = np.nan # p_value for 
meta_data["p_value_high_rich"] = np.nan # p_value for correlation for high richness
meta_data["rand_corr"] = np.nan # correlation of perturbed data

meta_data["include"] = True

###############################################################################
# compute correlation between invasion and extinction

def inv_ext_correlation(pres):
    # compute correlation between invasion and extinction
    
    # total species richness
    richness = np.sum(pres, axis = 1)
    
    # not present last year, but present this year
    inv = np.sum(pres[1:] & ~pres[:-1], axis = 1)
    # how many of the present species are invaders?
    inv = inv/richness[1:]
    
    # present last year, but not present this year
    ext = np.sum(~pres[1:] & pres[:-1], axis = 1)
    # how many of the present species are invaders?
    ext = ext/richness[:-1]
    
    #d = 1
    # return correlation, extinctions cause invasions
    #return np.corrcoef(inv[d:], ext[:-d])[0,1]
    # invasions cause extinctions
    #return np.corrcoef(inv[:-1], ext[1:])[0,1]
    # 
    return np.corrcoef(inv, ext)[0,1]


# compute correlation between invasion and extinction
itera = 100
p_distribution = np.empty((2, len(meta_data), itera))
for j, study_id in enumerate(meta_data.index):
    pres = presences[str(study_id)] == 1 # conversion to boolean
    richness = np.sum(pres, axis = 1)
    
    # for exclusion of communities with to strong fluctuations
    if min(richness)/max(richness)<0.25:
        meta_data.loc[study_id, "include"] = False
    
    # for exclusion of timepoints with very low species richness
    ind_high_rich = richness>=5
    
    meta_data.loc[study_id, "corr"] = inv_ext_correlation(pres)
    meta_data.loc[study_id, "corr_high_rich"] = inv_ext_correlation(pres[ind_high_rich])
    
    
    for i in range(itera):
        order = np.argsort(np.random.rand(len(pres)))
        p_distribution[0, j, i] = inv_ext_correlation(pres[order])
    
    pres_high_rich = pres[ind_high_rich]
    for i in range(itera):
        order = np.argsort(np.random.rand(len(pres_high_rich)))
        p_distribution[1, j, i] = inv_ext_correlation(pres_high_rich[order])
        
    meta_data.loc[study_id, "p_value"] = np.sum(meta_data.loc[study_id, "corr"]>p_distribution[0, j])/itera
    meta_data.loc[study_id, "rand_corr"] = np.mean(p_distribution[0, j])
    meta_data.loc[study_id, "p_value_high_rich"] = np.sum(meta_data.loc[study_id, "corr_high_rich"]>p_distribution[1, j])/itera
    if np.isnan(meta_data.loc[study_id, "p_value"]):
        raise
        
meta_data["p_value"] = 1 - meta_data["p_value"]
meta_data["p_value_high_rich"] = 1 - meta_data["p_value_high_rich"]

    
###############################################################################
# correlation between invasion and extinction over time

meta_data["p_inv"] = np.nan
meta_data["p_ext"] = np.nan
meta_data["R2_inv"] = np.nan
meta_data["R2_ext"] = np.nan
meta_data["p_rich"] = np.nan
meta_data["R2_rich"] = np.nan
meta_data["rich_max"] = np.nan
meta_data["rich_min"] = np.nan

for study_id in study_ids:
    # load dataset, shape (years, species)
    pres = presences[study_id] ==1
    
    years = presences["year_{}".format(study_id)][:-1]
    
    # proportion of invaders
    invasion = np.sum(pres[1:] & ~pres[:-1], axis = 1)/np.sum(pres[1:], axis = 1)
    
    meta_data.loc[int(study_id), ["R2_inv", "p_inv"]] = pearsonr(years, invasion)

    # proportion of extinctions
    extinction = np.sum(pres[:-1] & ~pres[1:], axis = 1)/np.sum(pres[:-1], axis = 1)
    meta_data.loc[int(study_id), ["R2_ext", "p_ext"]] = pearsonr(years, extinction)
    
    # species richness
    rich = np.sum(pres, axis = 1)
    meta_data.loc[int(study_id), ["R2_rich", "p_rich"]] = pearsonr(presences["year_{}".format(study_id)],
                                                                   rich)
    meta_data.loc[int(study_id), "rich_max"] = np.amax(rich)
    meta_data.loc[int(study_id), "rich_min"] = np.amin(rich)
    

meta_data["R2_ext"] = meta_data["R2_ext"]**2
meta_data["R2_inv"] = meta_data["R2_inv"]**2
meta_data["R2_rich"] = meta_data["R2_rich"]**2
#"""
    
meta_data = meta_data
presences = presences
# condition to be a useful candidate
# not driven by climate change, i.e. no correlation between time and
# invasion, extinctions or richness
p_tresh = 0.05
R2_tresh = 0.1
ind_climate = (((meta_data.p_inv>p_tresh) | (meta_data.R2_inv < R2_tresh))
               & ((meta_data.p_ext>p_tresh) | (meta_data.R2_ext < R2_tresh))
               & ((meta_data.p_rich>p_tresh) | (meta_data.R2_rich < R2_tresh)))

# strong correlation between invasion and extinctions
ind = ind_climate & (meta_data.p_value_high_rich<0.05)

# not too strong fluctuations in richness
#ind = ind & (meta_data.rich_max < 4*meta_data.rich_min)

print(sum(ind))


def plot_pres(study_id):
    pres = presences[str(study_id)] == 1
    pres = pres[:, np.any(~pres, axis = 0)]
    if np.any(np.all(~pres, axis = 1)):
        raise
    pres = np.append(np.full((1, pres.shape[1]), False), pres, axis = 0)
    pres = np.append(pres, np.full((1, pres.shape[1]), False), axis = 0)
    years = presences["year_{}".format(study_id)]
    trait_starts, starts = np.where((pres[1:] & ~pres[:-1]).T)
    trait_ends, ends = np.where((pres[:-1] & ~pres[1:]).T)
    
    for i in range(len(starts)):
        plt.plot([years[starts[i]], years[ends[i]-1]],
                [trait_starts[i], trait_ends[i]], 'b-', alpha = 0.5)
        if trait_starts[i] != trait_ends[i]:
            raise
        if starts[i] > ends[i]:
            raise

    
for study_id in study_ids:
    if meta_data.loc[int(study_id), "p_value"]>0.05:
        continue
    if meta_data.loc[int(study_id), "rich_max"] > 100:
        continue
    plt.figure()
    plot_pres(study_id)
    plt.title(study_id)
    plt.show()
    