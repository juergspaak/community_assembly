import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import warnings

from scipy.stats import spearmanr, pearsonr


path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
presences = np.load(path + "biotime_converted.npz")

study_ids = presences["study_ids"]
length = np.array([len(presences["year_{}".format(study_id)])
                   for study_id in study_ids])

long_studies = study_ids[length>=35]

fig, ax = plt.subplots(5,5, figsize = (12,12))
for i in range(5):
    ax[i,0].set_ylabel("Frequency")
    ax[-1,i].set_xlabel("Persistence length")
ax = ax.flatten()

p_all = []

for i, s_id in enumerate(long_studies):
    
    # compute persistence length
    persistence = presences[s_id].sum(axis = 0)
    p_all.extend(persistence/len(presences[s_id]))
    if np.amax(persistence>100):
        raise
    ax[i].hist(persistence, bins = np.amax(persistence))
    
    
ax[-1].hist(p_all, bins = 50, color = "red")
fig.tight_layout()

