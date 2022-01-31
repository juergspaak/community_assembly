import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress

"""
path = "C:/Users/Juerg Spaak/Documents/Science backup/P14_community_assembly/"
data = pd.read_csv(path + "TRY_plant_dataset.txt", sep = "\t",
                   encoding = 'unicode_escape')
data.loc[np.isnan(data.TraitID), "TraitID"] = 0
"""

traits = sorted(list(set(data.TraitID)))[1:] # remove 0 trait, corresponds to nan

for trait in traits:
    data_sub = data[data.TraitID == trait]
    data_sub.to_csv(path + "TRY_subdata_{}.csv".format(int(trait)), index = False)
    print(trait, data_sub.shape)
    