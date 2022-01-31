import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import linregress
import warnings
import pandas as pd

from scipy.stats import theilslopes, kendalltau

def plotting(x, y, ax, presences, color = None):

    keys = list(y.keys())
    metrics = pd.DataFrame(data = np.nan, index = keys,
                           columns = ["s", "intercept", "r",
                                      "p", "std", "corr",
                                   "n_spec", "n_time"])
    for key in keys:
        
        
        ind = np.isfinite(x[key]*y[key])
        if np.sum(ind) <=1:
            ax[0].plot(x[key], y[key], ',',
                   linestyle = '--', color = color)
            continue
        with warnings.catch_warnings(record = True):
            metrics.loc[key, metrics.columns[:5]] = linregress(
                                        x[key][ind],y[key][ind])
            metrics.loc[key, ["corr", "p"]] = kendalltau(x[key][ind],y[key][ind])
            
        
        p = metrics.loc[key, "p"]
        ls = '-' if p<0.05 else '--'
        if color == 'k':
            ls = ''
        ax[0].plot(x[key], y[key], ',',
                   linestyle = ls, color = color, label = key)
        with warnings.catch_warnings(record = True):
            metrics.loc[key, "corr"] = np.corrcoef(
                                                x[key][ind], y[key][ind])[0,1]
        metrics.loc[key, ["n_time", "n_spec"]] = presences[key].shape
        ax[1].plot(presences[key].shape[0], metrics.loc[key, "corr"],
                   "bo" if p<0.05 else 'r^')
        ax[2].plot(presences[key].shape[1], metrics.loc[key, "corr"],
                   "bo" if p<0.05 else 'r^')
        
    ax[1].axhline(0, color = "k", linestyle = "--")
    ax[2].axhline(0, color = "k", linestyle = "--")
    ax[3].axvline(0, color = "k", linestyle = "--")
    
    ax[3].hist(metrics["corr"], bins = np.linspace(-1,1,31), color = "grey")
    ax[3].hist(metrics["corr"][metrics["p"]<0.05],
               bins = np.linspace(-1,1,31), color = "blue")
    ax[3].set_xlabel("Correlation")
    ax[3].set_ylabel("Frequence")
    
    
    ax[1].text(.9,.9,"Datasets sig: {}".format(np.sum(metrics["p"]<0.05)),
                transform=ax[1].transAxes, ha = "right")
    ax[1].text(.9,.8,"Datasets tot: {}".format(len(metrics)),
                transform=ax[1].transAxes, ha = "right")
    ax[1].set_xlabel("Number of years")
    ax[2].set_xlabel("Number of species")
    ax[1].set_ylabel("Correlation")
    ax[2].set_ylabel("Correlation")
    return metrics
