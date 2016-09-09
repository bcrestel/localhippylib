'''
Created on Sep 9, 2016

@author: uvilla
'''
import numpy as np
import matplotlib.pyplot as plt

def acorr(mean_free_samples, lag, norm = 1):
    #http://stackoverflow.com/questions/14297012/estimate-autocorrelation-using-python    
    return (mean_free_samples[:mean_free_samples.size-lag]*mean_free_samples[lag:]).ravel().mean() / norm

def acorr_vs_lag(samples, max_lag = 500):
    mean = samples.mean()
    mean_free_samples = samples - mean
    
    norm = acorr(mean_free_samples, 0)
    
    lags = np.arange(0,max_lag+1)
    acorrs = np.ones(max_lag+1)
    
    for lag in lags[1:]:
        acorrs[lag] = acorr(mean_free_samples, lag, norm)
        
    return lags, acorrs

names = ['gpCN', 'pCN', 'inf-MALA']
style = ['-b', '-r', '-g']

q = [np.loadtxt(name+'.txt') for name in names]
plt.figure()
[plt.plot(q[i], style[i], label=names[i]) for i in range(len(names)) ]
plt.legend()
plt.figure()
for i in range(len(names)):
    lags, acorrs = acorr_vs_lag(q[i],1000)
    plt.plot(lags, acorrs, style[i], label=names[i]) 
plt.legend()

plt.show()
