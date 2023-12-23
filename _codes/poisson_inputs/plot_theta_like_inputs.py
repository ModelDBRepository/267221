#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 18:57:36 2021

@author: spiros
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

plt.rcParams.update(
    {
        "text.usetex": False,
        "font.size": 18,
    }
)

mpl.rcParams["axes.spines.right"] = False
mpl.rcParams["axes.spines.top"] = False
mpl.rcParams["axes.linewidth"] = 2
mpl.rcParams["ytick.major.width"] = 2
mpl.rcParams["xtick.major.width"] = 2
mpl.rcParams["ytick.major.size"] = 7
mpl.rcParams["xtick.major.size"] = 7

N = 1000  # num neurons

spiketimes = []

dirname = 'rate50/run_1/'
fname = 'noise_'
for i in range(N):
    filename = f'{dirname}{fname}{i}.txt'
    sp = np.loadtxt(filename, dtype=float)
    spiketimes.append(sp)
    
    
plt.figure()
plt.xlabel('time [ms]')
plt.ylabel('pre-synaptic neuron id')
random_idx = np.random.choice(range(N), size=80, replace=False)
for i, idx in enumerate(random_idx):
    plt.scatter(spiketimes[idx], i*np.ones_like(spiketimes[idx]), color='k')
plt.yticks(range(0, 80, 10), [str(x) for x in range(0, 80, 10)])
plt.xlim([500, 800])

mlist = []
for j in spiketimes:
    mlist += (j-400).tolist()
    

plt.hist(np.sort(mlist), density=True, bins=250)
plt.xlabel('time [ms]')
plt.ylabel('probability density')
theta_freq = 8
theta_phase = 0
spike = np.arange(0, 4000, 0.1)
probability = 0.001*(np.sin(2.0*np.pi*theta_freq*spike/1000. + theta_phase) + 1.0)/2.0
# plt.plot(spike, probability/0.001, linewidth=3)
plt.xlabel('time [ms]')
plt.ylabel('normalized probability')
plt.xlim([0, 1500])