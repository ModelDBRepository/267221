#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  3 09:53:01 2022

@author: spiros
"""

import pickle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from opt import mystyle

# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())


dirname = "data_analysis_poisson/"

trials = 10
repo = 20
    
# Variable to store the results
case = ["control",  "vip_del", "vipcr_del", "vipcck_del", "cck_del"]
cols = ["10Hz", "20Hz", "40Hz", "100Hz", "case"]
df_all = pd.DataFrame(columns=cols, index=range(len(case)*trials))

for repo in [10, 20, 40, 100]:
    counter = 0
    for c in case:
        with open(dirname + f'rate_{repo}.pkl', 'rb') as handle:
            DATA = pickle.load(handle)
        probs = DATA['probabilities'][c]
        for i in range(trials):
            df_all.loc[counter][f'{repo}Hz'] = probs[i]
            df_all.loc[counter]['case'] = c
            counter += 1

# Variable to store the results
case = ["control",  "vip_del", "vipcr_del", "vipcck_del", "cck_del"]
cols = ["probability", "case", "rate"]
repos = [10, 20, 40, 100]
df_all2 = pd.DataFrame(columns=cols, index=range(len(case)*trials*len(repos)))

counter = 0
for repo in repos:
    for c in case:
        with open(dirname + f'rate_{repo}.pkl', 'rb') as handle:
            DATA = pickle.load(handle)
        probs = DATA['probabilities'][c]
        for i in range(trials):
            df_all2.loc[counter]['probability'] = probs[i]
            df_all2.loc[counter]['case'] = c
            df_all2.loc[counter]['rate'] = repo
            counter += 1

plt.figure()
sns.boxplot(data=df_all2, x="rate", y="probability", hue="case")
plt.xlabel('input rate (Hz)')
plt.ylabel('somatic AP prob.')
plt.show()