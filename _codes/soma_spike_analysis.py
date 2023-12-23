#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 14 11:45:32 2021

@author: spiros
"""
import os, pickle
import numpy as np
import pandas as pd
# import seaborn as sns


repo = 150

dirname = "data_analysis_poisson/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + f'rate_{repo}.pkl', 'rb') as file:
    data = pickle.load(file)
    
cases = ['control', 'vip_del', 'vipcr_del', 'vipcck_del', 'cck_del']

df = data['spikes']

idx_ = np.where(np.array(df.index) == 0)[0]
df = df.reset_index()

probs = {c: [] for c in cases}

for i in range(len(idx_)):
    if i < len(idx_)-1:
        i1 = idx_[i]
        i2 = idx_[i+1]
        
        df_ = df.loc[i1:i2-1]
    else:
        df_ = df.loc[idx_[i]:]
    
    for c in cases:
        if len(df_) != 0:
            probs[c].append(len(df_[df_[c]!=0])/len(df_))
        else:
            probs[c].append(np.nan)
            

df_probs = pd.DataFrame.from_dict(probs)


df_sp = data['spikes_filtered']
df_sp = df_sp.reset_index()
rnd = np.random.choice(range(len(df_sp)), size=200, replace=False)

df_sp = df_sp.loc[rnd]
df_sp = df_sp.reset_index()

df_isi = data['isi_filtered']
df_isi = df_isi.reset_index()
df_isi = df_isi.loc[~(df_isi==0).any(axis=1)]
df_isi = df_isi.reset_index()

rnd = np.random.choice(range(len(df_isi)), size=200, replace=False)

df_isi = df_isi.loc[rnd]
