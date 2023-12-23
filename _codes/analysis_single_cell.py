#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 12 21:12:13 2021

@author: spiros
"""

import numpy as np
import pandas as pd
import pickle, os

dirname = "data_psp_single/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

conditions = ["Control", "-Inh", "-Somatic Inh", "-Prox Dend Inh", "-Dist Dend Inh"]
probs_dspike = {c: [] for c in conditions}

for rep in range(1, 51):
    print(f"Replication: {rep}")
    with open(dirname + f'single_cell_final_new_replication_{rep}.pkl', 'rb') as f:
        DATA = pickle.load(f)
    
    df_soma = DATA["df_soma"]
    df_dend = DATA["df_dend"]
    
    cols = ["peak", "time_rise", "time_decay", "dvdt", "latency",
            "thalf", "time_peak", "case"]

    trials = 1000
    
    # Remove all dSpikes and Somatic spike data from all conditions
    conditions = ["Control", "-Inh", "-Somatic Inh", "-Prox Dend Inh", "-Dist Dend Inh"]
    df_dend_sep = {}
    df_soma_sep = {}
    for c in conditions:
        df_ = df_dend[df_dend["case"]==c]
        df_ = df_.reset_index()
        df_dend_sep[c] = df_
    
        df_ = df_soma[df_soma["case"]==c]
        df_ = df_.reset_index()
        df_soma_sep[c] = df_
    
    # Indices of spike observation and removal from all conditions
    idx_ = []
    for c in conditions:
        df_ = df_dend_sep[c]
        idx_ += df_.index[df_["dvdt"] > 10].tolist()
        
        df_ = df_soma_sep[c]
        idx_ += df_.index[df_["peak"] > 50].tolist()
    
    idx_ = list(set(idx_))
    keep_idx = [i for i in range(trials) if i not in idx_]
    for c in conditions:
        df_soma_sep[c] = df_soma_sep[c].iloc[keep_idx].reset_index()
        df_dend_sep[c] = df_dend_sep[c].iloc[keep_idx].reset_index()
    
    idx_psp_all = df_dend.index[df_dend["dvdt"] <= 10].tolist()
    idx_dspike_all = df_dend.index[df_dend["dvdt"] > 10].tolist()
    
    
    df_dend3 = pd.DataFrame(columns = cols)
    df_soma3 = pd.DataFrame(columns = cols)
    for c in conditions:
        df_dend3 = pd.concat((df_dend3, df_dend_sep[c]), ignore_index=True)
        df_soma3 = pd.concat((df_soma3, df_soma_sep[c]), ignore_index=True)
    
    df_dend3.drop("level_0", inplace=True, axis=1)
    df_dend3.drop("index", inplace=True, axis=1)
    df_soma3.drop("level_0", inplace=True, axis=1)
    df_soma3.drop("index", inplace=True, axis=1)
    
    df_dend4 = df_dend.iloc[idx_dspike_all]
    
    df_dend5 = df_dend4
    df_dend5["case2"] = "dspike"
    df_dend6 = df_dend3
    df_dend6["case2"] = "dPSP"
    df_dend7 = pd.concat([df_dend6, df_dend5])
    
    
    for m in cols:
        df_dummy_s = pd.DataFrame(columns=conditions)
        df_dummy_d = pd.DataFrame(columns=conditions)
        if m != "case":
    
            for c in conditions:
                df_ = df_soma3[[m, "case"]]
                df_2 = df_[df_["case"]==c]
                df_dummy_s[c] = list(np.abs(df_2[m]))
                
                df_ = df_dend3[[m, "case"]]
                df_2 = df_[df_["case"]==c]
                df_dummy_d[c] = list(np.abs(df_2[m]))
                
            
            df_dummy_s.to_csv(dirname+f"{m}_soma.csv", index=False)
            df_dummy_d.to_csv(dirname+f"{m}_dend.csv", index=False)
            
    
    conditions = ["Control", "-Inh", "-Somatic Inh", "-Prox Dend Inh", "-Dist Dend Inh"]
    df_dend_sep = {}
    for c in conditions:
        df_ = df_dend[df_dend["case"]==c]
        df_ = df_.reset_index()
        df_dend_sep[c] = df_
            
    idx_ = []
    for c in conditions:
        df_ = df_dend_sep[c]
        idx_.append(df_.index[df_["dvdt"] > 10].tolist())
    
    # Keep the common elements across all lists/conditions
    keep_idx = list(set.intersection(*[set(list) for list in idx_]))
    
    df_dend_spikes_sep = {}
    for c in conditions:
        df_dend_spikes_sep[c] = df_dend_sep[c].iloc[keep_idx].reset_index()
    
    df_dend_spikes = pd.DataFrame(columns = cols)
    for c in conditions:
        df_dend_spikes = pd.concat((df_dend_spikes, df_dend_spikes_sep[c]), ignore_index=True)
    
    for m in cols:
        df_dummy_d = pd.DataFrame(columns=conditions)
        if m != "case":
    
            for c in conditions:
    
                df_ = df_dend_spikes[[m, "case"]]
                df_2 = df_[df_["case"]==c]
                df_dummy_d[c] = list(np.abs(df_2[m]))
    
            df_dummy_d.to_csv(dirname+f"{m}_dend_spikes_rep_{rep}.csv", index=False)
            
            
    for c in conditions:
        df_ = df_dend[df_dend["case"]==c]
        probs_dspike[c].append(len(df_[df_["dvdt"]>10])/len(df_)*100)
        
        
df_dspikes_ = pd.DataFrame.from_dict(probs_dspike)
df_dspikes_.to_csv(dirname+f"{m}_dend_spikes_probabilities.csv", index=False)