#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 13:06:56 2021

@author: spiros
"""

import sys
import pickle, os
import numpy as np
import pandas as pd

from scipy.signal import find_peaks

trials = 500
replications = 10

# Variable to store the results
case = ["control", "vipcr_del", "vipcck_del", "vip_del", "cck_del"]

probs = {}

# repo = str(10)
repo = sys.argv[1]

for reps in range(1, replications + 1):
    num_peaks = {}
    isi_total = {}
    print(f"\nRepetition: {reps}")
    idx_cck = []
    for trial in range(1, trials + 1):
        print(trial)
        fname = f"data_synapses_poisson_rate_{repo}/circuitry_trial{trial}_rep{reps}_vip_del.pkl"
        if not os.path.isfile(fname):
            continue
        with open(fname, 'rb') as handle:
            data = pickle.load(handle)
        vipcr_flag = np.max(data['vip_del_voltage_vipcr'])
        vipcck_flag = np.max(data['vip_del_voltage_vipcck'])
        cck_flag = np.max(data['vip_del_voltage_cck'])
        if cck_flag > 0 and vipcr_flag > 0 and vipcck_flag > 0:
            idx_cck.append(trial)
    
    for c in case:
        print(f'Case analyzed: {c}')
        peak_nums = []
        isis = []
        # for reps in range(1, replications + 1):
        print(f"\nRepetition: {reps}")
        for trial in range(1, trials + 1):
            # print(trial)
    
            fname = f"data_synapses_poisson_rate_{repo}/circuitry_trial{trial}_rep{reps}_{c}.pkl"
            if not os.path.isfile(fname):
                continue
            with open(fname, 'rb') as handle:
                data = pickle.load(handle)
            
            fname = f"{c}_voltage_soma"
            pc_cells = data[fname]
            time = data[f"{c}_time"]
            
            peak_num = len(find_peaks(pc_cells, height=0)[0])/3.6
            
            if peak_num > 1:
                isi = np.min(np.diff(time[find_peaks(pc_cells, height=0)[0]]))
            else:
                isi = 0
    
            if trial in idx_cck:
                peak_nums.append(peak_num)
                isis.append(isi)
        num_peaks[c] = peak_nums
        isi_total[c] = isis
    
    if reps == 1:
        df = pd.DataFrame.from_dict(num_peaks)
        df_filtered = df.loc[(df!=0).any(axis=1)]
        
        df2 = pd.DataFrame.from_dict(isi_total)
        df_filtered2 = df2.loc[(df2!=0).any(axis=1)]
    else:
        df_d = pd.DataFrame.from_dict(num_peaks)
        df_filtered_d = df_d.loc[(df_d!=0).any(axis=1)]
        
        df2_d = pd.DataFrame.from_dict(isi_total)
        df_filtered2_d = df2_d.loc[(df2_d!=0).any(axis=1)]
        
        df = pd.concat([df, df_d])
        df_filtered = pd.concat([df_filtered, df_filtered_d])
        df2 = pd.concat([df2, df2_d])
        df_filtered2 = pd.concat([df_filtered2, df_filtered2_d])
            
    
    for c in case:
        if reps == 1:
            if c in probs.keys():
                if len(df) != 0:
                    probs[c].append(len(df[df[c]!=0])/len(df))
                else:
                    probs[c].append(0)
            else:
                if len(df) != 0:
                    probs[c] = [len(df[df[c]!=0])/len(df)]
                else:
                    probs[c] = [0]
        else:
            if c in probs.keys():
                if len(df_d) != 0:
                    probs[c].append(len(df_d[df_d[c]!=0])/len(df_d))
                else:
                    probs[c].append(0)
            else:
                if len(df_d) != 0:
                    probs[c] = [len(df_d[df_d[c]!=0])/len(df_d)]
                else:
                    probs[c] = [0]
            

final_results = {}
final_results['spikes'] = df
final_results['spikes_filtered'] = df_filtered
final_results['isi'] = df2
final_results['isi_filtered'] = df_filtered2
final_results['probabilities'] = probs

dirname = "data_analysis_poisson/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + f'rate_{repo}.pkl', 'wb') as handle:
    pickle.dump(final_results, handle, protocol=pickle.HIGHEST_PROTOCOL)

