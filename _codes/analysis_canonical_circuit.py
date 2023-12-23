#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 18 00:27:48 2021

@author: spiros
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pickle
import os
import pandas as pd
from opt import mystyle

# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())


trials = 200
replications = 50

# Variable to store the results
case = ["control", "vipcr_del", "vipcck_del", "vip_del", "cck_del"]
cols = ["num_spikes", "peak_dpsp", "dvdt_dpsp", "peak_dspike", "dvdt_dspike", "psp_soma", "case"]
df = pd.DataFrame(columns=cols, index=range(len(case)*replications))
counter = 0

repo = str(3) + '_noise'
for reps in range(1, replications + 1):
    print(f"\nRepetition: {reps}")
    idx_cck = []
    dSpike_ctr = []
    dpsp_ctr = []
    for trial in range(1, trials + 1):
        print(f'trial...{trial}')
        c = "vip_del"
        fname = f"data_synapses_new{repo}/circuitry_trial{trial}_rep{reps}_{c}.pkl"
        with open(fname, 'rb') as handle:
            data = pickle.load(handle)
        
        vipcr_flag = np.max(data['results'][f'{c}_voltage_vipcr'])
        vipcck_flag = np.max(data['results'][f'{c}_voltage_vipcck'])
        cck_flag = np.max(data['results'][f'{c}_voltage_cck'])
        if cck_flag > 100 and vipcr_flag > 100 and vipcck_flag > 100:
            idx_cck.append(trial)

        fname2 = f"data_synapses_new{repo}/circuitry_trial{trial}_rep{reps}_control.pkl"
        with open(fname2, 'rb') as handle:
            data = pickle.load(handle)
            
        if np.array(data["dend"]["dvdt"]).item() > 10:
            dSpike_ctr.append(trial)
        else:
            dpsp_ctr.append(trial)
            
    cck_activity = []
    for c in case:
        counts = 0
        for trial in range(1, trials + 1):
            print(f"Trial: {trial}")

            fname = f"data_synapses_new{repo}/circuitry_trial{trial}_rep{reps}_{c}.pkl"

            with open(fname, 'rb') as handle:
                data = pickle.load(handle)
            
            if trial in dpsp_ctr:
                if cck_flag > 100:
                    counts += 1
                is_local = "df_soma_psp" in locals() and "df_dend_psp" in locals()
                if is_local:
                    df_soma_psp = pd.concat([df_soma_psp, data['soma']], ignore_index=True)
                    df_dend_psp = pd.concat([df_dend_psp, data['dend']], ignore_index=True)
                else:
                    df_soma_psp = data['soma']
                    df_dend_psp = data['dend']

            if trial in dSpike_ctr:
                if cck_flag > 100:
                    counts += 1
                is_local = "df_soma_dSpike" in locals() and "df_dend_dSpike" in locals()
                if is_local:
                    df_soma_dSpike = pd.concat([df_soma_dSpike, data['soma']], ignore_index=True)
                    df_dend_dSpike = pd.concat([df_dend_dSpike, data['dend']], ignore_index=True)
                else:
                    df_soma_dSpike = data['soma']
                    df_dend_dSpike = data['dend']

            if trial in idx_cck:
                if cck_flag > 100:
                    counts += 1
                is_local = "df_soma" in locals() and "df_dend" in locals()
                if is_local:
                    df_soma = pd.concat([df_soma, data['soma']], ignore_index=True)
                    df_dend = pd.concat([df_dend, data['dend']], ignore_index=True)
                else:
                    df_soma = data['soma']
                    df_dend = data['dend']
                    
        cck_activity.append(counts)
    # remove data with somatic spike
    df_dend2 = df_dend[df_soma["peak"] <= 60]
    df_soma2 = df_soma[df_soma["peak"] <= 60]
    
    
    ## Filter using the Control condition
    # idx_psp = df_dend.index[df_dend["dvdt"] <= 10.0].tolist()
    # idx_dspike = df_dend.index[df_dend["dvdt"] > 10.0].tolist() 
    
    df_dend3 = df_dend2[df_dend2["dvdt"] <= 10.0]
    df_soma3 = df_soma2[df_dend2["dvdt"] <= 10.0]
    
    df_dend4 = df_dend2[df_dend2["dvdt"] > 10.0]
    
    # plt.figure()
    # plt.subplot(2, 2, 1)
    # sns.boxplot(data=df_dend3, y="peak", x="case")
    # plt.subplot(2, 2, 2)
    # sns.boxplot(data=df_dend3, y="time_rise", x="case")
    # plt.subplot(2, 2, 3)
    # sns.boxplot(data=df_dend3, y="thalf", x="case")
    # plt.subplot(2, 2, 4)
    # sns.boxplot(data=df_dend3, y="dvdt", x="case")
    
    # plt.figure()
    # plt.subplot(2, 2, 1)
    # sns.boxplot(data=df_soma3, y="peak", x="case")
    # plt.subplot(2, 2, 2)
    # sns.boxplot(data=df_soma3, y="time_rise", x="case")
    # plt.subplot(2, 2, 3)
    # sns.boxplot(data=df_soma3, y="thalf", x="case")
    # plt.subplot(2, 2, 4)
    # sns.boxplot(data=df_soma3, y="dvdt", x="case")
    
    
    df_dend5 = df_dend4
    df_dend5["case2"] = "dspike"
    df_dend6 = df_dend3
    df_dend6["case2"] = "dPSP"
    df_dend7 = pd.concat([df_dend6, df_dend5])
    
    # plt.figure()
    # plt.subplot(2, 2, 1)
    # sns.boxplot(data=df_dend7[df_dend7["case2"]=="dspike"], y="peak", x="case")
    # plt.subplot(2, 2, 2)
    # sns.boxplot(data=df_dend7[df_dend7["case2"]=="dspike"], y="time_rise", x="case")
    # plt.subplot(2, 2, 3)
    # sns.boxplot(data=df_dend7[df_dend7["case2"]=="dspike"], y="thalf", x="case")
    # plt.subplot(2, 2, 4)
    # sns.boxplot(data=df_dend7[df_dend7["case2"]=="dspike"], y="dvdt", x="case")

    df_control1 = df_dend[(df_dend["case"]=="control") & (df_dend["dvdt"]>10)]
    df_vipdel1 = df_dend[(df_dend["case"]=="vip_del") & (df_dend["dvdt"]>10)]
    df_vipcck1 = df_dend[(df_dend["case"]=="vipcck_del") & (df_dend["dvdt"]>10)]
    df_vipcr1 = df_dend[(df_dend["case"]=="vipcr_del") & (df_dend["dvdt"]>10)]
    df_cck1 = df_dend[(df_dend["case"]=="cck_del") & (df_dend["dvdt"]>10)]

    df_control2 = df_dend[(df_dend["case"]=="control") & (df_dend["dvdt"]<=10)]
    df_vipdel2= df_dend[(df_dend["case"]=="vip_del") & (df_dend["dvdt"]<=10)]
    df_vipcck2 = df_dend[(df_dend["case"]=="vipcck_del") & (df_dend["dvdt"]<=10)]
    df_vipcr2 = df_dend[(df_dend["case"]=="vipcr_del") & (df_dend["dvdt"]<=10)]
    df_cck2 = df_dend[(df_dend["case"]=="cck_del") & (df_dend["dvdt"]<=10)]
    
    print(f'Control: {len(df_control1)/(len(df_control1)+len(df_control2))}')
    print(f'VIP del: {len(df_vipdel1)/(len(df_vipdel1)+len(df_vipdel2))}')
    print(f'VIP/CCK del: {len(df_vipcck1)/(len(df_vipcck1)+len(df_vipcck2))}')
    print(f'VIP/CR del: {len(df_vipcr1)/(len(df_vipcr1)+len(df_vipcr2))}')
    print(f'CCK del: {len(df_cck1)/(len(df_cck1)+len(df_cck2))}')
    
    case2 = ["control", "vipcr_del", "vipcck_del", "vip_del", "cck_del"]
    for c in case2:
        
        dfxx1 = df_dend[(df_dend["case"]==c) & (df_dend["dvdt"]>10)]
        dfxx2 = df_dend[(df_dend["case"]==c) & (df_dend["dvdt"]<=10)]
        dfxx3 = df_soma[df_soma["case"]==c]
        nums = len(dfxx1)/(len(dfxx1) + len(dfxx2))
        
        df.loc[counter].num_spikes = nums
        df.loc[counter].peak_dspike = np.mean(dfxx2["peak"])
        df.loc[counter].dvdt_dspike = np.mean(dfxx2["dvdt"])
        df.loc[counter].peak_dpsp = np.mean(dfxx1["peak"])
        df.loc[counter].dvdt_dpsp = np.mean(dfxx1["dvdt"])
        df.loc[counter].psp_soma = np.mean(dfxx3["peak"])
        df.loc[counter].case = c
        counter += 1


plt.figure()
sns.barplot(data=df, y="num_spikes", x="case", ci="sd", capsize=.2,
            palette="Set2")
plt.ylim([0, 0.5])
plt.ylabel("probability")
plt.xlabel("")
plt.title("LEC induced dSpikes")
plt.show()

# Remove all dSpikes and Somatic spike data from all conditions
conditions = ["control", "vip_del", "vipcck_del", "vipcr_del"]
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
df_ = df_dend_sep["control"]
idx_ += df_.index[df_["dvdt"] > 10].tolist()

df_ = df_soma_sep[c]
idx_ += df_.index[df_["peak"] > 50].tolist()

idx_ = list(set(idx_))
keep_idx = [i for i in range(len(df_soma_sep['control'])) if i not in idx_]
for c in conditions:
    df_soma_sep[c] = df_soma_sep[c].iloc[keep_idx].reset_index()
    df_dend_sep[c] = df_dend_sep[c].iloc[keep_idx].reset_index()

idx_psp_all = df_dend.index[df_dend["dvdt"] <= 10].tolist()
idx_dspike_all = df_dend.index[df_dend["dvdt"] > 10].tolist()


cols = ["peak", "time_rise", "time_decay", "dvdt", "latency",
        "thalf", "time_peak", "case"]
df_dend3 = pd.DataFrame(columns=cols)
df_soma3 = pd.DataFrame(columns=cols)
for c in conditions:
    df_dend3 = pd.concat((df_dend3, df_dend_sep[c]), ignore_index=True)
    df_soma3 = pd.concat((df_soma3, df_soma_sep[c]), ignore_index=True)

df_dend3.drop("level_0", inplace=True, axis=1)
df_dend3.drop("index", inplace=True, axis=1)
df_soma3.drop("level_0", inplace=True, axis=1)
df_soma3.drop("index", inplace=True, axis=1)

df_dend4 = df_dend.iloc[idx_dspike_all]


conditions = ["control", "vip_del", "vipcr_del", "vipcck_del", "cck_del"]
dirname = "canonical_circuit_noise/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

cols = ["peak", "time_rise", "time_decay", "dvdt", "latency",
        "thalf", "case"]
for m in cols:
    df_dummy_dspike = pd.DataFrame(columns=conditions)
    df_dummy_dpsp = pd.DataFrame(columns=conditions)
    df_dummy_soma = pd.DataFrame(columns=conditions)
    if m != "case":
        for c in conditions:
            df_ = df_dend_dSpike[[m, "case"]]
            df_2 = df_[df_["case"]==c]
            df_dummy_dspike[c] = list(df_2[m])

            df_ = df_dend_psp[[m, "case"]]
            df_2 = df_[df_["case"]==c]
            df_dummy_dpsp[c] = list(df_2[m])

            df_ = df_soma_psp[[m, "case"]]
            df_2 = df_[df_["case"]==c]
            df_dummy_soma[c] = list(df_2[m])

        df_dummy_dspike.to_csv(dirname+f"{m}_dend_spikes.csv", index=False)
        df_dummy_dpsp.to_csv(dirname+f"{m}_dend_psp.csv", index=False)
        df_dummy_soma.to_csv(dirname+f"{m}_soma_psp.csv", index=False)


# Save in csv the dSpike probabilities per case
df_dummy_probs = pd.DataFrame(columns=conditions)

for c in conditions:

    df_ = df[["num_spikes", "case"]]
    df_2 = df_[df_["case"]==c]
    df_dummy_probs[c] = list(df_2["num_spikes"])


df_dummy_probs.to_csv(dirname+"num_spikes_dend_probs.csv", index=False)


# Load the csv files and remove all zero rows from all conditions
# Find the indices to keep

idx_zeros_dends_spikes = []
idx_zeros_dends = []
idx_zeros_soma = []

for m in cols:
    if m != "case":
        df_zeros = pd.read_csv(dirname+f"{m}_dend_spikes.csv")
        # indices to keep ('false' in the previous list).
        zeros = df_zeros[(df_zeros.values == 0).any(axis=1)].index.values.tolist()
        # append the new indices from case m to be removed from all metrics.
        for j in zeros:
            if j not in idx_zeros_dends_spikes:
                idx_zeros_dends_spikes.append(j)

        df_zeros = pd.read_csv(dirname+f"{m}_dend_psp.csv")
        # indices to keep ('false' in the previous list).
        zeros = df_zeros[(df_zeros.values == 0).any(axis=1)].index.values.tolist()
        # append the new indices from case m to be removed from all metrics.
        for j in zeros:
            if j not in idx_zeros_dends:
                idx_zeros_dends.append(j)

        df_zeros = pd.read_csv(dirname+f"{m}_soma_psp.csv")
        # indices to keep ('false' in the previous list).
        zeros = df_zeros[(df_zeros.values == 0).any(axis=1)].index.values.tolist()
        # append the new indices from case m to be removed from all metrics.
        for j in zeros:
            if j not in idx_zeros_soma:
                idx_zeros_soma.append(j)


idx_zeros_dends_spikes = sorted(idx_zeros_dends_spikes)
idx_zeros_dends = sorted(idx_zeros_dends)
idx_zeros_soma = sorted(idx_zeros_soma)

for m in cols:
    if m != "case":
        df_keep = pd.read_csv(dirname+f"{m}_dend_spikes.csv")
        ix = [i for i in df_keep.index if i not in idx_zeros_dends_spikes]
        df_keep1 = df_keep.iloc[ix]
        df_keep2 = df_keep.sample(500, random_state=0)
        df_keep2.to_csv(dirname+f"{m}_dend_spikes.csv", index=False)

        df_keep = pd.read_csv(dirname+f"{m}_dend_psp.csv")
        ix = [i for i in df_keep.index if i not in idx_zeros_dends]
        df_keep1 = df_keep.iloc[ix]
        df_keep2 = df_keep.sample(500, random_state=0)
        df_keep2.to_csv(dirname+f"{m}_dend_psp.csv", index=False)

        df_keep = pd.read_csv(dirname+f"{m}_soma_psp.csv")
        ix = [i for i in df_keep.index if i not in idx_zeros_soma]
        df_keep1 = df_keep.iloc[ix]
        df_keep2 = df_keep.sample(500, random_state=0)
        df_keep2.to_csv(dirname+f"{m}_soma_psp.csv", index=False)

analyzed_data = {}
analyzed_data['d_Spike'] = df_dend_dSpike
analyzed_data['d_psp'] = df_dend_psp
analyzed_data['s_psp'] = df_soma_psp
analyzed_data['dend'] = df_dend
analyzed_data['soma'] = df_soma
analyzed_data['all'] = df

dirname = f"data_synapses_new{repo}/"
with open(dirname + 'circuitry_analyzed.pkl', 'wb') as handle:
    pickle.dump(analyzed_data, handle, protocol=pickle.HIGHEST_PROTOCOL)


