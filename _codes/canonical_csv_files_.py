#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 17:45:40 2022

@author: spiros
"""

import os
import pickle
import pandas as pd


def zero_entries(df):
    return df[(df.values == 0).any(axis=1)].index.values.tolist()


dirname = "data_synapses_new2/"
with open(dirname + 'circuitry_analyzed.pkl', 'rb') as handle:
    DATA = pickle.load(handle)
    
df = DATA['all']
df_dend = DATA['dend']
df_soma = DATA['soma']
df_dend_psp = DATA['d_psp']
df_soma_psp = DATA['s_psp']
df_dend_dSpike = DATA['d_Spike']


conditions = ["control", "vip_del", "vipcr_del", "vipcck_del", "cck_del"]
dirname = "canonical_circuit/"
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
        zeros = zero_entries(df_zeros)
        # append the new indices from case m to be removed from all metrics.
        for j in zeros:
            if j not in idx_zeros_dends_spikes:
                idx_zeros_dends_spikes.append(j)

        df_zeros = pd.read_csv(dirname+f"{m}_dend_psp.csv")
        # indices to keep ('false' in the previous list).
        zeros = zero_entries(df_zeros)
        # append the new indices from case m to be removed from all metrics.
        for j in zeros:
            if j not in idx_zeros_dends:
                idx_zeros_dends.append(j)

        df_zeros = pd.read_csv(dirname+f"{m}_soma_psp.csv")
        # indices to keep ('false' in the previous list).
        zeros = zero_entries(df_zeros)
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
        zeros = zero_entries(df_keep)
        print(m, zeros)
        ix = [i for i in df_keep.index if i not in idx_zeros_dends_spikes]
        df_keep1 = df_keep.iloc[ix]
        zeros = zero_entries(df_keep1)
        print(m, zeros)
        df_keep2 = df_keep1.sample(500, random_state=0)
        zeros = zero_entries(df_keep2)
        print(m, zeros)
        df_keep2.to_csv(dirname+f"{m}_dend_spikes.csv", index=False)

        df_keep = pd.read_csv(dirname+f"{m}_dend_psp.csv")
        ix = [i for i in df_keep.index if i not in idx_zeros_dends]
        df_keep1 = df_keep.iloc[ix]
        df_keep2 = df_keep1.sample(500, random_state=0)
        df_keep2.to_csv(dirname+f"{m}_dend_psp.csv", index=False)

        df_keep = pd.read_csv(dirname+f"{m}_soma_psp.csv")
        ix = [i for i in df_keep.index if i not in idx_zeros_soma]
        df_keep1 = df_keep.iloc[ix]
        df_keep2 = df_keep1.sample(500, random_state=0)
        df_keep2.to_csv(dirname+f"{m}_soma_psp.csv", index=False)
