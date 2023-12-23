#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 14 09:37:20 2021.

@author: spiros
"""
import time
import sys
import pickle
import os
from neuron import h
import numpy as np
import pandas as pd
from synaptic_metrics import synaptic_metrics
from opt import myinit
from cell_models import CCKCell, PyramidalCell, OLMCell, VIPCCKCell, VIPCRCell


t1sim = time.time()

clear = True

# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")

# Options given from console
reps = int(sys.argv[1])
trials = int(sys.argv[2])
deletion = int(sys.argv[3])

# Random seed
np.random.seed(1000*reps + trials)

# Construct the DataFrame
cols = ["peak", "time_rise", "time_decay", "dvdt", "latency", "thalf", "case"]
df_soma = pd.DataFrame(columns=cols, index=range(1))
df_dend = pd.DataFrame(columns=cols, index=range(1))

# Preallocate variables
results = {}
counter = 0
synsE, synsI = [], []
voltages_s, voltages_d, times = [], [], []
voltages_cck, voltages_olm, voltages_vipcck, voltages_vipcr = [], [], [], []

if deletion == 0:
    factorVIPCCK = 1
    factorVIPCR = 1
    factorCCK = 1
    case_name = "control"
elif deletion == 1:
    factorVIPCCK = 1
    factorVIPCR = 0
    factorCCK = 1
    case_name = "vipcr_del"
elif deletion == 2:
    factorVIPCCK = 0
    factorVIPCR = 1
    factorCCK = 1
    case_name = "vipcck_del"
elif deletion == 3:
    factorVIPCCK = 0
    factorVIPCR = 0
    factorCCK = 1
    case_name = "vip_del"
elif deletion == 4:
    factorVIPCCK = 1
    factorVIPCR = 1
    factorCCK = 0
    case_name = "cck_del"

# plot_all, plot_single = True, False
plot_single = False

# ========================================================================
# LOAD THE CELLS HOC FILES
# ========================================================================
# Load the Pyramidal Cell template
pyramidal = PyramidalCell(0)
pyramidal.soma.ic_constant *= 3.5
for sec in pyramidal.trunk:
    sec.ic_constant *= 0.8
for sec in pyramidal.rad:
    sec.ic_constant *= 0.8
for sec in pyramidal.slm:
    sec.ic_constant *= 0.8
for sec in pyramidal.ori:
    sec.ic_constant *= 0.8
# Load the CCK cell template
cck = CCKCell(0)
# Load the VIP/CCK cell template
vipcck = VIPCCKCell(0)
# Load the VIP/CR cell template
vipcr = VIPCRCell(0)
# Load the OLM cell template
olm = OLMCell(0)

# Stimulus time -- onset
spiketimes = [700.0]

# ========================================================================
# SYNAPSES ONTO PYRAMIDAL CELLS
# ========================================================================
slms = []
for sec in pyramidal.slm:
    if "thin" not in str(sec):
        slms.append(sec(0.5))

rads = []
for sec in pyramidal.rad:
    if "thin" not in str(sec):
        rads.append(sec(0.5))

trunks_distal, trunks_proximal = [], []
for sec in pyramidal.trunk:
    if "dist" in str(sec):
        trunks_distal.append(sec(0.5))
    else:
        trunks_proximal.append(sec(0.5))

trunks = trunks_proximal[1:]

proximal_dends = trunks + rads + trunks_distal[:-1]
distal_dends = slms + [trunks_distal[-1]]

proximal_dends_inh = []
for sec in proximal_dends:
    if "thin" not in str(sec):
        proximal_dends_inh.append(sec)

# Basal dendrites -- not used
oris = []
for sec in pyramidal.ori:
    oris.append(sec(0.5))

r1 = int(np.random.exponential(scale=12))
if r1 == 0:
    r1 = 3

gAMPA_LECtoPC = 0.3 * 8.20 * 6.0e-4
gNMDA_LECtoPC = 0.075 * 8.20 * 6.0e-4

# Synapses lists
synAMPA, vsAMPA, locLECtoPN, lecPN = [], [], [], []
synNMDA, vsNMDA = [], []
synGABA, vsGABA, locINtoPN, lec2 = [], [], [], []

# Synapses from LEC to PC
for i in range(r1):
    # Stimulus
    lecPN.append(h.NetStim())
    lecPN[-1].number = 1
    lecPN[-1].start = spiketimes[0]

    # New Synapse and store in a list
    # Choose a location at random @ SLM
    locLECtoPN.append(np.random.randint(low=0, high=len(distal_dends)))
    ######################################################################
    # AMPA synapse
    synAMPA.append(h.Exp2Syn(distal_dends[locLECtoPN[-1]]))
    synAMPA[-1].e = 0  # reversal potential
    synAMPA[-1].tau1 = 1.0  # rise time
    synAMPA[-1].tau2 = 10.0  # decay time
    # New VecStim object
    vsAMPA.append(h.NetCon(lecPN[-1], synAMPA[-1]))
    vsAMPA[-1].delay = 10.0  # delay in ms
    vsAMPA[-1].weight[0] = gAMPA_LECtoPC
    ######################################################################
    # NMDA synapse
    synNMDA.append(h.NMDA(distal_dends[locLECtoPN[-1]]))
    synNMDA[-1].tcon = 2.3  # rise time
    synNMDA[-1].tcoff = 100  # decay time
    # New VecStim object
    vsNMDA.append(h.NetCon(lecPN[-1], synNMDA[-1]))
    vsNMDA[-1].delay = 10.0  # delay in ms
    vsNMDA[-1].weight[0] = gNMDA_LECtoPC
    ######################################################################

# Inhibitory synapses on Pyramidal Cells
gGABA_INstoPCs = 0.4 * 0.80 * 2.0e-4*3
r2 = int(np.random.randn()*5 + 40)
if r2 == 0:
    r2 = 1

for jj in range(r2):
    r_loc = np.random.rand()
    if r_loc < 0.2:
        # New Synapse
        # print('New GABA @ trunk')
        locINtoPN.append(np.random.randint(low=0,
                                           high=len(proximal_dends_inh)))
        synGABA.append(h.Exp2Syn(proximal_dends_inh[locINtoPN[-1]]))
        if np.random.rand() < 0.5:
            # GABA-A
            synGABA[-1].tau1 = 0.11  # rise time
            synGABA[-1].tau2 = 9.70  # decay time
        else:
            # GABA-B
            synGABA[-1].tau1 = 1.0  # rise time
            synGABA[-1].tau2 = 50.0  # decay time
        synGABA[-1].e = -75  # reversal potential
        # New VecStim object
        vsGABA.append(h.NetCon(cck.soma(0.5)._ref_v, synGABA[-1],
                               sec=cck.soma))
        vsGABA[-1].threshold = -10
        vsGABA[-1].delay = 1.0 + np.random.rand() * 0.2  # delay in ms
        vsGABA[-1].weight[0] = factorCCK*gGABA_INstoPCs
    elif 0.2 <= r_loc < 0.8:
        # New Synapse
        # print('New GABA @ slm')
        locINtoPN.append(np.random.randint(low=0, high=len(distal_dends)))
        synGABA.append(h.Exp2Syn(distal_dends[locINtoPN[-1]]))
        if np.random.rand() < 0.5:
            # GABA-A
            synGABA[-1].tau1 = 1.0  # rise time
            synGABA[-1].tau2 = 11.0  # decay time
        else:
            # GABA-B
            synGABA[-1].tau1 = 1.0  # rise time
            synGABA[-1].tau2 = 50.0  # decay time
        synGABA[-1].e = -75  # reversal potential
        # New VecStim object
        vsGABA.append(h.NetCon(cck.soma(0.5)._ref_v, synGABA[-1],
                               sec=cck.soma))
        vsGABA[-1].threshold = -10
        vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
        vsGABA[-1].weight[0] = factorCCK*gGABA_INstoPCs*8
    else:
        rpick = np.random.rand()
        if rpick < 0.2:
            # print('New GABA @ soma')
            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(cck.soma(0.5)._ref_v, synGABA[-1],
                                   sec=cck.soma))
            vsGABA[-1].threshold = -10
            vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = factorCCK*gGABA_INstoPCs*2
        elif 0.2 <= rpick < 0.7:
            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(vipcck.soma(0.5)._ref_v, synGABA[-1],
                                   sec=vipcck.soma))
            vsGABA[-1].threshold = -10
            vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = factorVIPCCK*gGABA_INstoPCs*2
        else:
            # Stimulus
            lec2.append(h.NetStim())
            lec2[-1].number = 1
            lec2[-1].start = spiketimes[0]

            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(lec2[-1], synGABA[-1]))
            vsGABA[-1].delay = 10 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = gGABA_INstoPCs*2
# ========================================================================
# SYNAPSES ONTO CCK CELLS
# ========================================================================
eps = np.random.rand()
if eps < 19/23:
    rCCK = np.random.poisson(lam=10, size=1).item()
else:
    rCCK = 0

rCCK = 10
dendsCCK = []
for sec in cck.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsCCK.append(sec(0.5))

gAMPA_LECtoCCK = 0.0013038*1.4

lecCCK, locCCK = [], []
synAMPACCK, vsAMPACCK = [], []
synGABACCK, vsGABACCK = [], []
for i in range(rCCK):
    # Stimulus
    lecCCK.append(h.NetStim())
    lecCCK[-1].number = 1
    lecCCK[-1].start = spiketimes[0]

    # Choose a location at random @ distal dends
    locCCK.append(np.random.randint(low=0, high=len(dendsCCK)))
    # New Synapse and store in a list
    # AMPA synapse
    synAMPACCK.append(h.Exp2Syn(dendsCCK[locCCK[-1]]))
    synAMPACCK[-1].e = 0  # reversal potential
    synAMPACCK[-1].tau1 = 0.5  # rise time
    synAMPACCK[-1].tau2 = 3.0  # decay time
    # New VecStim object
    vsAMPACCK.append(h.NetCon(lecCCK[-1], synAMPACCK[-1]))
    vsAMPACCK[-1].delay = 2  # delay in ms
    vsAMPACCK[-1].weight[0] = gAMPA_LECtoCCK

rVIPtoCCK = int(np.random.exponential(scale=60))
for i in range(rVIPtoCCK):
    synGABACCK.append(h.Exp2Syn(cck.soma(0.5)))
    synGABACCK[-1].tau1 = 0.73  # rise time
    synGABACCK[-1].tau2 = 5.1  # decay time
    synGABACCK[-1].e = -75  # reversal potential
    # New VecStim object
    vsGABACCK.append(h.NetCon(vipcr.soma(0.5)._ref_v, synGABACCK[-1],
                              sec=vipcr.soma))
    vsGABACCK[-1].threshold = -10
    vsGABACCK[-1].delay = 1.7  # delay in ms
    vsGABACCK[-1].weight[0] = factorVIPCR*gGABA_INstoPCs*1000
# ========================================================================
# SYNAPSES ONTO VIP/CCK CELLS
# ========================================================================
eps = np.random.rand()
if eps < 21/25:
    rVIPCCK = np.random.poisson(lam=21, size=1).item()
else:
    rVIPCCK = 0

dendsVIPCCK = []
for sec in vipcck.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsVIPCCK.append(sec(0.5))

gAMPA_LECtoVIPCCK = 0.00055151*1.4

lecVIPCCK, locVIPCCK = [], []
synAMPAVIPCCK, vsAMPAVIPCCK = [], []
for i in range(rVIPCCK):
    # Stimulus
    lecVIPCCK.append(h.NetStim())
    lecVIPCCK[-1].number = 1
    lecVIPCCK[-1].start = spiketimes[0]

    # Choose a location at random @ distal dends
    locVIPCCK.append(np.random.randint(low=0, high=len(dendsVIPCCK)))
    # New Synapse and store in a list
    # AMPA synapse
    synAMPAVIPCCK.append(h.Exp2Syn(dendsVIPCCK[locVIPCCK[-1]]))
    synAMPAVIPCCK[-1].e = 0  # reversal potential
    synAMPAVIPCCK[-1].tau1 = 0.5  # rise time
    synAMPAVIPCCK[-1].tau2 = 3.0  # decay time
    # New VecStim object
    vsAMPAVIPCCK.append(h.NetCon(lecVIPCCK[-1], synAMPAVIPCCK[-1]))
    vsAMPAVIPCCK[-1].delay = 1.0  # delay in ms
    vsAMPAVIPCCK[-1].weight[0] = gAMPA_LECtoVIPCCK  #

# ========================================================================
# SYNAPSES ONTO VIP/CR CELLS
# ========================================================================
eps = np.random.rand()
if eps < 21/23:
    rVIPCR = np.random.poisson(lam=5, size=1).item()
else:
    rVIPCR = 0

rVIPCR = 10
dendsVIPCR = []
for sec in vipcr.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsVIPCR.append(sec(0.5))

gAMPA_LECtoVIPCR = 0.0008533371*1.4

lecVIPCR, locVIPCR = [], []
synAMPAVIPCR, vsAMPAVIPCR = [], []

for i in range(rVIPCR):
    # Stimulus
    lecVIPCR.append(h.NetStim())
    lecVIPCR[-1].number = 1
    lecVIPCR[-1].start = spiketimes[0]

    # Choose a location at random @ distal dends
    locVIPCR.append(np.random.randint(low=0, high=len(dendsVIPCR)))
    # New Synapse and store in a list
    # AMPA synapse
    synAMPAVIPCR.append(h.Exp2Syn(dendsVIPCR[locVIPCR[-1]]))
    synAMPAVIPCR[-1].e = 0  # reversal potential
    synAMPAVIPCR[-1].tau1 = 0.5  # rise time
    synAMPAVIPCR[-1].tau2 = 3.0  # decay time
    # New VecStim object
    vsAMPAVIPCR.append(h.NetCon(lecVIPCR[-1], synAMPAVIPCR[-1]))
    vsAMPAVIPCR[-1].delay = 0.5  # delay in ms
    vsAMPAVIPCR[-1].weight[0] = gAMPA_LECtoVIPCR  #

# ========================================================================
# SYNAPSES ONTO OLM CELLS
# ========================================================================
rOLM = np.random.poisson(lam=0.3, size=1).item()
dendsOLM = []
for sec in olm.all:
    if "dend" in str(sec):
        dendsOLM.append(sec(0.5))

gAMPA_LECtoOLM = 0.00091266*1.4

lecOLM, locOLM = [], []
synAMPAOLM, vsAMPAOLM = [], []
synGABAOLM, vsGABAOLM = [], []
for i in range(rOLM):
    # Stimulus
    lecOLM.append(h.NetStim())
    lecOLM[-1].number = 1
    lecOLM[-1].start = spiketimes[0]

    # Choose a location at random @ distal dends
    locOLM.append(np.random.randint(low=0, high=len(dendsOLM)))
    # New Synapse and store in a list
    # AMPA synapse
    synAMPA.append(h.Exp2Syn(dendsOLM[locOLM[-1]]))
    synAMPA[-1].e = 0  # reversal potential
    synAMPA[-1].tau1 = 0.5  # rise time
    synAMPA[-1].tau2 = 8.0  # decay time
    # New VecStim object
    vsAMPA.append(h.NetCon(lecOLM[-1], synAMPA[-1]))
    vsAMPA[-1].delay = 2.0  # delay in ms
    vsAMPA[-1].weight[0] = gAMPA_LECtoOLM

rVIPtoOLM = int(np.random.exponential(scale=60))
for i in range(rVIPtoOLM):
    synGABAOLM.append(h.Exp2Syn(olm.soma(0.5)))
    synGABAOLM[-1].tau1 = 1.3  # rise time
    synGABAOLM[-1].tau2 = 9.0  # decay time
    synGABAOLM[-1].e = -75  # reversal potential
    # New VecStim object
    vsGABAOLM.append(h.NetCon(vipcr.soma(0.5)._ref_v, synGABAOLM[-1],
                              sec=vipcr.soma))
    vsGABAOLM[-1].threshold = -10
    vsGABAOLM[-1].delay = 1.9 + np.random.randn() * 0.2  # delay in ms
    vsGABAOLM[-1].weight[0] = factorVIPCR*gGABA_INstoPCs*100
# =========================================================================
# SAVE VECTORS
# =========================================================================
# Membrane potential vector
soma_v_vec = h.Vector()
radTdist3_v_vec = h.Vector()
cck_v_vec = h.Vector()
vipcck_v_vec = h.Vector()
vipcr_v_vec = h.Vector()
olm_v_vec = h.Vector()

# Recording sites
soma_v_vec.record(pyramidal.soma(0.5)._ref_v)
radTdist3_v_vec.record(pyramidal.radTdist3(0.1)._ref_v)
cck_v_vec.record(cck.soma(0.5)._ref_v)
vipcck_v_vec.record(vipcck.soma(0.5)._ref_v)
vipcr_v_vec.record(vipcr.soma(0.5)._ref_v)
olm_v_vec.record(olm.soma(0.5)._ref_v)

# Time stamp vector
t_vec = h.Vector()
t_vec.record(h._ref_t)

# =========================================================================
# STIMULATION
# =========================================================================
# stim2 = h.IClamp(cck.soma(0.5))
# stim2.delay = 0
# stim2.dur = 2000
# stim2.amp = 0#-0.02  # 0.25 #nA

# Simulation parameters
simdur = 2000.0
myinit(modify=True)
h.continuerun(simdur)

# =========================================================================
# SAVE and PLOT the results
# =========================================================================
t_vec = np.array(t_vec)
soma_v_vec = np.array(soma_v_vec)
radTdist3_v_vec = np.array(radTdist3_v_vec)
cck_v_vec = np.array(cck_v_vec)
vipcck_v_vec = np.array(vipcck_v_vec)
vipcr_v_vec = np.array(vipcr_v_vec)
olm_v_vec = np.array(olm_v_vec)

n1 = np.abs(t_vec - 690).argmin()
n2 = np.abs(t_vec - 900).argmin()

v_soma = soma_v_vec[n1:n2] - soma_v_vec[n1]
v_dend = radTdist3_v_vec[n1:n2] - radTdist3_v_vec[n1]
time_ = t_vec[n1:n2]

peak_s, trise_s, thalf_s, tdecay_s, dvdt_s, latency_s, time_peak_s = synaptic_metrics(
    v_soma, time_
)
peak_d, trise_d, thalf_d, tdecay_d, dvdt_d, latency_d, time_peak_d = synaptic_metrics(
    v_dend, time_
)

df_soma.loc[counter].peak = peak_s
df_soma.loc[counter].time_rise = trise_s
df_soma.loc[counter].time_decay = tdecay_s
df_soma.loc[counter].dvdt = dvdt_s
df_soma.loc[counter].thalf = thalf_s
df_soma.loc[counter].latency = latency_s

df_dend.loc[counter].peak = peak_d
df_dend.loc[counter].time_rise = trise_d
df_dend.loc[counter].time_decay = tdecay_d
df_dend.loc[counter].dvdt = dvdt_d
df_dend.loc[counter].thalf = thalf_d
df_dend.loc[counter].latency = latency_d

# Add a column with the deletion
df_soma.loc[counter].case = case_name
df_dend.loc[counter].case = case_name

counter += 1

v_cck = cck_v_vec[n1:n2] - cck_v_vec[n1]
v_olm = olm_v_vec[n1:n2] - olm_v_vec[n1]
v_vipcck = vipcck_v_vec[n1:n2] - vipcck_v_vec[n1]
v_vipcr = vipcr_v_vec[n1:n2] - vipcr_v_vec[n1]

results[f"{case_name}_voltage_soma"] = v_soma
results[f"{case_name}_voltage_dend"] = v_dend
results[f"{case_name}_time"] = time_

results[f"{case_name}_voltage_cck"] = v_cck
results[f"{case_name}_voltage_olm"] = v_olm
results[f"{case_name}_voltage_vipcck"] = v_vipcck
results[f"{case_name}_voltage_vipcr"] = v_vipcr

all_dict = {}
all_dict['results'] = results
all_dict['soma'] = df_soma
all_dict['dend'] = df_dend


dirname = "data_synapses_new2/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + f'circuitry_trial{trials}_rep{reps}_{case_name}.pkl', 'wb') as handle:
    pickle.dump(all_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)

print(f'Simulation time: {np.round(time.time() - t1sim, 2)} secs')
