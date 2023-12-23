# -*- coding: utf-8 -*-
"""
Spyder Editor.

"""
import time, sys
from neuron import h
import numpy as np
from opt import myinit
from cell_models import CCKCell, PyramidalCell, OLMCell, VIPCCKCell, VIPCRCell

t1sim = time.time()

# Load the console parameters
reps = int(sys.argv[1])  # repetitions
trials = int(sys.argv[2])  # trials
deletion = int(sys.argv[3])  # simulated condition

# Set the random seed for reproducibility
np.random.seed(1000*reps + trials)

if deletion == 0: 
    factorVIPCCK = 1
    factorVIPCR = 1
    factorCCK = 1
    factorOLM = 1
    case_name = "control"
elif deletion == 1:
    factorVIPCCK = 1
    factorVIPCR = 0
    factorCCK = 1
    factorOLM = 1
    case_name = "vipcr_del"
elif deletion == 2:
    factorVIPCCK = 0
    factorVIPCR = 1
    factorCCK = 1
    factorOLM = 1
    case_name = "vipcck_del"
elif deletion == 3:
    factorVIPCCK = 0
    factorVIPCR = 0
    factorCCK = 1
    factorOLM = 1
    case_name = "vip_del"
elif deletion == 4:
    factorVIPCCK = 1
    factorVIPCR = 1
    factorCCK = 0
    factorOLM = 1
    case_name = "cck_del"
elif deletion == 5:
    factorVIPCCK = 1
    factorVIPCR = 1
    factorCCK = 1
    factorOLM = 0
    case_name = "olm_del"

# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")

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

# Location of poisson inputs, theta cycle (8 Hz)
rate_l = int(sys.argv[4])
fname_ = f"background_noise/rate{rate_l}/run_1/noise_"

# Shuffled inputs

inputs_ids = np.random.choice(range(10000), size=10, replace=False)

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

rLECtoPN = int(np.random.exponential(scale=12))
if rLECtoPN == 0:
    rLECtoPN = 3

gAMPA_LECtoPC = 0.3 * 8.20 * 6.0e-4
gNMDA_LECtoPC = 0.075 * 8.20 * 6.0e-4

# Synapses lists
synAMPA, vsAMPA, locLECtoPN, lecPN = [], [], [], []
synNMDA, vsNMDA = [], []
synGABA, vsGABA, locINtoPN, lec2 = [], [], [], []
inputsPN = []

# Synapses from LEC to PC
inputs_idxs = np.random.choice(inputs_ids, rLECtoPN, replace=True)

for i in range(rLECtoPN):
    # Stimulus
    spiketimes = np.loadtxt(f'{fname_}{inputs_idxs[i]}.txt')
    if len(spiketimes) == 0:
        spiketimes = [1]
    inputsPN.append(h.Vector(spiketimes))
    
    lecPN.append(h.VecStim())
    lecPN[-1].play(inputsPN[-1])

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
lecIntoPN, inputsINtoPN = [], []
gGABA_INstoPCs = 0.4 * 0.80 * 2.0e-4*3
rINtoPN = int(np.random.randn()*5 + 40)
if rINtoPN == 0:
    rINtoPN = 1

for jj in range(rINtoPN):
    r_loc = np.random.rand()
    if r_loc < 0.2:
        # New Synapse
        # print('New GABA @ trunk')
        locINtoPN.append(np.random.randint(low=0, high=len(proximal_dends_inh)))
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
        r_dinh = np.random.rand()
        if r_dinh < 0.4:
            vsGABA.append(h.NetCon(cck.soma(0.5)._ref_v, synGABA[-1],
                                   sec=cck.soma))
            vsGABA[-1].weight[0] = factorCCK*gGABA_INstoPCs*8
        else:
            vsGABA.append(h.NetCon(olm.soma(0.5)._ref_v, synGABA[-1],
                                   sec=olm.soma))
            vsGABA[-1].weight[0] = factorOLM*gGABA_INstoPCs*8
        vsGABA[-1].threshold = -10
        vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
    else:
        rpick = np.random.rand()
        if rpick < 0.2:
            # print('New GABA @ soma')
            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(cck.soma(0.5)._ref_v, synGABA[-1], sec=cck.soma))
            vsGABA[-1].threshold = -10
            vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = factorCCK*gGABA_INstoPCs*2
        elif 0.2 <= rpick < 0.7:
            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(vipcck.soma(0.5)._ref_v, synGABA[-1], sec=vipcck.soma))
            vsGABA[-1].threshold = -10
            vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = factorVIPCCK*gGABA_INstoPCs*2
        else:
            # Stimulus
            jIN = np.random.choice(inputs_ids, 1).item()
            spiketimes = np.loadtxt(f'{fname_}{jIN}.txt')
            inputsINtoPN.append(h.Vector(spiketimes))

            lecIntoPN.append(h.VecStim())
            lecIntoPN[-1].play(inputsINtoPN[-1])

            synGABA.append(h.Exp2Syn(pyramidal.soma(0.5)))
            synGABA[-1].tau1 = 0.3  # rise time
            synGABA[-1].tau2 = 8.0  # decay time
            synGABA[-1].e = -75  # reversal potential
            # New VecStim object
            vsGABA.append(h.NetCon(lecIntoPN[-1], synGABA[-1]))
            vsGABA[-1].delay = 10 + np.random.rand() * 0.1  # delay in ms
            vsGABA[-1].weight[0] = gGABA_INstoPCs*2
# ========================================================================
# SYNAPSES ONTO CCK CELLS
# ========================================================================
eps = np.random.rand()
if eps < 19/23:
    rLECtoCCK = np.random.poisson(lam=9, size=1).item()
else:
    rLECtoCCK = 0

dendsCCK = []
for sec in cck.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsCCK.append(sec(0.5))

gAMPA_LECtoCCK = 0.0013038*1.4

lecCCK, locCCK = [], []
synAMPACCK, vsAMPACCK = [], []
synGABACCK, vsGABACCK = [], []
inputs_CCK = []
# Synapses from LEC to CCK
inputs_idxs = np.random.choice(inputs_ids, rLECtoCCK, replace=True)

for i in range(rLECtoCCK):
    # Stimulus
    spiketimes = np.loadtxt(f'{fname_}{inputs_idxs[i]}.txt')
    inputs_CCK.append(h.Vector(spiketimes))

    lecCCK.append(h.VecStim())
    lecCCK[-1].play(inputs_CCK[-1])

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
    vsAMPACCK[-1].weight[0] = gAMPA_LECtoCCK  #

rVIPtoCCK = int(np.random.exponential(scale=60))
for i in range(rVIPtoCCK):
    synGABACCK.append(h.Exp2Syn(cck.soma(0.5)))
    synGABACCK[-1].tau1 = 0.73 # rise time
    synGABACCK[-1].tau2 = 5.1  # decay time
    synGABACCK[-1].e = -75  # reversal potential
    # New VecStim object
    vsGABACCK.append(h.NetCon(vipcr.soma(0.5)._ref_v, synGABACCK[-1], sec=vipcr.soma))
    vsGABACCK[-1].threshold = -10
    vsGABACCK[-1].delay = 1.7  # delay in ms
    vsGABACCK[-1].weight[0] = factorVIPCR*gGABA_INstoPCs*1000
# ========================================================================
# SYNAPSES ONTO VIP/CCK CELLS
# ========================================================================
eps = np.random.rand()
if eps < 21/25:
    rLECtoVIPCCK = np.random.poisson(lam=21, size=1).item()
else:
    rLECtoVIPCCK = 0

dendsVIPCCK = []
for sec in vipcck.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsVIPCCK.append(sec(0.5))

gAMPA_LECtoVIPCCK = 0.00055151*1.4

lecVIPCCK, locVIPCCK = [], []
synAMPAVIPCCK, vsAMPAVIPCCK = [], []
inputs_VIPCCK = []
# Synapses from LEC to VIP/CCK
inputs_idxs = np.random.choice(inputs_ids, rLECtoVIPCCK, replace=True)

for i in range(rLECtoVIPCCK):
    # Stimulus
    spiketimes = np.loadtxt(f'{fname_}{inputs_idxs[i]}.txt')
    inputs_VIPCCK.append(h.Vector(spiketimes))

    lecVIPCCK.append(h.VecStim())
    lecVIPCCK[-1].play(inputs_VIPCCK[-1])

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
    rLECtoVIPCR = np.random.poisson(lam=5, size=1).item()
else:
    rLECtoVIPCR = 0

dendsVIPCR = []
for sec in vipcr.all:
    if "lmM" in str(sec) or "radDist" in str(sec):
        dendsVIPCR.append(sec(0.5))

gAMPA_LECtoVIPCR = 0.0008533371*1.4

lecVIPCR, locVIPCR = [], []
synAMPAVIPCR, vsAMPAVIPCR = [], []
inputs_VIPCR = []
# Synapses from LEC to VIP/CR
inputs_idxs = np.random.choice(inputs_ids, rLECtoVIPCR, replace=True)

for i in range(rLECtoVIPCR):
    # Stimulus
    spiketimes = np.loadtxt(f'{fname_}{inputs_idxs[i]}.txt')
    inputs_VIPCR.append(h.Vector(spiketimes))

    lecVIPCR.append(h.VecStim())
    lecVIPCR[-1].play(inputs_VIPCR[-1])

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
rLECtoOLM = np.random.poisson(lam=0.3, size=1).item()
dendsOLM = []
for sec in olm.all:
    if "dend" in str(sec):
        dendsOLM.append(sec(0.5))

gAMPA_LECtoOLM = 0.00091266*1.4
gAMPA_PNtoOLM = 0.0025

lecOLM, locOLM = [], []
synAMPAOLM, vsAMPAOLM = [], []
synGABAOLM, vsGABAOLM = [], []
synAMPA_PN, vsAMPA_PN = [], []
inputs_OLM = []
# Synapses from LEC to OLM
inputs_idxs = np.random.choice(inputs_ids, rLECtoOLM, replace=True)

for i in range(rLECtoOLM):
    # Stimulus
    spiketimes = np.loadtxt(f'{fname_}{inputs_idxs[i]}.txt')
    inputs_OLM.append(h.Vector(spiketimes))

    lecOLM.append(h.VecStim())
    lecOLM[-1].play(inputs_OLM[-1])

    # Choose a location at random @ distal dends
    locOLM.append(np.random.randint(low=0, high=len(dendsOLM)))
    # New Synapse and store in a list
    # AMPA synapse
    synAMPAOLM.append(h.Exp2Syn(dendsOLM[locOLM[-1]]))
    synAMPAOLM[-1].e = 0  # reversal potential
    synAMPAOLM[-1].tau1 = 0.5  # rise time
    synAMPAOLM[-1].tau2 = 8.0  # decay time
    # New VecStim object
    vsAMPAOLM.append(h.NetCon(lecVIPCR[-1], synAMPAOLM[-1]))
    vsAMPAOLM[-1].delay = 2.0  # delay in ms
    vsAMPAOLM[-1].weight[0] = gAMPA_LECtoOLM

# Synapses from PN
rPNtoOLM = int(np.random.exponential(scale=10))
for i in range(rPNtoOLM):
    locOLM.append(np.random.randint(low=0, high=len(dendsOLM)))
    # AMPA synapse
    synAMPA_PN.append(h.Exp2Syn(dendsOLM[locOLM[-1]]))
    synAMPA_PN[-1].e = 0  # reversal potential
    synAMPA_PN[-1].tau1 = 0.5  # rise time
    synAMPA_PN[-1].tau2 = 8.0  # decay time
    # New VecStim object
    vsAMPA_PN.append(h.NetCon(pyramidal.axon(0.5)._ref_v, synAMPA_PN[-1], sec=pyramidal.axon))
    vsAMPA_PN[-1].threshold = -10
    vsAMPA_PN[-1].delay = 2.0  # delay in ms
    vsAMPA_PN[-1].weight[0] = gAMPA_PNtoOLM

# Synapses from VIP/CR
rVIPtoOLM = int(np.random.exponential(scale=60))
for i in range(rVIPtoOLM):
    synGABAOLM.append(h.Exp2Syn(olm.soma(0.5)))
    synGABAOLM[-1].tau1 = 1.3  # rise time
    synGABAOLM[-1].tau2 = 9.0  # decay time
    synGABAOLM[-1].e = -75  # reversal potential
    # New VecStim object
    vsGABAOLM.append(h.NetCon(vipcr.soma(0.5)._ref_v, synGABAOLM[-1], sec=vipcr.soma))
    vsGABAOLM[-1].threshold = -10
    vsGABAOLM[-1].delay = 1.9 + np.random.randn() * 0.2  # delay in ms
    vsGABAOLM[-1].weight[0] = factorVIPCR*gGABA_INstoPCs*100
# =========================================================================
# SAVE VECTORS
# =========================================================================
# Membrane potential vectors and recording locations
soma_v_vec = h.Vector().record(pyramidal.soma(0.5)._ref_v)
radTdist3_v_vec = h.Vector().record(pyramidal.radTdist3(0.1)._ref_v)
cck_v_vec = h.Vector().record(cck.soma(0.5)._ref_v)
vipcck_v_vec = h.Vector().record(vipcck.soma(0.5)._ref_v)
vipcr_v_vec = h.Vector().record(vipcr.soma(0.5)._ref_v)
olm_v_vec = h.Vector().record(olm.soma(0.5)._ref_v)

# Time stamp vector
t_vec = h.Vector().record(h._ref_t)
# =========================================================================
# STIMULATION
# =========================================================================
stim2 = h.IClamp(pyramidal.soma(0.5))
stim2.delay = 400
stim2.dur = 3600
stim2.amp = 0.12 #nA

# Simulation parameters
simdur = 4000.0
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

n1 = np.abs(t_vec - 390).argmin()
n2 = np.abs(t_vec - 3990).argmin()


v_soma = soma_v_vec[n1:n2]
v_dend = radTdist3_v_vec[n1:n2]
v_cck = cck_v_vec[n1:n2]
v_olm = olm_v_vec[n1:n2]
v_vipcck = vipcck_v_vec[n1:n2]
v_vipcr = vipcr_v_vec[n1:n2]
time_ = t_vec[n1:n2]

results = {} 
results[f"{case_name}_voltage_soma"] = v_soma
results[f"{case_name}_voltage_dend"] = v_dend
results[f"{case_name}_voltage_cck"] = v_cck
results[f"{case_name}_voltage_olm"] = v_olm
results[f"{case_name}_voltage_vipcck"] = v_vipcck
results[f"{case_name}_voltage_vipcr"] = v_vipcr
results[f"{case_name}_time"] = time_

import pickle, os
dirname = f"data_synapses_poisson_rate_{rate_l}/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + f'circuitry_trial{trials}_rep{reps}_{case_name}.pkl', 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

print(f'Simulation time: {np.round(time.time() - t1sim, 2)} secs')
