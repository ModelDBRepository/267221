from neuron import h
import numpy as np
import pandas as pd
import sys
from cell_models import PyramidalCell
from synaptic_metrics import synaptic_metrics
from opt import myinit
import pickle, os


# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")

# Console argument -- number of repetition
rep = 1

trials = 1000
inh_deletions = [2, 1, 3, 4, 5]

cols = ["peak", "time_rise", "time_decay", "dvdt", "latency",
        "thalf", "time_peak", "case"]
df_soma = pd.DataFrame(columns=cols, index=range(trials * len(inh_deletions)))
df_dend = pd.DataFrame(columns=cols, index=range(trials * len(inh_deletions)))
results = {}
counter = 0
synsE = []
synsI = []

# for inh in inh_deletions:
for inh in inh_deletions:
    if inh == 1:
        print("Inhibition OFF...\n")
    elif inh == 2:
        print("Inhibition ON...\n")
    elif inh == 3:
        print("Somatic Inhibition OFF...\n")
    elif inh == 4:
        print("Proximal Dendritic Inhibition OFF...\n")
    elif inh == 5:
        print("Distal Dendritic Inhibition OFF...\n")

    voltages_s, voltages_d, times = [], [], []
    for ntrial in range(trials):
        np.random.seed(1000*rep + ntrial)

        # plot_all, plot_single = True, False
        plot_single = False

        # Create a Pyramidal Cell instance
        newpyr = PyramidalCell(0)

        slms = []
        for sec in newpyr.slm:
            if "thin" not in str(sec):
                slms.append(sec(0.5))

        rads = []
        for sec in newpyr.rad:
            if "thin" not in str(sec):
                rads.append(sec(0.5))

        trunks_distal, trunks_proximal = [], []
        for sec in newpyr.trunk:
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
        for sec in newpyr.ori:
            oris.append(sec(0.5))

        # Synapses lists
        synAMPA1, vsAMPA1, loc1, lec = [], [], [], []
        synNMDA1, vsNMDA1, inputs = [], [], []

        synGABA, vsGABA, loc4 = [], [], []

        gAMPA = 0.3 * 8.20 * 6.0e-4
        gNMDA = 0.075 * 8.20 * 6.0e-4
        gGABA = 0.4 * 0.80 * 2.0e-4*3

        r1 = int(np.random.exponential(scale=12))
        if r1 == 0:
            r1 = 3
        r2 = int(np.random.randn()*5 + 40)
        if r2 == 0:
            r2 = 1

        synsE.append(r1)
        synsI.append(r2)
        print(f"Trial {ntrial}...Esyn: {r1}, Isyn: {r2}")
        for i in range(r1):
            spiketimes = [700.0]
            inputs.append(h.Vector(spiketimes))
            lec.append(h.VecStim())
            lec[-1].play(inputs[-1])

            # Choose a location randomly @ distal_dends
            loc1.append(np.random.randint(low=0, high=len(distal_dends)))

            # New Synapse and store in a list
            # AMPA synapse
            synAMPA1.append(h.Exp2Syn(distal_dends[loc1[-1]]))
            synAMPA1[-1].e = 0  # reversal potential
            synAMPA1[-1].tau1 = 1.0  # rise time
            synAMPA1[-1].tau2 = 10.0  # decay time

            # NMDA synapse
            synNMDA1.append(h.NMDA(distal_dends[loc1[-1]]))
            synNMDA1[-1].tcon = 2.3  # rise time
            synNMDA1[-1].tcoff = 100  # decay time

            # New VecStim object
            vsAMPA1.append(h.NetCon(lec[-1], synAMPA1[-1]))
            vsAMPA1[-1].delay = 0.0  # delay in ms
            vsAMPA1[-1].weight[0] = gAMPA  #

            vsNMDA1.append(h.NetCon(lec[-1], synNMDA1[-1]))
            vsNMDA1[-1].delay = 0.0  # delay in ms
            vsNMDA1[-1].weight[0] = gNMDA

        for jj in range(r2):
            inputs.append(h.Vector(spiketimes))
            lec.append(h.VecStim())
            lec[-1].play(inputs[-1])
            r_loc = np.random.rand()
            if r_loc < 0.2:
                # New Synapse
                #print('New GABA @ proximal_dends')
                loc4.append(np.random.randint(low=0, high=len(proximal_dends_inh)))
                synGABA.append(h.Exp2Syn(proximal_dends_inh[loc4[-1]]))
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
                vsGABA.append(h.NetCon(lec[-1], synGABA[-1]))
                vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
                if inh == 4 or inh == 1:
                    vsGABA[-1].weight[0] = 0
                else:
                    vsGABA[-1].weight[0] = gGABA
            elif 0.2 <= r_loc < 0.8:
                # New Synapse
                #print('New GABA @ distal_dends')
                loc4.append(np.random.randint(low=0, high=len(distal_dends)))
                synGABA.append(h.Exp2Syn(distal_dends[loc4[-1]]))
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
                vsGABA.append(h.NetCon(lec[-1], synGABA[-1]))
                vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
                if inh == 5 or inh == 1:
                    vsGABA[-1].weight[0] = 0
                else:
                    vsGABA[-1].weight[0] = gGABA*2
            else:
                #print('New GABA @ soma')
                synGABA.append(h.Exp2Syn(newpyr.soma(0.5)))
                synGABA[-1].tau1 = 0.5  # rise time
                synGABA[-1].tau2 = 8.0  # decay time
                synGABA[-1].e = -75  # reversal potential
                # New VecStim object
                vsGABA.append(h.NetCon(lec[-1], synGABA[-1]))
                vsGABA[-1].delay = 1.0 + np.random.rand() * 0.1  # delay in ms
                if inh == 3 or inh == 1:
                    vsGABA[-1].weight[0] = 0
                else:
                    vsGABA[-1].weight[0] = gGABA*2

        # =========================================================================
        # SAVE VECTORS
        # =========================================================================
        soma_v_vec = h.Vector()  # Membrane potential vector
        radTdist3_v_vec = h.Vector()
        radTdist2_v_vec = h.Vector()
        radTdist1_v_vec = h.Vector()
        radTmed_v_vec = h.Vector()
        radTprox_v_vec = h.Vector()
        t_vec = h.Vector()  # Time stamp vector

        soma_v_vec.record(newpyr.soma(0.5)._ref_v)
        radTdist3_v_vec.record(newpyr.radTdist3(0.1)._ref_v)
        radTdist2_v_vec.record(newpyr.radTdist2(0.5)._ref_v)
        radTdist1_v_vec.record(newpyr.radTdist1(0.5)._ref_v)
        radTmed_v_vec.record(newpyr.radTmed(0.1)._ref_v)
        radTprox_v_vec.record(newpyr.radTprox(0.5)._ref_v)
        t_vec.record(h._ref_t)

        # =========================================================================
        # STIMULATION
        # =========================================================================
        stim2 = h.IClamp(newpyr.soma(0.5))
        stim2.delay = 0
        stim2.dur = 2000
        stim2.amp = 0.00  # 0.25 #nA
        
        simdur = 2000.0
        myinit()
        h.continuerun(simdur)

        # =========================================================================
        # SAVE and PLOT the results
        # =========================================================================
        t_vec = np.array(t_vec)
        soma_v_vec = np.array(soma_v_vec)
        radTprox_v_vec = np.array(radTprox_v_vec)
        radTmed_v_vec = np.array(radTmed_v_vec)
        radTdist1_v_vec = np.array(radTdist1_v_vec)
        radTdist2_v_vec = np.array(radTdist2_v_vec)
        radTdist3_v_vec = np.array(radTdist3_v_vec)
        n1 = np.abs(t_vec - 690).argmin()
        n2 = np.abs(t_vec - 900).argmin()

        # Normalize voltage at zero by subtracting the resting potential.
        v_soma = soma_v_vec[n1:n2] - soma_v_vec[n1]
        v_dend = radTdist3_v_vec[n1:n2] - radTdist3_v_vec[n1]
        time = t_vec[n1:n2]

        peak_s, trise_s, thalf_s, tdecay_s, dvdt_s, latency_s, time_peak_s = synaptic_metrics(
            v_soma, time
        )
        peak_d, trise_d, thalf_d, tdecay_d, dvdt_d, latency_d, time_peak_d = synaptic_metrics(
            v_dend, time
        )

        # Save the data in DataFrame
        # Somatic recordings
        df_soma.loc[counter].peak = peak_s
        df_soma.loc[counter].time_rise = trise_s
        df_soma.loc[counter].time_decay = tdecay_s
        df_soma.loc[counter].dvdt = dvdt_s
        df_soma.loc[counter].thalf = thalf_s
        df_soma.loc[counter].latency = latency_s
        df_soma.loc[counter].time_peak = time_peak_s

        # Dendritic recordings
        df_dend.loc[counter].peak = peak_d
        df_dend.loc[counter].time_rise = trise_d
        df_dend.loc[counter].time_decay = tdecay_d
        df_dend.loc[counter].dvdt = dvdt_d
        df_dend.loc[counter].thalf = thalf_d
        df_dend.loc[counter].latency = latency_d
        df_dend.loc[counter].time_peak = time_peak_d

        if inh == 1:
            df_soma.loc[counter].case = "-Inh"
            df_dend.loc[counter].case = "-Inh"
        elif inh == 2:
            df_soma.loc[counter].case = "Control"
            df_dend.loc[counter].case = "Control"
        elif inh == 3:
            df_soma.loc[counter].case = "-Somatic Inh"
            df_dend.loc[counter].case = "-Somatic Inh"
        elif inh == 4:
            df_soma.loc[counter].case = "-Prox Dend Inh"
            df_dend.loc[counter].case = "-Prox Dend Inh"
        elif inh == 5:
            df_soma.loc[counter].case = "-Dist Dend Inh"
            df_dend.loc[counter].case = "-Dist Dend Inh"
            
        counter += 1

        voltages_s.append(v_soma)
        voltages_d.append(v_dend)
        times.append(time)

    if inh == 1:
        name = "no_inh"
    elif inh == 2:
        name = "control"
    elif inh == 3:
        name = "no_soma_inh"
    elif inh == 4:
        name = "no_prox_dend_inh"
    elif inh == 5:
        name = "no_dist_dend_inh"
        
    results[f"{name}_voltage_soma"] = voltages_s
    results[f"{name}_voltage_dend"] = voltages_d
    results[f"{name}_time"] = times


## Filter out data with dSpike using the Control condition
# ctrl = df_dend[df_dend['case']=='Control']
# idx_psp = ctrl.index[ctrl["thalf"] >= 3].tolist()
# idx_dspike = ctrl.index[ctrl["thalf"] < 3].tolist()

# # Expand in all cases
# idx_psp_all = []
# idx_dspike_all = []
# for i in range(len(inh_deletions)):
#     idx_psp_all += [j + trials*i for j in idx_psp]
#     idx_dspike_all += [j + trials*i for j in idx_dspike]

simdur
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


all_dict = {}
all_dict['results'] = results
all_dict['df_dend'] = df_dend
all_dict['df_dend3'] = df_dend3
all_dict['df_dend4'] = df_dend4
all_dict['df_dend5'] = df_dend5
all_dict['df_dend6'] = df_dend6
all_dict['df_dend7'] = df_dend7
all_dict['df_soma'] = df_soma
all_dict['df_soma3'] = df_soma3

# Save the results
dirname = "data_psp_single_cell/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + f'single_cell_final_new_replication_{rep}.pkl', 'wb') as handle:
    pickle.dump(all_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
