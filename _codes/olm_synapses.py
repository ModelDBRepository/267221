import os
import pickle
from neuron import h
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from synaptic_metrics import synaptic_metrics

from cell_models import OLMCell
from opt import myinit, mystyle


# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())

# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")
replications = 10
trials = 200

# Variable to store the results
cols = ["peak", "time_rise", "time_decay", "dvdt", "latency", "thalf",
        "replication"]
df_soma = pd.DataFrame(columns=cols, index=range(trials*replications))
counter = 0
synsE = []
results = {}

for nrep in range(replications):
    print(f"\nReplication of the experiment: {nrep}\n")
    voltages_s, times = [], []
    for ntrial in range(trials):

        np.random.seed(1000*nrep + ntrial)
        r1 = np.random.poisson(lam=0.3, size=1).item()
        synsE.append(r1)
        
        print(f"Trial {ntrial}...Esyn: {r1}")

        # plot_all, plot_single = True, False
        plot_single = False
        soma_stim = True

        # Create an OLM Cell instance
        cell = OLMCell(0)

        dends = []
        for sec in cell.all:
            if "dend" in str(sec):
                dends.append(sec(0.5))

        # Synapses lists
        synAMPA1, vsAMPA1 = [], []
        loc1, lec, inputs = [], [], []

        gAMPA = 0.00091266

        for i in range(r1):

            spiketimes = [700.0]
            inputs.append(h.Vector(spiketimes))
            lec.append(h.VecStim())
            lec[-1].play(inputs[-1])

            # Choose a location at random @ SLM
            loc1.append(np.random.randint(low=0, high=len(dends)))

            # New Synapse and store in a list
            # AMPA synapse
            synAMPA1.append(h.Exp2Syn(dends[loc1[-1]]))
            synAMPA1[-1].e = 0  # reversal potential
            synAMPA1[-1].tau1 = 0.5  # rise time
            synAMPA1[-1].tau2 = 8.0  # decay time

            # New VecStim object
            vsAMPA1.append(h.NetCon(lec[-1], synAMPA1[-1]))
            vsAMPA1[-1].delay = 0.0  # delay in ms
            vsAMPA1[-1].weight[0] = gAMPA  #

        # =========================================================================
        # SAVE VECTORS
        # =========================================================================
        soma_v_vec = h.Vector()  # Membrane potential vector
        soma_v_vec.record(cell.soma(0.5)._ref_v)

        t_vec = h.Vector()  # Time stamp vector
        t_vec.record(h._ref_t)
        # =========================================================================
        # STIMULATION
        # =========================================================================
        stim2 = h.IClamp(cell.soma(0.5))
        stim2.delay = 0
        stim2.dur = 2000
        stim2.amp = -0.0315
        simdur = 2000.0
        myinit(vinit=-70)
        h.continuerun(simdur)

        # =========================================================================
        # SAVE and PLOT the results
        # =========================================================================
        t_vec = np.array(t_vec)
        soma_v_vec = np.array(soma_v_vec)

        n1 = np.abs(t_vec - 690).argmin()
        n2 = np.abs(t_vec - 900).argmin()

        v_soma = soma_v_vec[n1:n2] - soma_v_vec[n1]
        time = t_vec[n1:n2]

        if r1 > 0:
            peak_s, trise_s, thalf_s, tdecay_s, dvdt_s, latency_s, time_peak = synaptic_metrics(
                v_soma, time
            )
        else:
            peak_s, trise_s, tdecay_s, dvdt_s, thalf_s, latency_s, time_peak = 0, 0, 0, 0, 0, 0, 0

        df_soma.loc[counter].peak = peak_s
        df_soma.loc[counter].time_rise = trise_s
        df_soma.loc[counter].time_decay = tdecay_s
        df_soma.loc[counter].dvdt = dvdt_s
        df_soma.loc[counter].thalf = thalf_s
        df_soma.loc[counter].latency = latency_s

        df_soma.loc[counter].replication = nrep

        counter += 1

        voltages_s.append(v_soma)
        times.append(time)

    results[f"Replication_{nrep}_voltage_soma"] = voltages_s
    results[f"Replication_{nrep}_time"] = times


# plt.figure()
# sns.boxplot(data=df_soma, y="peak")
# # sns.swarmplot(data=df_soma, y="peak")
# plt.ylabel('voltage [ms]')
# plt.xticks([])
# plt.title("LEC driven response on SST+ cells")
# plt.show()

plt.figure()
samples = np.random.choice(range(trials), size=10, replace=False)
nreps = np.random.choice(range(replications), size=10, replace=True)
for i in range(10):
    ntr = samples[i]
    nrep = nreps[i]
    v1s = results[f"Replication_{nrep}_voltage_soma"][ntr]
    t1 = results[f"Replication_{nrep}_time"][ntr]

    plt.plot(t1, v1s, label="control", color="black", linewidth=1.5)
    plt.xlabel("time (ms)")
    plt.ylabel("voltage (mV)")
    plt.title("SST+ cells")


if plot_single:
    plt.figure(figsize=(8, 6))
    plt.plot(t_vec[n1:n2],
             soma_v_vec[n1:n2] - soma_v_vec[n1],
             linewidth=3, label="soma")
    plt.ylabel("Membrane Voltage (mV)")
    ymin = np.min(soma_v_vec[n1:n2] - soma_v_vec[n1]) - 2
    ymax = np.max(soma_v_vec[n1:n2] - soma_v_vec[n1]) + 2
    plt.ylim([ymin, ymax])
    plt.xlabel("Time (ms)")
    plt.legend()
    plt.show()

dirname = "data_interneurons_synapses/"
if not os.path.exists(dirname):
    os.mkdir(dirname)

with open(dirname + 'olm.pkl', 'wb') as handle:
    pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)
