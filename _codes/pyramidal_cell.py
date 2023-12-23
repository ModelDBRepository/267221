# -*- coding: utf-8 -*-
"""
Spyder Editor.

This is a temporary script file.
"""

from neuron import h
import matplotlib.pyplot as plt
import numpy as np
from cell_models import PyramidalCell
from scipy.signal import find_peaks
from opt import myinit, mystyle


# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())


# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")

# Store variables
peakss = []
sag_all = []
Rin_all = []

# plot_all, plot_single = True, False
plot_all, plot_single = False, True
soma_stim = True

if soma_stim:
    start = 0.0
    stop = 0.299
    step = 0.025
else:
    start = 0.0
    stop = 1.001
    step = 0.025

# for icur in [-0.4, -0.2, 0.0, 0.400, 0.800]:
# for icur in np.flip(np.arange(start=-0.400, stop=-0.020, step=0.025)):
# for icur in np.arange(start=start, stop=stop, step=step):
for icur in [0.2]:
    # Create a Pyramidal Cell instance
    newpyr = PyramidalCell(0)
    # =========================================================================
    # SAVE VECTORS
    # =========================================================================
    soma_v_vec = h.Vector()  # Membrane potential vector
    radTdist2_v_vec = h.Vector()
    radTdist1_v_vec = h.Vector()
    radTmed_v_vec = h.Vector()
    radTprox_v_vec = h.Vector()
    axon_v_vec = h.Vector()
    t_vec = h.Vector()  # Time stamp vector
    soma_ik_kca_vec = h.Vector()
    soma_ik_vec = h.Vector()
    soma_ina_vec = h.Vector()

    soma_cai_vec = h.Vector()
    soma_cao_vec = h.Vector()
    soma_ca_vec = h.Vector()

    soma_ik_kca_vec.record(newpyr.soma(0.5)._ref_ik_kca)
    soma_ik_vec.record(newpyr.soma(0.5)._ref_ik)
    soma_ina_vec.record(newpyr.soma(0.5)._ref_ina)

    soma_cai_vec.record(newpyr.soma(0.5)._ref_cai)
    soma_ca_vec.record(newpyr.soma(0.5)._ref_ica)

    soma_v_vec.record(newpyr.soma(0.5)._ref_v)
    radTdist2_v_vec.record(newpyr.radTdist2(0.5)._ref_v)
    radTdist1_v_vec.record(newpyr.radTdist1(0.5)._ref_v)
    radTmed_v_vec.record(newpyr.radTmed(0.5)._ref_v)
    radTprox_v_vec.record(newpyr.radTprox(0.5)._ref_v)
    axon_v_vec.record(newpyr.axon(0.5)._ref_v)
    t_vec.record(h._ref_t)

    # =========================================================================
    # STIMULATION
    # =========================================================================
    simdur = 2000.0
    if soma_stim:
        stim = h.IClamp(newpyr.soma(0.5))
    else:
        stim = h.IClamp(newpyr.radTdist1(0.1))

    stim.delay = 400
    stim.dur = 1000
    stim.amp = icur  # 0.25 #nA

    stim2 = h.IClamp(newpyr.soma(0.5))
    stim2.delay = 0
    stim2.dur = 2000
    stim2.amp = -0.044  # 0.25 #nA

    # h.v_init = -70
    myinit()
    h.tstop = simdur
    h.dt = 0.025
    h.run()

    t_vec = np.array(t_vec)
    soma_v_vec = np.array(soma_v_vec)
    radTprox_v_vec = np.array(radTprox_v_vec)
    radTmed_v_vec = np.array(radTmed_v_vec)
    radTdist1_v_vec = np.array(radTdist1_v_vec)
    radTdist2_v_vec = np.array(radTdist2_v_vec)
    axon_v_vec = np.array(axon_v_vec)
    n1 = np.abs(t_vec - 300).argmin()
    n2 = np.abs(t_vec - 1900).argmin()

    if plot_single:

        if soma_stim:
            plt.figure(num=10, figsize=(10, 8))
            plt.plot(
                t_vec[n1:n2], soma_v_vec[n1:n2], linewidth=3, label=f"I={1e3*icur} pA"
            )
            plt.legend()
        else:
            plt.figure(num=10, figsize=(10, 8))
            plt.plot(
                t_vec[n1:n2],
                radTdist1_v_vec[n1:n2],
                linewidth=3,
                label=f"I={1e3*icur} pA",
            )
            plt.legend()

    if stim.amp >= 0:
        if soma_stim:
            V = np.array(soma_v_vec)
            maxima = -10
        else:
            V = np.array(radTdist1_v_vec)
            maxima = -25

        peaks, _ = find_peaks(V, height=maxima)
        print(f"Number of spikes: {len(peaks)}")
        peakss.append(len(peaks))

        if len(peaks) == 0:
            if soma_stim:
                V = np.array(soma_v_vec)
            else:
                V = np.array(radTdist1_v_vec)
            time_vec = np.array(t_vec)
            n = np.abs(t_vec - 1400).argmin()
            Vtau = V[n:]
            ttau = time_vec[n:]
            ttau -= ttau[0]
            Vtau = (Vtau - np.min(Vtau)) / (np.max(Vtau) - np.min(Vtau))
            Vt = 1 / np.exp(1)

            tau = ttau[np.abs(Vtau - Vt).argmin()]
            print(f"tau: {tau}ms")

    elif stim.amp < 0:
        t_vec = np.array(t_vec)
        n1 = np.abs(t_vec - 350).argmin()
        n2 = np.abs(t_vec - 1450).argmin()
        if soma_stim:
            V = np.array(soma_v_vec)
        else:
            V = np.array(radTdist1_v_vec)
        Vopt = V[n1:n2]

        ninit = np.abs(t_vec - 399).argmin()
        Vrest = V[ninit]
        Vpeak = np.min(Vopt)
        nss = np.abs(t_vec - 1390).argmin()
        Vss = V[nss]

        sag = (Vpeak - Vss) / Vpeak * 100
        sag_all.append(sag)
        Rin = (Vss - Vrest) / ((stim.amp))
        Rin_all.append(Rin)

        time_vec = np.array(t_vec)
        n = np.abs(t_vec - 1400).argmin()
        Vtau = V[n:]
        ttau = time_vec[n:]
        ttau -= ttau[0]
        Vtau = (Vtau - np.min(Vtau)) / (np.max(Vtau) - np.min(Vtau))
        Vt = 1 - 1 / np.exp(1)

        tau = ttau[np.abs(Vtau - Vt).argmin()]

        print(f"Rin: {Rin}MOhms; sag ratio: {sag} (%); tau: {tau}ms")

plt.xlabel('time (ms)')
plt.ylabel('voltage (mV)')
if plot_all:
    if stim.amp > 0:
        # Visualize model history
        plt.figure(figsize=(16, 10))
        Inj = np.arange(start=start, stop=stop, step=step) * 1e3
        plt.plot(Inj, peakss, linewidth=3, markersize=15,
                 marker=".", linestyle="solid")
        if soma_stim:
            plt.xlabel("Somatic current injection (pA)")
        else:
            plt.xlabel("Dendritic current injection (pA)")
        plt.ylabel("Number of spikes")
        plt.title("Pyramidal neurons")
        plt.ylim([0, 20])
        plt.yticks([0, 10, 20, 30, 40], ["0", "10", "20", "30", "40"])

    elif stim.amp < 0:
        # Visualize model history
        plt.figure(figsize=(16, 10))
        Inj = np.flip(np.arange(start=-0.400, stop=-0.020, step=0.025)) * 1e3
        plt.plot(
            Inj, sag_all, linewidth=3, markersize=15, marker=".", linestyle="solid"
        )
        if soma_stim:
            plt.xlabel("Somatic current injection (pA)")
        else:
            plt.xlabel("Dendritic current injection (pA)")
        plt.ylabel("Sag Ratio (%)")
        plt.title("Pyramidal neurons")
        plt.ylim([0, 15])
        plt.yticks([0, 5, 10, 15], ["0", "5", "10", "15"])
        plt.xlim(0, min(Inj))

        # Rin_all = [1000*x for x in Rin_all] # in MOhms

        plt.figure(figsize=(16, 10))
        plt.boxplot(Rin_all)
        plt.xlabel("input resistance")
        plt.ylabel("MegaOhms")
