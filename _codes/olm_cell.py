from neuron import h
import matplotlib.pyplot as plt
import numpy as np

from scipy.signal import find_peaks
from cell_models import OLMCell
from opt import myinit, mystyle


# Make plots nicer!
plt.style.use("seaborn-colorblind")
plt.rcParams.update(mystyle())

cellname = 'olm'
# If we had not included gui in the list of things to import
h.load_file("stdrun.hoc")

# Store variables
peakss = []
sag_all = []
Rin_all = []

# plot_all, plot_single = True, False
plot_all, plot_single = False, True

# positive currents!
start = 0.0
stop = 0.451
step = 0.025

# for icur in np.flip(np.arange(start=-0.250, stop=-0.020, step=0.025)):
# for icur in np.arange(start=start, stop=stop, step=step):
# for icur in [-0.2, -0.10, 0.1, 0.250]:
for icur in [-.25]:

    # Create an OLM Cell instance
    cell = OLMCell(0)

    # =========================================================================
    # SAVE VECTORS
    # =========================================================================
    soma_v_vec = h.Vector()  # Membrane potential vector

    dend_v_vec = h.Vector()
    t_vec = h.Vector()  # Time stamp vector

    soma_v_vec.record(cell.soma(0.5)._ref_v)
    t_vec.record(h._ref_t)

    # =========================================================================
    # STIMULATION
    # =========================================================================
    simdur = 2000.0
    stim = h.IClamp(cell.soma(0.5))

    stim.delay = 400
    stim.dur = 1000
    stim.amp = icur  # 0.25 # nA

    stim2 = h.IClamp(cell.soma(0.5))
    stim2.delay = 0
    stim2.dur = 2000
    stim2.amp = -0.041  # 0.25 #nA

    myinit()
    h.continuerun(simdur)

    t_vec = np.array(t_vec)
    soma_v_vec = np.array(soma_v_vec)
    n1 = np.abs(t_vec - 300).argmin()
    n2 = np.abs(t_vec - 1900).argmin()

    if plot_single:

        plt.figure(num=1, figsize=(10, 8))
        plt.plot(t_vec[n1:n2], soma_v_vec[n1:n2],
                 linewidth=3, label=f"I={1e3*icur} pA")
        plt.legend()

    if stim.amp >= 0:
        V = np.array(soma_v_vec)
        maxima = -10

        peaks, _ = find_peaks(V, height=maxima)
        print(f"Number of spikes: {len(peaks)}")
        peakss.append(len(peaks))

        if len(peaks) == 0:
            V = np.array(soma_v_vec)

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
        V = np.array(soma_v_vec)

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
    plt.show()

if plot_all:
    if stim.amp > 0:

        plt.figure(figsize=(16, 10))
        Inj = np.arange(start=start, stop=stop, step=step) * 1e3
        plt.plot(Inj, peakss, linewidth=3, markersize=15,
                 marker=".", linestyle="solid")
        plt.xlabel("Electrical stimulation (pA)")
        plt.ylabel("Number of spikes")
        plt.title(f"{cellname} neurons")
        plt.ylim([0, 100])
        plt.yticks([0, 25, 50, 75, 100], ["0", "25", "50", "75", "100"])
        plt.show()

    elif stim.amp < 0:

        plt.figure(figsize=(16, 10))
        Inj = np.flip(np.arange(start=-0.250, stop=-0.020, step=0.025)) * 1e3
        plt.plot(Inj, sag_all, linewidth=3, markersize=15,
                 marker=".", linestyle="solid")
        plt.xlabel("Electrical stimulation (pA)")

        plt.ylabel("Sag Ratio (%)")
        plt.title(f"{cellname} neurons")
        plt.ylim([0, 15])
        plt.yticks([0, 5, 10, 15], ["0", "5", "10", "15"])
        plt.xlim(0, min(Inj))
        plt.show()

        # Rin_all = [1000*x for x in Rin_all] # in MOhms

        plt.figure(figsize=(16, 10))
        plt.boxplot(Rin_all)
        plt.xlabel("input resistance")
        plt.ylabel("MegaOhms")
        plt.xlabel("input resistance")
        plt.ylabel("MegaOhms")
        plt.show()

