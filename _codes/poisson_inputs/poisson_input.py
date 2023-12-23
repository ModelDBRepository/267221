import numpy as np
import os, sys, brian2

nrun = int(sys.argv[1])
rate = int(sys.argv[2])

brian2.seed(nrun)
print ('RUN: ' + str(nrun))

foldername = 'rate'+str(rate)+'/run_'+str(nrun)
os.system('mkdir -p -v '+foldername)

N = 10000  # number of neurons
time_input = 4000 * brian2.ms  # simulation time
P = brian2.PoissonGroup(N, rates=rate * brian2.Hz)
S = brian2.SpikeMonitor(P)

brian2.run(time_input, report='text', report_period = 10 * brian2.second)

fname = 'noise_'
delay = 400  # ms

theta_freq = 8  # Hz
theta_phase = 0
noise = 0.2
for s in range(len(S.spike_trains().keys())):
    spiketimes = np.round(1000*S.spike_trains()[s]/brian2.second, 1)
    # Make the spike train theta modulated
    spikes = []
    for spike in spiketimes:
        probability = (np.sin(2.0*np.pi*theta_freq*spike/1000. + theta_phase) + 1.0)/2.0
        r_ = np.random.rand()
        if ((probability > 0.7) and (r_ < probability / 2.0)) or np.random.rand() < noise:
            spikes.append(spike)

    spikes = delay + np.array(spikes)
    np.savetxt(f"{foldername}/{fname}{s}.txt", spikes,
               fmt='%10.1f', newline='\n')