#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 08:44:29 2021.

@author: spiros
"""

import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d


def synaptic_metrics(voltage, time, onset=700):
    """
    .

    Parameters
    ----------
    voltage : np.ndarray
        Membrane voltage in mV.
    time : np.ndarray
        Time vector in ms.

    Returns
    -------
    peak : float
        Peak amplitude in mV.
    trise : float
        Rise time in ms.
    thalf : float
        Hlaf-width in ms.
    tdecay : float
        Decay time in ms.
    derivative : float
        Derivative of the membrane voltage w.r.t time (mV/ms).
    latency : float
        Time from onset to 5% of voltage in ms.
    time_peak : float
        Time to rach the peak in ms.

    """
    # Find the peak voltage
    voltage2 = voltage.copy()
    # voltage2[voltage2 < 0] = 0
    peak = np.max(voltage2)

    # Rise time
    # Find the index of the peak
    idx_max = np.argmax(voltage2)
    if idx_max == 0:
        peak, trise, thalf, tdecay, derivative, latency, time_peak = 0, 0, 0, 0, 0, 0 ,0
    else:
        # Find the index of the 10% peak value (rise phase)
        idx_10 = (np.abs(voltage2[:idx_max] - 0.1 * peak)).argmin()
        # Find the index of the 90% peak value (rise phase)
        idx_90 = (np.abs(voltage2[:idx_max] - 0.9 * peak)).argmin()
        # Calculate the tau rise
        trise = time[idx_90] - time[idx_10]

        # Half-width
        # Find the index of the 50% of peak (rise phase)
        idx_50minus = (np.abs(voltage2[:idx_max] - 0.5 * peak)).argmin()
        # Find the index of the 50% of peak (decay phase)
        idx_50plus = (np.abs(voltage2[idx_max:] - 0.5 * peak)).argmin() + idx_max
        # Calculate the half-width
        thalf = time[idx_50plus] - time[idx_50minus]

        # Decay time
        tdecay = time[idx_50plus] - time[idx_max]

        # Latency
        # Find the index of the 5% of the peak (rise phase)
        idx_5 = (np.abs(voltage2[:idx_max] - 0.05 * peak)).argmin()
        idx_onset = (np.abs(time - onset)).argmin()
        latency = time[idx_5] - time[idx_onset] 

        # Calculate the derivative of the signal (dv/dt -- mV/ms)
        if (idx_max - idx_5) > 2:
            time_diff = np.gradient(time[idx_5:idx_max])
            # time_diff = gaussian_filter1d(time_diff, sigma=5)
            voltage_diff = np.gradient(voltage2[idx_5:idx_max])
            # voltage_diff = gaussian_filter1d(voltage_diff, sigma=5)
            voltage_diff[time_diff < 1e-5] = np.nan
            time_diff[time_diff < 1e-5] = np.nan
            derivatives = voltage_diff / time_diff
            derivative = np.nanmax(derivatives)
        else:
            derivative = 0.0

        # Time of peak
        time_peak = time[idx_max] - time[idx_onset]
    
    # check for negative values and remove them from the analysis.
    my_list = np.array([peak, trise, thalf, tdecay, derivative, latency, time_peak])
    flag = sum(my_list<0)
    if flag > 0:
         peak, trise, thalf, tdecay, derivative, latency, time_peak = 0, 0, 0, 0, 0, 0 ,0
    

    return peak, trise, thalf, tdecay, derivative, latency, time_peak
