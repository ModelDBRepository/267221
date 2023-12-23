#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 10 16:01:45 2022

@author: spiros
"""

from neuron import h


def myinit(vinit=-65, temp=34, cao=2, cai=5e-5, dt=0.025,
           enable=True, modify=False):
    """
    Initialize function.

    Parameters
    ----------
    vinit: float, optional
        Voltage initialization in mV. The default is -65.
    temp: float, optional
        Tmperature in degC. The default is 34.
    cao: float, optional
        Extracellular Ca2+ concentration in mM. The default is 2.
    cai: float, optional
        Intracellular Ca2+ concentration in mM. The default is 5e-5.
    dt : float, optional
        Integration time-step in ms. The default is 0.1.
    enable: bool, optional
        Enable CVODE for faster simulations. The default is `True`.
    modify: bool, optional
        Enable `long_double` and `daspk`. The default is `False`.

    Returns
    -------
    None.

    """
    # Enable CVODE for faster Simulations
    # For faster Simulations
    cvode = h.CVode()
    cvode.active(enable)
    if modify:
        # x = cvode.atol(1e-3)
        cvode.use_long_double(True)
        cvode.use_daspk(True)
    # new code to happen after initialization here
    print('Initializing NEURON...')
    # only need the following if states have been changed
    h.celsius = temp
    h.cao0_ca_ion = cao
    h.cai0_ca_ion = cai
    h.dt = dt
    h.stdinit()
    h.init()
    h.finitialize()
    if h.cvode.active():
        h.cvode.re_init()
    else:
        h.fcurrent()
    h.frecord_init()
    print(f'Initialization has been completed. Temperature is: {temp} degC,\n'
          f'[Ca2+] extracellular is: {cao} mM, and intracellular is: {cai} mM')


def mystyle():
    """
    Create custom plotting style.

    Returns
    -------
    my_style : dict
        Dictionary with matplotlib parameters.

    """
    # color pallette
    my_style = {
        # Use LaTeX to write all text
        "text.usetex": False,
        "font.family": "Times New Roman",
        "font.weight": "bold",
        # Use 16pt font in plots, to match 16pt font in document
        "axes.labelsize": 16,
        "axes.titlesize": 20,
        "font.size": 16,
        # Make the legend/label fonts a little smaller
        "legend.fontsize": 14,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "axes.linewidth": 2.5,
        "lines.markersize": 10.0,
        "lines.linewidth": 2.5,
        "xtick.major.width": 2.2,
        "ytick.major.width": 2.2,
        "axes.labelweight": "bold",
        "axes.spines.right": False,
        "axes.spines.top": False
    }

    return my_style


def set_size(width, fraction=1):
    """
    Set figure dimensions to avoid scaling in LaTeX.

    Parameters
    ----------
    width: float
        Document textwidth or columnwidth in pts
    fraction: float, optional
        Fraction of the width which you wish the figure to occupy
    Returns
    -------
    fig_dim: tuple
        Dimensions of figure in inches
    """
    # Width of figure (in pts)
    fig_width_pt = width * fraction
    # Convert from pt to inches
    inches_per_pt = 1 / 72.27
    # Golden ratio to set aesthetic figure height
    # https://disq.us/p/2940ij3
    golden_ratio = (5**.5 - 1) / 2
    # Figure width in inches
    fig_width_in = fig_width_pt * inches_per_pt
    # Figure height in inches
    fig_height_in = fig_width_in * golden_ratio
    fig_dim = (fig_width_in, fig_height_in)
    return fig_dim