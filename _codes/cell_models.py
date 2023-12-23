#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 09:19:08 2022

@author: spiros
"""

from neuron import h
from neuron.units import um

h.load_file("stdrun.hoc")

class PyramidalCell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_topology()
        self._create_lists()
        self._setup_geometry()
        self._setup_passive()
        self._setup_segments()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.radTprox = h.Section(name='radTprox', cell=self)
        self.radTmed = h.Section(name='radTmed', cell=self)
        self.radTdist1 = h.Section(name='radTdist1', cell=self)
        self.radTdist2 = h.Section(name='radTdist2', cell=self)
        self.radTdist3 = h.Section(name='radTdist3', cell=self)
        self.lm_thick1 = h.Section(name='lm_thick1', cell=self)
        self.lm_medium1 = h.Section(name='lm_medium1', cell=self)
        self.lm_thin1a = h.Section(name='lm_thin1a', cell=self)
        self.lm_thin1b = h.Section(name='lm_thin1b', cell=self)
        self.lm_thick2 = h.Section(name='lm_thick2', cell=self)
        self.lm_medium2 = h.Section(name='lm_medium2', cell=self)
        self.lm_thin2a = h.Section(name='lm_thin2a', cell=self)
        self.lm_thin2b = h.Section(name='lm_thin2b', cell=self)
        self.rad_thick1 = h.Section(name='rad_thick1', cell=self)
        self.rad_medium1 = h.Section(name='rad_medium1', cell=self)
        self.rad_thin1a = h.Section(name='rad_thin1a', cell=self)
        self.rad_thin1b = h.Section(name='rad_thin1b', cell=self)
        self.rad_thick2 = h.Section(name='rad_thick2', cell=self)
        self.rad_medium2 = h.Section(name='rad_medium2', cell=self)
        self.rad_thin2a = h.Section(name='rad_thin2a', cell=self)
        self.rad_thin2b = h.Section(name='rad_thin2b', cell=self)
        self.oriprox1 = h.Section(name='oriprox1', cell=self)
        self.oridist1a = h.Section(name='oridist1a', cell=self)
        self.oridist1b = h.Section(name='oridist1b', cell=self)
        self.oriprox2 = h.Section(name='oriprox2', cell=self)
        self.oridist2a = h.Section(name='oridist2a', cell=self)
        self.oridist2b = h.Section(name='oridist2b', cell=self)
        self.axon = h.Section(name='axon', cell=self)

    def _create_lists(self):
        self.all = self.soma.wholetree()
        self.ori = [sec for sec in self.all if sec.name().__contains__('ori')]
        self.rad = [sec for sec in self.all if sec.name().__contains__('rad_')]
        self.slm = [sec for sec in self.all if sec.name().__contains__('lm_')]
        self.trunk = [sec for sec in self.all if sec.name().__contains__('radT')]

    def _setup_topology(self):
        # Connect sections
        # Apical trunk
        self.radTprox.connect(self.soma(1))
        self.radTmed.connect(self.radTprox(1))
        self.radTdist1.connect(self.radTmed(1))
        self.radTdist2.connect(self.radTdist1(1))
        self.radTdist3.connect(self.radTdist2(1))
        # Apical oblique tree
        # Right
        self.rad_thick1.connect(self.radTmed(1))
        self.rad_medium1.connect(self.rad_thick1(1))
        self.rad_thin1a.connect(self.rad_medium1(1))
        self.rad_thin1b.connect(self.rad_medium1(1))
        # Left
        self.rad_thick2.connect(self.radTmed(1))
        self.rad_medium2.connect(self.rad_thick2(1))
        self.rad_thin2a.connect(self.rad_medium2(1))
        self.rad_thin2b.connect(self.rad_medium2(1))
        # Apical tuft tree
        # Right
        self.lm_thick1.connect(self.radTdist3(1))
        self.lm_medium1.connect(self.lm_thick1(1))
        self.lm_thin1a.connect(self.lm_medium1(1))
        self.lm_thin1b.connect(self.lm_medium1(1))
        # Left
        self.lm_thick2.connect(self.radTdist3(1))
        self.lm_medium2.connect(self.lm_thick2(1))
        self.lm_thin2a.connect(self.lm_medium2(1))
        self.lm_thin2b.connect(self.lm_medium2(1))
        # Basal tree
        # Right
        self.oriprox1.connect(self.soma(0))
        self.oridist1a.connect(self.oriprox1(1))
        self.oridist1b.connect(self.oriprox1(1))
        # Left
        self.oriprox2.connect(self.soma(0))
        self.oridist2a.connect(self.oriprox2(1))
        self.oridist2b.connect(self.oriprox2(1))
        # Axon
        self.axon.connect(self.soma(0))

    def _setup_geometry(self):
        self.soma.L = self.soma.diam = 10 * um
        self.radTprox.L, self.radTprox.diam = 80 * um, 3 * um
        self.radTmed.L, self.radTmed.diam = 70 * um, 2.5 * um
        self.radTdist1.L, self.radTdist1.diam = 50 * um, 2 * um
        self.radTdist2.L, self.radTdist2.diam = 50 * um, 1.75 * um
        self.radTdist3.L, self.radTdist3.diam = 50 * um, 1.5 * um
        self.lm_thick1.L, self.lm_thick1.diam = 50 * um, 0.95 * um
        self.lm_medium1.L, self.lm_medium1.diam = 50 * um, 0.9 * um
        self.lm_thin1a.L, self.lm_thin1a.diam = 75 * um, 0.6 * um
        self.lm_thin1b.L, self.lm_thin1b.diam = 75 * um, 0.6 * um
        self.lm_thick2.L, self.lm_thick2.diam = 50 * um, 0.95 * um
        self.lm_medium2.L, self.lm_medium2.diam = 50 * um, 0.90 * um
        self.lm_thin2a.L, self.lm_thin2a.diam = 75 * um, 0.6 * um
        self.lm_thin2b.L, self.lm_thin2b.diam = 75 * um, 0.6 * um
        self.rad_thick1.L, self.rad_thick1.diam = 50 * um, 1.9 * um
        self.rad_medium1.L, self.rad_medium1.diam = 50 * um, 1.5 * um
        self.rad_thin1a.L, self.rad_thin1a.diam = 75 * um, 0.95 * um
        self.rad_thin1b.L, self.rad_thin1b.diam = 75 * um, 0.95 * um
        self.rad_thick2.L, self.rad_thick2.diam = 50 * um, 1.9 * um
        self.rad_medium2.L, self.rad_medium2.diam = 50 * um, 1.5 * um
        self.rad_thin2a.L, self.rad_thin2a.diam = 75 * um, 0.95 * um
        self.rad_thin2b.L, self.rad_thin2b.diam = 75 * um, 0.95 * um
        self.oriprox1.L, self.oriprox1.diam = 75 * um, 1.6 * um
        self.oridist1a.L, self.oridist1a.diam = 75 * um, 1.1 * um
        self.oridist1b.L, self.oridist1b.diam = 75 * um, 1.1 * um
        self.oriprox2.L, self.oriprox2.diam = 75 * um, 1.6 * um
        self.oridist2a.L, self.oridist2a.diam = 75 * um, 1.1 * um
        self.oridist2b.L, self.oridist2b.diam = 75 * um, 1.1 * um
        self.axon.L, self.axon.diam = 150 * um, 1 * um

    def _setup_passive(self):
        Ra_soma = 300.0  # Axial resistance in Ohm * cm
        # Soma
        self.soma.Ra = Ra_soma
        self.soma.cm = 1  # Membrane capacitance in uF/cm2
        # Axon
        self.axon.Ra = Ra_soma
        self.axon.cm = 1  # Membrane capacitance in uF/cm2
        # Apical trunk
        Ra_trunk = [j*Ra_soma for j in [0.95, 0.90, 0.85, 0.75, 0.72]]
        for i, sec in enumerate(self.trunk):
            sec.Ra = Ra_trunk[i]
            sec.cm = 2.2  # Membrane capacitance in uF/cm2
        # Apical tuft
        for i, sec in enumerate(self.slm):
            sec.Ra = 0.5*Ra_soma
            sec.cm = 3.0  # Membrane capacitance in uF/cm2
        # Oblique dendrtites
        for i, sec in enumerate(self.rad):
            sec.Ra = 0.5*Ra_soma
            sec.cm = 3.0  # Membrane capacitance in uF/cm2
        # Basal dendrites
        for i, sec in enumerate(self.ori):
            sec.Ra = Ra_soma
            sec.cm = 3.0  # Membrane capacitance in uF/cm2


    def _setup_segments(self):
        # Create segments based on `lambda_f`
        d_lambda = 0.1
        frequency = 100
        for sec in self.all:
            sec.nseg = int((sec.L/(d_lambda*h.lambda_f(frequency, sec=sec)) + 0.9) / 2)*2 + 1


    def _setup_biophysics(self):

        # =====================================================================
        # General Parameters
        gna_hh_dend = 0.007  # sodium (Na+) dend channel
        gk_hh_dend  = 0.007/8.065  # Delayed rectified potassium (K+) channel
        gka_soma = 0.0075  # A-type K @ soma
        gka_dend = 0.048672  # A-type @ dends
        gh_soma = 0.000015  # for sag ratio
        gkm_trunk = 0.0006  # K m-type @ trunk
        gCaL_trunk = .1*0.000316  # Ca L-type @ trunk
        gkca_trunk = 0.0005  # slow AHP K+ current
        Rm_soma = 200000  # ohm cm2
        Ra_soma = 300.0  # ohm cm
        # =====================================================================

        # =====================================================================
        # Somatic compartment
        self.soma.Ra = Ra_soma
        self.soma.cm = 1

        # Na+, Kdr, and leak channels
        self.soma.insert('hha2')
        for seg in self.soma:
            seg.hha2.gnabar = 0.035
            seg.hha2.gkbar = 0.0014
            seg.hha2.gl = 0.01*1/Rm_soma
            seg.hha2.el = -70

        # h channels
        self.soma.insert('h')
        for seg in self.soma:
            seg.h.ghbar = gh_soma
            seg.h.vhalf = -90
            seg.h.eh = -10

        # A-type K+ channels
        self.soma.insert('kap')
        for seg in self.soma:
            seg.kap.gkabar = gka_soma

        # m-type K+ channels
        self.soma.insert('km')
        for seg in self.soma:
            seg.km.gbar = 0.0166

        # medium Ca2+-dependent K+ channels (mAHP)
        self.soma.insert('mykca')
        for seg in self.soma:
            seg.mykca.gkbar = 0.09075

        # slow Ca2+-dependent K+ channels (sAHP)
        self.soma.insert('kca')
        for seg in self.soma:
            seg.kca.gbar = gkca_trunk

        # T-type Ca2+ channels
        self.soma.insert('cat')
        for seg in self.soma:
            seg.cat.gcatbar = 0.00005

        # R-type Ca2+ channels
        self.soma.insert('somacar')
        for seg in self.soma:
            seg.somacar.gcabar = 0.003

        # L-type Ca2+ channels
        self.soma.insert('cal')
        for seg in self.soma:
            seg.cal.gcalbar = 0.007

        # Constant current for Vrest calibration
        self.soma.insert('constant')
        for seg in self.soma:
            seg.constant.ic = 0.008
        # =====================================================================

        # =====================================================================
        # Axonal compartment(s)
        Rm_axon = 200000

        # Na+, Kdr, and leak channels
        self.axon.insert('hha2')
        for seg in self.axon:
            seg.hha2.gnabar = 1.5
            seg.hha2.gkbar = 0.1
            seg.hha2.gl = 0.01*1/Rm_axon
            seg.hha2.el = -70

        # m-type K+ channels
        self.axon.insert('km')
        for seg in self.axon:
            seg.km.gbar = 0.000003

        # Constant current for Vrest calibration
        self.axon.insert('constant')
        for seg in self.axon:
            seg.constant.ic = -0.001
        # =====================================================================


        # =====================================================================
        # Apical trunk -- all compartments
        Rm_trunk = [190000, 150000, 35000, 35000, 35000]
        Ra_trunk = [j*Ra_soma for j in [0.95, 0.9, 0.85, 0.75, 0.72]]
        gka_trunk = [2.4*gka_soma] + [j*gka_dend for j in [200/350, 270/350, 340/350, 1]]
        gCaL_scale = [0.1, 0.2, 0.4, 0.6, 0.8]
        gh_trunk = [j*gh_soma for j in [1.2, 4, 6, 6.5, 7]]
        gkca_ = [j*gkca_trunk for j in [1, 1, 0.1, 0.1, 0.1]]
        gcagk = [0.066, 0.066, 0.004125, 0.004125, 0.004125]
        ic_ = [-0.0016, -0.0022, -0.0028, -0.0034, -0.004]
        for i, sec in enumerate(self.trunk):
            sec.Ra = Ra_trunk[i]
            sec.cm = 2.2

            # Na+ and Kdr channels
            sec.insert('hha_old')
            for seg in sec:
                seg.hha_old.gl = 0.01*1/Rm_trunk[i]
                seg.hha_old.gnabar = 2.0*gna_hh_dend
                seg.hha_old.gkbar = 2.0*gk_hh_dend
                seg.hha_old.el = -70
                seg.hha_old.ar2 = 0.9

            # h channels
            sec.insert('h')
            for seg in sec:
                seg.h.ghbar = gh_trunk[i]
                seg.h.vhalf = -90
                seg.h.eh = -10

            # A-type K+ channels
            if sec.name().__contains__('radTprox'):
                sec.insert('kap')
                for seg in sec:
                    seg.kap.gkabar = gka_trunk[i]
            else:
                sec.insert('kad')
                for seg in sec:
                    seg.kad.gkabar = gka_trunk[i]

            # m-type K+ channels
            sec.insert('km')
            for seg in sec:
                seg.km.gbar = gkm_trunk

            # L-type Ca2+ channels
            sec.insert('calH')
            for seg in sec:
                seg.calH.gcalbar = gCaL_scale[i]*gCaL_trunk

            # slow Ca2+-dependent K+ channels (sAHP)
            sec.insert('kca')
            for seg in sec:
                seg.kca.gbar = gkca_[i]

            # medium Ca2+-dependent K+ channels (mAHP)
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = gcagk[i]

            # T-type Ca2+ channel
            sec.insert('cat')
            for seg in sec:
                seg.cat.gcatbar = 0.0004

            # R-type Ca2+ channel
            sec.insert('car')
            for seg in sec:
                seg.car.gcabar = 0.1*0.003

            # Constant current for Vrest calibration
            sec.insert('constant')
            for seg in sec:
                seg.constant.ic = ic_[i]
        # =====================================================================


        # =====================================================================
        # Apical tuft (s.lm.) -- all compartments
        Rm_slm = [12000]*8
        for i, sec in enumerate(self.slm):
            sec.Ra = 0.5*Ra_soma
            sec.cm = 3.0

            # Na+ and Kdr channels
            sec.insert('hha_old')
            for seg in sec:
                seg.hha_old.gl = 0.01*1/Rm_slm[i]
                seg.hha_old.gnabar = 2.2*gna_hh_dend
                seg.hha_old.gkbar = 2.2*gk_hh_dend
                seg.hha_old.el = -70
                seg.hha_old.ar2 = 0.95

            # h channels
            sec.insert('h')
            for seg in sec:
                seg.h.ghbar = 8*gh_soma
                seg.h.vhalf = -90
                seg.h.eh = -10

            # A-type K+ channels
            sec.insert('kad')
            for seg in sec:
                # K channels
                seg.kad.gkabar = gka_dend

            # m-type K+ channel
            sec.insert('km')
            for seg in sec:
                seg.km.gbar = 2*gkm_trunk

            # L-type Ca2+ channels
            sec.insert('calH')
            for seg in sec:
                seg.calH.gcalbar = 0.1*gCaL_trunk

            # medium Ca2+-dependent K+ channel (mAHP)
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = 0.004125

            # slow Ca2+-dependent K+ channel (sAHP)
            sec.insert('kca')
            for seg in sec:
                seg.kca.gbar = 0.1*gkca_trunk

            # Constant current for Vrest calibration
            sec.insert('constant')
            for seg in sec:
                seg.constant.ic = -0.004
        # =====================================================================

        # =====================================================================
        # Apical oblique (s.rad.) -- all compartments
        Rm_rad = [13000, 12000, 12000, 12000]*2
        for i, sec in enumerate(self.rad):
            sec.Ra = 0.5*Ra_soma
            sec.cm = 3.0

            # Na+, Kdr, and leak channels
            sec.insert('hha_old')
            for seg in sec:
                seg.hha_old.gl = 0.01*1/Rm_rad[i]
                seg.hha_old.gnabar = 1.2*gna_hh_dend
                seg.hha_old.gkbar = 1.2*gk_hh_dend
                seg.hha_old.el = -70

            # h channels
            sec.insert('h')
            for seg in sec:
                seg.h.ghbar = 4*gh_soma
                seg.h.vhalf = -80
                seg.h.eh = -10

            # A-type K+ channels
            sec.insert('kad')
            for seg in sec:
                seg.kad.gkabar = gka_dend * 300/350

            # m-type K+ channels
            sec.insert('km')
            for seg in sec:
                seg.km.gbar = 2*gkm_trunk

            # L-type Ca2+ channels
            sec.insert('calH')
            for seg in sec:
                seg.calH.gcalbar = 0.4*gCaL_trunk

            # medium Ca2+-dependent K+ channel (mAHP)
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = 0.004125

            # slow Ca2+-dependent K+ channel (sAHP)
            sec.insert('kca')
            for seg in sec:
                seg.kca.gbar = 0.1*gkca_trunk

            # Constant current for Vrest calibration
            sec.insert('constant')
            for seg in sec:
                seg.constant.ic = -0.0024
        # =====================================================================

        # =====================================================================
        # Basal dendrites (s.o.) -- all compartments
        Rm_ori = [18000]*6
        for i, sec in enumerate(self.ori):
            sec.Ra = Ra_soma
            sec.cm = 3.0

            # Na+, Kdr, and leak channels
            sec.insert('hha_old')
            for seg in sec:
                seg.hha_old.gl = 0.01*1/Rm_ori[i]
                seg.hha_old.gnabar = gna_hh_dend
                seg.hha_old.gkbar = gk_hh_dend
                seg.hha_old.el = -70

            # h channels
            gh_basal = [j*gh_soma for j in [1, 1.5, 1.5]*2]
            sec.insert('h')
            for seg in sec:
                seg.h.ghbar = gh_basal[i]
                seg.h.vhalf = -81
                seg.h.eh = -10

            # A-type K+ channels
            sec.insert('kap')
            for seg in sec:
                seg.kap.gkabar = gka_soma

            # m-type K+ channels
            sec.insert('km')
            for seg in sec:
                seg.km.gbar = gkm_trunk

            # slow Ca2+-dependent K+ channel (sAHP)
            sec.insert('kca')
            for seg in sec:
                seg.kca.gbar = gkca_trunk

            # medium Ca2+-dependent K+ channel (mAHP)
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = 0.0165

            # Constant current for Vrest calibration
            sec.insert('constant')
            for seg in sec:
                seg.constant.ic = -0.001
        # =====================================================================

        # =====================================================================
        # All compartments
        for sec in self.all:
            sec.insert('cad')
            sec.ena = 50  # Na reversal potential in mV
            sec.ek = -80  # K reversal potential in mV
            sec.eca = 140  # Ca reversal potential in mV
        # =====================================================================


    def __repr__(self):
        return 'PyramidalCell[{}]'.format(self._gid)


class CCKCell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_topology()
        self._create_lists()
        self._setup_geometry()
        self._setup_passive()
        self._setup_segments()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)

        self.radProx1 = h.Section(name='radProx1', cell=self)
        self.radMed1 = h.Section(name='radMed1', cell=self)
        self.radDist1 = h.Section(name='radDist1', cell=self)
        self.lmM1 = h.Section(name='lmM1', cell=self)
        self.lmt1 = h.Section(name='lmt1', cell=self)

        self.radProx2 = h.Section(name='radProx2', cell=self)
        self.radMed2 = h.Section(name='radMed2', cell=self)
        self.radDist2 = h.Section(name='radDist2', cell=self)
        self.lmM2 = h.Section(name='lmM2', cell=self)
        self.lmt2 = h.Section(name='lmt2', cell=self)

        self.oriProx1 = h.Section(name='oriProx1', cell=self)
        self.oriMed1 = h.Section(name='oriMed1', cell=self)
        self.oriDist1 = h.Section(name='oriDist1', cell=self)

        self.oriProx2 = h.Section(name='oriProx2', cell=self)
        self.oriMed2 = h.Section(name='oriMed2', cell=self)
        self.oriDist2 = h.Section(name='oriDist2', cell=self)

    def _create_lists(self):
        self.all = self.soma.wholetree()
        self.dends = [sec for sec in self.all if not sec.name().__contains__('soma')]

    def _setup_topology(self):
        # Connect sections
        self.radProx1.connect(self.soma(0))
        self.radMed1.connect(self.radProx1(1))
        self.radDist1.connect(self.radMed1(1))
        self.lmM1.connect(self.radDist1(1))
        self.lmt1.connect(self.lmM1(1))
        self.radProx2.connect(self.soma(1))
        self.radMed2.connect(self.radProx2(1))
        self.radDist2.connect(self.radMed2(1))
        self.lmM2.connect(self.radDist2(1))
        self.lmt2.connect(self.lmM2(1))
        self.oriProx1.connect(self.soma(0))
        self.oriMed1.connect(self.oriProx1(1))
        self.oriDist1.connect(self.oriMed1(1))
        self.oriProx2.connect(self.soma(1))
        self.oriMed2.connect(self.oriProx2(1))
        self.oriDist2.connect(self.oriMed2(1))


    def _setup_geometry(self):
        self.soma.L, self.soma.diam = 20 * um, 10 * um

        self.radProx1.L, self.radProx1.diam = 100 * um, 4 * um
        self.radMed1.L, self.radMed1.diam = 100 * um, 3 * um
        self.radDist1.L, self.radDist1.diam = 200 * um, 2 * um
        self.lmM1.L, self.lmM1.diam = 100 * um, 1.5 * um
        self.lmt1.L, self.lmt1.diam = 100 * um, 1 * um

        self.radProx2.L, self.radProx2.diam = 100 * um, 4 * um
        self.radMed2.L, self.radMed2.diam = 100 * um, 3 * um
        self.radDist2.L, self.radDist2.diam = 200 * um, 2 * um
        self.lmM2.L, self.lmM2.diam = 100 * um, 1.5 * um
        self.lmt2.L, self.lmt2.diam = 100 * um, 1 * um

        self.oriProx1.L, self.oriProx1.diam = 100 * um, 2 * um
        self.oriMed1.L, self.oriMed1.diam = 100 * um, 1.5 * um
        self.oriDist1.L, self.oriDist1.diam = 100 * um, 1 * um

        self.oriProx2.L, self.oriProx2.diam = 100 * um, 2 * um
        self.oriMed2.L, self.oriMed2.diam = 100 * um, 1.5 * um
        self.oriDist2.L, self.oriDist2.diam = 100 * um, 1 * um


    def _setup_passive(self):
        for sec in self.all:
            sec.Ra = 100  # Axial resistance
            sec.cm = 1.1  # membrane capacitance


    def _setup_segments(self):
        # Create segments based on `lambda_f`
        d_lambda = 0.1
        frequency = 100
        for sec in self.all:
            sec.nseg = int((sec.L/(d_lambda*h.lambda_f(frequency, sec=sec)) + 0.9) / 2)*2 + 1


    def _setup_biophysics(self):

        # =====================================================================
        # All compartments
        for sec in self.all:

            # HH channels (Na, Kdr, and leak)
            sec.insert('ichan2cck')
            for seg in sec:
                seg.ichan2cck.gnabar = 0.188*0.9*1.2
                seg.ichan2cck.gkbar = 0.013*0.85*1.2
                seg.ichan2cck.gl = 0.00006*0.61
                seg.ichan2cck.el = -70

            # Ca buffering mechanism
            sec.insert('ccanl')

            # A-type K+ channel
            sec.insert('borgka')
            for seg in sec:
                seg.borgka.gkabar = 0.0006

            # N-type Ca2+ channel
            sec.insert('nca')
            for seg in sec:
                seg.nca.gncabar = 0.0000016

            # L-type Ca2+ channel
            sec.insert('lca')
            for seg in sec:
                seg.lca.glcabar = 0.000025

            # Ca2+-dependent K+ (SK) channel
            sec.insert('gskch')
            for seg in sec:
                seg.gskch.gskbar = 0.0004

            # Ca2+ and Voltage-dependent K+ (BK) channel
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = 0.072

            # Ih
            sec.insert('Ih')
            for seg in sec:
                seg.Ih.gkhbar = 0.000025
                seg.Ih.alpha = 100
                seg.Ih.slope = 10
                seg.Ih.amp = 0.01
                seg.Ih.taumin = 20

            # Properties
            sec.eca = 130  # Ca reversal potential in mV
            sec.ena = 55  # Na reversal potential in mV
            sec.ek = -85  # K reversal potential in mV
            sec.eh = -40  # h-channel reversal potential in mV
        # =====================================================================

        # =====================================================================
        self.soma.insert('constant')
        for seg in self.soma:
            seg.constant.ic = -0.00455
        # =====================================================================


    def __repr__(self):
        return 'CCKCell[{}]'.format(self._gid)


class OLMCell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_topology()
        self._create_lists()
        self._setup_geometry()
        self._setup_passive()
        self._setup_segments()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)
        self.dend1 = h.Section(name='dend1', cell=self)
        self.dend2 = h.Section(name='dend2', cell=self)
        self.axon = h.Section(name='axon', cell=self)

    def _create_lists(self):
        self.all = self.soma.wholetree()
        self.dends = [sec for sec in self.all if sec.name().__contains__('dend')]

    def _setup_topology(self):
        # Connect sections
        self.dend1.connect(self.soma(0))
        self.dend2.connect(self.soma(1))
        self.axon.connect(self.soma(0.5))

    def _setup_geometry(self):
        self.soma.L, self.soma.diam = 25 * um, 10 * um

        self.dend1.L, self.dend1.diam = 250 * um, 3 * um
        self.dend2.L, self.dend2.diam = 250 * um, 3 * um
        self.axon.L, self.axon.diam = 150 * um, 1.5 * um


    def _setup_passive(self):
        for sec in self.all:
            sec.Ra = 150  # Axial resistance
            sec.cm = 1.6  # membrane capacitance


    def _setup_segments(self):
        # Create segments based on `lambda_f`
        d_lambda = 0.1
        frequency = 100
        for sec in self.all:
            sec.nseg = int((sec.L/(d_lambda*h.lambda_f(frequency, sec=sec)) + 0.9) / 2)*2 + 1


    def _setup_biophysics(self):

        Rm = 20000*2
        # =====================================================================
        # Soma
        self.soma.insert('IA')
        for seg in self.soma:
            seg.IA.gkAbar = 0.008*1.5

        self.soma.insert('Ih')
        for seg in self.soma:
            seg.Ih.gkhbar = 0.00035*0.1*0.6*0.025*26
            seg.Ih.alpha = 125
            seg.Ih.slope = 12

        self.soma.insert('k_olm')
        for seg in self.soma:
            seg.k_olm.gkbar = 0.0319*1.5

        self.soma.insert('na_olm')
        for seg in self.soma:
            seg.na_olm.gnabar = 0.0107*1.2*2
            seg.na_olm.gl = 1/Rm*2
            seg.na_olm.el = -67

        self.soma.insert('constant')
        for seg in self.soma:
            seg.constant.ic = -0.012
        # =====================================================================

        # =====================================================================
        # Dendrites
        for sec in self.dends:

            sec.insert('IA')
            for seg in sec:
                seg.IA.gkAbar = 0.008*0.9

            sec.insert('Ih')
            for seg in sec:
                seg.Ih.gkhbar = 0.00035*0.06*0.025*26
                seg.Ih.alpha = 125
                seg.Ih.slope = 12

            sec.insert('k_olm')
            for seg in sec:
                seg.k_olm.gkbar = 20*0.023
                seg.k_olm.vhalf_a = 20
                seg.k_olm.qa = 21
                seg.k_olm.vhalf_b = 30
                seg.k_olm.qb = 12


            sec.insert('na_olm')
            for seg in sec:
                seg.na_olm.gnabar = 2*0.0117*2
                seg.na_olm.gl = 1/Rm
                seg.na_olm.el = -65
                seg.na_olm.vhalfa_m = 39
                seg.na_olm.vhalfb_m = 64
                seg.na_olm.vhalfa_h = 65
                seg.na_olm.vhalfb_h = 35
        # =====================================================================

        # =====================================================================
        # Axon
        self.axon.insert('k_olm')
        for seg in self.axon:
            seg.k_olm.gkbar = 0.05104

        self.axon.insert('na_olm')
        for seg in self.axon:
            seg.na_olm.gnabar = 0.01712*2
            seg.na_olm.gl = 1/Rm
            seg.na_olm.el = -67
        # =====================================================================

    def __repr__(self):
        return 'OLMCell[{}]'.format(self._gid)


class VIPCCKCell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_topology()
        self._create_lists()
        self._setup_geometry()
        self._setup_passive()
        self._setup_segments()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)

        self.radProx1 = h.Section(name='radProx1', cell=self)
        self.radMed1 = h.Section(name='radMed1', cell=self)
        self.radDist1 = h.Section(name='radDist1', cell=self)
        self.lmM1 = h.Section(name='lmM1', cell=self)
        self.lmt1 = h.Section(name='lmt1', cell=self)

        self.radProx2 = h.Section(name='radProx2', cell=self)
        self.radMed2 = h.Section(name='radMed2', cell=self)
        self.radDist2 = h.Section(name='radDist2', cell=self)
        self.lmM2 = h.Section(name='lmM2', cell=self)
        self.lmt2 = h.Section(name='lmt2', cell=self)

        self.oriProx1 = h.Section(name='oriProx1', cell=self)
        self.oriMed1 = h.Section(name='oriMed1', cell=self)
        self.oriDist1 = h.Section(name='oriDist1', cell=self)

        self.oriProx2 = h.Section(name='oriProx2', cell=self)
        self.oriMed2 = h.Section(name='oriMed2', cell=self)
        self.oriDist2 = h.Section(name='oriDist2', cell=self)

    def _create_lists(self):
        self.all = self.soma.wholetree()
        self.dends = [sec for sec in self.all if not sec.name().__contains__('soma')]

    def _setup_topology(self):
        # Connect sections
        self.radProx1.connect(self.soma(0))
        self.radMed1.connect(self.radProx1(1))
        self.radDist1.connect(self.radMed1(1))
        self.lmM1.connect(self.radDist1(1))
        self.lmt1.connect(self.lmM1(1))
        self.radProx2.connect(self.soma(1))
        self.radMed2.connect(self.radProx2(1))
        self.radDist2.connect(self.radMed2(1))
        self.lmM2.connect(self.radDist2(1))
        self.lmt2.connect(self.lmM2(1))
        self.oriProx1.connect(self.soma(0))
        self.oriMed1.connect(self.oriProx1(1))
        self.oriDist1.connect(self.oriMed1(1))
        self.oriProx2.connect(self.soma(1))
        self.oriMed2.connect(self.oriProx2(1))
        self.oriDist2.connect(self.oriMed2(1))


    def _setup_geometry(self):
        self.soma.L, self.soma.diam = 20 * um, 10 * um

        self.radProx1.L, self.radProx1.diam = 50 * um, 4 * um
        self.radMed1.L, self.radMed1.diam = 50 * um, 3 * um
        self.radDist1.L, self.radDist1.diam = 100 * um, 2 * um
        self.lmM1.L, self.lmM1.diam = 50 * um, 1.5 * um
        self.lmt1.L, self.lmt1.diam = 50 * um, 1 * um

        self.radProx2.L, self.radProx2.diam = 50 * um, 4 * um
        self.radMed2.L, self.radMed2.diam = 50 * um, 3 * um
        self.radDist2.L, self.radDist2.diam = 100 * um, 2 * um
        self.lmM2.L, self.lmM2.diam = 50 * um, 1.5 * um
        self.lmt2.L, self.lmt2.diam = 50 * um, 1 * um


        self.oriProx1.L, self.oriProx1.diam = 50 * um, 2 * um
        self.oriMed1.L, self.oriMed1.diam = 50 * um, 1.5 * um
        self.oriDist1.L, self.oriDist1.diam = 50 * um, 1 * um

        self.oriProx2.L, self.oriProx2.diam = 50 * um, 2 * um
        self.oriMed2.L, self.oriMed2.diam = 50 * um, 1.5 * um
        self.oriDist2.L, self.oriDist2.diam = 50 * um, 1 * um


    def _setup_passive(self):
        for sec in self.all:
            sec.Ra = 100  # Axial resistance
            sec.cm = 1.3  # membrane capacitance


    def _setup_segments(self):
        # Create segments based on `lambda_f`
        d_lambda = 0.1
        frequency = 100
        for sec in self.all:
            sec.nseg = int((sec.L/(d_lambda*h.lambda_f(frequency, sec=sec)) + 0.9) / 2)*2 + 1


    def _setup_biophysics(self):

        # =====================================================================
        # All compartments
        for sec in self.all:

            # HH channels (Na, Kdr, and leak)
            sec.insert('ichan2vip')
            for seg in sec:
                seg.ichan2vip.gnabar = 0.18/0.1*2
                seg.ichan2vip.gkbar = 0.013/0.1*2
                seg.ichan2vip.gl = 0.00018/2*0.1*2
                seg.ichan2vip.el = -65

            # Ca buffering mechanism
            sec.insert('ccanl')
            for seg in sec:
                seg.ccanl.catau = 70

            # A-type K+ channel
            sec.insert('borgka')
            for seg in sec:
                seg.borgka.gkabar = 0.00015

            # N-type Ca2+ channel
            sec.insert('nca')
            for seg in sec:
                seg.nca.gncabar = 0.0008*0.15

            # L-type Ca2+ channel
            sec.insert('lca')
            for seg in sec:
                seg.lca.glcabar = 0.005*0.15

            # Ca2+-dependent K+ (SK) channel
            sec.insert('gskch')
            for seg in sec:
                seg.gskch.gskbar = 0.000002*2*10000
                seg.gskch.beta = 0.00025*0.5

            # Ca2+ and Voltage-dependent K+ (BK) channel
            sec.insert('mykca')
            for seg in sec:
                seg.mykca.gkbar = 0.0002*10

            # Ih
            sec.insert('Ih')
            for seg in sec:
                seg.Ih.gkhbar = 0.000025*1.1
                seg.Ih.alpha = 100
                seg.Ih.slope = 8
                seg.Ih.amp = 0.001
                seg.Ih.taumin = 30


            sec.ena = 55  # Na reversal potential in mV
            sec.ek = -85  # K reversal potential in mV
            sec.eca = 130 # Ca reversal potential in mV
            sec.eh = -40 # h-channel reversal potential in mV
        # =====================================================================

        # =====================================================================
        self.soma.insert('constant')
        for seg in self.soma:
            seg.constant.ic = -0.0003
        # =====================================================================


    def __repr__(self):
        return 'VIPCCKCell[{}]'.format(self._gid)


class VIPCRCell:
    def __init__(self, gid):
        self._gid = gid
        self._setup_morphology()
        self._setup_topology()
        self._create_lists()
        self._setup_geometry()
        self._setup_passive()
        self._setup_segments()
        self._setup_biophysics()

    def _setup_morphology(self):
        self.soma = h.Section(name='soma', cell=self)

        self.radProx1 = h.Section(name='radProx1', cell=self)
        self.radMed1 = h.Section(name='radMed1', cell=self)
        self.radDist1 = h.Section(name='radDist1', cell=self)
        self.lmM1 = h.Section(name='lmM1', cell=self)
        self.lmt1 = h.Section(name='lmt1', cell=self)

        self.radProx2 = h.Section(name='radProx2', cell=self)
        self.radMed2 = h.Section(name='radMed2', cell=self)
        self.radDist2 = h.Section(name='radDist2', cell=self)
        self.lmM2 = h.Section(name='lmM2', cell=self)
        self.lmt2 = h.Section(name='lmt2', cell=self)

        self.oriProx1 = h.Section(name='oriProx1', cell=self)
        self.oriMed1 = h.Section(name='oriMed1', cell=self)
        self.oriDist1 = h.Section(name='oriDist1', cell=self)

        self.oriProx2 = h.Section(name='oriProx2', cell=self)
        self.oriMed2 = h.Section(name='oriMed2', cell=self)
        self.oriDist2 = h.Section(name='oriDist2', cell=self)

    def _create_lists(self):
        self.all = self.soma.wholetree()
        self.dends = [sec for sec in self.all if not sec.name().__contains__('soma')]

    def _setup_topology(self):
        # Connect sections
        self.radProx1.connect(self.soma(0))
        self.radMed1.connect(self.radProx1(1))
        self.radDist1.connect(self.radMed1(1))
        self.lmM1.connect(self.radDist1(1))
        self.lmt1.connect(self.lmM1(1))
        self.radProx2.connect(self.soma(1))
        self.radMed2.connect(self.radProx2(1))
        self.radDist2.connect(self.radMed2(1))
        self.lmM2.connect(self.radDist2(1))
        self.lmt2.connect(self.lmM2(1))
        self.oriProx1.connect(self.soma(0))
        self.oriMed1.connect(self.oriProx1(1))
        self.oriDist1.connect(self.oriMed1(1))
        self.oriProx2.connect(self.soma(1))
        self.oriMed2.connect(self.oriProx2(1))
        self.oriDist2.connect(self.oriMed2(1))


    def _setup_passive(self):
        for sec in self.all:
            sec.Ra = 150  # Axial resistance
            sec.cm = 1.4  # membrane capacitance


    def _setup_segments(self):
        # Create segments based on `lambda_f`
        d_lambda = 0.1
        frequency = 100
        for sec in self.all:
            sec.nseg = int((sec.L/(d_lambda*h.lambda_f(frequency, sec=sec)) + 0.9) / 2)*2 + 1


    def _setup_geometry(self):
        self.soma.L, self.soma.diam = 20 * um, 10 * um

        self.radProx1.L, self.radProx1.diam = 100 * um, 4 * um
        self.radMed1.L, self.radMed1.diam = 100 * um, 3 * um
        self.radDist1.L, self.radDist1.diam = 200 * um, 2 * um
        self.lmM1.L, self.lmM1.diam = 100 * um, 1.5 * um
        self.lmt1.L, self.lmt1.diam = 100 * um, 1 * um

        self.radProx2.L, self.radProx2.diam = 100 * um, 4 * um
        self.radMed2.L, self.radMed2.diam = 100 * um, 3 * um
        self.radDist2.L, self.radDist2.diam = 200 * um, 2 * um
        self.lmM2.L, self.lmM2.diam = 100 * um, 1.5 * um
        self.lmt2.L, self.lmt2.diam = 100 * um, 1 * um


        self.oriProx1.L, self.oriProx1.diam = 100 * um, 2 * um
        self.oriMed1.L, self.oriMed1.diam = 100 * um, 1.5 * um
        self.oriDist1.L, self.oriDist1.diam = 100 * um, 1 * um

        self.oriProx2.L, self.oriProx2.diam = 100 * um, 2 * um
        self.oriMed2.L, self.oriMed2.diam = 100 * um, 1.5 * um
        self.oriDist2.L, self.oriDist2.diam = 100 * um, 1 * um

    def _setup_biophysics(self):

        # =====================================================================
        Rm = 20000*1.6
        # Somatic compartment
        self.soma.insert('pas')
        for seg in self.soma:
            seg.pas.g = 1/Rm
            seg.pas.e = -70

        # Na channels
        self.soma.insert('Nafcr')
        for seg in self.soma:
            seg.Nafcr.gnafbar = 0.015*2

        # Kdr channels
        self.soma.insert('kdrcr')
        for seg in self.soma:
            seg.kdrcr.gkdrbar = 0.0018*2

        # IKscr
        self.soma.insert('IKscr')
        for seg in self.soma:
            seg.IKscr.gKsbar = 0.000725*5*1.8

        # iCcr
        self.soma.insert('iCcr')
        for seg in self.soma:
            seg.iCcr.gkcbar = 0.00003*10*1.5

        # cancr
        self.soma.insert('cancr')
        for seg in self.soma:
            seg.cancr.gcabar = 0.001*0.1*10

        # gskch
        self.soma.insert('gskch')
        for seg in self.soma:
            seg.gskch.gskbar = 0.000002*1000*2

        # Constant current for Vrest calibration
        self.soma.insert('constant')
        for seg in self.soma:
            seg.constant.ic = -0.003

        # ccanl
        self.soma.insert('ccanl')
        for seg in self.soma:
            seg.ccanl.catau = 20

        # Ih
        self.soma.insert('Ih')
        for seg in self.soma:
            seg.Ih.gkhbar = 0.000025*8
            seg.Ih.alpha = 100
            seg.Ih.slope = 10
            seg.Ih.amp = 0.001
            seg.Ih.taumin = 5

        self.soma.eca = 130
        self.soma.eh = -40
        # =====================================================================

        # =====================================================================
        # Dendritic compartments
        for sec in self.dends:
            sec.insert('pas')
            for seg in sec:
                seg.pas.g = 1/Rm
                seg.pas.e = -70

            sec.insert('Nafcr')
            for seg in sec:
                seg.Nafcr.gnafbar = 0.018*5

            sec.insert('kdrcr')
            for seg in sec:
                seg.kdrcr.gkdrbar = 0.018*0.5

            sec.insert('Ih')
            for seg in sec:
                seg.Ih.gkhbar = 0.000025*1.1
                seg.Ih.alpha = 100
                seg.Ih.slope = 8
                seg.Ih.amp = 0.001
                seg.Ih.taumin = 30

        # =====================================================================

        # =====================================================================
        # All compartments
        for sec in self.all:
            sec.ena = 55  # Na reversal potential in mV
            sec.ek = -85  # K reversal potential in mV
        # =====================================================================


    def __repr__(self):
        return 'VIPCRCell[{}]'.format(self._gid)
