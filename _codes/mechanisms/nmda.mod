TITLE NMDA synapse 
: from Durstewitz & Gabriel (2006), Cerebral Cortex

NEURON {
	POINT_PROCESS NMDA
	NONSPECIFIC_CURRENT i
	RANGE g, a, b, gNMDAmax, i
	RANGE tauD, tauF, util, tcon, tcoff, enmda
}

UNITS {
	(uS) = (microsiemens)
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tcon = 2.3 (ms)
	tcoff = 95.0 (ms)
	enmda = 0 (mV)
	gNMDAmax = 1 (uS)
	tauD = 800 (ms)
	tauF = 800 (ms)
	util = .3
}

ASSIGNED {
	v (mV)
	i (nA)
	g (1)
	factor
}

STATE {
	a
	b
}

BREAKPOINT {
	LOCAL s
	SOLVE states METHOD cnexp
	s = 1.50265/(1+0.33*exp(-0.0625(/mV)*v))
	g = b - a
	i = gNMDAmax*g*s*(v-enmda)
}

DERIVATIVE states {
	a' = -a/tcon
	b' = -b/tcoff
}

INITIAL { 
	a = 0  
	b = 0 
	factor = tcon*tcoff*1(/ms) / (tcoff - tcon)
}

NET_RECEIVE(wgt, R, u, tlast (ms), nspike) {
	LOCAL x

	if (nspike==0) { R=1  u=util }
	else {
		if (tauF>0) { u = util+(1-util)*u*exp(-(t-tlast)/tauF) }
		if (tauD>0) { R = 1+(R*(1-u)-1)*exp(-(t-tlast)/tauD) }
	}
	x = wgt*factor*R*u
	state_discontinuity(a, a+x)
	state_discontinuity(b, b+x)
	tlast = t
	nspike = nspike+1
}
