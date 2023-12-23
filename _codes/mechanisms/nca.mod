TITLE N-type calcium channel
 
COMMENT
Voltage-dependent calcium channel
Used in Basket, Axoaxonic, Bistartfified, CCK+, VIP+/CCK+ cells.
Based on Cutsuridis et al., 2010, Hippocampus, 20:423:446, doi: 10.1002/hipo.20661
ENDCOMMENT
 
 NEURON { 
	SUFFIX nca
	USEION ca READ eca WRITE ica
	RANGE gncabar, eca, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}
 
PARAMETER {
	gncabar = 0 (S/cm2)
}

ASSIGNED {
	v (mV)
	ica (mA/cm2)
	eca (mV)
	cinf (1)
	dinf (1)
	tauc (ms)
	taud (ms) 
}

STATE {
	c
	d
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gncabar*pow(c, 2)*d*(v - eca)
}

DERIVATIVE states {
	rates(v)
	c' = (cinf - c)/tauc
	d' = (dinf - d)/taud
}

INITIAL {
	rates(v)
	c = cinf
	d = dinf
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns. Taylor expansion is used.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

PROCEDURE rates(v (mV)) {
	:Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL alpha, beta
	
	:"c" NCa activation system
	alpha = -0.19(/ms)*vtrap(v - 19.88(mV), -10(mV))
	beta = 0.046(/ms)*exp(-v/20.73(mV))
	tauc = 1/(alpha + beta)
	cinf = alpha/(alpha + beta)
	
	:"d" NCa inactivation system
	alpha = 0.00016(/ms)/exp(-v/48.4(mV))
	beta = 1(/ms)/(exp((-v + 39(mV))/10(mV)) + 1)
	taud = 1/(alpha + beta)
	dinf = alpha/(alpha + beta)
}
