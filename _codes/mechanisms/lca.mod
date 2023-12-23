TITLE L-type calcium channel

COMMENT
L-type Ca2+ current.
Voltage and Ca2+-dependent current.
Used in Basket, Axoaxonic, Bistratified, CCK+, VIP+/CCK+ cells.
ENDCOMMENT

NEURON {
	SUFFIX lca
	USEION ca READ eca, cai, cao WRITE ica
	RANGE glcabar, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	glcabar	= 0.0 (S/cm2)
	ki = .001 (mM)
	tfa = 1 (1)
}

ASSIGNED {
	v (mV)
	eca (mV)
	celsius (degC)
	ica (mA/cm2)
	glca (S/cm2)
	minf
	taum (ms)
	cai (mM)
	cao (mM)
}

STATE {
	m
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	glca = glcabar*pow(m, 2)*h2(cai)
	ica = glca*ghk(v, cai, cao)
}

DERIVATIVE state {
	rates(v)
	m' = (minf - m)/taum
}

INITIAL {
	rates(v)
	m = minf
}

FUNCTION h2(cai (mM)) {
	h2 = ki/(ki + cai)
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (mV) {
	LOCAL nu, f
	f = KTF(celsius)/2 : in mV
	nu = v/f : unitless
	ghk = -f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
	KTF = (25.(mV)/293.15(degC))*(celsius + 273.15(degC))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns. Taylor expansion is used.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION alp(v( mV)) (/ms) {
	alp = 15.69(/ms)*vtrap(-(v - 81.5(mV)), 10.0(mV))
}

FUNCTION bet(v (mV)) (/ms) {
	bet = 0.29(/ms)*exp(-v/10.86(mV))
}

PROCEDURE rates(v (mV)) {
	:callable from hoc
	taum = 1/(tfa*(alp(v) + bet(v)))
	minf = tfa*alp(v)*taum
}
