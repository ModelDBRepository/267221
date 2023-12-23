TITLE Calcium-activated potassium channel

COMMENT
Ca2+-dependent K+ current.
non-voltage-dependent.
Used in Basket, Axoaxonic, Bistratified, CCK+, VIP+/CCK+, and VIP+/CR+ cells.
ENDCOMMENT

NEURON {
	SUFFIX gskch
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gskbar, ek, ik
	RANGE beta, factor, qinf, tauq
}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gskbar = 0 (S/cm2)
	beta = 0.00025 (/ms)
	factor = 1.25e1 (/mM/mM/ms) : unit conversion factor
}

ASSIGNED {
	v (mV)
	ek (mV)
	cai (mM)
	ik (mA/cm2)
	qinf (1)
	tauq (ms)
}

STATE {
	q
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gskbar*pow(q, 2)*(v - ek)
}

DERIVATIVE states {
	rates(cai)
	q' = (qinf - q)/tauq
}

INITIAL {
	rates(cai)
	q = qinf
}

PROCEDURE rates(ca (mM)) {
	:Computes rate and other constants at current cai.
	LOCAL alpha
	
	:"q" activation system
	alpha = factor * ca * ca
	tauq = 1 / (alpha + beta)
	qinf = alpha * tauq
}
