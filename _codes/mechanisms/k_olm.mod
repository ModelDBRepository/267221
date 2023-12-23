TITLE Potassium current

COMMENT
Potassium current for the soma of OLM cells.
ENDCOMMENT

NEURON {
	SUFFIX k_olm
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	RANGE Ra, Rb, vhalf_a, vhalf_b, qa, qb
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gkbar = 0 (S/cm2)
	Ra = 0.018(/ms)
	Rb = 0.0036(/ms)
	vhalf_a = 25(mV)
	vhalf_b = 35(mV)
	qa = 25 (mV)
	qb = 12 (mV)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	ninf 
	taun (ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*pow(n, 4)*(v - ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/taun
}

INITIAL {
	rates(v)
	n = ninf
}

FUNCTION Exp(x) {
	if (x < -100) {
		Exp = 0
	} else {
		Exp = exp(x)
	}
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {
		vtrap = 1(/mV)*x/(Exp(x/y) - 1)
	}
}

PROCEDURE rates(v (mV) ) {  
	:Computes rate and other constants at current v.
	LOCAL alpha, beta

	alpha = Ra*vtrap(-(v - vhalf_a), qa)
	beta = Rb*vtrap(v - vhalf_b, qb)
	
	taun = 1/(alpha + beta)
	ninf = alpha*taun
}
