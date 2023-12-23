TITLE Delayed rectifier potassium channel

COMMENT
Voltage-dependent K+ current.
Used in VIP+/CR+ cells.
ENDCOMMENT

NEURON {
	SUFFIX kdrcr
	USEION k READ ek WRITE ik
	RANGE gkdrbar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gkdrbar = 0.0338 (S/cm2)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	ninf
	taun (ms)
	gk (mho/cm2)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkdrbar*pow(n, 4)*(v-ek)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf-n)/taun
}

INITIAL {
	rates(v)
	n = ninf
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns. Taylor expansion is used.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION alf(v (mV)) (/ms){
	alf = 0.018(/ms)*vtrap(-(v-13(mV)), 25(mV))
}

FUNCTION bet(v (mV)) (/ms) {
	bet = 0.0054(/ms)*vtrap(v-23(mV), 12(mV))
}

PROCEDURE rates(v (mV)) {
	taun = 1/(alf(v) + bet(v))
	ninf = alf(v)*taun
	
}
