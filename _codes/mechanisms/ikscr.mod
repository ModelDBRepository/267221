TITLE Slowly inactivating potassium channel

COMMENT
K+ current, used in VIP+/CR+ cell.
ENDCOMMENT

NEURON {
	SUFFIX IKscr
	USEION k READ ek WRITE ik
	RANGE gKsbar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gKsbar = 0 (S/cm2)
}

ASSIGNED {
	v (mV)
	ek (mV)
	ik (mA/cm2)
	ainf binf
	taua (ms)
	taub (ms)
	gk (mho/cm2)
}

STATE {
	a
	b
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gKsbar*a*b*(v - ek)
}

DERIVATIVE states {
	rates(v)
	a' = (ainf - a)/taua
	b' = (binf - b)/taub
}

INITIAL {
	rates(v)
	a = ainf
	b = binf
}

PROCEDURE rates(v (mV)) {
	ainf = 1/(1 + exp(-(v + 34(mV))/6.5(mV)))
	taua = 10

	binf = 1/(1 + exp((v + 65(mV))/6.6(mV)))
	taub = 200(ms) + 3200(ms) / (1 + exp(-(v + 63.6(mV))/4(mV)))
}
