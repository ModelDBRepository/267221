TITLE Fast calcium- and voltage-dependent potassium channel

COMMENT
Voltage-dependent and Ca2+-dependent K+ current.
Used in VIP+/CR+ cells.
ENDCOMMENT

NEURON {
	SUFFIX iCcr
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkcbar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gkcbar = 0.0 (S/cm2)
	taumin = 1.1 (ms)
	K = 1000 (/mM) : calcium scaling factor
	vsc = 40(mV) : voltage scaling factor (shifting)
	vth = -250 (mV) : voltage threshold to avoid numerical overflows
}

ASSIGNED {
	v (mV)
	cai	(mM)
	ik (mA/cm2)
	cinf (1)
	tauc (ms)
	ek (mV)
}

STATE {
	c
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkcbar*pow(c, 2)*(v-ek)
}

DERIVATIVE states {
	rates(v, cai)
	c' = (cinf - c)/tauc
}

INITIAL {
	rates(v, cai)
	c = cinf
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns. Taylor expansion is used.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION ca2vol(v (mV), ca (mM)) (mV) {
	ca2vol = v + vsc*log10(K*ca)
	: To avoid numerical overflow in exponentials
	if (ca2vol < vth) { ca2vol = vth } 
}

FUNCTION calf(v (mV), ca (mM)) (/ms) {
	calf = 0.00642(/ms)*vtrap(-(ca2vol(v, ca) + 18(mV)), 12(mV))
}

FUNCTION cbet(v (mV), ca (mM)) (/ms) {
	cbet = 1.7(/ms)*exp(-(ca2vol(v, ca) + 152(mV))/30(mV))
}

PROCEDURE rates(v (mV), ca(mM) ) {
	cinf = calf(v, ca)/(calf(v, ca) + cbet(v, ca))
	tauc = 1/(calf(v, ca) + cbet(v, ca))

	if (tauc < taumin) {
		tauc = taumin
	}
}
