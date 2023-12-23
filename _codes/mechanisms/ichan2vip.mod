TITLE HH channel that includes a sodium, a delayed rectifier and a leak channel

COMMENT
Hodgin-Huxley style channels for Na+, delayed rectifier K+, and leak currents.
This channel specifically used in vip/cck cells.
ENDCOMMENT

NEURON {
	SUFFIX ichan2vip
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gna, gk
	RANGE gnabar, gkbar, gl, el
	RANGE ina, ik, il
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gnabar = 0 (S/cm2)
	gkbar = 0 (S/cm2)
	gl = 0 (S/cm2)
}

ASSIGNED {
	v (mV)

	ena (mV)
	ek (mV)
	el (mV)

	ina (mA/cm2)
	ik (mA/cm2)
	il (mA/cm2)

	minf (1)
	taum (ms)

	hinf (1)
	tauh (ms)

	ninf (1)
	taun (ms)

	gna (S/cm2) 
	gk (S/cm2)
}

STATE {
	m
	h
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gna = gnabar*pow(m, 3)*h
	ina = gna*(v - ena)  :Sodium current

	gk = gkbar*pow(n, 4)
	ik = gk*(v - ek)   :fast potassium current

	il = gl*(v - el)
}

DERIVATIVE states {
	rates(v)
	:Na activation variable
	m' = (minf - m)/taum
	:Na inactivation variable
	h' = (hinf - h)/tauh
	:K fast activation variable
	n' = (ninf - n)/taun
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
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

PROCEDURE rates(v (mV)) {
	:Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL  alpha, beta
	
	:"m" sodium activation system - act and inact cross at -40
	alpha = 0.3(/ms)*vtrap(-(v + 26(mV)), 5(mV))
	beta = 0.3(/ms)*vtrap((v - 2(mV)), 5(mV))
	taum = 1/(alpha + beta)
	minf = alpha*taum

	:"h" sodium inactivation system
	alpha = 0.23(/ms)/exp((v + 48(mV))/10(mV))
	beta = 3.33(/ms)/(1 + exp(-(v - 4.5(mV))/5(mV)))
	tauh = 1/(alpha + beta)
	hinf = alpha*tauh

	:"nf" fKDR activation system
	alpha = 0.07(/ms)*vtrap(-(v - 39(mV)), 6(mV))
	beta = 0.264(/ms)/exp((v - 14(mV))/40(mV))
	taun = 1/(alpha + beta)
	ninf = alpha*taun
}
