TITLE Fast sodium channel

COMMENT
Fast activating Na+ current
Used in VIP+/CR+ cell.
ENDCOMMENT

NEURON {
	SUFFIX Nafcr
	USEION na READ ena WRITE ina
	RANGE gnafbar, ina
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gnafbar= 0.086 (S/cm2)
}

ASSIGNED {
	v (mV)
	ena (mV)
	ina (mA/cm2)
	minf (1)
	hinf (1)
	taum (ms)
	tauh (ms)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnafbar*pow(m, 3)*h*(v-ena)
}

DERIVATIVE states {
	rates(v)
	m' = (minf-m)/taum
	h' = (hinf-h)/tauh
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	:Traps for 0 in denominator of rate eqns. Taylor expansion is used.
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION malf(v (mV)) (/ms){
	malf = 0.2816(/ms)*vtrap(-(v + 33(mV)), 9.3(mV))
}

FUNCTION mbet(v (mV)) (/ms) {
	mbet = 0.2464(/ms)*vtrap((v + 6(mV)), 6(mV))
}	

FUNCTION half(v (mV)) (/ms) {
	half = 0.098(/ms)/exp((v + 23.1(mV))/20(mV))
}

FUNCTION hbet(v (mV)) (/ms) {
	hbet = 1.4(/ms)/(1 + exp(-(v + 25.1(mV))/10(mV)))
}

PROCEDURE rates(v (mV)) {
	LOCAL ma, mb, ha, hb
	
	ma = malf(v)
	mb = mbet(v)
	minf = ma/(ma+mb)
	taum = 1/(ma+mb)

	ha = half(v)
	hb = hbet(v)
	hinf = ha/(ha+hb)
	tauh = 1 / (ha+hb)
}
