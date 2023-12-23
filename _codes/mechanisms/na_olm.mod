TITLE Sodium current

COMMENT
Used in dendrites of OLM cells.

References:

References:

1.	Martina et al., 2000, Science, 287:295-300, doi: 10.1126/science.287.5451.295.

				soma			axon-lacking dend	axon-bearing dend
Na+	gmax	    107 ps/um2		117 ps/um2			107 ps/um2
	slope 	    10.9 mV/e		11.2 mV/e		   	11.2 mV/e
	V1/2        -37.8 mV       -45.6 mV				-45.6 mV

2.	Marina and Jonas, 1997, J Physiol, 505:593-603, doi: 10.1111/j.1469-7793.1997.593ba.x

	*Note* The interneurons here are basket cells from the dentate gyrus.

Na+	Activation V1/2				-25.1 mV
	slope			 		11.5
	Activation t (-20 mV)	 		0.16 ms
	Deactivation t (-40 mV)	 		0.13 ms
 	Inactivation V1/2			-58.3 mV
	slope			 		6.7
	onset of inactivation t (-20 mV)	1.34 ms
	onset of inactivation t (-55 mV)	18.6 ms
	recovery from inactivation t		2.0 ms
	(30 ms conditioning pulse)
	recovery from inactivation t		2.7 ms
	(300 ms conditioning pulse)
ENDCOMMENT

NEURON {
	SUFFIX na_olm
	USEION na READ ena WRITE ina
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gl, el, ina
	RANGE Ra_m, vhalfa_m, qa_m, Rb_m, vhalfb_m, qb_m
	RANGE Ra_h, vhalfa_h, qa_h, Rb_h, vhalfb_h, qb_h
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gnabar = 0 (S/cm2)
	gl = 0 (S/cm2)

	: rate coefficient params for soma, and axon
	Ra_m = 0.1(/ms)
	vhalfa_m = 32(mV)
	qa_m = 10(mV)

	Rb_m = 4(/ms)
	vhalfb_m = 57(mV)
	qb_m = 18(mV)

	Ra_h = 0.07(/ms)
	vhalfa_h = 58(mV)
	qa_h = 20(mV)

	Rb_h = 1(/ms)
	vhalfb_h = 28(mV)
	qb_h = 10(mV)
}

ASSIGNED {
	v (mV)
	el (mV)
	ena (mV)
	ina (mA/cm2)
	il (mA/cm2)
	minf
	hinf
	taum (ms)
	tauh (ms)
}

STATE {
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*pow(m, 3)*h*(v - ena)
	il = gl*(v - el)
}

DERIVATIVE states {
	:exact when v held constant
	rates(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
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

PROCEDURE rates(v(mV)) {
	:Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL alpha, beta

	alpha = Ra_m*vtrap(-(v + vhalfa_m), qa_m)
	beta  = Rb_m*Exp(-(v + vhalfb_m)/qb_m)
	taum  = 1/(alpha + beta)
	minf  = alpha*taum

	alpha = Ra_h*Exp(-(v + vhalfa_h)/qa_h)
	beta  = Rb_h/(1 + Exp(-(v + vhalfb_h)/qb_h))
	tauh  = 1/(alpha + beta)
	hinf  = alpha*tauh
}