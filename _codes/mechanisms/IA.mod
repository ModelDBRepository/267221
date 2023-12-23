TITLE A-type potassium channel

COMMENT
Used in OLM cells.

Reference:

1.	Zhang and McBain, 1995, J Physiol, 488:647-660, doi: 10.1113/jphysiol.1995.sp020997

	Activation V1/2 = -14 mV
	slope = 16.6
	activation t = 5 ms
	Inactivation V1/2 = -71 mV
	slope = 7.3
	inactivation t = 15 ms
	recovery from inactivation = 142 ms

2.	Martina et al., 1998, J Neurosci 18:8111-8125, doi: 10.1523/JNEUROSCI.18-20-08111.1998
	(only the gkAbar is from this paper)

	gkabar = 0.0175 mho/cm2
	Activation V1/2 = -6.2 +/- 3.3 mV
	slope = 23.0 +/- 0.7 mV
	Inactivation V1/2 = -75.5 +/- 2.5 mV
	slope = 8.5 +/- 0.8 mV
	recovery from inactivation t = 165 +/- 49 ms  

3.	Warman et al., 1994, J Neurophysiol, 71:2033-2045, doi:10.1152/jn.1994.71.6.2033

	gkabar = 0.01 mho/cm2
	
	**(number taken from the work by Numann et al., 1987, J Physiol, 393:331-353, doi: 10.1113/jphysiol.1987.sp016826 in guinea pig CA1 neurons)
ENDCOMMENT

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

NEURON {
	SUFFIX IA
	USEION k READ ek WRITE ik
	RANGE gkAbar, ik
}

PARAMETER {
	gkAbar = 0.0165 (S/cm2)	:from Martina et al.
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	ainf
	binf
	tau_b (ms)
	tau_a (ms)
}

STATE {
	a
	b
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkAbar*a*b*(v - ek)
}

DERIVATIVE states { 
	: Computes state variables m, h, and n rates(v) at the current v and dt.
	: call of rates(v) required to update inf and tau values
	rates(v)
	a' = (ainf - a)/(tau_a)
	b' = (binf - b)/(tau_b)
}

INITIAL {
	rates(v)
	a = ainf
	b = binf
}

PROCEDURE rates(v (mV)) {
	:Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL alpha_b, beta_b

	alpha_b = 0.000009(/ms) / exp((v - 26(mV))/18.5(mV))
	beta_b = 0.014 (/ms) / (exp(-(v + 70(mV))/11(mV)) + 0.2)
	tau_b = 1/(alpha_b + beta_b)
	binf = 1/(1 + exp((v + 71(mV))/7.3(mV)))

	tau_a = 5
	ainf = 1/(1 + exp(-(v + 14(mV))/20.6(mV)))
}
