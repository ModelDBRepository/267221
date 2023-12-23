TITLE hyperpolarization-activated channel - h current

COMMENT
Ih current - hyperpolarization-activated nonspecific Na and K channel
- contributes to the resting membrane potential
- controls the afterhyperpolarization

Used in CCK+, OLM, VIP+/CCK+, VIP+/CR+ cells.

Reference:

	1.	Maccaferri and McBain, 1996, J Physiol, 497:119-130, doi: 10.1113/jphysiol.1996.sp021754
		.

		V1/2 = -84.1 mV
		k = 10.2
		reversal potential = -32.9 +/- 1.1 mV

		at -70 mV, currents were fitted by a single exponetial of t = 2.8+/- 0.76 s
		at -120 mV, two exponentials were required, t1 = 186.3+/-33.6 ms 
		t2 = 1.04+/-0.16 s

	2.	Maccaferri et al., 1993, J Neurophysiol, 69:2129-2136, doi: 10.1152/jn.1993.69.6.2129

		V1/2 = -97.9 mV
		k = 13.4
		reversal potential = -18.3 mV

	3.	Pape, 1996, Annu Rev Physiol, 58:299-327, doi: 10.1146/annurev.ph.58.030196.001503

			single channel conductance is around 1 pS
			average channel density is below 0.5 um-2
			0.5 pS/um2 = 0.00005 mho/cm2 = 0.05 umho/cm2

	4.	Magee, 1998, J Neurosci, 18:7613-7624, doi: 10.1523/JNEUROSCI.18-19-07613.1998

Deals with Ih in CA1 pyramidal cells. Finds that conductance density increases with distance from the soma.
	soma g = 0.0013846 mho/cm2
	dendrite g (300-350 um away) = 0.0125 mho/cm2
	see Table 1 in the paper
ENDCOMMENT
 
NEURON {
	SUFFIX Ih
	USEION h READ eh WRITE ih VALENCE 1
	RANGE gkhbar, ih
	RANGE slope, alpha, taumin, amp
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gkhbar = 0.0 (S/cm2)
	eh = -32.9 (mV)
	slope = 10.2 (mV)
	alpha = 84.1 (mV)
	taumin = 100 (ms)
	amp = 1 (ms)
}

ASSIGNED {
	v (mV)
	ih (mA/cm2)
	rinf (1)
	taur (ms)
}

STATE {
	r
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	ih = gkhbar*r*(v - eh)
}

DERIVATIVE state {
	:Computes state variable h at current v and dt.
	rates(v)
	r' = (rinf - r)/taur
}

INITIAL {
	rates(v)
	r = rinf
}

PROCEDURE rates(v (mV)) {
	:Computes rate and other constants at current v.
	:Call once from HOC to initialize inf at resting v.
	LOCAL ar, br

	rinf = 1/(1 + exp((v + alpha)/slope))

	ar = exp(-(0.116(/mV)*v + 0.116))
	br = exp(0.09(/mV)*v - 1.84)
	taur = taumin + amp/(ar + br)
}
