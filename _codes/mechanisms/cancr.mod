TITLE N-type calcium channel

COMMENT
Used in somatic and dendritic regions in VIP+/CR+ cell.
Voltage and Ca2+-dependent channel.
After Borg, Updated by Maria Markaki 03/12/03
ENDCOMMENT

NEURON {
	SUFFIX cancr
	USEION ca READ cai, eca WRITE ica
	RANGE gcabar, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER { :parameters that can be entered when function is called in cell-setup 
	gcabar = 0 (S/cm2)	: initialized conductance
	ki = 0.025 (mM)		:test middle point of inactivation fct
	zetam = -3.4
	zetah = 2
	vhalfm = -21 (mV)
	vhalfh = -40 (mV)
	cst = 1(/ms) : scale factor
	tm0 = 1.5 (ms)
	th0 = 75 (ms)
	taumin = 2 (ms)		: minimal value of the time cst
}

ASSIGNED {: parameters needed to solve DE
	v (mV)
	celsius (degC)
	ica (mA/cm2)
	cai (mM)       :5e-5 initial internal Ca2+ concentration
	eca (mV)
	minf (1)
	hinf (1)
	sinf (1)
}

STATE {
	m
	h
	s
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar*pow(m, 2)*h*h2(cai) * (v - eca)
}

DERIVATIVE states {
	rates(v, cai)
	m' = (minf -m)/tm0
	h'= (hinf - h)/th0
	s' = (sinf - s)/taumin
}

INITIAL {
	rates(v, cai)
	m = minf
	h = hinf
	s = sinf
}

FUNCTION h2(ci (mM)) {
	h2 = ki/(ki + ci)
}

FUNCTION alpm(v (mV)) (/ms) {
	alpm = 1(/ms)*exp(zetam*(v-vhalfm)*FARADAY/(R*(273.16(degC) + celsius))) 
}

FUNCTION alph(v (mV)) (/ms) {
	alph = 1(/ms)*exp(zetah*(v-vhalfh)*FARADAY/(R*(273.16(degC) + celsius))) 
}

PROCEDURE rates(v (mV), cai (mM)) {
   	minf = cst/(cst + alpm(v))
	hinf = cst/(cst + alph(v))
	sinf = (ki/cai)^2 / ((ki/cai)^2 + 1)
}
