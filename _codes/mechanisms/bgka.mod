TITLE Borg-Graham type generic K-A channel

NEURON {
	SUFFIX borgka
	USEION k READ ek WRITE ik
    RANGE gkabar, ek, ik, vhalfl, vhalfl
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	FARADAY = (faraday) (kilocoulombs)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	gkabar = .0 (S/cm2)
	vhalfn = -33.6 (mV)
	vhalfl = -83 (mV)
	a0l = 0.08 (/ms)
	a0n = 0.02 (/ms)
	zetan = -3
	zetal = 4
	gmn = 0.6
	gml = 1
	cst = 1(/ms) : scale factor
}

ASSIGNED {
	v (mV)
	ek (mV)
	celsius (degC)
	ik (mA/cm2)
    ninf
    linf      
    taul (ms)
    taun (ms)
    gka (S/cm2)
}

STATE {
	n
	l
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkabar*n*l*(v - ek)
}

DERIVATIVE states { 
	rates(v)
	n' = (ninf - n)/taun
	l' = (linf - l)/taul
}

INITIAL {
	rates(v)
	n = ninf
	l = linf
}

FUNCTION alpn(v (mV)) (/ms) {
	alpn = 1(/ms)*exp(zetan*(v - vhalfn)*FARADAY/(R*(273.16(degC) + celsius)))
}

FUNCTION betn(v (mV)) (/ms) {
	betn = 1(/ms)*exp(zetan*gmn*(v - vhalfn)*FARADAY/(R*(273.16(degC) + celsius)))
}

FUNCTION alpl(v (mV)) (/ms) {
	alpl = 1(/ms)*exp(zetal*(v - vhalfl)*FARADAY/(R*(273.16(degC) + celsius)))
}

FUNCTION betl(v (mV)) (/ms) {
	betl = 1(/ms)*exp(zetal*gml*(v - vhalfl)*FARADAY/(R*(273.16(degC) + celsius)))
}

PROCEDURE rates(v (mV)) { :callable from hoc
	LOCAL q10

	q10 = 3^((celsius - 30(degC))/10(degC))

	ninf = cst/(cst + alpn(v))
	taun = betn(v)/(q10*a0n*(1 + alpn(v)))

	linf = cst/(cst + alpl(v))
	taul = betl(v)/(q10*a0l*(1 + alpl(v)))
}
