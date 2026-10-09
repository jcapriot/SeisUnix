/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUOP: $Revision: 1.32 $ ; $Date: 2013/06/03 18:23:54 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 								",
" SUOP - do unary arithmetic operation on segys 		",
" 								",
" suop <stdin >stdout op=abs					",
" 								",
" Required parameters:						",
"	none							",
"								",
" Optional parameter:						",
"	op=abs		operation flag				",
"			abs   : absolute value			",
"			avg   : remove average value		",
"			ssqrt : signed square root		",
"			sqr   : square				",
"			ssqr  : signed square			",
"			sgn   : signum function			",
"			exp   : exponentiate			",
"			sexp  : signed exponentiate		",
"			slog  : signed natural log		",
"			slog2 : signed log base 2		",
"			slog10: signed common log		",
"			cos   : cosine				",
"			sin   : sine				",
"			tan   : tangent				",
"			cosh  : hyperbolic cosine		",
"			sinh  : hyperbolic sine			",
"			tanh  : hyperbolic tangent		",
"			cnorm : norm complex samples by modulus ", 
"			norm  : divide trace by Max. Value	",
"			db    : 20 * slog10 (data)		",
"			neg   : negate value			",
"			posonly : pass only positive values	",
"			negonly : pass only negative values	",
"                       sum   : running sum trace integration   ",
"                       diff  : running diff trace differentiation",
"                       refl  : (v[i+1] - v[i])/(v[i+1] + v[i]) ",
"			mod2pi : modulo 2 pi			",
"			inv   : inverse				",
"			rmsamp : rms amplitude			",
"                       s2v   : sonic to velocity (ft/s) conversion     ",
"                       s2vm  : sonic to velocity (m/s) conversion     ",
"                       d2m   : density (g/cc) to metric (kg/m^3) conversion ",
"                       drv2  : 2nd order vertical derivative ",
"                       drv4  : 4th order vertical derivative ",
"                       integ : top-down integration            ",
"                       spike : local extrema to spikes         ",
"                       saf   : spike and fill to next spike    ",
"                       freq  : local dominant freqeuncy        ",
"                       lnza  : preserve least non-zero amps    ",
"                       --------- window operations ----------- ",
"                       mean  : arithmetic mean                 ",
"                       despike  : despiking based on median filter",
"                       std   : standard deviation              ",
"                       var   : variance                        ",
"       nw=21           number of time samples in window        ",
"                       --------------------------------------- ",
"			nop   : no operation			",
"								",
" Note:	Binary ops are provided by suop2.			",
" Operations inv, slog, slog2, and slog10 are \"punctuated\",	", 
" meaning that if, the input contains 0 values,			",
" 0 values are returned.					",	
"								",
" For file operations on non-SU format binary files use:  farith",
NULL};
#endif


/* Credits:
 *
 * CWP: Shuki Ronen, Jack K Cohen (c. 1987)
 *  Toralf Foerster: norm and db operations, 10/95.
 *  Additions by Reg Beardsley, Chris Liner, and others.
 *
 * Notes:
 *	If efficiency becomes important consider inverting main loop
 *      and repeating operation code within the branches of the switch.
 *
 *	Note on db option.  The following are equivalent:
 *	... | sufft | suamp | suop op=norm | suop op=slog10 |\
 *		sugain scale=20| suxgraph style=normal
 *
 *	... | sufft | suamp | suop op=db | suxgraph style=normal
 */
/**************** end self doc ***********************************/


/* Library version of SUOP
 *
 * The main program is not part of the library. suop is one long switch on the operation, and most of the cases
 * are plain array arithmetic (abs, sqr, exp, a cumulative sum, ...), which is left to the caller to do with
 * whatever it has for arrays. The operations that are more than that are here, each as a function of its own.
 *
 * Differences from the program:
 *  - The scratch arrays are arguments, and so owned by the caller. (The program has tmp and tmp1 with nt samples,
 *    but some operations write one sample past the end of them.)
 *  - Operations saf and freq finish by clearing tmp[nt], which is one past the end of the trace, where the last
 *    sample tmp[nt-1] is meant. The program never sets that sample (so it is whatever malloc left there).
 *  - Operation freq clears the part of its output that the program leaves unset.
 *  - Operation despike uses the last window instead of the first for the sample i == (nw-1)/2.
 */

static float mymedian(float *a, int nw);
static int floatcomp(const void* elem1, const void* elem2);

/* Operation saf: local extrema to spikes, and fill to the next spike
 *
 * data[nt] is the trace, tmp[nt] is scratch space
 */
void su_op_saf(float *data, int nt, float *tmp)
{
	float x1, x2, x3;
	int i;
	int iold;    /* index of last spike */
	int j;	     /* fill counter	    */
	float vold;  /* value of last spike */

	iold = 0;
	vold = 0.0;
	for (i = 1; i < nt-1; ++i) {
		x1 = data[i-1];
		x2 = data[i];
		x3 = data[i+1];

		/* local min or max */
		if ( ( (x1 < x2) && (x3 < x2) ) ||
		     ( (x1 > x2) && (x3 > x2) ) ){
			tmp[i] = x2;

			/* fill from last spike */
			for (j = iold; j < i; ++j ) {
				tmp[j] = vold;
			}

			/* reset old values */
			iold = i;
			vold = x2;

		} else {
			/* neither */
			tmp[i] = 0.0;
		}
	}
	/* edge effects */
	tmp[0] = 0.0;
	tmp[nt-1] = 0.0;	/* (the program has tmp[nt] here) */

	/* reload trace for output */
	for (i = 0; i < nt; ++i) {
		data[i] = tmp[i];
	}
}

/* Operation freq: local dominant frequency
 *
 * data[nt] is the trace, dt is the sample interval (in seconds), tmp[nt] and tmp1[nt] are scratch space
 */
void su_op_freq(float *data, int nt, float dt, float *tmp, float *tmp1)
{
	float x1, x2, x3, delay, freq;
	int i, iold, j;

	/* (the program leaves the part of tmp1 after the last minimum unset) */
	for (i = 0; i < nt; ++i) {
		tmp1[i] = 0.0;
	}

	/*
	  scan for extrema
	  local max:  (x2 > x1) and (x2 > x3)... assign +1
	  local min:  (x2 < x1) and (x2 < x3)... assign -1
	*/
	for (i = 1; i < nt-1; ++i) {
		x1 = data[i-1];
		x2 = data[i];
		x3 = data[i+1];

		/* local min or max */
		if ( (x2 > x1) && (x2 > x3) ) {
			tmp[i] =  1.0;
		} else if ( (x2 < x1) && (x2 < x3) ) {
			tmp[i] = -1.0;
		} else {
			/* neither */
			tmp[i] =  0.0;
		}
	}
	/* edge effects */
	tmp[0] = 0.0;
	tmp[nt-1] = 0.0;	/* (the program has tmp[nt] here) */

	/*
	  delay between consecutive local minima
	  is local dominant period, which is inverse
	  of local dom frequency
	*/

	/* assume first sample is a minima... kludge */
	iold = 0;

	/* skip a few samples to avoid false high freq at start */
	for (i = 5; i < nt; ++i) {

		if ( tmp[i] == -1.0 ) {

			/* calc delay from last min */
			delay = (i - iold) * dt;

			/* calc freq Hz
			   This delay is the period */
			freq  = 1 / ( delay);

			/* fill with this freq from last min to here */
			for (j = iold; j < i; ++j) {
				tmp1[j] = freq;
			}

			/* hold loc of this minima */
			iold = i;
		}
	}

	/* reload trace for output */
	for (i = 0; i < nt; ++i) {
		data[i] = tmp1[i];
	}
}

/* Operation despike: despiking based on a median filter
 *
 * data[nt] is the trace, nw is the number of samples in the window (odd, and at most nt), and tmp[nt] and
 * tomed[nw] are scratch space
 */
void su_op_despike(float *data, int nt, int nw, float *tmp, float *tomed)
{
	int i, j, mt;

	/* copy trace into temp array */
	for (i = 0; i < nt; ++i) {
		tmp[i] = data[i];
	}

	/* half width of window */
	mt = (nw-1)/2;

	/* loop over output times */
	for (i = 0; i < nt; ++i) {
		/* check we are not off the data ends (the program has i-mt > 0 here) */
		if (i-mt >= 0 && i+mt < nt)
		{
			for(j=0; j<nw; j++)
				tomed[j] = tmp[i-mt+j];
		}
		else
		{
			for(j=0; j<nw; j++)
				tomed[j] = (i-mt < 0) ? tmp[j] : tmp[nt-nw+j];
		}
		data[i] = mymedian(tomed,  nw);
	}
}

static int floatcomp(const void* elem1, const void* elem2)
{
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

static float mymedian(float *buff, int nw)
{
	qsort(buff, nw, FSIZE, floatcomp);
	if(nw %2)
		return  buff[nw/2];
 	else
		return  0.5*(buff[nw/2]+ buff[nw/2 + 1]);
}
