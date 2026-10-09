/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUICLOGFFT: $Revision: 1.10 $ ; $Date: 2015/11/17 22:35:40 $	*/

#include "su.h"
#include "segy.h"

#if 0 /* the program, which is not part of the library */
/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SUICLOGFFT - fft of complex log frequency traces to real time traces",
" 								",
"  suiclogftt <stdin >sdout sign=-1				",
" 								",
" Required parameters:						",
" 	none							",
" 								",
" Optional parameter:						",
" 	sign=-1		sign in exponent of inverse fft		",
"	sym=0		=1 center  output 			",
" Output traces are normalized by 1/N where N is the fft size.	",
" 								",
" Note:								",
" Nominally this is the inverse to the complex log fft, but	",
" suclogfft | suiclogfft is not quite a no-op since the trace	",
" 	length will usually be longer due to fft padding.	",
" 								",
" 								",
" Example:  Homomorphic wavelet estimation                              ",
" suclogfft < shotgather | suamp mode=real | sustack key=dt > real.su   ",
" suclogfft < shotgather | suamp mode=imag | sustack key=dt > imag.su   ",
" suop2 real.su imag.su op=zipper | suiclogfft | suminphase > wavelet.su",
" 								",
" 								",
" 								",
NULL};

/* Credits:
 * 
 *   CWP: John Stockwell, Dec 2010 based on
 *     suifft.c by:
 *	CWP: Shuki Ronen, Chris Liner, Jack K. Cohen,  c. 1989
 *
 * Trace header fields accessed: ns, trid
 * Trace header fields modified: ns, trid
 */
/**************** end self doc ***********************************/


#define PFA_MAX	720720		/* Largest allowed fft	*/

segy tr;

int
main(int argc, char **argv)
{
	complex *ct;	/* complex input trace			*/
	float *rt;	/* real output trace			*/
	int nfft;		/* fft size 				*/
	int nf;			/* number of frequencies		*/
	int sign;		/* sign in exponent of transform	*/
	float onfft;		/* 1.0/nfft				*/
	float newd1;		/* reset time domain sampling		*/

	int sym=0;		/* =1 symmetric output			*/



	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Get info from first trace */ 
	if (!gettr(&tr))  err("can't get first trace");
	if (tr.trid != FUNPACKNYQ)
		err("input not complex freq data, trid=%d", tr.trid);
	nfft = tr.ns - 2; /* see sufft */
	nf = nfft/2 + 1;
	onfft = 1.0/nfft;
	if(tr.d1) {
	    newd1 = onfft/tr.d1;
	} else {
	   if(tr.dt) newd1 = (float) (((double) tr.dt)/1000000.0);
	   else newd1 = 0.0f;
	}

	/* Set sign in exponent of transform */
	if (!getparint   ("sym", &sym)) sym = 0;
	if (!getparint   ("sign", &sign)) sign = -1;

        checkpars();

	if (sign != 1 && sign != -1)  err("sign = %d must be 1 or -1", sign);


	/* Allocate fft arrays */
	ct   = ealloc1complex(nf);
	rt   = ealloc1float(nfft);


	/* Main loop over traces */
	do {
		int i;

		/* zero out array */
		memset( (void *) ct, 0, nf*sizeof(complex));
		memset( (void *) rt, 0, nfft*FSIZE);
		
		/* Load traces into ct (pfa fills in negative freqs) */
		/* exponentiate prior to inverse transform */
		for (i = 0; i < nf; ++i) {
			if (tr.data[2*i]) {
				ct[i].r = (float) exp(tr.data[2*i]) * cos(tr.data[2*i+1]);
				ct[i].i = (float) exp(tr.data[2*i]) * sin(tr.data[2*i+1]);

				/* if sym=1 then flip the sign on every */
				/* other value to make output symmetric */
                        	if (sym) {
                                	if (ISODD(i)) {
                                        	ct[i].r = -ct[i].r;
                                        	ct[i].i = -ct[i].i;
                                	}
                        	}
     
					
			} else {
			
					ct[i].r = 0.0;
					ct[i].i = 0.0;
				
			}
		}


		/* Inverse FFT */
		pfacr(sign, nfft, ct, rt);

		/* Load back and scale for inverse fft */
		for (i = 0; i < nfft; ++i)  
			tr.data[i] = rt[i] * onfft;

		tr.trid = TREAL;
		tr.ns = nfft;
		tr.f1 = 0.0f;
		tr.d1 = newd1;
		if (sym) tr.delrt=-(tr.delrt + (nfft*newd1*1000.0)/2.0);

		puttr(&tr);

	} while (gettr(&tr));


	return(CWP_Exit());
}
#endif

/* Library version of SUICLOGFFT (and SUICEPSTRUM)
 *
 * The exponential of the complex log spectrum of nf frequencies: exp(re) (cos im, sin im), 0 where the real part is 0, and
 * with the sign of every other frequency changed if sym (which centers the transformed trace).
 * The main program is not part of the library, and its Fourier transforms are not either: what it does between the
 * transforms is. The spectra are arrays of floats, the real and the imaginary part of each frequency next to each
 * other (the layout of a numpy complex64 array). The transforms are the caller's (the SU ones have exp(+i w t) in their
 * forward kernel, which is the conjugate of numpy's).
 */
void su_iclogfft_spectrum(int nf, const float *clog, int sym, float *ct)
{
	int i;
	float e;

	for (i = 0; i < nf; ++i) {
		if (clog[2*i]) {
			e = exp(clog[2*i]);
			ct[2*i] = e * cos(clog[2*i+1]);
			ct[2*i+1] = e * sin(clog[2*i+1]);
			if (sym && ISODD(i)) {
				ct[2*i] = -ct[2*i];
				ct[2*i+1] = -ct[2*i+1];
			}
		} else {
			ct[2*i] = 0.0;
			ct[2*i+1] = 0.0;
		}
	}
}
