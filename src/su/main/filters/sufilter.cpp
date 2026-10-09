/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFILTER: $Revision: 1.23 $ ; $Date: 2011/11/12 00:09:00 $        */


#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUFILTER - applies a zero-phase, sine-squared tapered filter		",
"									",
" sufilter <stdin >stdout [optional parameters]         		",
"									",
" Required parameters:                                         		",
"       if dt is not set in header, then dt is mandatory        	",
"									",
" Optional parameters:							",
"       f=f1,f2,...             array of filter frequencies(HZ) 	",
"       amps=a1,a2,...          array of filter amplitudes		",
"       dt = (from header)      time sampling interval (sec)        	",
"	verbose=0		=1 for advisory messages		",
"									",
" Defaults:f=.10*(nyquist),.15*(nyquist),.45*(nyquist),.50*(nyquist)	",
"                        (nyquist calculated internally)		",
"          amps=0.,1.,...,1.,0.  trapezoid-like bandpass filter		",
"									",
" Examples of filters:							",
" Bandpass:   sufilter <data f=10,20,40,50 | ...			",
" Bandreject: sufilter <data f=10,20,30,40 amps=1.,0.,0.,1. | ..	",
" Lowpass:    sufilter <data f=10,20,40,50 amps=1.,1.,0.,0. | ...	",
" Highpass:   sufilter <data f=10,20,40,50 amps=0.,0.,1.,1. | ...	",
" Notch:      sufilter <data f=10,12.5,35,50,60 amps=1.,.5,0.,.5,1. |..	",
NULL};
#endif


/* Credits:
 *      CWP: John Stockwell, Jack Cohen
 *	CENPET: Werner M. Heigl - added well log support
 *
 * Possible optimization: Do assignments instead of crmuls where
 * filter is 0.0.
 *
 * Trace header fields accessed: ns, dt, d1
 */
/**************** end self doc ***********************************/


/* Library version of SUFILTER
 *
 * The main program is not part of the library. polygonalFilter() is the one in sufilter.c. The program builds its
 * filter from the first trace, here the caller builds it (for whatever number of points it transforms) and does the
 * filtering, which is a Fourier transform, a multiplication with the filter, and an inverse Fourier transform. Note
 * that the filter has the 1/nfft of the (unnormalized) pfa transforms in it.
 */

/* Prototype of function used internally */
void polygonalFilter(float *f, float *amps,
			int npoly, int nfft, float dt, float *filter, int *intfr);

void polygonalFilter(float *f, float *amps, int npoly,
				int nfft, float dt, float *filter, int *intfr)
/*************************************************************************
polygonalFilter -- polygonal filter with sin^2 tapering
**************************************************************************
Input:
f		array[npoly] of frequencies defining the filter
amps		array[npoly] of amplitude values
npoly		size of input f and amps arrays
dt		time sampling interval
nfft		number of points in the fft
intfr		scratch array of npoly ints (the original allocates it)

Output:
filter		array[nfft/2+1] filter values (the original has this as nfft, but it is only
		used for the nfft/2+1 frequencies up to the Nyquist)
**************************************************************************
Notes: Filter is to be applied in the frequency domain
**************************************************************************
Author:  CWP: John Stockwell   1992
*************************************************************************/
#define PIBY2   1.57079632679490
{
        int icount,ifs;		/* loop counting variables              */
	int taper=0;		/* flag counter				*/
        int nf;                 /* number of frequencies (incl Nyq)     */
        int nfm1;               /* nf-1                                 */
        float onfft;            /* reciprocal of nfft                   */
        float df;               /* frequency spacing (from dt)          */

        nf = nfft/2 + 1;
        nfm1 = nf - 1;
        onfft = 1.0 / nfft;

        /* Compute array of integerized frequencies that define the filter*/
        df = onfft / dt;
        for(ifs=0; ifs < npoly ; ++ifs) {
                intfr[ifs] = NINT(f[ifs]/df);
                if (intfr[ifs] > nfm1) intfr[ifs] = nfm1;
        }

	/* Build filter, with scale, and taper specified by amps[] values*/

	/* Do low frequency end first*/
	for(icount=0; icount < intfr[0] ; ++icount)
		filter[icount] = amps[0] * onfft;

	/* now do the middle frequencies */
	for(ifs=0 ; ifs<npoly-1 ; ++ifs){

	   if(amps[ifs] < amps[ifs+1]) {
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
		    float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
		    float s = sin(c*(icount - intfr[ifs] + 1));
		    float adiff = amps[ifs+1] - amps[ifs];
		    filter[icount] = (amps[ifs] + adiff*s*s) * onfft;
		}

	   } else if (amps[ifs] > amps[ifs+1]) {
		++taper;
		for(icount=intfr[ifs]; icount<=intfr[ifs+1]; ++icount) {
			   float c = PIBY2 / (intfr[ifs+1] - intfr[ifs] + 2);
                	   float s = sin(c*(intfr[ifs+1] - icount + 1));
			   float adiff = amps[ifs] - amps[ifs+1];
                	   filter[icount] = (amps[ifs+1] + adiff*s*s) * onfft;
		  }

	   } else
		if(!(taper)){
		for(icount=intfr[ifs]; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		} else {
		for(icount=intfr[ifs]+1; icount <= intfr[ifs+1]; ++icount)
		   	   filter[icount] = amps[ifs] * onfft;
		}
	}

	/* finally do the high frequency end */
	for(icount=intfr[npoly-1]+1; icount<nf; ++icount){
		filter[icount] = amps[npoly-1] * onfft;
	}
}
