/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFXDECON: $Revision: 1.10 $ ; $Date: 2011/11/16 17:47:47 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>
#include <algorithm>

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUFXDECON - random noise attenuation by FX-DECONvolution              ",
"									",
" sufxdecon <stdin >stdout [...]	                                ",
"									",
" Required Parameters:							",
"									",
" Optional Parameters:							",
" taper=.1	length of taper                                         ",
" fmin=6.       minimum frequency to process in Hz  (accord to twlen)   ",
" fmax=.6/(2*dt)  maximum frequency to process in Hz                    ",
" twlen=entire trace  time window length (minimum .3 for lower freqs)   ",
" ntrw=10       number of traces in window                              ",
" ntrf=4        number of traces for filter (smaller than ntrw)         ",
" verbose=0	=1 for diagnostic print					",
" tmpdir=	if non-empty, use the value as a directory path	prefix	",
"		for storing temporary files; else, if the CWP_TMPDIR	",
"		environment variable is set, use its value for the path;",
"		else use tmpfile()					",
"									",
" Notes: Each trace is transformed to the frequency domain.             ",
"        For each frequency, Wiener filtering, with unity prediction in ",
"        space, is used to predict the next sample.                     ",
"        At the end of the process, data is mapped back to t-x domain.  ",
"									",
"                                                                       ",
NULL};
#endif

/* Credits:
 *
 *	CWP: Carlos E. Theodoro (10/07/97)
 *
 * References:
 *		Canales(1984):'Random noise reduction' 54th. SEGM
 *		Gulunay(1986):'FXDECON and complex Wiener Predicition
 *                             filter' 56th. SEGM
 *		Galbraith(1991):'Random noise attenuation by F-X
 *                             prediction: a tutorial' 61th. SEGM
 *
 * Algorithm:
 *	- read data
 *	- loop over time windows
 *		- select data
 *		- FFT (t -> f)
 *		- loop over space windows
 *			- select data
 *			- loop over frequencies
 *				- autocorelation
 *				- matrix problem
 *				- construct filter
 *				- filter data
 *			- loop along space window
 *				- FFT (f -> t)
 *				- reconstruct data
 * 	- output data
 *
 * Trace header fields accessed: ns, dt, d1
 * Trace header fields modified:
 */

/**************** end self doc ***********************************/

#define LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define PFA_MAX	720720	/* Largest allowed nfft	          */

/* complex convolution */
static void cconv (int lx, int ifx, complex *x,
       	    int ly, int ify, complex *y,
	    int lz, int ifz, complex *z)
{
        int ilx=ifx+lx-1,ily=ify+ly-1,ilz=ifz+lz-1,i,j,jlow,jhigh;
	complex sum;

	x -= ifx;  y -= ify;  z -= ifz;
	for (i=ifz; i<=ilz; ++i) {
	        jlow = i-ily;  if (jlow<ifx) jlow = ifx;
		jhigh = i-ify;  if (jhigh>ilx) jhigh = ilx;
		for (j=jlow,sum=complex(0.0f,0.0f); j<=jhigh; ++j){
		        sum += x[j]*y[i-j];
		}
		z[i] = sum;
	}
}

/* complex correlation */
static void cxcor (int lx, int ifx, complex *x,
            int ly, int ify, complex *y,
            int lz, int ifz, complex *z)
{
        int i,j;
	std::vector<complex> xr(lx);
	for (i=0,j=lx-1; i<lx; ++i,--j)
	        xr[i] = std::conj(x[j]);
	cconv(lx,1-ifx-lx,xr.data(),ly,ify,y,lz,ifz,z);
}

/* The algorithm of the program, on a gather that is in memory.

   nsp, ntrc	samples per trace, number of traces
   tdatain	[ntrc][nsp] the traces
   tdataout	[ntrc][nsp] the filtered traces (all of the samples are written)
   dt		sample interval (s)
   taper	length of the taper (s), more than 0
   fmin, fmax	range of frequencies (Hz) that are filtered
   twlen	time window length (s)
   ntrw		number of traces in a window (at least 1)
   ntrf		number of traces of the filter (at least 1)

   Returns 0, or 1 if the transform is too long.

   Differences from the program, which are bugs of it: at least 1 time window and 1 space window are used (the program
   gives all zeros if rounding makes either of them 0), and the traces that a window takes from outside of itself are
   the nearest ones of the gather (for a single window of more than ntrw traces the program takes the last trace of the
   gather for the traces after the first ntrw). */
int su_fxdecon(int nsp, int ntrc, const float *tdatain, float *tdataout, float dt, float taper, float fmin,
	float fmax, float twlen, int ntrw, int ntrf)
{
	int nspw, nspws, nspwf, nspwi, nstaper, ntrwu, ntw, nwx, nfft, nf;
	int jx, jw, jr, itt, it, iw, ifq, ir, ig, ifv;
	float d1, dummy;

	/* setting taper and spatial and temporal windows */
	nstaper = NINT(taper/dt);
	nstaper = (nstaper%2 ? nstaper+1 : nstaper);

	ntw = NINT((nsp-1)*dt/twlen);
	if (ntw < 1) ntw = 1;
	if (ntw==1) taper = 0.0;

	nspws = NINT(twlen/dt) + 1;
	nspwf = nspws + nstaper/2;
	nspwi = nspws + nstaper;

	nwx = ntrc/ntrw;
	if (nwx < 1) nwx = 1;

	/* Set up pfa fft */
	nfft = npfaro(2*nspwf, LOOKFAC * 2*nspwf);
	if (nfft < nsp) nfft = npfaro(nsp, LOOKFAC * nsp);
	if (nfft >= PFA_MAX) return 1;

	nf = nfft/2 + 1;
	d1 = 1.0/(nfft*dt);

	/* space allocation */
	std::vector<complex> ct(nfft), sfreq(2*ntrw+2*ntrf), autocorr(ntrf+1), fvector(2*ntrf+1), sfreqout(2*ntrw),
		freqv(nfft);
	std::vector<float> tidataw(nfft), ttidataw(nfft), todataw(nfft), ttodataw(nfft), gvector(2*ntrf),
		rmat((size_t) 4*ntrf*ntrf);
	std::vector<complex> fdata((size_t) nf*ntrc), fdataw((size_t) nf*(2*ntrw+2*ntrf)), ffreq((size_t) nf*2*ntrw);
	std::vector<float*> rmatrix(2*ntrf);
	std::vector<int> ipvt(2*ntrf);
	const complex zero(0.0f, 0.0f);
	for (ir=0;ir<2*ntrf;ir++) rmatrix[ir] = &rmat[(size_t) ir*2*ntrf];

	/* zero output */
	memset((void *) tdataout, 0, (size_t) ntrc*nsp*FSIZE);

	/* loop over time windows */
	for (iw=0;iw<ntw;iw++) {

	        if (iw>0 && iw<ntw-1) nspw=nspwi;
		else if (iw==0) {
		        if (ntw>1) nspw = nspwf;
			else        nspw = nsp;
		} else
		        nspw = nsp - nspws*iw + nstaper/2;

		/* zero fdata */
		std::fill(fdata.begin(), fdata.end(), zero);

		/* select data */
		for (jw=0;jw<ntrc;jw++) {
		        const float *trace = tdatain + (size_t) jw*nsp;

		        if (iw>0)
			   for (it=0;it<nspw;it++)
			     tidataw[it]=trace[it + iw*nspws - nstaper/2];
			else
			        for (it=0;it<nspw;it++)
				         tidataw[it]=trace[it];

			memset((void *) (tidataw.data() + nspw), 0, (nfft-nspw)*FSIZE);
			std::fill(ct.begin(), ct.end(), zero);

			/* FFT from t to f */
			for (it=0;it<nfft;it++)
			        ttidataw[it]=(it%2 ? -tidataw[it] : tidataw[it]);
			pfarc(1, nfft, ttidataw.data(), ct.data());

			/* Store values */
			for (ifq = 0; ifq < nf; ifq++) {
			        fdata[(size_t) jw*nf + ifq] = ct[nf-1-ifq];
			}
		}

		/* Loop over space windows */
		for (jx=0;jx<nwx;jx++){

		        /* to take care of a possible incomplete last window */
		        if (ntrc<jx*ntrw+2*ntrw)
			        ntrwu = ntrc - jx*ntrw;
			else
			        ntrwu = ntrw;

			/* select data: the traces of the window, and ntrf on each side (the nearest ones, for those
			   that are not there) */
			for (jw=0;jw<ntrwu+2*ntrf;jw++) {
			        int jt = jw + jx*ntrw - ntrf;
				if (jt < 0) jt = 0;
				if (jt > ntrc-1) jt = ntrc-1;
				memcpy((void *) &fdataw[(size_t) jw*nf], (const void *) &fdata[(size_t) jt*nf],
					nf*sizeof(complex));
			}

			/* loop over frequencies */
			for (ifq=0;ifq<nf;ifq++) {

			        if ((float)ifq*d1>=fmin && (float)ifq*d1<=fmax) {

				       for (jw=0;jw<ntrwu+2*ntrf;jw++)
					       sfreq[jw]=fdataw[(size_t) jw*nf + ifq];

				       std::fill(autocorr.begin(), autocorr.end(), zero);

				       /* complex autocorrelation */
				       cxcor(ntrwu,0,sfreq.data(),ntrwu,0,sfreq.data(),ntrf+1,0,autocorr.data());

				       if (autocorr[0].real() == 0.0) {
					       /* all of the data is 0: so is the prediction (the matrix is singular) */
					       for (jw=0;jw<ntrwu;jw++) ffreq[(size_t) jw*nf + ifq]=zero;
					       continue;
				       }

				       /* zeroing */
				       std::fill(gvector.begin(), gvector.end(), 0.0f);
				       std::fill(fvector.begin(), fvector.end(), zero);
				       std::fill(rmat.begin(), rmat.end(), 0.0f);

				       /* matrix problem */
				       for (ir=0;ir<ntrf;ir++)
					       for (jr=0;jr<ntrf;jr++) {
						       if (ir>=jr) rmatrix[ir][jr]=autocorr[ir-jr].real();
						       else        rmatrix[ir][jr]=autocorr[jr-ir].real();
					       }

				       for (ir=ntrf;ir<2*ntrf;ir++)
					       for (jr=0;jr<ntrf;jr++) {
						       if (ir-ntrf<jr) rmatrix[ir][jr]=-autocorr[jr-ir+ntrf].imag();
						       else            rmatrix[ir][jr]= autocorr[ir-jr-ntrf].imag();
					       }

				       for (ir=ntrf;ir<2*ntrf;ir++)
					       for (jr=ntrf;jr<2*ntrf;jr++)
						       rmatrix[ir][jr]=rmatrix[ir-ntrf][jr-ntrf];

				       for (ir=0;ir<ntrf;ir++)
					       for (jr=ntrf;jr<2*ntrf;jr++)
						       rmatrix[ir][jr]=-rmatrix[ir+ntrf][jr-ntrf];

				       for (ig=0;ig<2*ntrf;ig++) {
					       if (ig<ntrf) gvector[ig]=autocorr[ig+1].real();
					       else gvector[ig]=autocorr[ig-ntrf+1].imag();
				       }

				       LU_decomposition(2*ntrf,rmatrix.data(),ipvt.data(),&dummy);
				       backward_substitution(2*ntrf,rmatrix.data(),ipvt.data(),gvector.data());

				       /* construct filter */
				       for (ifv=0,ig=ntrf-1;ifv<ntrf;ifv++,ig--)
					       fvector[ifv]=std::conj(complex(gvector[ig]/2.,gvector[ig+ntrf]/2.));

				       for (ifv=ntrf+1,ig=0;ifv<2*ntrf+1;ifv++,ig++)
					       fvector[ifv]=complex(gvector[ig]/2.,gvector[ig+ntrf]/2.);

				       std::fill(sfreqout.begin(), sfreqout.end(), zero);

				       /* convolution of data with filter */
				       /* output is one sample ahead */
				       cconv(ntrwu+2*ntrf,-ntrf,sfreq.data(),2*ntrf+1,-ntrf,fvector.data(),ntrwu,0,sfreqout.data());

				       /* store filtered values */
				       for (jw=0;jw<ntrwu;jw++) ffreq[(size_t) jw*nf + ifq]=sfreqout[jw];
				}
			} /* end of frequencies loop */

			/* loop along space windows */
			for (jw=0;jw<ntrwu;jw++) {
			        float *out = tdataout + (size_t) (jx*ntrw+jw)*nsp;

			        /* select data */
			        for (ifq=0,itt=nf-1;ifq<nf;ifq++,itt--)
				        freqv[ifq] = ffreq[(size_t) jw*nf + itt];

				std::fill(freqv.begin()+nf, freqv.end(), zero);
				std::fill(todataw.begin(), todataw.end(), 0.0f);

				/* FFT back from f to t and scaling */
				pfacr(-1, nfft, freqv.data(), todataw.data());
				for (it=0;it<MIN(nsp,nfft);it++)
					todataw[it]/=nfft;
				for (it=0;it<MIN(nsp,nfft);it++)
				        ttodataw[it]=(it%2 ? -todataw[it] : todataw[it]);

				/*loop along time */
				if (ntw>1) {
				        /* first portion of time window */
				        if (iw>0)
					        for (it=0;it<nstaper;it++)
						        out[it+iw*nspws-nstaper/2]+=
							        ttodataw[it]*((float)(it)*dt/taper);
					else
					        for (it=0;it<nstaper;it++)
						        out[it]=ttodataw[it];

					/* intermediate portion of time window */
					if (iw>0)
					        for (it=nstaper;it<nspw-nstaper;it++)
						        out[it+iw*nspws-nstaper/2]=ttodataw[it];
					else
					        for (it=nstaper;it<nspw-nstaper;it++)
						        out[it]=ttodataw[it];

					/* last portion of time window */
					if (iw>0 && iw<ntw-1)
					        for (it=nspw-nstaper;it<nspw;it++)
						        out[it+iw*nspws-nstaper/2]+=
							        ttodataw[it]*(1.-((float)(it-nspw+nstaper))*dt/taper);
					else if (iw==ntw-1)
					        for (it=nspw-nstaper;it<nspw;it++)
						        out[it+iw*nspws-nstaper/2]=ttodataw[it];
					else
					        for (it=nspw-nstaper;it<nspw;it++)
						        out[it]+=ttodataw[it]*(1.-((float)(it-nspw+nstaper))*dt/taper);
				}
				else {
				        for (it=0;it<nsp;it++)
					        out[it]=ttodataw[it];
				}

			} /* end loop along the traces of the window */

		} /* end loop over space windows */

	} /* end of time windows loop */

	return 0;
}
