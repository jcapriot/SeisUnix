/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUIMP2D: $Revision: 1.22 $ ; $Date: 2015/06/02 20:15:23 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"							",
" SUIMP2D - generate shot records for a line scatterer	",
"           embedded in three dimensions using the Born	",
"	    integral equation				",							
"							",
" suimp2d [optional parameters] >stdout			",
"							",
" Optional parameters					",
"	nshot=1		number of shots			",
"	nrec=1		number of receivers		",
"	c=5000		speed				",
"	dt=.004		sampling rate			",
"	nt=256		number of samples		",
"	x0=1000		point scatterer location	",
"	z0=1000		point scatterer location	",
"	sxmin=0		first shot location		",
"	szmin=0		first shot location		",
"	gxmin=0		first receiver location		",
"	gzmin=0		first receiver location		",
"	dsx=100		x-step in shot location		",
"	dsz=0	 	z-step in shot location		",
"	dgx=100		x-step in receiver location	",
"	dgz=0		z-step in receiver location	",
"							",
" Example:						",
"	suimp2d nrec=32 | sufilter | supswigp | ...	",
"							",
NULL};
#endif


/* Credits:
 *	CWP: Norm Bleistein, Jack K. Cohen
 */

/* Theory: Use the 3D Born integral equation (e.g., Geophysics,
 * v51, n8, p1554(7)). Use 2-D delta function for alpha and do
 * remaining y-integral by stationary phase.
 *
 * Note: Setting a 2D offset in a single offset field beats the
 *       hell out of us.  We did _something_.
 *
 * Trace header fields set: ns, dt, tracl, tracr, fldr, sx, selev,
 *                          gx, gelev, offset
 */
/**************** end self doc ***************************/

/* Library version of SUIMP2D
 *
 * The main program is not part of the library. It loops over the shots and the receivers, and writes the traces with their
 * headers. Doing that is left to the caller, and what is here is the response at one receiver of a shot:
 *
 *	su_imp_nfft:		the length of the FFT that the traces are padded to (the same for SUIMP3D)
 *	su_imp2d_trace:		the response of the line scatterer
 *
 * All the arrays are the caller's, none are allocated here.
 */

#define IMP_LOOKFAC	2	/* Look ahead factor for npfaro	  */
#define IMP_PFA_MAX	720720	/* Largest allowed nfft	          */

/* The length of the FFT for traces of nt samples, or -1 if it would be too big */
int su_imp_nfft(int nt)
{
	int nfft = npfaro(nt, IMP_LOOKFAC * nt);

	if (nfft >= SU_NFLTS || nfft >= IMP_PFA_MAX) return -1;
	return nfft;
}

/* The response of a line scatterer embedded in three dimensions at a receiver, by the Born integral equation:
 * a spike at the time (rs + rg)/c, with the amplitude of the geometrical spreading, that is multiplied by omega^(3/2) in the
 * frequency domain.
 *
 * nt, dt: number of samples and the sample interval
 * nfft: the length of the FFT (su_imp_nfft)
 * c: speed. rs, rg: distances from the scatterer to the shot and to the receiver
 * tout[nt]: the times of the samples, it*dt
 * rt[nfft], ct[nfft/2+1]: scratch. data[nt]: the trace
 */
void su_imp2d_trace(int nt, float dt, int nfft, float c, float rs, float rg, const float *tout, float *rt, complex *ct,
	float *data)
{
	size_t ntsize = nt * FSIZE;
	size_t nzeros = (nfft - nt) * FSIZE;
	int nfby2p1 = nfft/2 + 1;
	float k, d, t, spread, i32;
	float amplitude[1];
	int i;

	/* Set the constant in the response amplitude including scale for inverse fft below */
	k = 1.0 / (4.0 * sqrt(2.0*c*dt*nfft) * c * dt * dt * nfft * nfft);

	memset((void *) data, 0, ntsize);
	d = rs + rg;
	t = d/c;
	spread = sqrt(rs*rg*d);
	amplitude[0] = k/spread;

	/* Distribute response over full trace*/
	ints8r(1,dt,t,amplitude,0,0,nt,(float *) tout,data);

	/* Load trace into rt (zero-padded) */
	memcpy((void *) rt, (const void *) data, ntsize);
	memset((void *)(rt + nt), 0, nzeros);

	/* FFT */
	pfarc(1, nfft, rt, ct);

	/* Formula requires multiplication by abs(omega) to 3/2 power */
	for (i = 0; i < nfby2p1; ++i) {
		i32 = i * sqrt((double) i);
		ct[i] *= i32;
	}

	/* Invert and take real part */
	pfacr(-1, nfft, ct, rt);
	memcpy((void *) data, (const void *) rt, ntsize);
}
