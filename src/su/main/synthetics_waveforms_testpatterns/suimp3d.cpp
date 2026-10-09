/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUIMP3D: $Revision: 1.25 $ ; $Date: 2015/06/02 20:15:23 $	*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"							",
"SUIMP3D - generate inplane shot records for a point 	",
"          scatterer embedded in three dimensions using	",
"          the Born integral equation			",							
"							",
"suimp3d [optional parameters] >stdout 			",
"							",
"Optional parameters					",
"	nshot=1		number of shots			",
"	nrec=1		number of receivers		",
"	c=5000		speed				",
"	dt=.004		sampling rate			",
"	nt=256		number of samples		",
"	x0=1000		point scatterer location	",
"	y0=0		point scatterer location	",
"	z0=1000		point scatterer location	",
"   dir=0		do not include direct arrival	",
"	            =1 include direct arrival	",
"	sxmin=0		first shot location		",
"	symin=0		first shot location		",
"	szmin=0		first shot location		",
"	gxmin=0		first receiver location		",
"	gymin=0		first receiver location		",
"	gzmin=0		first receiver location		",
"	dsx=100		x-step in shot location		",
"	dsy=0	 	y-step in shot location		",
"	dsz=0	 	z-step in shot location		",
"	dgx=100		x-step in receiver location	",
"	dgy=0		y-step in receiver location	",
"	dgz=0		z-step in receiver location	",
"							",
" Example:                                              ",
"       suimp3d nrec=32 | sufilter | supswigp | ...     ",
"							",
NULL};
#endif


/* Credits:
 *	CWP: Norm Bleistein, Jack K. Cohen
 *  UHouston: Chris Liner 2010 (added direct arrival option)
 *
 */
 
/* Theory: Use the 3D Born integral equation (e.g., Geophysics,
 * v51, n8, p1554(7)). Use 3-D delta function for alpha.
 *
 * Note: Setting a 3D offset in a single offset field beats the
 *       hell out of us.  We did _something_.
 *
 * Trace header fields set: ns, dt, tracl, tracr, fldr, tracf,
 *                          sx, sy, selev, gx, gy, gelev, offset
 */
/**************** end self doc ***************************/

/* Library version of SUIMP3D
 *
 * The main program is not part of the library. It loops over the shots and the receivers, and writes the traces with their
 * headers. Doing that is left to the caller, and what is here is the response at one receiver of a shot (see also su_imp_nfft in
 * SUIMP2D, which gives the length of the FFT).
 *
 * All the arrays are the caller's, none are allocated here.
 */

/* The response of a point scatterer in three dimensions at a receiver, by the Born integral equation: a spike at the time
 * (rs + rg)/c, with the amplitude of the geometrical spreading, that is multiplied by omega^2 in the frequency domain; with
 * dir the direct arrival (at the time rd/c, rd the horizontal distance from the shot to the receiver) is added.
 *
 * nt, dt: number of samples and the sample interval
 * nfft: the length of the FFT (su_imp_nfft)
 * c: speed. rs, rg: distances from the scatterer to the shot and to the receiver. rd: the distance for the direct arrival
 * tout[nt]: the times of the samples, it*dt
 * temp[nt], rt[nfft], ct[nfft/2+1]: scratch. data[nt]: the trace
 */
void su_imp3d_trace(int nt, float dt, int nfft, float c, float rs, float rg, int dir, float rd, const float *tout, float *temp,
	float *rt, complex *ct, float *data)
{
	size_t ntsize = nt * FSIZE;
	size_t nzeros = (nfft - nt) * FSIZE;
	int nfby2p1 = nfft/2 + 1;
	float k, d, t, spread;
	float amplitude[1];
	int i;

	/* Set the constant in the response amplitude including scale for inverse fft below */
	k = 1.0 / (4.0 * c * c * dt * dt * dt * nfft * nfft * nfft);

	memset((void *) data, 0, ntsize);
	d = rs + rg;
	t = d/c;
	spread = rs*rg;
	amplitude[0] = k/spread;

	/* Distribute diffraction response over full trace */
	ints8r(1,dt,t,amplitude,0,0,nt,(float *) tout,data);

	if (dir == 1) {
		/* direct arrival distance, time, and spreading */
		t = rd/c;
		spread = rd;
		amplitude[0] = k/spread;

		/* Distribute direct response over temp trace, and add it in */
		ints8r(1,dt,t,amplitude,0,0,nt,(float *) tout,temp);
		for (i=0; i<nt; i++) data[i] = data[i] + temp[i];
	}

	/* Load trace into rt (zero-padded) */
	memcpy((void *) rt, (const void *) data, ntsize);
	memset((void *)(rt + nt), 0, nzeros);

	/* FFT */
	pfarc(1, nfft, rt, ct);

	/* Multiply by omega^2 */
	for (i = 0; i < nfby2p1; ++i)
		ct[i] *= (float) (i*i);

	/* Invert and take real part */
	pfacr(-1, nfft, ct, rt);
	memcpy((void *) data, (const void *) rt, ntsize);
}
