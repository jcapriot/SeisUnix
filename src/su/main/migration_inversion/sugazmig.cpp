/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUGAZMIGQ: $Revision: 1.21 $ ; $Date: 2011/11/16 22:14:43 $	    */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUGAZMIGQ - SU version of Jeno GAZDAG's phase-shift migration 	",
"	     for zero-offset data, with attenuation Q.			",
"									",
" sugazmig <infile >outfile vfile= [optional parameters]		",
"									",
" Optional Parameters:							",
" dt=from header(dt) or	.004	time sampling interval			",
" dx=from header(d2) or 1.0	midpoint sampling interval		",
" ft=0.0			first time sample			",
" ntau=nt(from data)	number of migrated time samples			",
" dtau=dt(from header)	migrated time sampling interval			",
" ftau=ft		first migrated time sample			",
" tmig=0.0		times corresponding to interval velocities in vmig",
" vmig=1500.0	interval velocities corresponding to times in tmig	",
" vfile=		name of file containing velocities		",
" Q=1e6			quality factor					",
" ceil=1e6		gain ceiling beyond which migration ceases	",
"									",
" verbose=0	verbose = 1 echoes information				",
"									",
" tmpdir= 	 if non-empty, use the value as a directory path	",
"		 prefix for storing temporary files; else if the	",
"	         the CWP_TMPDIR environment variable is set use		",
"	         its value for the path; else use tmpfile()		",
"									",
" Note: ray bending effects not accounted for in this version.		",
"									",
" The tmig and vmig arrays specify an interval velocity function of time.",
" Linear interpolation and constant extrapolation is used to determine	",
" interval velocities at times not specified.  Values specified in tmig	",
" must increase monotonically.						",
"									",
" Alternatively, interval velocities may be stored in a binary file	",
" containing one velocity for every time sample in the data that is to be",
" migrated.  If vfile is specified, then the tmig and vmig arrays are ignored.",
"									",
" Caveat: Adding Q is a first attempt to address GPR issues.		",
NULL};
#endif

/* 
 * Credits: 
 *  Constant Q attenuation correction by Chuck Oden 5 May 2004
 *  CWP John Stockwell 12 Oct 1992
 * 	Based on a constant v version by Dave Hale.
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */ 
/**************** end self doc ********************************/


/* Library version of SUGAZMIG (SUGAZMIGQ)
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_gazmig:	migrates the zero-offset section pin[nx][nt] (the traces in order) to qout[nx][ntau], by Gazdag's
 *			phase shift for the interval velocities vt[ntau] of the migrated times, with the quality factor
 *			qual (1e6: none) and the gain ceiling gainceil
 *
 * Returns 0, or -3 for parameters that are not allowed.
 */

static void gazdagvt (float k,
	int nt, float dt, float ft,
	int ntau, float dtau, float ftau,
	const float *vt, complex *p, complex *q, float qual, float gainceil)
{
	int ntfft,nw,it,itau,iw;
	float dw,fw,tmax,w,tau,phase,coss, gain, alpha;
	complex cshift;
	std::vector<complex> pp;
	std::vector<float> cumgain;

	/* determine frequency sampling */
	ntfft = npfa(nt);
	nw = ntfft;
	dw = 2.0*PI/(ntfft*dt);
	fw = -PI/dt;

	/* determine maximum time */
	tmax = ft+(nt-1)*dt;

	/* allocate workspace */
	pp.assign(nw, complex(0.0f,0.0f));
	cumgain.assign(nw, 1.0f);

	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it++)
		pp[it] = (it%2 ? -p[it] : p[it]);
	for (it=nt; it<ntfft; it++)
		pp[it] = complex(0.0f,0.0f);
	pfacc(1,ntfft,pp.data());

	/* account for non-zero ft and non-zero ftau */
	for (itau=0 ; itau < ftau && itau < ntau ; itau++) {
		for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
			if (w==0.0) w = 1e-10/dt;
			coss = 1.0-pow(0.5 * vt[itau] * k/w,2.0);
			if (coss>=pow(ftau/tmax,2.0)) {
				phase = w*(ft-ftau*sqrt(coss));
				cshift = complex(cos(phase),sin(phase));
				pp[iw] = pp[iw]*cshift;
			} else {
				pp[iw] = complex(0.0f,0.0f);
			}
		}
	}

	/* loop over migrated times tau */
	for (itau=0,tau=ftau; itau<ntau; itau++,tau+=dtau) {

		/* initialize migrated sample */
		q[itau] = complex(0.0f,0.0f);

		/* loop over frequencies w */
		for (iw=0,w=fw; iw<nw; iw++,w+=dw) {

			/* accumulate image (summed over frequency) */
			q[itau] = q[itau]+pp[iw];

			/* compute cosine squared of propagation angle */
			if (w==0.0) w = 1e-10/dt;
			coss = 1.0-pow(0.5 * vt[itau] * k/w,2.0);

			/* if wave could have been recorded in time */
			if (coss>=pow(tau/tmax,2.0)) {

				/* extrapolate down one migrated time step */
				phase = -w*dtau*sqrt(coss);
				cshift = complex(cos(phase),sin(phase));

				/* apply gain until gain ceiling is reached */
				if (cumgain[iw] < gainceil) {
					alpha = w/(2.0*vt[itau]*qual);
					gain = exp(fabs(0.5*vt[itau]*dtau*alpha));
					pp[iw] = pp[iw]*(cshift*gain);
					cumgain[iw] *= gain;
				} else {
					pp[iw] = complex(0.0f,0.0f);
				}

			/* else, if wave couldn't have been recorded in time */
			} else {

				/* zero the wave */
				pp[iw] = complex(0.0f,0.0f);
			}
		}

		/* scale accumulated image just as we would for an FFT */
		q[itau] = q[itau]*(1.0f/nw);
	}
}

int su_gazmig(int nx, int nt, float dt, float ft, float dx, int ntau, float dtau, float ftau, const float *vt,
	float qual, float gainceil, const float *pin, float *qout)
{
	int ix,it,itau,ik,nxfft,nk;
	float dk,k,fk;
	float **p,**q;
	complex **cp,**cq;

	if (nx<1 || nt<2 || ntau<1 || dt<=0.0 || dx<=0.0 || dtau<=0.0 || qual==0.0) return -3;
	for (itau=0; itau<ntau; itau++) if (vt[itau]<=0.0) return -3;

	/* determine wavenumber sampling (for real to complex FFT) */
	nxfft = npfar(nx);
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);
	fk = 0.0;

	/* allocate space */
	p = alloc2float(nt,nxfft);
	q = alloc2float(ntau,nxfft);
	cp = alloc2complex(nt,nk);
	cq = alloc2complex(ntau,nk);

	for (ix=0; ix<nx; ix++)
		for (it=0; it<nt; it++) p[ix][it] = pin[(size_t)ix*nt+it];

	/* pad with zeros and Fourier transform x to k */
	for (ix=nx; ix<nxfft; ix++)
		for (it=0; it<nt; it++)
			p[ix][it] = 0.0;
	pfa2rc(-1,2,nt,nxfft,p[0],cp[0]);

	/* migrate each wavenumber */
	for (ik=0,k=fk; ik<nk; ik++,k+=dk)
		gazdagvt(k,nt,dt,ft,ntau,dtau,ftau,vt,cp[ik],cq[ik],qual,gainceil);

	/* Fourier transform k to x (including FFT scaling) */
	pfa2cr(1,2,ntau,nxfft,cq[0],q[0]);
	for (ix=0; ix<nx; ix++)
		for (itau=0; itau<ntau; itau++)
			qout[(size_t)ix*ntau+itau] = q[ix][itau]/nxfft;

	free2float(p);
	free2float(q);
	free2complex(cp);
	free2complex(cq);
	return 0;
}
