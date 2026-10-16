/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGTK: $Revision: 1.17 $ ; $Date: 2011/11/16 22:14:43 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUMIGTK - MIGration via T-K domain method for common-midpoint stacked data",
"									",
" sumigtk <stdin >stdout dxcdp= [optional parms]			",
"									",
" Required Parameters:							",
" dxcdp                   distance between successive cdps		",
"									",
" Optional Parameters:							",
" fmax=Nyquist            maximum frequency				",
" tmig=0.0                times corresponding to interval velocities in vmig",
" vmig=1500.0             interval velocities corresponding to times in tmig",
" vfile=                  binary (non-ascii) file containing velocities v(t)",
" nxpad=0                 number of cdps to pad with zeros before FFT	",
" ltaper=0                length of linear taper for left and right edges", 
" verbose=0               =1 for diagnostic print			",
"									",
" tmpdir= 	 if non-empty, use the value as a directory path	",
"		 prefix for storing temporary files; else if the	",
"	         the CWP_TMPDIR environment variable is set use		",
"	         its value for the path; else use tmpfile()		",
" 									",
" Notes:								",
" Input traces must be sorted by either increasing or decreasing cdp.	",
" 									",
" The tmig and vmig arrays specify an interval velocity function of time.",
" Linear interpolation and constant extrapolation is used to determine	",
" interval velocities at times not specified.  Values specified in tmig	",
" must increase monotonically.						",
"									",
" Alternatively, interval velocities may be stored in a binary file	",
" containing one velocity for every time sample.  If vfile is specified,",
" then the tmig and vmig arrays are ignored.				",
"									",
" The time of first sample is assumed to be zero, regardless of the value",
" of the trace header field delrt.					",
" 									",
" The migration is a reverse time migration in the (t,k) domain. In the	",
" first step, the data g(t,x) are Fourier transformed x->k into the	",	
" the time-wavenumber domain g(t,k).					",
"									",
" Then looping over wavenumbers, the data are then reverse-time		",
" finite-difference migrated, wavenumber by wavenumber.  The resulting	",
" migrated data m(tau,k), now in the tau (migrated time) and k domain,	",
" are inverse fourier transformed back into m(tau,xout) and written out.",	
"									",
NULL};
#endif

/* Credits:
 *	CWP: Dave Hale
 *
 * Trace header fields accessed:  ns and dt
 */
/**************** end self doc *******************************************/

static void mig1k (float k, float fmax, float speed, int nt, float dt, 
	float *v, complex *p, complex *q)
/*****************************************************************************
migration in t-k domain for one wavenumber
******************************************************************************
Input:
k		wavenumber
fmax		maximum frequency
speed		speed parameter - >>1.0 for lots of dispersion
nt		number of time samples
dt		time sampling interval
v		array[nt] containing interval velocities v(t)
p		array[nt] containing input p(t;k)

Output:
q		array[nt] containing migrated q(t;k)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/05/90
*****************************************************************************/
{
	int nfac=5,ns,is,it,i1,i2,i1stop;
	float fac=0.75,ds,g00r,g00i,g01r,g01i,g10r,g10i,g11r,g11i,
		temp,vmin,vmax,ct,
		*vs,*cs,*cp,*grp,*gip;
	complex czero=complex(0.0f,0.0f),*gs;
	
	/* determine time sampling to avoid excessive grid dispersion */
	for (it=1,vmin=vmax=v[0]; it<nt; ++it) {
		if (v[it]<vmin) vmin = v[it];
		if (v[it]>vmax) vmax = v[it];
	}
	if (k!=0.0) {
		ds = fac*MAX(0.4*PI/ABS(vmax*k),0.1*vmin/(vmax*fmax));
		ds *= speed;
		ds = MIN(ds,dt);
	} else {
		ds = dt;
	}
	ns = 1+(nt-1)*dt/ds;
	ns = MIN(ns,nfac*nt);
	ds = (nt-1)*dt/(ns-1);
	
	/* allocate workspace */
	vs = ealloc1float(ns);
	cs = ealloc1float(ns*2+1) + 1;
	gs = ealloc1complex(ns);
	
	/* resample v(t) and p(t) */
	ress8r(nt,dt,0.0,v,v[0],v[nt-1],ns,ds,0.0,vs);
	ress8c(nt,dt,0.0,p,czero,czero,ns,ds,0.0,gs);

        /* compute finite-difference coefficients */
	for (is=0; is<ns; is++) {
		temp = 0.125*vs[is]*k*ds;
		temp = temp*temp;
		temp = (1.0-temp)/(1.0+temp);
		cs[2*is] = cs[2*is+1] = temp;
	}

	/* loop over t2 = (tau-t)/sqrt(2) */
	for (i2=2-ns; i2<ns; i2++) {

		/* determine t1 stop index */
		i1stop = (i2<=0)?1-i2:i2;

		/* initialize finite-difference star and coefficient */
		g00r = g00i = g01r = g01i = 0.0;
		grp = (float*)(&gs[ns-1]);
		gip = grp+1;
		cp = &cs[i2+ns-2];
		ct = *cp--;

		/* loop over t1 = (tau+t)/sqrt(2) */
		for (i1=ns-1; i1>=i1stop; i1--) {

			/* update real part */
			g10r = g00r;
			g11r = g01r;
			g00r = *grp;
			*grp = g01r = ct*(g11r+g00r)-g10r;

			/* update imaginary part */
			g10i = g00i;
			g11i = g01i;
			g00i = *gip;
			*gip = g01i = ct*(g11i+g00i)-g10i;

			/* update pointers and finite-difference coefficient */
			grp -= 2; gip -= 2;
			ct = *cp--;
		}
	}
	
	/* resample q(t) */
	ress8c(ns,ds,0.0,gs,czero,czero,nt,dt,0.0,q);
	
	/* free workspace */
	free1float(vs);
	free1float(cs-1);
	free1complex(gs);
}



/* Library version of SUMIGTK
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migtk:	migrates the stacked section data[nx][nt] (the traces in order of the cdps) to out[nx][nt] by the
 *			t-k method, with the interval velocities vt[nt] of the times of the samples
 *
 * The coefficients of the finite differences were read one before their array by the last step of each sweep (the value was
 * not used): they have a spare element in front. The program printed the number of steps of every wavenumber, which is not
 * here. Returns 0, or -3 for parameters that are not allowed.
 */
int su_migtk(int nx, int nt, float dt, float dx, float fmax, int nxpad, int ltaper, float speed, const float *vt,
	const float *data, float *out)
{
	int nxfft,nk,ix,it,ik;
	float dk,taper,k,fftscl;
	float **gtx;
	complex **gtk;
	std::vector<float> v(vt, vt+nt);

	if (nx<1 || nt<2 || dt<=0.0 || dx<=0.0 || nxpad<0 || ltaper<0 || speed<=0.0 || fmax<=0.0) return -3;
	for (it=0; it<nt; ++it) if (vt[it]<=0.0) return -3;

	/* determine wavenumber sampling */
	nxfft = npfaro(nx+nxpad,2*(nx+nxpad));
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);

	/* allocate space for Fourier transform */
	gtk = ealloc2complex(nt,nk);
	gtx = (float**)ealloc1(nxfft,sizeof(float*));
	for (ix=0; ix<nxfft; ++ix)
		gtx[ix] = (float*)gtk[0]+ix*nt;

	/* apply fft scaling to traces and pad with zeros */
	fftscl = 1.0/nxfft;
	for (ix=0; ix<nx; ++ix) {
		for (it=0; it<nt; ++it)
			gtx[ix][it] = data[(size_t)ix*nt+it]*fftscl;
		if (ix<ltaper) {
			taper = (float)(ix+1)/(float)(ltaper+1);
			for (it=0; it<nt; ++it)
				gtx[ix][it] *= taper;
		} else if (ix>=nx-ltaper) {
			taper = (float)(nx-ix)/(float)(ltaper+1);
			for (it=0; it<nt; ++it)
				gtx[ix][it] *= taper;
		}
	}
	for (ix=nx; ix<nxfft; ++ix)
		for (it=0; it<nt; ++it)
			gtx[ix][it] = 0.0;

	/* Fourier transform g(t,x) to g(t,k) */
	pfa2rc(-1,2,nt,nxfft,gtx[0],gtk[0]);

	/* loop over wavenumbers */
	for (ik=0,k=0.0; ik<nk; ++ik,k+=dk)
		mig1k(k,fmax,speed,nt,dt,v.data(),gtk[ik],gtk[ik]);

	/* Fourier transform g(t,k) to g(t,x) */
	pfa2cr(1,2,nt,nxfft,gtk[0],gtx[0]);

	for (ix=0; ix<nx; ++ix)
		for (it=0; it<nt; ++it)
			out[(size_t)ix*nt+it] = gtx[ix][it];

	free1(gtx);
	free2complex(gtk);
	return 0;
}
