/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUOCEXT: $Revision: 1.6 $ ; $Date: 2011/11/16 22:12:22 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUOCEXT - smaller Offset EXTrapolation via Offset Continuation        ",
"           method for common-offset gathers                            ",
"									",
" suocext <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= offextr= [...]	",
"									",
" Required Parameters:							",
" cdpmin=	minimum cdp (integer number) for which to apply DMO	",
" cdpmax=	maximum cdp (integer number) for which to apply DMO	",
" dxcdp=	distance between adjacent cdp bins (m)			",
" noffmix=	number of offsets to mix (see notes)			",
" offextr=	offset to extrapolate					",
"									",
" Optional Parameters:							",
" tdmo=0.0	times corresponding to rms velocities in vdmo (s)	",
" vdmo=1500.0	rms velocities corresponding to times in tdmo (m/s)	",
" sdmo=1.0	DMO stretch factor; try 0.6 for typical v(z)		",
" fmax=0.5/dt	maximum frequency in input traces (Hz)			",
" verbose=0	=1 for diagnostic print					",
" tmpdir=	if non-empty, use the value as a directory path	prefix	",
"		for storing temporary files; else if the CWP_TMPDIR	",
"		environment variable is set use	its value for the path;	",
"		else use tmpfile()					",
"									",
" Notes:								",
" Input traces should be sorted into common-offset gathers.  One common- ",
" offset gather ends and another begins when the offset field of the trace",
" headers changes. One common-offset gather usually is enough.		",
"									",
" The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
" the cdp location expressed in units of meters or feet.		",
"									",
" The number of offsets to mix (noffmix) should typically equal the ratio of",
" the shotpoint spacing to the cdp spacing.  This choice ensures that every",
" cdp will be represented in each offset mix.  Traces in each mix will	",
" contribute through DMO to other traces in adjacent cdps within that mix.",
"									",
" The tdmo and vdmo arrays specify a velocity function of time that is	",
" used to implement a first-order correction for depth-variable velocity.",
" The times in tdmo must be monotonically increasing.			",
"									",
" For each offset, the minimum time at which a non-zero sample exists is ",
" used to determine a mute time.  Output samples for times earlier than this", 
" mute time will be zeroed.  Computation time may be significantly reduced",
" if the input traces are zeroed (muted) for early times at large offsets.",
"                                                                       ",
" A term for better amplitude reconstruction was added to Hale's formulation.",
NULL};
#endif

/* Credits: Carlos E. Theodoro (modification of Hale's SUDMOFK program)
 *
 * Technical Reference:
 *	C. Theodoro & K. Larner, 1998
 *      Extrapolation of seismic data to small offsets (CWP-276). 
 *
 *	Dip-Moveout Processing - SEG Course Notes
 *	Dave Hale, 1988
 *
 *	Bleistein, Cohen & Jaramillo, 1997
 *      True amplitude transformation to zero offset of data from 
 *      curved reflectors (CWP-262). 
 *
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp.
 * Trace header fields modified:  offset.
 */
/**************** end self doc *******************************************/


static void maketu (float offset, int itmin, float fmax,
	int nt, float dt, float ft, float *vrms, float **uoftp,
	int *nup, float *dup, float *fup, float **tofup, float *tconp);

/*****************************************************************************
extrapolate one common offset gather to the offset offextr, in the f-k domain (dmooff of the program)
******************************************************************************
offset		source receiver offset
offextr		offset to extrapolate to
fmax		maximum frequency
sdmo		DMO stretch factor
nx		number of midpoints
dx		midpoint sampling interval
nt		number of time samples
dt		time sampling interval
ft		first time
vrms_in		array[nt] of rms velocities
ptx		array[nxfft+2][nt] for p(t,x), zero-padded for fft after the nx traces; the first nx traces are
		replaced by the extrapolated ones (nxfft is su_ocext_nxfft)
******************************************************************************
Author: Carlos E. Theodoro (modification of the dmooff of sudmofk by Dave Hale)
*****************************************************************************/
void su_ocext_offset(float offset, float offextr, float fmax, float sdmo, int nx, float dx, int nt, float dt, float ft,
	const float *vrms_in, float *ptx)
{
	int nxfft,itmin,nu,nufft,nw,nk,ix,iu,iw,ik,it,iwn,
		iwmin,iwmax,nupad,ikmax;
	float dw,dk,tcon,wwscl,scale,scales,kmax,
		amp,phase,fr,fi,pwr,pwi,deltoff,numerator,
		wmin,wmax,fftscl,du,fu,w,k,osdmo,*uoft,*tofu;
	float *vrms = (float*)ealloc1float(nt);
	complex czero(0.0f,0.0f),**ptk,*pu,*pw;

	memcpy(vrms,vrms_in,sizeof(float)*nt);

	/* number of cdps after padding for fft */
	nxfft = npfar(nx+(int)(0.5*ABS(offset/dx)));

	/* get minimum time of first non-zero sample */
	for (ix=0,itmin=nt; ix<nx; ++ix) {
		for (it=0; it<itmin && ptx[(size_t)ix*nt+it]==0.0; ++it);
		itmin = it;
	}

	/* if all zeros, simply return */
	if (itmin>=nt) { free1float(vrms); return; }

	/* make stretch and compress functions t(u) and u(t) */
	maketu(offset,itmin,fmax,nt,dt,ft,vrms,&uoft,&nu,&du,&fu,&tofu,&tcon);

	/* adjust DMO stretch factor for nominal error in log stretch; */
	/* solve sdmo*(sqrt(1-a/sdmo)-1) = 0.5*log(1-a), where a=0.5 */
	sdmo *= .62;

	/* inverse of dmo stretch factor */
	osdmo = 1.0/sdmo;

	/* maximum DMO shift (in samples) for any wavenumber k */
	nupad = 1.5*sdmo*tcon/du;

	/* frequency sampling */
	nufft = npfa(nu+nupad);
	nw = nufft;
	dw = 2.0*PI/(nufft*du);

	/* allocate workspace */
	pu = pw = ealloc1complex(nufft);

	/* wavenumber sampling and maximum wavenumber to apply dmo */
	nk = nxfft/2+1;
	dk = 2.0*PI/ABS(nxfft*dx);
	kmax = PI/ABS(dx);
	ikmax = NINT(kmax/dk);

	/* pointers to complex p(t,k) */
	ptk = (complex**)ealloc1(nk,sizeof(complex*));
	for (ik=0; ik<nk; ++ik)
		ptk[ik] = (complex*)ptx+(size_t)ik*nt;

	/* fft scale factor */
	fftscl = (float)nk/(float)(ikmax+1)/(nufft*nxfft);

	/* Fourier transform p(t,x) to p(t,k) */
	pfa2rc(-1,2,nt,nxfft,ptx,ptk[0]);

	/* loop over wavenumbers less than maximum */
	for (ik=0,k=0.0; ik<=ikmax; ++ik,k+=dk) {

		/* stretch p(t;k) to p(u) */
		ints8c(nt,dt,ft,ptk[ik],czero,czero,nu,tofu,pu);

		/* pad with zeros and Fourier transform p(u) to p(w) */
		for (iu=nu; iu<nufft; ++iu)
			pu[iu] = czero;
		pfacc(1,nufft,pu);

		/* minimum and maximum frequencies to process */
		wmin = ABS(0.5*vrms[0]*k);
		wmax = ABS(PI/du);
		iwmin = MAX(1,MIN((nw-1)/2,NINT(wmin/dw)));
		iwmax = MAX(0,MIN((nw-1)/2,NINT(wmax/dw)));

		/* difference between squared offsets */
		deltoff = pow(offset*.5,2.0)-pow(offextr*.5,2.0);

		/* constant independent of w */
		wwscl = osdmo*deltoff*pow(k/tcon,2.0);

		/* zero dc (should be zero anyway) */
		pw[0] = czero;

		/* zero frequencies below minimum */
		for (iw=1,iwn=nw-iw; iw<iwmin; ++iw,--iwn)
			pw[iw] = pw[iwn] = czero;

		/* do dmo between minimum and maximum frequencies */
		for (iw=iwmin,iwn=nw-iwmin,w=iwmin*dw;
			iw<=iwmax; ++iw,--iwn,w+=dw) {
			scales = 1.0+wwscl/(w*w);
			scale = sqrt(scales);
			numerator = 1+2*(osdmo*deltoff*pow(k/tcon,2.)/(w*w));
			phase = sdmo*w*tcon*(scale-1.0);
			amp = fftscl*(1.0-sdmo+sdmo*numerator/scale);
			fr = amp*cos(phase);
			fi = amp*sin(phase);
			pwr = pw[iw].real();
			pwi = pw[iw].imag();
			pw[iw] = complex(pwr*fr-pwi*fi, pwr*fi+pwi*fr);
			pwr = pw[iwn].real();
			pwi = pw[iwn].imag();
			pw[iwn] = complex(pwr*fr+pwi*fi, pwi*fr-pwr*fi);
		}

		/* zero frequencies above maximum to Nyquist (if present) */
		for (iw=iwmax+1,iwn=nw-iw; iw<=nw/2; ++iw,--iwn)
			pw[iw] = pw[iwn] = czero;

		/* Fourier transform p(w) to p(u) */
		pfacc(-1,nufft,pu);

		/* compress p(u) to p(t;k) */
		ints8c(nu,du,fu,pu,czero,czero,nt,uoft,ptk[ik]);
	}

	/* zero wavenumber between maximum and Nyquist */
	for (; ik<nk; ++ik)
		for (it=0; it<nt; ++it)
			ptk[ik][it] = czero;

	/* Fourier transform p(t,k) to p(t,x) */
	pfa2cr(1,2,nt,nxfft,ptk[0],ptx);

	/* free workspace */
	free1float(tofu);
	free1float(uoft);
	free1complex(pu);
	free1(ptk);
	free1float(vrms);
}

static void maketu (float offset, int itmin, float fmax,
	int nt, float dt, float ft, float *vrms, float **uoftp,
	int *nup, float *dup, float *fup, float **tofup, float *tconp)
/*****************************************************************************
make stretch and compress functions t(u) and u(t)
******************************************************************************
Input:
offset		source receiver offset
itmin		index of minimum first non-zero sample for this offset
fmax		maximum frequency
nt		number of time samples
dt		time sampling interval
ft		first time
vrms		array[nt] of rms velocities

Output:
uoftp		array[nt] of u(t)
nup		number of u (stretched t) samples
dup		u sampling interval
fup		first u
tofup		array[nu] of t(u)
tconp		time constant relating t(u) and u(t)
******************************************************************************
Author:	 Dave Hale, Colorado School of Mines, 08/08/91
*****************************************************************************/
{
	int it,numax,nu;
	float tmin,dumin,et,eu,t1,t2,
		v2,v22,v44,gamma,
		v2m,v22m,v44m,gammam,t,dv2,vi2,vi4,v24,
		*uoft,du,fu,*tofu,tcon;
	
	/* determine maximum number of u */
	numax = 500+log((float)nt)*(float)(nt-1);
		
	/* allocate space for u(t) */
	uoft = ealloc1float(nt);
	
	/* determine t1 and t2, rounded to nearest sampled times */
	tmin = ft+itmin*dt;
	et = ft+(nt-1)*dt;
	t1 = MIN(et,MAX(ft+dt,tmin));
	if (offset!=0.0)
		t2 = MAX(t1,1.0/(1.0/et+0.2*dt*vrms[0]*vrms[0]/
			(offset*offset)));
	else
		t2 = t1;
	t1 = ft+NINT(t1/dt)*dt;
	t2 = ft+NINT(t2/dt)*dt;
	
	/* compute u(t) */
	v2 = vrms[0];
	v22 = v2*v2;
	v44 = v22*v22;
	gamma = 1.0;
	for (it=0,t=ft; it<nt; ++it,t+=dt) {
		v2m = v2;
		v22m = v22;
		v44m = v44;
		gammam = gamma;
		if (t>0.0) {
			v2 = vrms[it];
			v22 = v2*v2;
			vi2 = (t*v22-(t-dt)*v22m)/dt;
			vi4 = vi2*vi2;
			v44 = (dt*vi4+(t-dt)*v44m)/t;
		} else {
			v2 = v2m;
			v22 = v22m;
			v44 = v44m;
		}
		dv2 = (v2-v2m)/dt;
		v24 = v22*v22;
		gamma = 1.5*v44/v24-t*dv2/v2-0.5;
		if (t<=t1) {
			uoft[it] = t-t1;
		} else if (t>t1 && t<=t2) {
			du = t1*(gamma*log(t/(t-0.5*dt)) -
				gammam*log((t-dt)/(t-0.5*dt)));
			dumin = 0.1*dt*t1/t;
			uoft[it] = uoft[it-1]+MAX(dumin,du);
		} else if (t>t2) {
			uoft[it] = 2.0*uoft[it-1]-uoft[it-2];
		}
	}
	
	/* determine minimum u(t)-u(t-dt) */
	dumin = uoft[1]-uoft[0];
	for (it=1; it<nt; ++it)
		dumin = MIN(dumin,uoft[it]-uoft[it-1]);
	
	/* determine u sampling for t(u) to avoid aliasing */
	fu = 0.0;
	eu = uoft[nt-1];
	du = dumin/MIN(1.0,2.0*fmax*dt);
	nu = 1+NINT((eu-fu)/du);
	if (nu>numax) {
		nu = numax;
		du = (eu-fu)/(nu-1);
	}
	
	/* allocate space for t(u) */
	tofu = ealloc1float(nu);
	
	/* compute t(u) by inverse linear interpolation of u(t) */
	yxtoxy(nt,dt,ft,uoft,nu,du,fu,ft,et,tofu);
	
	/* set time constant */
	tcon = t1;
	
	/* set returned values */
	*uoftp = uoft;
	*nup = nu;
	*dup = du;
	*fup = fu;
	*tofup = tofu;
	*tconp = tcon;
}
