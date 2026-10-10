/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSTOLT: $Revision: 1.22 $ ; $Date: 2011/11/16 22:14:43 $		*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUSTOLT - Stolt migration for stacked data or common-offset gathers	",
"									",
" sustolt <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]		",
"									",
" Required Parameters:							",
" cdpmin=		  minimum cdp (integer number) in dataset	",
" cdpmax=		  maximum cdp (integer number) in dataset	",
" dxcdp=		  distance between adjacent cdp bins (m)	",
"									",
" Optional Parameters:							",
" noffmix=1		number of offsets to mix (for unstacked data only)",
" tmig=0.0		times corresponding to rms velocities in vmig (s)",
" vmig=1500.0		rms velocities corresponding to times in tmig (m/s)",
" smig=1.0		stretch factor (0.6 typical if vrms increasing)",
" vscale=1.0		scale factor to apply to velocities		",
" fmax=Nyquist		maximum frequency in input data (Hz)		",
" lstaper=0		length of side tapers (# of traces)		",
" lbtaper=0		length of bottom taper (# of samples)		",
" verbose=0		=1 for diagnostic print				",
" tmpdir=		if non-empty, use the value as a directory path	",
"			prefix for storing temporary files; else if the	",
"			the CWP_TMPDIR environment variable is set use	",
"			its value for the path; else use tmpfile()	",
" Notes:								",
" If unstacked traces are input, they should be NMO-corrected and sorted",
" into common-offset  gathers.  One common-offset gather ends and another",
" begins when the offset field of the trace headers changes. If both	",
" NMO and DMO are applied, then this is equivalent to prestack time 	",
" migration (though the velocity profile is assumed v(t), only).	",
"									",
" The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
" the cdp location expressed in units of meters or feet.		",
"									",
" The number of offsets to mix (noffmix) should be specified for	",
" unstacked data only.	noffmix should typically equal the ratio of the	",
" shotpoint spacing to the cdp spacing.	 This choice ensures that every	",
" cdp will be represented in each offset mix.  Traces in each mix will	",
" contribute through migration to other traces in adjacent cdps within	",
" that mix.								",
"									",
" The tmig and vmig arrays specify a velocity function of time that is	",
" used to implement Stolt's stretch for depth-variable velocity.  The	",
" stretch factor smig is often referred to as the \"W\" factor.		",
" The times in tmig must be monotonically increasing.			",
NULL};
#endif


/* Credits:
 *	CWP: Dave Hale c. 1990
 *
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp
 */
/**************** end self doc *******************************************/


/* Library version of SUSTOLT
 *
 * The main program (the parameters, the traces and their headers, the grouping of the traces into offsets and mixes of offsets,
 * the temporary file of headers) is left to the caller. What is here is the migration of one common-offset gather, which the
 * caller sums over the offsets of a mix:
 *
 *	su_stolt:	migrates the gather gather[nx][nt] (the traces of the cdps cdpmin..cdpmax, in order, with zeros where there
 *			is none) to out[nx][nt]
 *
 * The velocity function is given as the nmig pairs (tmig, vmig) of the program (with the last time at the end of the data and
 * the last velocity repeated), which are interpolated to v(t). The times of the first sample is 0.
 *
 * Returns 0, or -3 for parameters that are not allowed.
 */

static void makev (int nmig, float *tmig, float *vmig, float vscale,
	int nt, float dt, float ft, float **v, float *vmin, float *vmax)
{
	int it;
	float t,*vel=NULL,velmin=0.0,velmax=0.0,(*vmigd)[4];

	vmigd = (float(*)[4])ealloc1float(nmig*4);

	cmonot(nmig,tmig,vmig,vmigd);

	vel = ealloc1float(nt);
	memset((void *) vel, 0 , nt*FSIZE);

	for (it=0,t=ft; it<nt; ++it,t+=dt)
		intcub(0,nmig,tmig,vmigd,1,&t,&vel[it]);

	for (it=0; it<nt; ++it) {
		vel[it] *= vscale;
	}

	for (it=1,velmin=velmax=vel[0]; it<nt; ++it) {
		velmin = MIN(velmin,vel[it]);
		velmax = MAX(velmax,vel[it]);
	}
	free1float((float*)vmigd);
	*v = vel;
	*vmin = velmin;
	*vmax = velmax;
}

static void makeu (float vstolt, float *v, int nt, float dt, float *u)
{
	int it;
	float t,scale,sum;

	scale = 2.0/(vstolt*vstolt);
	u[0] = sum = 0.0;
	for (it=1,t=dt; it<nt; ++it,t+=dt) {
		sum += 0.5*dt*(t*v[it]*v[it]+(t-dt)*v[it-1]*v[it-1]);
		u[it] = sqrt(scale*sum);
	}
}

static void makeut (float vstolt, float fmax, float *vrms,
	int nt, float dt, float **ut, int *nu, float *du, float **tu)
{
	int it;

	/* check for constant velocity */
	for (it=1; it<nt; ++it)
		if (vrms[it]!=vrms[0]) break;

	/* if constant velocity */
	if (it==nt) {
		*ut = NULL;
		*tu = NULL;
		*nu = nt;
		*du = dt;

	/* else if velocity not constant */
	} else {
		int nuu;
		float duu=0.0;
		float delu=0.0;
		float umax=0.0;
		float *u=NULL;
		float *t=NULL;

		/* u(t) */
		u = alloc1float(nt);

		makeu(vstolt,vrms,nt,dt,u);

		/* smallest du and maximum u */
		duu = FLT_MAX;
		for (it=1; it<nt; ++it) {
			delu =	u[it]-u[it-1];
			if (delu<duu  ) duu = delu;
		}
		umax = u[nt-1];

		/* u sampling */
		duu = duu/(2.0*fmax*dt);
		nuu = 1+NINT(umax/duu);

		/* t(u) */
		t = alloc1float(nuu);
		yxtoxy(nt,dt,0.0,u,nuu,duu,0.0,0.0,(nt-1)*dt,t);

		/* set output parameters before returning */
		*ut = u;
		*tu = t;
		*nu = nuu;
		*du = duu;
	}
}

static void stolt1k (float k, float v, float s, float fmax, int nt, float dt,
	complex *p, complex *q)
{
	int nw,it,nwtau,iwtau,ntau,itau,iwtaul,iwtauh;
	float vko2s,wmax,dw,fw,dwtau,fwtau,wtau,dtau,
		wtauh,wtaul,scale,fftscl,a,b,*wwtau=NULL;
	complex czero(0.0f,0.0f),*pp=NULL,*qq=NULL;

	/* modify stolt stretch factor to simplify calculations below */
	if (s!=1.0) s = 2.0-s;

	/* (v*k/2)^2 */
	vko2s = 0.25*v*v*k*k;

	/* maximum frequency to migrate in radians per unit time */
	wmax = 2.0*PI*MIN(fmax,0.5/dt);

	/* frequency sampling - must pad to avoid interpolation error;
	 * pad by factor of 2 because time axis is not centered;
	 * pad by factor of 1/0.6 because 8-point sinc is valid
	 * only to about 0.6 Nyquist
	 */
	nw = nt*2/0.6;
	nw = npfao(nw,nw*2);
	dw = 2.0*PI/(nw*dt);
	fw = -PI/dt;

	/* migrated time */
	ntau = nt;
	dtau = dt;

	/* migrated frequency - no need to pad since no interpolation */
	nwtau = npfao(ntau,ntau*2);
	dwtau = 2.0*PI/(nwtau*dtau);
	fwtau = -PI/dtau;

	/* tweak first migrated frequency to avoid wtau==0.0 below */
	fwtau += 0.001*dwtau;

	/* high and low migrated frequencies - don't migrate evanescent */
	wtauh = sqrt(MAX(0.0,wmax*wmax-s*vko2s));
	iwtauh = MAX(0,MIN(nwtau-1,NINT((wtauh-fwtau)/dwtau)));
	iwtaul = MAX(0,MIN(nwtau-1,NINT((-wtauh-fwtau)/dwtau)));
	wtauh = fwtau+iwtauh*dwtau;
	wtaul = fwtau+iwtaul*dwtau;

	/* workspace */
	pp = alloc1complex(nw);
	qq = alloc1complex(nwtau);
	wwtau = alloc1float(nwtau);
	for (it=0; it<nwtau; ++it) qq[it] = czero;

	/* pad with zeros and Fourier transform t to w, with w centered */
	for (it=0; it<nt; it+=2)
		pp[it] = p[it];
	for (it=1; it<nt; it+=2)
		pp[it] = -p[it];
	for (it=nt; it<nw; ++it)
		pp[it] = czero;
	pfacc(1,nw,pp);

	/* zero -Nyquist frequency for symmetry */
	pp[0] = czero;

	/* frequencies at which to interpolate */
	if (s==1.0) {
		for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
			wwtau[iwtau] = -sqrt(wtau*wtau+vko2s);
		for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
			wwtau[iwtau] = sqrt(wtau*wtau+vko2s);
	} else {
		a = 1.0/s;
		b = 1.0-a;
		for (iwtau=iwtaul,wtau=wtaul; wtau<0.0; ++iwtau,wtau+=dwtau)
			wwtau[iwtau] = b*wtau-a*sqrt(wtau*wtau+s*vko2s);
		for (; iwtau<=iwtauh; ++iwtau,wtau+=dwtau)
			wwtau[iwtau] = b*wtau+a*sqrt(wtau*wtau+s*vko2s);
	}

	/* interpolate */
	ints8c(nw,dw,fw,pp,czero,czero,
		iwtauh-iwtaul+1,wwtau+iwtaul,qq+iwtaul);

	/* fft scaling and obliquity factor */
	fftscl = 1.0/nwtau;
	if (s==1.0) {
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh;
			++iwtau,wtau+=dwtau) {
			scale = fftscl*wtau/wwtau[iwtau];
			qq[iwtau] *= scale;
		}
	} else {
		a = 1.0/(s*s);
		b = 1.0-1.0/s;
		for (iwtau=iwtaul,wtau=wtaul; iwtau<=iwtauh;
			++iwtau,wtau+=dwtau) {
			scale = fftscl*(b+a*wtau/(wwtau[iwtau]-b*wtau));
			qq[iwtau] *= scale;
		}
	}

	/* zero evanescent frequencies */
	for (iwtau=0; iwtau<iwtaul; ++iwtau)
		qq[iwtau] = czero;
	for (iwtau=iwtauh+1; iwtau<nwtau; ++iwtau)
		qq[iwtau] = czero;

	/* Fourier transform wtau to tau, accounting for centered wtau */
	pfacc(-1,nwtau,qq);
	for (itau=0; itau<ntau; itau+=2)
		q[itau] = qq[itau];
	for (itau=1; itau<ntau; itau+=2)
		q[itau] = -qq[itau];

	/* free workspace */
	free1complex(pp);
	free1complex(qq);
	free1float(wwtau);
}

static void taper (int lxtaper, int lbtaper,
	int nx, int ix, int nt, float *trace)
{
	int it;
	float xtaper;

	/* if near left side */
	if (ix<lxtaper) {
		xtaper = 0.54+0.46*cos(PI*(lxtaper-ix)/lxtaper);

	/* else if near right side */
	} else if (ix>=nx-lxtaper) {
		xtaper = 0.54+0.46*cos(PI*(lxtaper+ix+1-nx)/lxtaper);

	/* else x tapering is unnecessary */
	} else {
		xtaper = 1.0;
	}

	/* if x tapering is necessary, apply it */
	if (xtaper!=1.0)
		for (it=0; it<nt; ++it)
			trace[it] *= xtaper;

	/* if requested, apply t tapering */
	for (it=MAX(0,nt-lbtaper); it<nt; ++it)
		trace[it] *= (0.54+0.46*cos(PI*(lbtaper+it+1-nt)/lbtaper));
}

int su_stolt(int nx, int nt, float dt, float dx, int ntmig, const float *tmig_in, const float *vmig_in, float vscale,
	float smig, float fmax, int lstaper, int lbtaper, const float *gather, float *out)
{
	int ix, it, iu, nu, nxpad, nxfft, nk, ik, jx;
	float vmin, vmax, vstolt, du, dk, scale, k;
	float *v = NULL, *ut = NULL, *tu = NULL, *tmig, *vmig;
	float **px, **pk_rows;
	complex **pk;

	if (nx<1 || nt<2 || dt<=0.0 || dx<=0.0 || ntmig<1 || smig<=0.0 || smig>=2.0 || fmax<=0.0 || lstaper<0 || lbtaper<0) return -3;
	tmig = ealloc1float(ntmig);
	vmig = ealloc1float(ntmig);
	for (it=0; it<ntmig; ++it) { tmig[it] = tmig_in[it]; vmig[it] = vmig_in[it]; }
	for (it=1; it<ntmig; ++it)
		if (tmig[it]<=tmig[it-1]) { free1float(tmig); free1float(vmig); return -3; }

	fmax = MIN(fmax,0.5/dt);

	/* make uniformly sampled rms velocity function of time */
	makev(ntmig,tmig,vmig,vscale,nt,dt,0.0,&v,&vmin,&vmax);
	free1float(tmig);
	free1float(vmig);
	if (vmin<=0.0) { free1float(v); return -3; }

	/* Stolt migration velocity is the minimum velocity */
	vstolt = vmin;

	/* make u(t) and t(u) for Stolt stretch */
	makeut(vstolt,fmax,v,nt,dt,&ut,&nu,&du,&tu);
	free1float(v);

	/* wavenumber (k) sampling */
	nxpad = 0.5*vmax*nt*dt/dx;
	nxfft = npfar(nx+nxpad);
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);

	/* allocate and zero common-offset gather p(t,x) */
	pk = alloc2complex(MAX(nu,nt),nk);
	pk_rows = (float**)ealloc1(nxfft,sizeof(float*));
	px = pk_rows;
	px[0] = (float*)pk[0];
	for (ix=1; ix<nxfft; ++ix)
		px[ix] = px[0]+ix*MAX(nu,nt);
	memset((void *) px[0],0,(size_t)nxfft*MAX(nu,nt)*sizeof(float));

	for (ix=0; ix<nx; ++ix)
		for (it=0; it<nt; ++it)
			px[ix][it] = gather[(size_t)ix*nt+it];

	/* apply side and bottom tapers */
	for (ix=0; ix<nx; ++ix)
		taper(lstaper,lbtaper,nx,ix,nt,px[ix]);

	/* if necessary, stretch */
	if (nu!=nt && tu!=NULL) {
		float *temp=ealloc1float(nu);
		for (ix=0; ix<nx; ++ix) {
			ints8r(nt,dt,0.0,px[ix],0.0,0.0,
				nu,tu,temp);
			for (iu=0; iu<nu; ++iu)
				px[ix][iu] = temp[iu];
		}
		free1float(temp);
	}

	/* Fourier transform p(u,x) to p(u,k) */
	pfa2rc(-1,2,MAX(nu,nt),nxfft,px[0],pk[0]);

	/* migrate each wavenumber */
	for (ik=1,k=dk; ik<nk; ++ik,k+=dk)
		stolt1k(k,vstolt,smig,fmax,
			MAX(nu,nt),du,pk[ik],pk[ik]);

	/* Fourier transform p(u,k) to p(u,x) and scale */
	scale = 1.0/nxfft;
	pfa2cr(1,2,MAX(nu,nt),nxfft,pk[0],px[0]);
	for (jx=0; jx<nx; ++jx)
		for (iu=0; iu<nu; ++iu)
			px[jx][iu] *= scale;

	/* if necessary, unstretch */
	if (nu!=nt && ut!=NULL) {
		float *temp=ealloc1float(nt);
		for (ix=0; ix<nx; ++ix) {
			ints8r(nu,du,0.0,px[ix],0.0,0.0,
				nt,ut,temp);
			for (it=0; it<nt; ++it)
				px[ix][it] = temp[it];
		}
		free1float(temp);
	}

	for (ix=0; ix<nx; ++ix)
		for (it=0; it<nt; ++it)
			out[(size_t)ix*nt+it] = px[ix][it];

	free1(pk_rows);
	free2complex(pk);
	if (ut) free1float(ut);
	if (tu) free1float(tu);
	return 0;
}
