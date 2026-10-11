/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDMOFKCW: $Revision: 1.8 $ ; $Date: 2011/11/16 17:51:02 $	*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUDMOFKCW - converted-wave DMO via F-K domain (log-stretch) method for",
" 		common-offset gathers					",
" 									",
" sudmofkcw <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]	",
" 									",
" Required Parameters:							",
" cdpmin		  minimum cdp (integer number) for which to apply DMO",
" cdpmax		  maximum cdp (integer number) for which to apply DMO",
" dxcdp		   distance between adjacent cdp bins (m)		",
" noffmix		 number of offsets to mix (see notes)		",
" 									",
" Optional Parameters:							",
" tdmo=0.0		times corresponding to rms velocities in vdmo (s)",
" vdmo=1500.0		rms velocities corresponding to times in tdmo (m/s)",
" gamma=0.5		 velocity ratio, upgoing/downgoing		",
" ntable=1000		 number of tabulated z/h and b/h (see notes)	",
" sdmo=1.0		DMO stretch factor; try 0.6 for typical v(z)	",
" flip=0		 =1 for negative shifts and exchanging s1 and s2",
" 			 (see notes below)				",
" fmax=0.5/dt		maximum frequency in input traces (Hz)		",
" verbose=0		=1 for diagnostic print				",
" 									",
" Notes:								",
" Input traces should be sorted into common-offset gathers.  One common-",
" offset gather ends and another begins when the offset field of the trace",
" headers changes.							",
" 									",
" The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
" the cdp location expressed in units of meters or feet.		",
" 									",
" The number of offsets to mix (noffmix) should typically equal the ratio of",
" the shotpoint spacing to the cdp spacing.  This choice ensures that every",
" cdp will be represented in each offset mix.  Traces in each mix will	",
" contribute through DMO to other traces in adjacent cdps within that mix.",
" 									",
" The tdmo and vdmo arrays specify a velocity function of time that is	",
" used to implement a first-order correction for depth-variable velocity.",
" The times in tdmo must be monotonically increasing. The velocity function",
" is assumed to have been gotten by traditional NMO. 			",
" 									",
" For each offset, the minimum time at which a non-zero sample exists is",
" used to determine a mute time.  Output samples for times earlier than this",
" mute time will be zeroed.  Computation time may be significantly reduced",
" if the input traces are zeroed (muted) for early times at large offsets.",
" 									",
" z/h is horizontal-reflector depth normalized to half source-reciver offset",
" h.  Normalized shift of conversion point is b/h.  The code now does not",
" support signed offsets, therefore it is recommended that only end-on data,",
" not split-spread, be used as input (of course after being sorted into	",
" common-offset gathers).  z/h vs b/h depends on gamma (see Alfaraj's Ph.D.",
" thesis, 1993).							",
" 									",
" Flip factor = 1 implies positive shift of traces (in the increasing CDP",
" bin number direction).  When processing split-spread data, for example,",
" if one side of the spread is processed with flip=0, then the other side",
" of the spread should be processed with flip=1.  The flip factor also	",
" determines the actions of the factors s1 and s2, i.e., stretching or	",
" squeezing.								",
" 									",
" Trace header fields accessed:  nt, dt, delrt, offset, cdp.		",
"									",
NULL};
#endif


/* Credits:
 *	CWP: Mohamed Alfaraj
 *		Dave Hale
 *
 * Technical Reference:
 *	Transformation to zero offset for mode-converted waves
 *	Mohammed Alfaraj, Ph.D. thesis, 1993, Colorado School of Mines
 *
 *	Dip-Moveout Processing - SEG Course Notes
 *	Dave Hale, 1988
 */
/**************** end self doc *******************************************/

/* Library version of SUDMOFKCW
 *
 * The main program (the parameters, the traces and their headers, the offsets and the mixes of them) is left to the caller;
 * the rms velocity function of time is made by su_dmofk_vrms of SUDMOFK, and the number of cdps of the padded gather by
 * su_dmofk_nxfft.
 *
 *	su_dmofkcw_table:	the table of z/h and b/h for the velocity ratio gamma (not made when gamma is 1)
 *	su_dmofkcw_stretch:	the stretch and squeeze factors s1 and s2 (exchanged by the caller for flip)
 *	su_dmofkcw_offset:	DMO of a common-offset gather p[nxfft + 2][nt] of converted waves, in place
 *
 * Differences from the program:
 *  - The gather is a flat array, p[(nxfft + 2) * nt].
 *  - The arrays of the shifts in time of getbbohitn (b/h as a function of time, and the sample indices, which the program
 *    does not use) are freed.
 */

static void maketu (float offset, int itmin, float fmax,
	int nt, float dt, float ft, float *vrms, float **uoftp,
	int *nup, float *dup, float *fup, float **tofup, float *tconp);
static void getbbohitn (float x, int itmin, int nt, float dt, float *v,
	int ntable, const float *boh, const float *zoh, float gamma,
	float **bbohp, int **t0p);

int su_dmofkcw_table(int ntable, float gamma, float *zoh, float *boh)
/*****************************************************************************
make table of z/h vs b/h, z is depth, b is lateral shift, and h is half offset
******************************************************************************
Input:
ntable		number of elements in table
gamma		upgoing-to-downgoing velocity ratio

Output:
zoh		array[ntable] of z/h
boh		array[ntable] of b/h
******************************************************************************
Author:  Mohammed Alfaraj, Colorado School of Mines, 01/08/92
*****************************************************************************/
{
	int	i;
	float 	f, alpha, beta, bstart, b, z, db;

	if (ntable<2 || gamma<=0.0) return -3;
	f=(1.0-gamma)/(1.0+gamma);
	alpha=1.0+1.0/(gamma*gamma);
	beta=1.0-1.0/(gamma*gamma);
	if(gamma<1.0)
		db = -((bstart=1.0)-f)/(ntable+2);
	else
		db = (f -(bstart=(-1.0)))/(ntable+2);

	/* generate table */
	for(i=0,b=bstart; i<ntable; i++,b+=(db)){
		z=(1+b)*(1-b*b)*(2.0-alpha-beta*b)/(alpha-2.0+
		b*(beta-alpha-2.0)-beta*b*b);
		if(z<0) return -4;	/* negative sqrt when generating table */
		zoh[i] = sqrt(z);
		boh[i] = b;
	}
	return 0;
}

void su_dmofkcw_stretch(float sdmo, float gamma, float *s1, float *s2)
/*****************************************************************************
calculate stretch/squeeze factors s1 and s2
******************************************************************************
Author:  Mohammed Alfaraj, Colorado School of Mines, 01/08/92
*****************************************************************************/
{
	/* solve for DMO stretch factors s1 and s2 for nominal error
	   in log stretch; for gamma=1, s1=s2=0.62 */
	if(gamma <= 1.0) {
		*s1 = 0.62+0.633*(pow(gamma,-3.0)-1.0);
		*s2 = 0.62+0.615*(pow(2.0-gamma,-7.0)-1.0);
	} else {
		*s1 = 0.62+0.615*(pow(2.0-1.0/gamma,-7.0)-1.0);
		*s2 = 0.62+0.633*(pow(1.0/gamma,-3.0)-1.0);
	}

	/* modify stretch/squeeze factors for v(z) */
	*s1 *= sdmo;
	*s2 *= sdmo;
}

void su_dmofkcw_offset(float offset, float fmax, int nx, float dx, int nt, float dt, float ft, const float *vrms_in,
	float *ptx, float gamma, const float *boh, const float *zoh, int ntable, float s1, float s2, float sign)
/*****************************************************************************
perform dmo in f-k domain for one offset (dmooff of the program)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 08/08/91
*****************************************************************************/
{
	int nxfft,itmin,nu,nufft,nw,nk,ix,iu,iw,ik,it,iwn,
		iwmin,iwmax,nupad,ikmax,*itn=NULL,lnt;
	float dw,dk,tcon,wwscl,wwscl2,scale,scales,kmax,lt,
		amp,phase,fr,fi,pwr,pwi,cphase,sphase,os1,os2,
		wmin,wmax,fftscl,du,fu,w,k,*uoft,*tofu,g1,h,vmin,hk,t,*bboh=NULL;
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

	/* constants depending on gamma, offset, vrms[0], and nt */
	g1 = 2.0*sqrt(gamma)/(1.0+gamma);
	h = offset/2.0;
	lnt = nt-1;
	lt = (nt-1)*dt;
	vmin = 0.5*MIN((1.0+1.0/gamma)*vrms[0],(1.0+gamma)*vrms[0]);

	/* if gamma != 1, get bboh[] and itn[] for this offset */
	if(gamma!=1.0)
		getbbohitn (offset,itmin,nt,dt,vrms,ntable,boh,zoh,
		gamma,&bboh,&itn);

	/* inverse of dmo stretch/squeeze factors */
	os1 = 1.0/s1;
	os2 = 1.0/s2;

	/* maximum DMO shift (in samples) for any wavenumber k */
	nupad = 1.5*((s1+s2)/2.0)*tcon/du;

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

		/* apply time vriant linear phase shift */
		hk=sign*h*k;
		for (it=lnt,t=lt; it>=itmin; --it,t-=dt) {

			/* calculate phase-shift=boh*h*k, h=offset/2 */
			if(gamma != 1.0)
				phase = bboh[it]*hk;
			else
				phase = 0.0;

			/* phase shift p(t,k) */
			cphase=cos(phase);
			sphase=sin(phase);
			fr = ptk[ik][it].real();
			fi = ptk[ik][it].imag();
			ptk[ik][it] = complex(fr*cphase + fi*sphase, -fr*sphase + fi*cphase);
		}

		/* stretch p(t;k) to p(u) */
		ints8c(nt,dt,ft,ptk[ik],czero,czero,nu,tofu,pu);

		/* pad with zeros and Fourier transform p(u) to p(w) */
		for (iu=nu; iu<nufft; ++iu)
			pu[iu] = czero;
		pfacc(1,nufft,pu);

		/* minimum and maximum frequencies to process */
		wmin = ABS(0.5*vmin*k);
		wmax = ABS(PI/du);
		iwmin = MAX(1,MIN((nw-1)/2,NINT(wmin/dw)));
		iwmax = MAX(0,MIN((nw-1)/2,NINT(wmax/dw)));

		/* constant independent of w */
		wwscl = os1*pow(g1*hk/tcon,2.0);
		wwscl2 = wwscl*os2/os1;

		/* zero dc (should be zero anyway) */
		pw[0] = czero;

		/* zero frequencies below minimum */
		for (iw=1,iwn=nw-iw; iw<iwmin; ++iw,--iwn)
			pw[iw] = pw[iwn] = czero;

		/* do dmo between minimum and maximum frequencies */
		for (iw=iwmin,iwn=nw-iwmin,w=iwmin*dw;
			iw<=iwmax; ++iw,--iwn,w+=dw) {

			/* positive w */
			scales = 1.0+wwscl/(w*w);
			scale = sqrt(scales);
			phase = s1*w*tcon*(scale-1.0);
			amp = fftscl*(1.0-s1+s1/scale);
			fr = amp*cos(phase);
			fi = amp*sin(phase);
			pwr = pw[iw].real();
			pwi = pw[iw].imag();
			pw[iw] = complex(pwr*fr-pwi*fi, pwr*fi+pwi*fr);

			/* negative w */
			scales = 1.0+wwscl2/(w*w);
			scale = sqrt(scales);
			phase = s2*w*tcon*(scale-1.0);
			amp = fftscl*(1.0-s2+s2/scale);
			fr = amp*cos(phase);
			fi = amp*sin(phase);
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
	if (bboh) free1float(bboh);
	if (itn) free1int(itn);
}

static void maketu (float offset, int itmin, float fmax,
	int nt, float dt, float ft, float *vrms, float **uoftp,
	int *nup, float *dup, float *fup, float **tofup, float *tconp)
/*****************************************************************************
make stretch and compress functions t(u) and u(t)
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 08/08/91
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
		if (t<=t1 || it<1) {
			uoft[it] = t-t1;
		} else if (t>t1 && t<=t2) {
			du = t1*(gamma*log(t/(t-0.5*dt)) -
				gammam*log((t-dt)/(t-0.5*dt)));
			dumin = 0.1*dt*t1/t;
			uoft[it] = uoft[it-1]+MAX(dumin,du);
		} else if (it>=2) {
			uoft[it] = 2.0*uoft[it-1]-uoft[it-2];
		} else {
			uoft[it] = uoft[it-1]+dt;
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

static void getbbohitn (float x, int itmin, int nt, float dt, float *v,
	int ntable, const float *boh, const float *zoh, float gamma,
	float **bbohp, int **itnp)
/*****************************************************************************
convert b/h from a function of depth to a function of NMO time
******************************************************************************
Input:
x		offset
itmin		mininimum NMO time index to process
nt		number of time samples
dt		time sampling interval
v		array[nt] of RMS velocities
ntable		number of tabulated zoh or boh
boh		array[ntable] of b/h
zoh		array[ntable] of z/h
gamma		velocity ratio (upgoing/downgoing)

Output:
bbohp		array[nt] of b/h
itnp		array[nt] of time-sample indexes
******************************************************************************
Author:  Mohammed Alfaraj, Colorado School of Mines, 01/08/92
*****************************************************************************/
{
	int 	i, j=0, it, gotit;
	float	alpha, beta, tossq, xov, nz, t, xsq, tog, temp;
	float	*bboh;
	int	*itn;

	/* allocate space */
	bboh = ealloc1float(nt);
	itn = ealloc1int(nt);
	memset(bboh,0,sizeof(float)*nt);
	memset(itn,0,sizeof(int)*nt);

	/* constants depending on gamma */
	alpha = 1.0+1.0/(gamma*gamma);
	beta = 1.0-1.0/(gamma*gamma);
	tossq = 2.0/((1.0+1.0/gamma)*(1.0+1.0/gamma));
	tog=2.0/gamma;
	xsq = x*x;

	/* loop over it */
	for (it=itmin,t=itmin*dt; it<nt; it++,t+=dt) {
		xov = x/v[it];
		nz = t/xov;
		gotit = 0;

		/* loop over table */
		for (i=j; i<ntable; i++)
			if(zoh[i]>=nz) {
				bboh[it] = boh[i];
				temp = (t*t/(1-boh[i]*boh[i])-xov*xov
				*(tog/(alpha+ beta*boh[i])-1))*
				(alpha+beta*boh[i])*tossq;

				/* zero time if evanscent */
				if (temp<0.0)
					itn[it] = 0;
				else{
					itn[it] = NINT(sqrt(temp)/dt);
				}
				j = i;
				gotit = 1;
				break;
			}
		if (gotit) continue;

		/* else set boh to asymtotic value then break */
		for (j=it; j<nt; j++,t+=dt) {
			bboh[j] = boh[ntable-1];
			temp = (t*t/(1-bboh[j]*bboh[j])-(xsq/(v[j]*
			v[j]))*(tog/(alpha+ beta*bboh[j])-1))*
			(alpha+beta*bboh[j])*tossq;
			itn[j] = temp<0.0 ? 0 : NINT(sqrt(temp)/dt);
		}
		break;
	}

	/* set returned values */
	*bbohp = bboh;
	*itnp = itn;
}
