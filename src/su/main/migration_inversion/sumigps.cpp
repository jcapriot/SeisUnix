/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGPS: $Revision: 1.15 $ ; $Date: 2011/11/16 22:14:43 $		*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUMIGPS - MIGration by Phase Shift with turning rays			",
"									",
" sumigps <stdin >stdout [optional parms]				",
"									",
" Required Parameters:							",
" 	None								",
"									",
" Optional Parameters:							",
" dt=from header(dt) or .004	time sampling interval			",
" dx=from header(d2) or 1.0	distance between sucessive cdp's	",
" ffil=0,0,0.5/dt,0.5/dt  trapezoidal window of frequencies to migrate	",
" tmig=0.0		times corresponding to interval velocities in vmig",
" vmig=1500.0		interval velocities corresponding to times in tmig",
" vfile=		binary (non-ascii) file containing velocities v(t)",
" nxpad=0		number of cdps to pad with zeros before FFT	",
" ltaper=0		length of linear taper for left and right edges", 
" verbose=0		=1 for diagnostic print				",
"									",
"									",
" tmpdir= 	 if non-empty, use the value as a directory path	",
"		 prefix for storing temporary files; else if the	",
"	         the CWP_TMPDIR environment variable is set use		",
"	         its value for the path; else use tmpfile()		",
" 									",
" Notes:								",
" Input traces must be sorted by either increasing or decreasing cdp.	",
"									",
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
NULL};
#endif


/* Credits:
 *	CWP: Dave Hale (originally called supsmig.c)
 *
 *  Trace header fields accessed:  ns, dt, d2
 */
/**************** end self doc *******************************************/


/* Library version of SUMIGPS
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migps:	migrates the section data[nx][nt] (the traces in order of the cdps) to out[nx][nt], with the interval
 *			velocities vt[nt] of the times of the samples
 *
 * The cosine and sine tables, which the program kept in static variables, are made for each wavenumber. The phase is wrapped
 * with floor, where the program truncated it, which indexed the table before its start for the negative phases of the turned image.
 * Returns 0, or -3 for parameters that are not allowed.
 */

typedef struct TATableStruct {
	int np;
	float *p;
	float *taumax;
	float *taumin;
	float *tturn;
	float **t;
	float **a;
} TATable;

static void mig1k (TATable *table, float k, float ffil[4], int ntflag,
	int nt, float dt, complex *cp,
	int ntau, float dtau, complex *cq)
{
	int ntfft,nw,iwnyq,iwmin,iwmax,iw,itau,it,np,
		ip1,ip2,iwfil0,iwfil1,iwfil2,iwfil3,
		itaumax,itaumin,itab;
	float dw,wnyq,wmin,wmax,w,ampf,p1,p2,a1,a2,wa1,wa2,
		taumaxi,taumini,kow,phst,pwsr,pwsi,pwdr,pwdi,
		ampi,phsi,phsn,
		*p,*taumax,*taumin,*tturn,**t,**a;
	complex *pt,*pw,*qn,*qt;
	int ntab=1025;
	float fntab=(float)ntab,*ctab,*stab,opi2=1.0/(PI*2.0);
	float angle,dangle=2.0*PI/(ntab-1);

	/* cosine/sine tables */
	ctab = alloc1float(ntab);
	stab = alloc1float(ntab);
	for (itab=0,angle=0.0; itab<ntab; ++itab,angle+=dangle) {
		ctab[itab] = cos(angle);
		stab[itab] = sin(angle);
	}

	/* parts of time/amplitude table */
	np = table->np;
	p = table->p;
	taumax = table->taumax;
	taumin = table->taumin;
	tturn = table->tturn;
	t = table->t;
	a = table->a;

	/* frequency sampling */
	ntfft = npfao(2*nt,4*nt);
	nw = ntfft;
	dw = 2.0*PI/(ntfft*dt);

	/* allocate workspace */
	pt = pw = alloc1complex(ntfft);
	qn = alloc1complex(ntau);
	qt = alloc1complex(ntau);

	/* nyquist frequency */
	wnyq = PI/dt;

	/* index of smallest frequency >= nyquist */
	iwnyq = (ntfft%2)?ntfft/2:ntfft/2+1;

	/* determine frequency filter sample indices */
	iwfil0 = MAX(0,MIN(iwnyq,NINT(2.0*PI*ffil[0]/dw)));
	iwfil1 = MAX(0,MIN(iwnyq,NINT(2.0*PI*ffil[1]/dw)));
	iwfil2 = MAX(0,MIN(iwnyq,NINT(2.0*PI*ffil[2]/dw)));
	iwfil3 = MAX(0,MIN(iwnyq,NINT(2.0*PI*ffil[3]/dw)));

	/* pad p(t) with zeros */
	for (it=0; it<nt; ++it)
		pt[it] = cp[it];
	for (it=nt; it<ntfft; ++it)
		pt[it] = complex(0.0f,0.0f);

	/* Fourier transform p(t) to p(w) (include scaling) */
	pfacc(1,ntfft,pt);
	for (iw=0; iw<nw; ++iw)
		pw[iw] *= 1.0f/ntfft;

	/* initially zero normal and turned images */
	for (itau=0; itau<ntau; ++itau)
		qn[itau] = qt[itau] = complex(0.0f,0.0f);

	/* minimum and maximum frequency indices */
	wmin = ABS(k)/p[np-1];
	iwmin = MAX(MAX(1,iwfil0),(int)(wmin/dw));
	if (p[0]<=ABS(k)/wnyq)
		wmax = wnyq;
	else
		wmax = ABS(k)/p[0];
	iwmax = MIN(MIN(iwnyq-1,iwfil3),(int)(wmax/dw));

	/* loop over frequencies */
	for (iw=iwmin,w=iwmin*dw; iw<=iwmax; ++iw,w+=dw) {

		/* if slope not within range of table, continue */
		kow = ABS(k)/w;
		if (kow>p[np-1]) continue;

		/* amplitude of frequency filter */
		if (iwfil0<=iw && iw<iwfil1)
			ampf = (float)(iw-iwfil0)/(float)(iwfil1-iwfil0);
		else if (iwfil2<iw && iw<=iwfil3)
			ampf = (float)(iwfil3-iw)/(float)(iwfil3-iwfil2);
		else
			ampf = 1.0;

		/* weights for interpolation in table */
		xindex(np,p,kow,&ip1);
		ip1 = MIN(ip1,np-2);
		ip2 = ip1+1;
		p1 = p[ip1];
		p2 = p[ip2];
		a1 = (p2*p2-kow*kow)/(p2*p2-p1*p1);
		a2 = (kow*kow-p1*p1)/(p2*p2-p1*p1);
		wa1 = w*a1;
		wa2 = w*a2;

		/* maximum tau index for normal image */
		taumaxi = a1*taumax[ip1]+a2*taumax[ip2];
		itaumax = MIN(ntau-1,NINT(taumaxi/dtau));

		/* minimum tau index for turned image */
		taumini = a1*taumin[ip1]+a2*taumin[ip2];
		itaumin = MAX(0,NINT(taumini/dtau));

		/* phase at turning point */
		phst = 2.0*(wa1*tturn[ip1]+wa2*tturn[ip2]);

		/* filtered sum and differences for positive and negative w */
		pwsr = ampf*(pw[iw].real()+pw[nw-iw].real());
		pwsi = ampf*(pw[iw].imag()+pw[nw-iw].imag());
		pwdr = ampf*(pw[iw].real()-pw[nw-iw].real());
		pwdi = ampf*(pw[iw].imag()-pw[nw-iw].imag());

		/* accumulate normal image */
		if (ntflag&1) {
		for (itau=0; itau<itaumax; ++itau) {
			ampi = a1*a[ip1][itau]+a2*a[ip2][itau];
			phsi = wa1*t[ip1][itau]+wa2*t[ip2][itau];
			phsn = phsi*opi2;
			itab = fntab*(phsn-floor(phsn));
			qn[itau] += complex(ampi*(pwsr*ctab[itab]+pwdi*stab[itab]),
				ampi*(pwsi*ctab[itab]-pwdr*stab[itab]));
		}
		}

		/* accumulate turned image */
		if (ntflag&2) {
		for (itau=itaumin+1; itau<itaumax; ++itau) {
			ampi = a1*a[ip1][itau]+a2*a[ip2][itau];
			phsi = phst-(wa1*t[ip1][itau]+wa2*t[ip2][itau]);
			phsn = phsi*opi2;
			itab = fntab*(phsn-floor(phsn));
			qt[itau] += complex(ampi*(pwsr*stab[itab]-pwdi*ctab[itab]),
				ampi*(pwsi*stab[itab]+pwdr*ctab[itab]));
		}
		}
	}

	/* sum normal and turned images */
	for (itau=0; itau<ntau; ++itau)
		cq[itau] = qn[itau]+qt[itau];

	/* free workspace */
	free1complex(pt);
	free1complex(qn);
	free1complex(qt);
	free1float(ctab);
	free1float(stab);
}

static TATable *tableta (int np, int ntau, float dtau, float vtau[],
	int nt, float dt)
{
	int jp,ktau,itau,jtau,ltau,it;
	float pj,taui,taul,tauj,ti,tl,ai,al,
		vel,dvel,angle,cosa,frac,
		*p,*taumax,*taumin,*tturn,**t,**a;
	TATable *table;

	/* allocate table */
	table = (TATable*)alloc1(1,sizeof(TATable));
	table->np = np;
	table->p = p = alloc1float(np);
	table->taumax = taumax = alloc1float(np);
	table->taumin = taumin =alloc1float(np);
	table->tturn = tturn = alloc1float(np);
	table->t = t = alloc2float(ntau,np);
	table->a = a = alloc2float(ntau,np);
	for (jp=0; jp<np; ++jp) tturn[jp] = taumax[jp] = taumin[jp] = 0.0;

	/* loop over slopes p */
	for (jp=0; jp<np; ++jp) {

		/* slope p */
		p[jp] = pj = 2.0/vtau[0]*sqrt((float)(jp)/(float)(np-1));

		/* time and amplitude at tau = 0 */
		t[jp][0] = 0.0;
		a[jp][0] = 1.0;

		/* tau index */
		ktau = 0;

		/* initial ray tracing parameters */
		taui = 0.0;
		ti = 0.0;
		ai = 1.0;
		vel = 0.5*vtau[0];
		dvel = 0.5*(vtau[1]-vtau[0]);
		angle = asin(MIN(1.0,pj*vel));

		/* loop over times t */
		for (it=1; it<nt; ++it) {

			/* remember last tau, t, and a */
			taul = taui;
			tl = ti;
			al = ai;

			/* update cosine of propagation angle */
			cosa = cos(angle);

			/* update tau, t, and a */
			taui += dt*cosa;
			ti += dt*cosa*cosa;
			ai = 1.0;

			/* if ray emerges at surface, break */
			if (taui<0.0) break;

			/* update turning time and max,min tau */
			if (taui>=taul) {
				tturn[jp] = ti;
				taumax[jp] = taui;
				taumin[jp] = taui;
			} else {
				taumin[jp] = taui;
			}

			/* compute tau sample indices */
			itau = (int)(taui/dtau);
			ltau = (int)(taul/dtau);

			/* loop over tau samples crossed by ray */
			for (jtau=ltau+1; jtau<=MIN(itau,ntau-1); ++jtau) {

				/* tau of sample crossed */
				tauj = jtau*dtau;

				/* time and amp via linear interpolation */
				frac = (tauj-taul)/(taui-taul);
				ktau++;
				t[jp][ktau] = (1.0-frac)*tl+frac*ti;
				a[jp][ktau] = (1.0-frac)*al+frac*ai;
			}

			/* update angle */
			angle += pj*dvel;

			/* update velocity and first difference */
			if (itau<ntau-1) {
				frac = (taui-(itau*dtau))/dtau;
				dvel = 0.5*(vtau[itau+1]-vtau[itau]);
				vel = 0.5*((1.0-frac)*vtau[itau] +
					frac*vtau[itau+1]);
			} else {
				dvel = 0.5*(vtau[ntau-1]-vtau[ntau-2]);
				vel = 0.5*vtau[ntau-1] +
					(taui-(ntau-1)*dtau)*dvel/dtau;
			}
		}

		/* extrapolate times and amplitudes for interpolation */
		for (jtau=ktau+1; jtau<ntau; ++jtau) {
			t[jp][jtau] = t[jp][ktau];
			a[jp][jtau] = a[jp][ktau];
		}
	}
	return table;
}

static void freetable(TATable *table)
{
	free1float(table->p);
	free1float(table->taumax);
	free1float(table->taumin);
	free1float(table->tturn);
	free2float(table->t);
	free2float(table->a);
	free1(table);
}

int su_migps(int nx, int nt, float dt, float dx, const float *ffil_in, int nxpad, int ltaper, int np, int ntflag,
	const float *vt, const float *data, float *out)
{
	int nxfft,nk,ix,it,ik;
	float dk,taper,k,fftscl,ffil[4];
	float **gtx;
	complex **gtk;
	TATable *table;

	if (nx<1 || nt<2 || dt<=0.0 || dx<=0.0 || np<2 || nxpad<0 || ltaper<0 || ntflag<1 || ntflag>3) return -3;
	for (it=0; it<nt; ++it) if (vt[it]<=0.0) return -3;
	for (it=0; it<4; ++it) ffil[it] = ffil_in[it];

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

	/* build time/amplitude table */
	table = tableta(np,nt,dt,(float*)vt,nt,dt);

	/* loop over wavenumbers */
	for (ik=0,k=0.0; ik<nk; ++ik,k+=dk)
		mig1k(table,k,ffil,ntflag,nt,dt,gtk[ik],nt,dt,gtk[ik]);

	/* Fourier transform g(t,k) to g(t,x) */
	pfa2cr(1,2,nt,nxfft,gtk[0],gtx[0]);

	for (ix=0; ix<nx; ++ix)
		for (it=0; it<nt; ++it)
			out[(size_t)ix*nt+it] = gtx[ix][it];

	freetable(table);
	free1(gtx);
	free2complex(gtk);
	return 0;
}
