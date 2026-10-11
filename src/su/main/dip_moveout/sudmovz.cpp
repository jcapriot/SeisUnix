/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDMOVZ: $Revision: 1.6 $ ; $Date: 2021/04/15 00:47:25 $	*/

#include "su.h"
#include "segy.h"
#define NPTAB 401
#define NITER 10
#define TOL 1.0e-4

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUDMOVZ - DMO for V(Z) media for common-offset gathers		",
"									",
" sudmovz <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [...]		",
"									",
" Required Parameters:							",
" cdpmin         minimum cdp (integer number) for which to apply DMO	",
" cdpmax         maximum cdp (integer number) for which to apply DMO	",
" dxcdp          distance between adjacent cdp bins (m)			",
" noffmix        number of offsets to mix (see notes)			",
"									",
" Optional Parameters:							",
" vfile=         binary (non-ascii) file containing interval velocities (m/s)",
" tdmo=0.0       times corresponding to interval velocities in vdmo (s)	",
" vdmo=1500.0    interval velocities corresponding to times in tdmo (m/s)",
" fmax=0.5/dt    maximum frequency in input traces (Hz)			",
" smute=1.5      stretch mute used for NMO correction			",
" speed=1.0      twist this knob for speed (and aliasing)		",
" verbose=0      =1 for diagnostic print				",
"									",
" Notes:								",
" Input traces should be sorted into common-offset gathers.  One common-",
" offset gather ends and another begins when the offset field of the trace",
" headers changes.							",
"									",
" The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
" the cdp location expressed in units of meters or feet.		",
"									",
" The number of offsets to mix (noffmix) should typically equal the ratio of",
" the shotpoint spacing to the cdp spacing.  This choice ensures that every",
" cdp will be represented in each offset mix.  Traces in each mix will	",
" contribute through DMO to other traces in adjacent cdps within that mix.",
"									",
" vfile should contain the regularly sampled interval velocities as a	",
" function of time.  If vfile is not supplied, the interval velocity	",
" function is defined by linear interpolation of the values in the tdmo	",
" and vdmo arrays.  The times in tdmo must be monotonically increasing.	",
"									",
" For each offset, the minimum time to process is determined using the	",
" smute parameter.  The DMO correction is not computed for samples that	",
" have experienced greater stretch during NMO.				",
"									",
" Trace header fields accessed:  nt, dt, delrt, offset, cdp.		",
NULL};
#endif

/**************** end self doc ********************************/

static int dmooff (float **x, float **a, float **tau, int nttab, float dt, float ft,
	int nptab, float dptab, float fp, float lp,
	int nx, float dx, float offset, int nt,
	float *vrms, float fmax, float smute, float speed, int verbose,
	float **q);
static void dmomap (float **x, float **a, float **tau, int nt, float dt, float ft,
	int nptab, float dptab, int np, float dp2, float fp, float lp,
	int itn, float tn, float h, float *vrms,
	int verbose, float **t0table);
static void dmoapply (int nx, float dx, int nt, float dt, int np, float dp2,
	float h, float **tntable, int *it0min, int *nt0,
	float tmin, float fmax, float **qx);
static void jakub1k (float k, float h, int np, float dp2, int nt, float dt, float ft,
	float **tntable, int *it0min, int *nt0, float fmax, complex *qq);
static float tmute(float *vrms, int nt, float dt, float ft, float h, float smute);
static void getX (float tsg, float h, float p0, float v, float *vars);
static void getFJ (float **x, float **a, float **tau,
	int nt, int np, float dt, float dp, float tsg, float h, float p0,
	float *vars, float *eqns, float **jacob);
static void intder (float **f, int nx, float dx, float fx,
	int ny, float dy, float fy, float x, float y,
	float *fxy, float *dfdx, float *dfdy);
static void rayvt (float p, int nt, float dt,
	float v[], float tau[], float x[], float a[]);

static int dmooff (float **x, float **a, float **tau, int nttab, float dt, float ft,
	int nptab, float dptab, float fp, float lp,
	int nx, float dx, float offset, int nt,
	float *vrms, float fmax, float smute, float speed, int verbose,
	float **q)
/***************************************************************************
Compute and apply DMO correction to one common offset section.
****************************************************************************
Input:
x	array[nptab][nttab] of horizontal position of ray 
a	array[nptab][nttab] of propagation angle along ray
tau	array[nptab][nttab] of time depth of ray
nttab	number of times in ray tables
dt	time sampling interval
ft	first time (must be zero for now)
nptab	number of ray parameters in tables
dptab	ray parameter sampling interval in tables
fp	first ray parameter
lp	last ray parameter
nx	number of traces in common offset section
dx	trace sampling interval
offset	source-receiver full offset
nt	number of time samples
vrms	array[nttab] of rms velocity as function of two-way vertical time
fmax	maximum frequency to process
smute	NMO stretch mute.  Samples experiencing greater NMO stretch are muted
speed	knob to increase speed and aliasing by coarsening the slope sampling
verbose	flag for diagnostic info
q	array[nx][nt] of input common offset section

Output:
q	array[nx][nt] of output common offset section (DMO corrected)
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 09/28/91
****************************************************************************/
{
	int it,ntn,itn,itnmin,itnmax,np,npdmo,ip,it0max,*it0min,*nt0;
	float h,lt,tn,ftn,dtn,dp2,tmin,**t0table,**tntable;

	/* compute minimum time based on stretch mute */
	tmin = tmute(vrms,nt,dt,ft,0.5*offset,smute);
	
	/* compute p sampling required to avoid aliasing */
	lt = (nt-1)*dt;
	if (tmin>=lt)
		return -4;  /* no data survive the mute */

	dtn = 10.0*dt;
	ftn = tmin;
	ntn = NINT((lt-ftn)/dtn)+1;
	dp2 = tmin/(fmax*0.25*offset*offset);
	dp2 *= speed;
	np = NINT(lp*lp/dp2)+1;

	/* allocate space for dmo mapping tables */
	it0min = ealloc1int(np);
	nt0 = ealloc1int(np);
	t0table = ealloc2float(ntn,np);
	tntable = ealloc2float(nt,np);

	/* zero the tables */
	for (ip=0; ip<np; ++ip)
		for (itn=0; itn<ntn; ++itn)
			t0table[ip][itn] = 0.0;
	for (ip=0; ip<np; ++ip)
		for (it=0; it<nt; ++it)
			tntable[ip][it] = 0.0;
	
	/* compute (unsigned) half-offset */
	h = ABS(offset/2.0);

	/* loop over nmo times tn */
	for (itn=0,tn=ftn; itn<ntn; ++itn,tn+=dtn) {

		/* compute dmo mapping t0(p) for this tn */
		dmomap(x,a,tau,nttab,dt,ft,nptab,dptab,np,dp2,fp,lp,
			itn,tn,h,vrms,verbose,t0table);
	}

	/* inverse interpolation to get tn(t0,p) */
	for (ip=0; ip<np; ++ip) {
		for (itn=0; itn<ntn; ++itn)
			if (t0table[ip][itn]>0.0) break;
		itnmin = itn;
		if (itnmin>=ntn) break;  /* no sample of this row is set (the program read after its end) */
		
		it0min[ip] = NINT(t0table[ip][itnmin]/dt);
		if (it0min[ip]>=nt) it0min[ip] = nt-1;

		for (itn=ntn-1; itn>=0; --itn)
			if (t0table[ip][itn]>0.0) break;
		itnmax = itn;
		if (itnmax <= itnmin) break;

		it0max = NINT(t0table[ip][itnmax]/dt);
		if (it0max>=nt) it0max = nt-1;

		nt0[ip] = it0max-it0min[ip]+1;

		/* check that table is monotonically increasing */
		for (itn=itnmin+1; itn<itnmax; ++itn)
			if (t0table[ip][itn]<=t0table[ip][itn-1] && verbose)
				fprintf(stderr,
					"t0table not monotonically "
					"increasing---extrapolating\n");

		yxtoxy(itnmax-itnmin+1,dtn,ftn+itnmin*dtn,t0table[ip]+itnmin,
			nt0[ip],dt,it0min[ip]*dt,
			ftn+itnmin*dtn,ftn+itnmax*dtn,tntable[ip]);

		/* extend tntable by linear extrapolation */
		for (it=nt0[ip]; it<nt; ++it) {
			tntable[ip][it] = 2.0*tntable[ip][it-1]
				-tntable[ip][it-2];
		}
		nt0[ip] = nt-it0min[ip];
	}

	/* reset the number of ray parameters for dmo */
	npdmo = MIN(ip,np);

	/* free space */
	free2float(t0table);

	/* dmo by dip decomposition */
	dmoapply(nx,dx,nt,dt,npdmo,dp2,h,tntable,it0min,nt0,tmin,fmax,q);

	/* free space */
	free1int(it0min);
	free1int(nt0);
	free2float(tntable);
	return 0;
}

static void dmomap (float **x, float **a, float **tau, int nt, float dt, float ft,
	int nptab, float dptab, int np, float dp2, float fp, float lp,
	int itn, float tn, float h, float *vrms,
	int verbose, float **t0table)
/***************************************************************************
Compute DMO correction for one NMO time, tn, using Newton's method to solve
system of equations.
****************************************************************************
Input:
x	array[nptab][nt] of horizontal position of ray 
a	array[nptab][nt] of propagation angle along ray
tau	array[nptab][nt] of time depth of ray
nt	number of times in ray tables
dt	time sampling interval
ft	first time (must be zero for now)
nptab	number of ray parameters in ray tables
dptab	ray parameter sampling interval in ray tables
np	number of ray parameters for which to find DMO mapping
dp2	sampling interval in (ray parameter)^2 in DMO mapping table
fp	first ray parameter
lp	last ray parameter
itn	index of NMO time for which to compute DMO mapping
tn	NMO time for which to compute DMO mapping
h	source-receiver half offset
vrms	array[nt] of rms velocity as function of two-way vertical time
verbose	flag for diagnostic info (does nothing here)

Output:
t0table	array[np][ntn] containing DMO mapping t0(tn,p0)
	One call to this function fills in the itn'th row of this table.
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 09/28/91
****************************************************************************/
{
	int ip,i,neqns=5,niter=NITER,iter,*ipvt;
	float p2,vtn,tsg,p0,tautp,dtaudt,dtaudp,error,rcond=1.0,
		lt,tol=TOL,*vars,*eqns,*work,t0,x0,**jacob;

	/* compute last time in ray tables */
	lt = ft+(nt-1)*dt + 0*verbose;

	/* allocate space for system of equations */
	vars = ealloc1float(neqns);		/* variables X */
	eqns = ealloc1float(neqns);		/* F(X) */
	jacob = ealloc2float(neqns,neqns);	/* Jacobian dF_i/dx_j */
	work = ealloc1float(neqns);		/* work array for sgeco() */
	ipvt = ealloc1int(neqns);	/* pivot index array for sgeco() */

	/* calculate tsg from tn via Dix approximation */
	vtn = vrms[(int)(tn/dt)];
	tsg = sqrt(tn*tn+4*h*h/(vtn*vtn));

	/* get initial values for parameters */
	getX(tsg,h,fp,vtn,vars);

	for (ip=0,p2=fp; ip<np; ++ip,p2+=dp2) {

		/* compute p0 */
		p0 = sqrt(p2);

		/* use Newton's method to improve the guess */
		for (iter=0; iter<niter; ++iter) {
	
			/* calculate F(X) and J(X) */
			getFJ(x,a,tau,nt,nptab,dt,dptab,tsg,h,p0,
				vars,eqns,jacob);

			/* compute total error */
			for (i=0,error=0.0; i<neqns; ++i)
				error += ABS(eqns[i]);

			/* if iteration has converged, break */
			if (error<=tol) break;

			/* solve system */
			sgeco(jacob,neqns,ipvt,&rcond,work);
			sgesl(jacob,neqns,ipvt,eqns,0);

			/* update X */
			if (rcond!=0.0) {
				for (i=0; i<neqns; ++i)
					vars[i] -= eqns[i];
			} else {
				/* system is ill-conditioned, break */
				break;
			}

			/* clip ps and pg values to pmax if necessary */
			if (ABS(vars[1])>lp) {
				vars[1] = SGN(vars[1])*lp;
			}
			if (ABS(vars[2])>lp) {
				vars[2] = SGN(vars[2])*lp;
			}

			/* clip tg to 0<tg<tsg */
			if (vars[3]<0.0) {
				vars[3] = 0.0;
			}
			if (vars[3]>tsg) {
				vars[3] = tsg;
			}

			/* clip t0 values to 0<t0<lt */
			if (vars[4]<0.0) {
				vars[4] = 0.0;
			}
			if (vars[4]>lt) {
				vars[4] = lt;
			}

		}

		/* vars[] now contains the solution to the system */

		/* compute depth of reflection point */
		intder(tau,nt,dt,ft,nptab,dptab,fp,vars[4],ABS(p0),
			&tautp,&dtaudt,&dtaudp);

		/* if reflector is at surface, break */
		if (tautp<dt) {
			break;
		}

		/* if system is ill-conditioned, break */
		if (rcond==0.0) {
			break;
		}

		/* now have t0(x0) for this tsg, p0, h -- store it */
		x0 = vars[0];
		t0 = vars[4];
		t0table[ip][itn] = t0+p0*x0;

	}

	/* free space */
	free1float(vars);
	free1float(eqns);
	free2float(jacob);
	free1float(work);
	free1int(ipvt);
}

static void dmoapply (int nx, float dx, int nt, float dt, int np, float dp2,
	float h, float **tntable, int *it0min, int *nt0,
	float tmin, float fmax, float **qx)
/*****************************************************************************
Apply DMO correction to common offset section using previously computed
mapping tables.
******************************************************************************
Input:
nx		number of traces in common offset section
dx		trace sampling interval
nt		number of time samples per trace
dt		time sampling interval
np		number of ray parameters in DMO mapping tables
dp2		ray parameter (squared) sampling interval
h		source-receiver half-offset
tntable		array[np][nt0] of DMO mapping table tn(t0,p0)
it0min		array[np] of index of first entry in DMO mapping table 
nt0		array[np] of number of t0's in DMO mapping table
fmax		maximum frequency of data (frequencies above fmax are muted)
qx		array[nx][nt] of input NMO corrected data

Output:
qx		array[nx][nt] of output DMO corrected data
******************************************************************************
Author: Craig Artley, Colorado School of Mines, 08/29/91
******************************************************************************/
{
	int ix,nxfft,it,itmin,ik,nk;
	float k,dk,fk,ft,scale;
	complex **qk;

	/* determine wavenumber sampling (for real to complex FFT) */
	nxfft = npfar(nx);
	nk = nxfft/2+1;
	dk = 2.0*PI/(nxfft*dx);
	fk = 0.0;

	ft = 0.0;

	/* pointers to complex q(t,k) */
	qk = (complex**)ealloc1(nk,sizeof(complex*));
	for (ik=0; ik<nk; ++ik)
		qk[ik] = (complex*)qx[0]+ik*nt;
	
	/* Fourier transform x to k */
	pfa2rc(-1,2,nt,nxfft,qx[0],qk[0]);

	/* loop over wavenumber */
	for (ik=0,k=fk; ik<nk; ++ik,k+=dk) {

		/* apply DMO correction */
		jakub1k(k,h,np,dp2,nt,dt,ft,tntable,
			it0min,nt0,fmax,qk[ik]);
	}

	/* Fourier transform k to x (including FFT scaling and mute) */
	pfa2cr(1,2,nt,nxfft,qk[0],qx[0]);
	scale = 1.0/nxfft;
	itmin = MIN(nt,NINT(tmin/dt));
	for (ix=0; ix<nx; ++ix) {
		for (it=0; it<itmin; ++it)
			qx[ix][it] = 0.0;
		for (it=itmin; it<nt; ++it)
			qx[ix][it] *= scale;
	}

	/* free space */
	free1(qk);
}

static void jakub1k (float k, float h, int np, float dp2, int nt, float dt, float ft,
	float **tntable, int *it0min, int *nt0, float fmax, complex *qq)
/*****************************************************************************
Jakubowicz's dip decomposition DMO for one wavenumber
******************************************************************************
Input:
k		wavenumber
h		half-offset (ignored)
np		number of ray paramaters
dp2		sampling interval in (ray parameter)^2
nt		number of time samples
dt		time sampling interval
ft		first time sample
qq		complex array[nt] of input NMO corrected data 

Output:
qq		complex array[nt] of output DMO corrected data
******************************************************************************
Author: Craig Artley, Colorado School of Mines, 08/29/91
******************************************************************************/
{
	int ntfft,nw,it,iw,iwlm,iwhm,iwlp,iwhp,ip;
	float p,ph,pl,ph2,dw,fw,wl,wh,scale;
	complex czero(0.0f,0.0f),*qn,*q0;

	/* determine frequency sampling */
	ntfft = npfa(nt);
	nw = ntfft;
	dw = 2.0*PI/(ntfft*dt);
	fw = -PI/dt + 0.0*h;
	
	/* allocate workspace */
	qn = ealloc1complex(ntfft);
	q0 = ealloc1complex(ntfft);

	/* zero w0 domain accumulator */
	for (iw=0; iw<nw; ++iw) {
		qn[iw] = czero;
		q0[iw] = czero;
	}

	/* loop over zero-offset ray parameter */
	for (ip=0,p=pl=0.0,ph2=dp2; ip<np; ++ip,ph2+=dp2) {

		/* update ray parameter */
		ph = sqrt(ph2);

		/* determine limits of w0 for this wavenumber and slope */
		wl = k/ph;
		if (pl==0.0) {
			wh = 2.0*PI*fmax;
		} else {
			wh = MIN(k/pl,2.0*PI*fmax);
		}

		/* compute limits of sum for -p */
		iwlm = (int) ((-wl-fw)/dw);
		iwhm = (int) ((-wh-fw)/dw);
		if (iwhm<0) iwhm = 0;
		++iwhm;

		/* compute limits of sum for +p */
		iwlp = (int) ((wl-fw)/dw);
		iwhp = (int) ((wh-fw)/dw);
		if (iwhp>=nw) iwhp = nw-1;
		++iwlp;

		/* if this range of frequencies contributes to the sum */
		if ((iwhm<=iwlm) || (iwlp<=iwhp)) {
		
			/* interpolate, padding with zeros */
			for (it=0; it<it0min[ip]; ++it)
				qn[it] = czero;
			ints8c(nt,dt,ft,qq,czero,czero,nt0[ip],tntable[ip],
				qn+it0min[ip]);
			for (it=it0min[ip]+nt0[ip]; it<ntfft; ++it)
				qn[it] = czero;

			/* Fourier transform t to w, with w centered */
			for (it=1; it<ntfft; it+=2) {
				qn[it] = -qn[it];
			}
			pfacc(1,ntfft,qn);

			/* accumulate into q0 for slopes near -p */
			for (iw=iwhm; iw<=iwlm; ++iw) {
				q0[iw] += qn[iw];
			}

			/* accumulate into q0 for slopes near +p */
			for (iw=iwlp; iw<=iwhp; ++iw) {
				q0[iw] += qn[iw];
			}
		}

		/* roll ray parameters */
		pl = p;
		p = ph;
	}

	/* Fourier transform w to t, with w centered */
	/* (including FFT scaling and division by 2 for plane filling) */
	pfacc(-1,ntfft,q0);
	scale = 0.5/ntfft;
	for (it=0; it<nt; it+=2) {
		qq[it] = q0[it]*scale;
	}
	for (it=1; it<nt; it+=2) {
		qq[it] = -q0[it]*scale;
	}

	/* free workspace */
	free1complex(qn);	
	free1complex(q0);	
}
	
static float tmute (float *vrms, int nt, float dt, float ft, float h, float smute)
/***************************************************************************
Returns MUTE Time for given rms velocity function and half-offset.
****************************************************************************
Input:
vrms	array[nt] of rms velocity as a function of NMO time
nt	number of times
dt	time sampling interval
ft	first time
h	source-receiver half-offset
smute	maximum allowable NMO stretch t/tn [(100+stretch_percentage)/100]

Return:
	gretest time that is subject to given NMO stretch mute
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 09/28/91
****************************************************************************/
{
	int it;
	float tn,*t;

	/* allocate space */
	t = ealloc1float(nt);

	for (it=0,tn=ft; it<nt; ++it,tn+=dt)
		t[it] = sqrt(tn*tn+4.0*h*h/(vrms[it]*vrms[it]));

	/* find greatest nmo time that is subject to stretch mute */
	for (it=nt-1; it>0; --it)
		if (dt/(t[it]-t[it-1])>smute) break;

	/* free space */
	free1float(t);

	return ft+it*dt;
}

static void getX (float tsg, float h, float p0, float v, float *vars)
/***************************************************************************
GET trial solution X to system of equations assuming constant velocity.
****************************************************************************
Input:
tsg	recording time (from source to receiver)
h	source-receiver half-offset
p0	ray parameter of zero-offset ray
v	constant velocity

Output:
vars	array[5] of variables that satisfies the system for const. v
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 09/21/91
****************************************************************************/
{
	float x,z,a2,b2,h2,x0,ps,pg,tg,t0,sine,tangent;

	/* compute constants that determine prestack migration ellipse */
	h2 = h*h;
	a2 = v*v*tsg*tsg/4.0;
	b2 = MAX(a2-h2,0.0);

	/* compute sine and tangent of propagation angle */
	sine = 0.5*p0*v;
	tangent = sine/sqrt(MAX(1.0-sine*sine,0.0));
	
	/* find reflection point (x,z) on prestack migration ellipse */
	x = tangent*sqrt(a2/(b2+tangent*tangent));
	z = sqrt(b2*MAX(1.0-x*x/a2,0.0));
	
	/* find zero-offset location */
	x0 = h2/a2*x;

	/* find times along receiver and zero-offset rays */
	tg = sqrt((h-x)*(h-x)+z*z)/v;
	t0 = sqrt((x0-x)*(x0-x)+z*z)*2.0/v;

	/* find ray parameters for source and receiver rays */
	ps = (x+h)/(v*v*(tsg-tg));
	pg = (x-h)/(v*v*tg);

	/* load trial solution into the variables array */
	vars[0] = x0;
	vars[1] = ps;
	vars[2] = pg;	
	vars[3] = tg;
	vars[4] = t0;
}

static void getFJ (float **x, float **a, float **tau, int nt, int np,
	float dt, float dp, float tsg, float h, float p0,
	float *vars, float *eqns, float **jacob)
/***************************************************************************
Compute equation array and Jacobian matrix for system of equations for
a given set of parameters.
****************************************************************************
Input:
x	array[np][nt] of horizontal location of tabled rays
a	array[np][nt] of propagation angle of tabled rays
tau	array[np][nt] of travel-time depth of tabled rays
nt	number of times in ray tables
np	number of ray parameters in ray tables
dt	time sampling interval in tables
dp	ray parameter sampling interval in tables
tsg	recording time (from source to receiver)
h	source-receiver half offset
p0	zero-offset ray parameter
vars	array[5] of variables of system of equations

Output:
eqns	array[5] of right-hand-sides of equations
jacob	pointer to array[5][5] of Jacobian partial derivatives
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 09/21/91
****************************************************************************/
{
	int i,j,neqns=5;
	float x0,ps,pg,ts,tg,t0,ft=0.0,fp=0.0,ftp,dfdt,dfdp,sign;

	/* extract parameters from variables array */
	x0 = vars[0];
	ps = vars[1];
	pg = vars[2];
	tg = vars[3];
	t0 = vars[4];
	ts = MAX(tsg-tg,0.0);

	/* zero out F and J */
	for (i=0; i<neqns; ++i)
		eqns[i] = 0.0;
	for (j=0; j<neqns; ++j)
		for (i=0; i<neqns; ++i)
			jacob[j][i] = 0.0;

	/* insert constant terms into equations */
	eqns[0] += 2.0*h;
	eqns[1] += h-x0;
	jacob[0][1] -= 1.0;

	/* get x(ps,2*ts) and partial derivatives, load into equations */
	intder(x,nt,dt,ft,np,dp,fp,2.0*ts,ABS(ps),&ftp,&dfdt,&dfdp);
	sign = SGN(ps);
	eqns[0] -= sign*ftp;
	jacob[1][0] -= dfdp;
	jacob[3][0] += 2.0*sign*dfdt;

	/* get x(pg,2*tg) and partial derivatives, load into equations  */
	intder(x,nt,dt,ft,np,dp,fp,2.0*tg,ABS(pg),&ftp,&dfdt,&dfdp);
	sign = SGN(pg);
	eqns[0] += sign*ftp;
	eqns[1] += sign*ftp;
	jacob[2][0] += dfdp;
	jacob[2][1] += dfdp;
	jacob[3][0] += 2.0*sign*dfdt;
	jacob[3][1] += 2.0*sign*dfdt;

	/* get x(p0,t0) and partial derivatives, load into equations */
	intder(x,nt,dt,ft,np,dp,fp,t0,ABS(p0),&ftp,&dfdt,&dfdp);
	sign = SGN(p0);
	eqns[1] -= sign*ftp;
	jacob[4][1] -= sign*dfdt;

	/* get tau(ps,2*ts) and partial derivatives, load into equations */
	intder(tau,nt,dt,ft,np,dp,fp,2.0*ts,ABS(ps),&ftp,&dfdt,&dfdp);
	sign = SGN(ps);
	eqns[2] -= ftp;
	jacob[1][2] -= sign*dfdp;
	jacob[3][2] += 2.0*dfdt;

	/* get tau(pg,2*tg) and partial derivatives, load into equations */
	intder(tau,nt,dt,ft,np,dp,fp,2.0*tg,ABS(pg),&ftp,&dfdt,&dfdp);
	sign = SGN(pg);
	eqns[2] += ftp;
	eqns[3] += ftp;
	jacob[2][2] += sign*dfdp;
	jacob[2][3] += sign*dfdp;
	jacob[3][2] += 2.0*dfdt;
	jacob[3][3] += 2.0*dfdt;

	/* get tau(p0,t0) and partial derivatives, load into equations */
	intder(tau,nt,dt,ft,np,dp,fp,t0,ABS(p0),&ftp,&dfdt,&dfdp);
	sign = SGN(p0);
	eqns[3] -= ftp;
	jacob[4][3] -= dfdt;

	/* get a(ps,2*ts) and partial derivatives, load into equations */
	intder(a,nt,dt,ft,np,dp,fp,2.0*ts,ABS(ps),&ftp,&dfdt,&dfdp);
	sign = SGN(ps);
	eqns[4] += sign*ftp;
	jacob[1][4] += dfdp;
	jacob[3][4] -= 2.0*sign*dfdt;

	/* get a(pg,2*tg) and partial derivatives, load into equations */
	intder(a,nt,dt,ft,np,dp,fp,2.0*tg,ABS(pg),&ftp,&dfdt,&dfdp);
	sign = SGN(pg);
	eqns[4] += sign*ftp;
	jacob[2][4] += dfdp;
	jacob[3][4] += 2.0*sign*dfdt;

	/* get a(p0,t0) and partial derivatives, load into equations */
	intder(a,nt,dt,ft,np,dp,fp,t0,ABS(p0),&ftp,&dfdt,&dfdp);
	sign = SGN(p0);
	eqns[4] -= 2.0*sign*ftp;
	jacob[4][4] -= 2.0*sign*dfdt;
}

static void intder (float **f, int nx, float dx, float fx,
	int ny, float dy, float fy, float x, float y,
	float *fxy, float *dfdx, float *dfdy)
/***************************************************************************
Finds f(x,y), df/dx, df/dy given function values on a mesh f[iy][ix],
via linear interpolation and finite differencing.
****************************************************************************
Input:
f	array[ny][nx] of sampled function values
nx	number of samples in x direction (fast dimension)
dx	x sampling interval
fx	first sample in x direction
ny	number of samples in y direction (slow dimension)
dy	y sampling interval
fy	first sample in y direction
x	x value at which to evaluate functions
y	y value at which to evaluate functions

Output:
fxy	interpolated function value f(x,y)
dfdx	partial derivative df/dx evaluated at (x,y)
dfdy	partial derivative df/dy evaluated at (x,y)
****************************************************************************
Note: Constant extrapolation of function and zero derivatives for points
	outside the table.  This is desirable, but somewhat awkward.
****************************************************************************
Author: Craig Artley, Colorado School of Mines, 08/31/91
****************************************************************************/
{
	int	ix0,ix1,iy0,iy1;
	float	x0i,y0i,fracx,fracy,lx,ly,f00,f10,f01,f11;

	/* compute x and y extent of table */
	lx = fx+(nx-1)*dx;
	ly = fy+(ny-1)*dy;

	if (x<=fx) {
		ix0 = ix1 = 0;
		fracx = 0.0;
	} else if (x>=lx) {
		ix0 = ix1 = nx-1;
		fracx = 1.0;
	} else {
		x0i = (x-fx)/dx;
		ix0 = x0i;
		ix1 = ix0+1;
		fracx = x0i-ix0;
	}

	if (y<=fy) {
		iy0 = iy1 = 0;
		fracy = 0;
	} else if (y>=ly) {
		iy0 = iy1 = ny-1;
		fracy = 1.0;
	} else {
		y0i = (y-fy)/dy;
		iy0 = y0i;
		iy1 = iy0+1;
		fracy = y0i-iy0;
	}

	/* get tabled function values nearest (x,y) */
	f00 = f[iy0][ix0];	f01 = f[iy1][ix0];
	f10 = f[iy0][ix1];	f11 = f[iy1][ix1];

	/* linear interpolation to find f(x,y) */
	*fxy = (1-fracx)*(1-fracy)*f00 + (1-fracx)*fracy*f01
		+ fracx*(1-fracy)*f10 + fracx*fracy*f11;

	/* interpolate and finite difference to find partial derivatives */
	*dfdx = ((1-fracy)*(f10-f00) + fracy*(f11-f01))/dx;
	*dfdy = ((1-fracx)*(f01-f00) + fracx*(f11-f10))/dy;
}

static void rayvt (float p, int nt, float dt,
	float v[], float tau[], float x[], float a[]) 
/*****************************************************************************
trace ray for v = v(tau), tau=ONE-way vertical time
******************************************************************************
Input:
p		ray parameter
nt		number of times (including t=0.0)
dt		time sampling interval
v		array[nt] of velocities (v = v(tau))

Output:
tau		array[nt] of vertical times tau(t)
x		array[nt] of horizontal distances x(t)
a		array[nt] of propagation angles a(t)
******************************************************************************
Notes:
Ray tracing begins at tau=x=t=0.0.  If the ray turns and reaches tau=0.0
at some time index it, then tau[it+1:nt-1] = tau[it], x[it+1:nt-1] = x[it],
and a[it+1:nt-1] = a[it].  In other words, constant extrapolation is used
to fill the remaining values in the output arrays.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 08/30/90
Modified:  Craig Artley, Colorado School of Mines, 09/25/90
	   Fixed bug in calculating propagation angles.
	   Added linear interpolation between velocity samples.
******************************************************************************/
{
	int it,jt,itau;
	float vt,dvt,at,taut,xt,cosat,taui,frac;

	/* initialize ray tracing parameters */
	vt = v[0];
	dvt = v[1]-v[0];
	at = asin(MIN(1.0,p*vt));
	taut = xt = 0.0;
	tau[0] = taut;
	x[0] = xt;
	a[0] = at;

	/* loop over times t */
	for (it=1; it<nt; ++it) {

		/* compute cosine of propagation angle */
		cosat = cos(at);

		/* update vertical time */
		tau[it] = tau[it-1]+dt*cosat;

		/* if ray has returned to surface, break */
		if (tau[it]<=0.0) break;

		/* update horizontal distance */
		x[it] = x[it-1]+dt*p*vt*vt;

		/* update angle */
		at += p*dvt;
		a[it] = at;

		/* update velocity and derivative */
		if (it<nt-1) {
			taui = tau[it]/dt;
			itau = taui;
			frac = taui-itau;
			vt = (1.0-frac)*v[itau] + frac*v[itau+1];
			dvt = v[itau+1]-v[itau];
		}
	}

	/* if ray returned to surface, extrapolate with surface values */
	for (jt=it; jt<nt; ++jt) {
		tau[jt] = tau[jt-1];
		x[jt] = x[jt-1];
		a[jt] = a[jt-1];
	}
}

/* Library version of SUDMOVZ
 *
 * The main program (the parameters, the traces and their headers, the offsets and the mixes of them, the files of
 * velocities) is left to the caller. The tables of the rays of the velocity are made once, and kept in a model:
 *
 *	su_dmovz_create:	the ray tables and the rms velocity from the interval velocities vt[nt] (as the program makes them)
 *	su_dmovz_nxfft:		the number of cdps of the common-offset gather after zero padding for the Fourier transform
 *	su_dmovz_offset:	DMO of a common-offset gather q[nxfft + 2][nt] (flat), in place
 *	su_dmovz_free:		frees the model
 *
 * Differences from the program:
 *  - The gather is a flat array, q[(nxfft + 2) * nt].
 *  - dmomap read rcond (the condition number that the solver of the system sets) before it was set when the first iteration
 *    had converged; it starts as 1.
 *  - A row of the table of times that has no sample set ended the loop over the rows only after the program had
 *    read the sample after the end of the row (the first of the next row, or past the table for the last).
 *  - A gather of zero offset is left as it is (the program divided by the offset to find the sampling of the ray parameter).
 *  - A gather that the stretch mute leaves nothing of gives a code (the program ended with an error).
 */

struct DmoVz {
	int nt, nttab, nptab;
	float dt, dptab, fp, lp;
	float *vrms;
	float **x, **a, **tau;
};

struct DmoVz *su_dmovz_create(int nt, float dt, const float *vt_in, int nptab, int *err)
{
	struct DmoVz *M;
	int it, ip, nttab = 2*nt-1;
	float *vt, sum, p;

	*err = 0;
	if (nt<2 || dt<=0.0 || nptab<2) { *err = -3; return NULL; }
	for (it=0; it<nt; ++it) if (!(vt_in[it]>0.0)) { *err = -3; return NULL; }
	M = (struct DmoVz*)calloc(1,sizeof(struct DmoVz));
	vt = ealloc1float(nttab);
	M->vrms = ealloc1float(nttab);
	memcpy(vt,vt_in,sizeof(float)*nt);

	/* constant extrapolation of interval velocities */
	for (it=nt; it<nttab; ++it)
		vt[it] = vt[nt-1];

	/* compute vrms */
	M->vrms[0] = vt[0];
	for (it=1,sum=0.0; it<nttab; ++it) {
		sum += vt[it-1]*vt[it-1];
		M->vrms[it] = sqrt(sum/it);
	}

	/* compute half-velocities */
	for (it=0; it<nttab; ++it)
		vt[it] *= 0.5;

	M->nt = nt; M->nttab = nttab; M->nptab = nptab; M->dt = dt;
	M->fp = 0.0;
	M->lp = 1.0/vt[0];
	M->dptab = M->lp/(nptab-1);

	/* compute ray tables tau(p,t), x(p,t), a(p,t) */
	M->tau = alloc2float(nttab,nptab);
	M->x = alloc2float(nttab,nptab);
	M->a = alloc2float(nttab,nptab);
	for (ip=0,p=M->fp; ip<nptab; ++ip,p+=M->dptab)
		rayvt(p,nttab,dt,vt,M->tau[ip],M->x[ip],M->a[ip]);
	free1float(vt);
	return M;
}

int su_dmovz_nxfft(int nx, float dx, float offset)
{
	return npfar(nx+(int)(0.5*ABS(offset/dx)));
}

int su_dmovz_offset(const struct DmoVz *M, int nx, float dx, float offset, float fmax, float smute, float speed, float *q)
{
	int nxfft = su_dmovz_nxfft(nx,dx,offset), i, status;
	float **rows = (float**)ealloc1(nxfft+2,sizeof(float*));

	if (dx==0.0 || fmax<=0.0 || smute<=1.0 || speed<=0.0) { free1(rows); return -3; }
	/* a gather of zero offset is not changed by DMO (the program made a table of no rays, and divided by the offset) */
	if (offset==0.0) { free1(rows); return 0; }
	for (i=0; i<nxfft+2; ++i) rows[i] = q+(size_t)i*M->nt;
	status = dmooff(M->x,M->a,M->tau,M->nttab,M->dt,0.0,M->nptab,M->dptab,M->fp,M->lp,nx,dx,offset,M->nt,M->vrms,fmax,smute,
		speed,0,rows);
	free1(rows);
	return status;
}

void su_dmovz_free(struct DmoVz *M)
{
	if (M==NULL) return;
	free1float(M->vrms);
	free2float(M->x); free2float(M->a); free2float(M->tau);
	free(M);
}
