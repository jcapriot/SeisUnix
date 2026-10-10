/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFDMOD2: $Revision: 1.29 $ ; $Date: 2015/06/02 20:15:23 $        */

/* sufdmod2 - finite difference modeling by simple 2nd order explict method */

#include "par.h"
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUFDMOD2 - Finite-Difference MODeling (2nd order) for acoustic wave equation",
" 									",
" sufdmod2 <vfile >wfile nx= nz= tmax= xs= zs= [optional parameters]	",
" 									",
" Required Parameters:							",
" <vfile		file containing velocity[nx][nz]		",
" >wfile		file containing waves[nx][nz] for time steps	",
" nx=			number of x samples (2nd dimension)		",
" nz=			number of z samples (1st dimension)		",
" xs=			x coordinates of source, or, alternatively, the name",
"			of a file that contains the x- and z-coordinates,",
"			with the number of pairs as the first record and",
"			the actual pairs of (x,z) locations following.  ",
" zs=			z coordinates of source				",
" tmax=			maximum time					",
" 									",
" Optional Parameters:							",
" sstrength=1.0		strength of source				",
" pw=0			use point or extended source geometry parameters",
" 			=1  use horizontal plane wave source 		",
" pwt=20		amp taper on ends of line src (in grid points)  ",
" mono=0		use ricker wavelet as source function 		",
" 			=1  use single frequency src (freq=2*fpeak)	",
" nt=1+tmax/dt		number of time samples (dt determined for stability)",
" mt=1			number of time steps (dt) per output time step	",
" 									",
" dx=1.0		x sampling interval				",
" fx=0.0		first x sample					",
" dz=1.0		z sampling interval				",
" fz=0.0		first z sample					",
" 									",
" fmax = vmin/(10.0*h)	maximum frequency in source 			",
" fpeak=0.5*fmax	peak frequency in ricker wavelet		",
" 									",
" dfile=		input file containing density[nx][nz]		",
" vsx=			x coordinate of vertical line of seismograms	",
" hsz=			z coordinate of horizontal line of seismograms	",
" vsfile=		output file for vertical line of seismograms[nz][nt]",
" hsfile=		output file for horizontal line of seismograms[nx][nt]",
" ssfile=		output file for source point seismograms[nt]	",
" verbose=0		=1 for diagnostic messages, =2 for more		",
" abs=1,1,1,1		absorbing boundary conditions on top,left,bottom,right",
"			sides of the model. 				",
" 			=0,1,1,1 for free surface condition on the top	",
" 									",
" Notes:								",
" 									",
" This program uses the traditional explicit second order differencing	",
" method. 								",
" 									",
NULL};
#endif


/*
 * Authors:  CWP:Dave Hale
 *           CWP:modified for SU by John Stockwell, 1993.
 *           U Houston: added plane wave and monochromatic wave 
 *                        source options.  Chris Liner, 2010
 *
 *
 * Trace header fields set: sx, gx, ns, delrt, tracl, tracr, offset, d1, d2,
 *                          sdepth, trid
 *
 * Modifications: Tony Kocurko (TK:)
 *                Memorial University in Newfoundland and Labrador
 *                - Allow user to supply the name of a file containing
 *                  shot point locations, rather than supplying them
 *                  as values to the xs= and zs= command line arguments.
 *                - Correct the calculation of izs[is].
 *
 * Technical reference:
 * Kelly, K. R., R. W. Ward, S. Treitel, and R. M. Alford (1976),
 * Synthetic Seismograms: A finite-difference approach, 
 * Geophysics, Vol. 41. No. I (February, 1976), p. 2-27.
 *
 */
/**************** end self doc ********************************/

static float ricker (float t, float fpeak, int mono);
static void exsrc (int ns, float *xs, float *zs, float *vs, float (*xsd)[4], float (*zsd)[4],
	int nx, float dx, float fx,
	int nz, float dz, float fz,
	float dt, float t, float fpeak, int pwt, int mono, float **s)
/*****************************************************************************
update source pressure function for an extended source
******************************************************************************
Input:
ns		number of x,z coordinates for extended source
xs		array[ns] of x coordinates of extended source
zs		array[ns] of z coordinates of extended source
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
dt		time step (ignored)
t		time at which to compute source function
fpeak		peak frequency (the program took 0.5*fmax, or fmax for the single frequency source, whatever fpeak was)
pwt		src taper in grid points
mono		=0 use ricker src... =1 use monofreq src (2*fpeak)

Output:
s		array[nx][nz] of source pressure at time t+dt
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/01/90
******************************************************************************/
{
	int ix,iz,izv,is;
	float ts,xn,zn,v,xv,zv,dxdv,dzdv,xvn,zvn;
	float amp,dv,dist,distprev;
	float a, pio2, opwt;
	float tdelay;

	/* zero source array */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			s[ix][iz] = 0.0 *dt ;
	
	/* compute time-dependent part of source function */
	tdelay = 1.0/fpeak;
	if (t>2.0*tdelay && mono==0) return;
	ts = ricker(t-tdelay,fpeak,mono);
	
	/* loop over extended source locations */
	for (v=vs[0],distprev=0.0,dv=1.0; dv!=0.0; distprev=dist,v+=dv) {
		
		/* determine x(v), z(v), dx/dv, and dz/dv along source */
		intcub(0,ns,vs,xsd,1,&v,&xv);
		intcub(0,ns,vs,zsd,1,&v,&zv);
		intcub(1,ns,vs,xsd,1,&v,&dxdv);
		intcub(1,ns,vs,zsd,1,&v,&dzdv);
		
		/* determine increment along extended source */
		if (dxdv==0.0)
			dv = dz/ABS(dzdv);
		else if (dzdv==0.0)
			dv = dx/ABS(dxdv);
		else
			dv = MIN(dz/ABS(dzdv),dx/ABS(dxdv));
		if (v+dv>vs[ns-1]) dv = vs[ns-1]-v;
		dist = dv*sqrt(dzdv*dzdv+dxdv*dxdv)/sqrt(dx*dx+dz*dz);
		
		/* determine source amplitude */
		amp = (dist+distprev)/2.0;
		
		/* let source contribute within limited distance */
		xvn = (xv-fx)/dx;
		zvn = (zv-fz)/dz;
		izv = NINT(zvn);
		pio2 = 3.1415926/2.0;
		opwt = 1.0 / (float) pwt;
		for (ix=0; ix<nx; ++ix) {
			for (iz=MAX(0,izv-3); iz<=MIN(nz-1,izv+3); ++iz) {
				xn = ix-xvn;
				zn = iz-zvn;
				a = 1.0;
				if (ix < pwt) {
				  a = sin( pio2*((float) ix * opwt ));
				} 
				if (ix > nx-pwt) {
				  a = cos((ix - (nx-pwt))*pio2*opwt); 
				}
				s[ix][iz] += a*ts*amp*exp(-xn*xn-zn*zn);
			}
		}
	}
}

static void ptsrc (float sstrength, float xs, float zs,
	int nx, float dx, float fx,
	int nz, float dz, float fz,
	float dt, float t, float fmax, float fpeak, float tdelay, int mono, float **s)
/*****************************************************************************
update source pressure function for a point source
******************************************************************************
Input:
xs		x coordinate of point source
zs		z coordinate of point source
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
dt		time step (ignored)
t		time at which to compute source function
fmax		maximum frequency (ignored)
fpeak		peak frequency

Output:
tdelay		time delay of beginning of source function
s		array[nx][nz] of source pressure at time t+dt
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/01/90
******************************************************************************/
{
	int ix,iz,ixs,izs;
	float ts,xn,zn,xsn,zsn;
	
	/* zero source array */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			s[ix][iz] = 0.0 * dt*fmax;
	
	/* compute time-dependent part of source function */
	/* fpeak = 0.5*fmax;  this is now getparred */
	tdelay = 1.0/fpeak;
	if (t>2.0*tdelay && mono==0) return;
	ts = ricker(t-tdelay,fpeak,mono);
	
	/* let source contribute within limited distance */
	xsn = (xs-fx)/dx;
	zsn = (zs-fz)/dz;
	ixs = NINT(xsn);
	izs = NINT(zsn);
	for (ix=MAX(0,ixs-3); ix<=MIN(nx-1,ixs+3); ++ix) {
		for (iz=MAX(0,izs-3); iz<=MIN(nz-1,izs+3); ++iz) {
			xn = ix-xsn;
			zn = iz-zsn;
			s[ix][iz] = sstrength*ts*exp(-xn*xn-zn*zn);
		}
	}
}

static float ricker (float t, float fpeak, int mono)
/*****************************************************************************
Compute Ricker wavelet as a function of time
******************************************************************************
Input:
t		time at which to evaluate Ricker wavelet
fpeak		peak (dominant) frequency of wavelet
mono		=0 use ricker... =1 use single frequency (2*fpeak)  CLL 11/27/06
******************************************************************************
Notes:
The amplitude of the Ricker wavelet at a frequency of 2.5*fpeak is 
approximately 4 percent of that at the dominant frequency fpeak.
The Ricker wavelet effectively begins at time t = -1.0/fpeak.  Therefore,
for practical purposes, a causal wavelet may be obtained by a time delay
of 1.0/fpeak.
The Ricker wavelet has the shape of the second derivative of a Gaussian.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 04/29/90
******************************************************************************/
{
	float x,xx;
	
	x = PI*fpeak*t;
	xx = x*x;
	/* return (-6.0+24.0*xx-8.0*xx*xx)*exp(-xx); */
	/* return PI*fpeak*(4.0*xx*x-6.0*x)*exp(-xx); */
	/* return exp(-xx)*(1.0-2.0*xx); */
	if (mono==0) {
		return exp(-xx)*(1.0-2.0*xx);
	} else {
		return sin(4*x);
	}
}

/* 2D finite differencing subroutine */

/* functions declared and used internally */
static void star1 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp);
static void star2 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp);
static void star3 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp);
static void star4 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp);
static void absorb (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **pm, float **p, float **pp,
	int *abs);

static void tstep2 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp, int *abs)
/*****************************************************************************
One time step of FD solution (2nd order in space) to acoustic wave equation
******************************************************************************
Input:
nx		number of x samples
dx		x sampling interval
nz		number of z samples
dz		z sampling interval
dt		time step
dvv		array[nx][nz] of density*velocity^2
od		array[nx][nz] of 1/density (NULL for constant density=1.0)
s		array[nx][nz] of source pressure at time t+dt
pm		array[nx][nz] of pressure at time t-dt
p		array[nx][nz] of pressure at time t

Output:
pp		array[nx][nz] of pressure at time t+dt
******************************************************************************
Notes:
This function is optimized for special cases of constant density=1 and/or
equal spatial sampling intervals dx=dz.  The slowest case is variable
density and dx!=dz.  The fastest case is density=1.0 (od==NULL) and dx==dz.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/13/90
******************************************************************************/
{
	/* convolve with finite-difference star (special cases for speed) */
	if (od!=NULL && dx!=dz) {
		star1(nx,dx,nz,dz,dt,dvv,od,s,pm,p,pp);
	} else if (od!=NULL && dx==dz) {
		star2(nx,dx,nz,dz,dt,dvv,od,s,pm,p,pp);
	} else if (od==NULL && dx!=dz) {
		star3(nx,dx,nz,dz,dt,dvv,od,s,pm,p,pp);
	} else {
		star4(nx,dx,nz,dz,dt,dvv,od,s,pm,p,pp);
	}
	
	/* absorb along boundaries */
	absorb(nx,dx,nz,dz,dt,dvv,od,pm,p,pp,abs);
}

/* convolve with finite-difference star for variable density and dx!=dz */
static void star1 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp)
{
	int ix,iz;
	float xscale1,zscale1,xscale2,zscale2;
		
	/* determine constants */
	xscale1 = (dt*dt)/(dx*dx);
	zscale1 = (dt*dt)/(dz*dz);
	xscale2 = 0.25*xscale1;
	zscale2 = 0.25*zscale1;
	
	/* do the finite-difference star */
	for (ix=1; ix<nx-1; ++ix) {
		for (iz=1; iz<nz-1; ++iz) {
			pp[ix][iz] = 2.0*p[ix][iz]-pm[ix][iz] +
				dvv[ix][iz]*(
					od[ix][iz]*(
						xscale1*(
							p[ix+1][iz]+
							p[ix-1][iz]-
							2.0*p[ix][iz]
						) +
						zscale1*(
							p[ix][iz+1]+
							p[ix][iz-1]-
							2.0*p[ix][iz]
						)
					) +
					(
						xscale2*(
							(od[ix+1][iz]-
							od[ix-1][iz]) *
							(p[ix+1][iz]-
							p[ix-1][iz])
						) +
						zscale2*(
							(od[ix][iz+1]-
							od[ix][iz-1])*
							(p[ix][iz+1]-
							p[ix][iz-1])
						)
					)
				) +
				s[ix][iz];
		}
	}
}

/* convolve with finite-difference star for variable density and dx==dz */
static void star2 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp)
{
	int ix,iz;
	float scale1,scale2;
	
	if ( dx != dz ) 
		fprintf(stderr,"ASSERT FAILED: dx != dz in star2\n");
	/* determine constants */
	scale1 = (dt*dt)/(dx*dx);
	scale2 = 0.25*scale1;
	
	/* do the finite-difference star */
	for (ix=1; ix<nx-1; ++ix) {
		for (iz=1; iz<nz-1; ++iz) {
			pp[ix][iz] = 2.0*p[ix][iz]-pm[ix][iz] +
				dvv[ix][iz]*(
					od[ix][iz]*(
						scale1*(
							p[ix+1][iz]+
							p[ix-1][iz]+
							p[ix][iz+1]+
							p[ix][iz-1]-
							4.0*p[ix][iz]
						)
					) +
					(
						scale2*(
							(od[ix+1][iz]-
							od[ix-1][iz]) *
							(p[ix+1][iz]-
							p[ix-1][iz]) +
							(od[ix][iz+1]-
							od[ix][iz-1]) *
							(p[ix][iz+1]-
							p[ix][iz-1])
						)
					)
				) +
				s[ix][iz];
		}
	}
}

/* convolve with finite-difference star for density==1.0 and dx!=dz */
static void star3 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp)
{
	int ix,iz;
	float xscale,zscale;
		
	if ( od != ((float **) NULL) ) 
		fprintf(stderr,"ASSERT FAILED: od !=  NULL in star3\n");
	/* determine constants */
	xscale = (dt*dt)/(dx*dx);
	zscale = (dt*dt)/(dz*dz);
	
	/* do the finite-difference star */
	for (ix=1; ix<nx-1; ++ix) {
		for (iz=1; iz<nz-1; ++iz) {
			pp[ix][iz] = 2.0*p[ix][iz]-pm[ix][iz] +
				dvv[ix][iz]*(
					xscale*(
						p[ix+1][iz]+
						p[ix-1][iz]-
						2.0*p[ix][iz]
					) +
					zscale*(
						p[ix][iz+1]+
						p[ix][iz-1]-
						2.0*p[ix][iz]
					)
				) +
				s[ix][iz];
		}
	}
}

/* convolve with finite-difference star for density==1.0 and dx==dz */
static void star4 (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **s,
	float **pm, float **p, float **pp)
{
	int ix,iz;
	float scale;
	
	/* determine constants */
	if ( od != ((float **) NULL) ) 
		fprintf(stderr,"ASSERT FAILED: od !=  NULL in star4\n");
	if ( dz != dx ) 
		fprintf(stderr,"ASSERT FAILED: dz !=  dx in star4\n");
	scale = (dt*dt)/(dx*dz);
	
	/* do the finite-difference star */
	for (ix=1; ix<nx-1; ++ix) {
		for (iz=1; iz<nz-1; ++iz) {
			pp[ix][iz] = 2.0*p[ix][iz]-pm[ix][iz] +
				scale*dvv[ix][iz]*(
					p[ix+1][iz]+
					p[ix-1][iz]+
					p[ix][iz+1]+
					p[ix][iz-1]-
					4.0*p[ix][iz]
				) +
				s[ix][iz];
		}
	}
}

static void absorb (int nx, float dx, int nz, float dz, float dt,
	float **dvv, float **od, float **pm, float **p, float **pp,
	int *abs)
{
	int ix,iz;
	float ov,ovs,cosa,beta,gamma,dpdx,dpdz,dpdt,dpdxs,dpdzs,dpdts;

	/* solve for upper boundary */
	iz = 1;
	for (ix=0; ix<nx; ++ix) {

		if (abs[0]!=0) {

			if (od!=NULL)
				ovs = 1.0/(od[ix][iz]*dvv[ix][iz]);
			else
				ovs = 1.0/dvv[ix][iz];
			ov = sqrt(ovs);
			if (ix==0)
				dpdx = (p[1][iz]-p[0][iz])/dx;
			else if (ix==nx-1)
				dpdx = (p[nx-1][iz]-p[nx-2][iz])/dx;
			else
				dpdx = (p[ix+1][iz]-p[ix-1][iz])/(2.0*dx);
			dpdt = (pp[ix][iz]-pm[ix][iz])/(2.0*dt);
			dpdxs = dpdx*dpdx;
			dpdts = dpdt*dpdt;
			if (ovs*dpdts>dpdxs)
				cosa = sqrt(1.0-dpdxs/(ovs*dpdts));
			else 
				cosa = 0.0;
			beta = ov*dz/dt*cosa;
			gamma = (1.0-beta)/(1.0+beta);

			pp[ix][iz-1] = gamma*(pp[ix][iz]-p[ix][iz-1])+p[ix][iz];
		} else {
			pp[ix][iz-1] = 0.0;
		}
	}


	/* extrapolate along left boundary */
	ix = 1;
	for (iz=0; iz<nz; ++iz) {
		if (abs[1]!=0) {
			if (od!=NULL)
				ovs = 1.0/(od[ix][iz]*dvv[ix][iz]);
			else
				ovs = 1.0/dvv[ix][iz];
			ov = sqrt(ovs);
			if (iz==0)
				dpdz = (p[ix][1]-p[ix][0])/dz;
			else if (iz==nz-1)
				dpdz = (p[ix][nz-1]-p[ix][nz-2])/dz;
			else
				dpdz = (p[ix][iz+1]-p[ix][iz-1])/(2.0*dz);
			dpdt = (pp[ix][iz]-pm[ix][iz])/(2.0*dt);
			dpdzs = dpdz*dpdz;
			dpdts = dpdt*dpdt;
			if (ovs*dpdts>dpdzs)
				cosa = sqrt(1.0-dpdzs/(ovs*dpdts));
			else
				cosa = 0.0;
			beta = ov*dx/dt*cosa;
			gamma = (1.0-beta)/(1.0+beta);
			pp[ix-1][iz] = gamma*(pp[ix][iz]-p[ix-1][iz])+p[ix][iz];
		} else {
			pp[ix-1][iz] = 0.0;
		}
	}

	/* extrapolate along lower boundary */
	iz = nz-2;
	for (ix=0; ix<nx; ++ix) {
		if (abs[2]!=0) {
			if (od!=NULL)
				ovs = 1.0/(od[ix][iz]*dvv[ix][iz]);
			else
				ovs = 1.0/dvv[ix][iz];
			ov = sqrt(ovs);
			if (ix==0)
				dpdx = (p[1][iz]-p[0][iz])/dx;
			else if (ix==nx-1)
				dpdx = (p[nx-1][iz]-p[nx-2][iz])/dx;
			else
				dpdx = (p[ix+1][iz]-p[ix-1][iz])/(2.0*dx);
			dpdt = (pp[ix][iz]-pm[ix][iz])/(2.0*dt);
			dpdxs = dpdx*dpdx;
			dpdts = dpdt*dpdt;
			if (ovs*dpdts>dpdxs)
				cosa = sqrt(1.0-dpdxs/(ovs*dpdts));
			else
				cosa = 0.0;
			beta = ov*dz/dt*cosa;
			gamma = (1.0-beta)/(1.0+beta);

			pp[ix][iz+1] = gamma*(pp[ix][iz]-p[ix][iz+1])+p[ix][iz];
		} else {
			pp[ix][iz+1] = 0.0;
		}
	}

	/* extrapolate along right boundary */
	ix = nx-2;
	for (iz=0; iz<nz; ++iz) {
		if (abs[3]!=0) {
			if (od!=NULL)
				ovs = 1.0/(od[ix][iz]*dvv[ix][iz]);
			else
				ovs = 1.0/dvv[ix][iz];
			ov = sqrt(ovs);
			if (iz==0)
				dpdz = (p[ix][1]-p[ix][0])/dz;
			else if (iz==nz-1)
				dpdz = (p[ix][nz-1]-p[ix][nz-2])/dz;
			else
				dpdz = (p[ix][iz+1]-p[ix][iz-1])/(2.0*dz);
			dpdt = (pp[ix][iz]-pm[ix][iz])/(2.0*dt);
			dpdzs = dpdz*dpdz;
			dpdts = dpdt*dpdt;
			if (ovs*dpdts>dpdzs)
				cosa = sqrt(1.0-dpdzs/(ovs*dpdts));
			else
				cosa = 0.0;
			beta = ov*dx/dt*cosa;
			gamma = (1.0-beta)/(1.0+beta);
			pp[ix+1][iz] =gamma*(pp[ix][iz]-p[ix+1][iz])+p[ix][iz];
		} else {
			pp[ix+1][iz] = 0.0;
		}
	}
}


/* Library version of SUFDMOD2
 *
 * The main program (the parameters, the files, the headers of the traces, the seismograms) is left to the caller. What is here are
 * the pieces of the computation, with no static state, and all of the arrays are owned by the caller. The 2-D arrays are flat, with z
 * the fast axis: [nx][nz].
 *
 *	su_fdmod2_exsrc_setup:	the cubic splines through the points of an extended source (the program made them on its first call)
 *	su_fdmod2_ptsrc:	the source pressure at a time step, for a point source
 *	su_fdmod2_exsrc:	the same for an extended source, or a plane wave
 *	su_fdmod2_tstep:	one time step of the finite-difference solution, with the absorbing boundaries
 *
 * Differences from the program:
 *  - The extended source takes the peak frequency fpeak that is given, where the program used 0.5*fmax (or fmax for the single
 *    frequency source, so that its frequency was 2*fmax and not the documented 2*fpeak).
 *  - The checks of the parameters are the caller's.
 */

static float **fd_rows(const float *base, int n, int m)
{
	int i;
	float **rows = (float**)ealloc1(n,sizeof(float*));
	for (i=0; i<n; ++i) rows[i] = (float*)base+(size_t)i*m;
	return rows;
}

/* vs[ns], xsd[ns][4], zsd[ns][4] are made from the ns points (xs,zs) */
void su_fdmod2_exsrc_setup(int ns, const float *xs, const float *zs, float *vs, float *xsd, float *zsd)
{
	int is;
	for (is=0; is<ns; ++is) vs[is] = is;
	cmonot(ns,vs,(float*)xs,(float(*)[4])xsd);
	cmonot(ns,vs,(float*)zs,(float(*)[4])zsd);
}

void su_fdmod2_ptsrc(float sstrength, float xs, float zs, int nx, float dx, float fx, int nz, float dz, float fz,
	float dt, float t, float fmax, float fpeak, int mono, float *s)
{
	float **rows = fd_rows(s,nx,nz);
	ptsrc(sstrength,xs,zs,nx,dx,fx,nz,dz,fz,dt,t,fmax,fpeak,0.0,mono,rows);
	free1(rows);
}

void su_fdmod2_exsrc(int ns, const float *xs, const float *zs, const float *vs, const float *xsd, const float *zsd,
	int nx, float dx, float fx, int nz, float dz, float fz,
	float dt, float t, float fpeak, int pwt, int mono, float *s)
{
	float **rows = fd_rows(s,nx,nz);
	exsrc(ns,(float*)xs,(float*)zs,(float*)vs,(float(*)[4])xsd,(float(*)[4])zsd,nx,dx,fx,nz,dz,fz,dt,t,fpeak,pwt,mono,rows);
	free1(rows);
}

/* od is NULL for a constant density of 1 */
void su_fdmod2_tstep(int nx, float dx, int nz, float dz, float dt, const float *dvv, const float *od, const float *s,
	const float *pm, const float *p, float *pp, const int *abs)
{
	float **rdvv = fd_rows(dvv,nx,nz), **rod = (od!=NULL) ? fd_rows(od,nx,nz) : NULL, **rs = fd_rows(s,nx,nz);
	float **rpm = fd_rows(pm,nx,nz), **rp = fd_rows(p,nx,nz), **rpp = fd_rows(pp,nx,nz);
	tstep2(nx,dx,nz,dz,dt,rdvv,rod,rs,rpm,rp,rpp,(int*)abs);
	free1(rdvv); if (rod!=NULL) free1(rod); free1(rs); free1(rpm); free1(rp); free1(rpp);
}
