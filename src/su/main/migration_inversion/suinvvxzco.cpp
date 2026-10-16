/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUINVXZCO: $Revision: 1.6 $ ; $Date: 2011/11/16 22:14:43 $	*/

#include "su.h"
#include <vector>
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" SUINVXZCO - Seismic INVersion of Common Offset data for a smooth 	",
"             velocity function V(X,Z) plus a slowness perturbation vp(x,z)",
" 									",
"     suinvvxzco <infile >outfile [optional parameters] 		",
"									",
" Required Parameters:							",
" vfile                  file containing velocity array v[nx][nz]	",
" nx=                    number of x samples (2nd dimension) in velocity",
" nz=                    number of z samples (1st dimension) in velocity",
" nxm=			number of midpoints of input traces		",
"									",
" Optional Parameters:							",
" dt= or from header (dt) 	time sampling interval of input data	",
" offs= or from header (offset) 	source-receiver offset	 	",
" dxm= or from header (d2) 	sampling interval of midpoints 		",
" fxm=0		first midpoint in input trace				",
" nxd=5		skipped number of midpoints (see note)			",
" dx=50.0	x sampling interval of velocity				",
" fx=0.0	first x sample of velocity				",
" dz=50.0	z sampling interval of velocity				",
" nxb=nx/2	band centered at midpoints (see note)			",
" nxc=0		hozizontal range in which velocity is changed		",
" nzc=0		vertical range in which velocity is changed		",
" fxo=0.0	x-coordinate of first output trace 			",
" dxo=15.0	horizontal spacing of output trace 			",
" nxo=101	number of output traces 				",	
" fzo=0.0	z-coordinate of first point in output trace 		",
" dzo=15.0	vertical spacing of output trace 			",
" nzo=101	number of points in output trace			",	
" fmax=0.25/dt	Maximum frequency set for operator antialiasing		",
" ang=180	Maximum dip angle allowed in the image			",
" ls=0		=1 for line source; =0 for point source			",
" pert=0	=1 calculate time correction from v_p[nx][nz]		",
" vpfile	file containing slowness perturbation array v_p[nx][nz]	",
" verbose=1              =1 to print some useful information		",
" 									",
" Notes:								",
" Traveltime and amplitude are calculated by finite difference which	",
" is done only in one of every NXD midpoints; in the skipped midpoint, 	",
" interpolation is used to calculate traveltime and amplitude.		", 
" For each midpoint, traveltime and amplitude are calculated in the 	",
" horizontal range of (xm-nxb*dx, xm+nxb*dx). Velocity is changed by 	",
" constant extropolation in two upper trianglar corners whose width is 	",
" nxc*dx and height is nzc*dz.						",
"									",
" Eikonal equation will fail to solve if there is a polar turned ray.	",
" In this case, the program shows the related geometric information. 	",
" There are three ways to remove the turned ray: smoothing velocity, 	",
" reducing nxb, and increaing nxc and nzc (if the turned ray occurs  	",
" in the shallow areas). To prevent traveltime distortion from a over	",
" smoothed velocity, traveltime is corrected based on the slowness 	",
" perturbation.								",
"									",
" Offsets are signed - may be positive or negative. 			",
"									",
NULL};
#endif

/*
 *
 * Author:  Zhenyue Liu, 08/28/93,  Colorado School of Mines 
 *
 * Reference:
 * Bleistein, N., Cohen, J. K., and Hagin, F., 1987,
 *  Two-and-one-half dimensional Born inversion with an arbitrary reference
 *         Geophysics Vol. 52, no.1, 26-36.
 *
 */
/**************** end self doc ***********************************/

static void delta_t (int na, float da, float r, float dr,
	float uc[], float wc[], float sc[], float bc[],
	float un[], float wn[], float sn[], float bn[]);
static void eiktam_pert (float xs, float zs, int nz, float dz, float fz, int nx,
	float dx, float fx, float **vel, float **vel_p, int pert,
	float **time, float **angle, float **sig, float **bet);


/* Library version of SUINVVXZCO
 *
 * The main program (the parameters, the traces and their headers, the velocity files) is left to the caller.
 *
 *	su_invvxzco:	the 2.5-D Born inversion (migration) of one common-offset gather data[nxm][nt] (the traces in order of
 *			the midpoints fxm + im * dxm, with the offset offs) to out[nxo][nzo] by the traveltimes and amplitudes of
 *			the velocity vel[nx][nz] (z fast; and the slowness perturbation vp[nx][nz] if pert is not 0)
 *
 * Differences from the program: the global variables (ierr, x_err, ...) with which the program would report that the eikonal
 * equation failed with a turned ray are gone (the library stopped the program with err(): it flags it, and the status is -5); the positions that the
 * program read outside of its arrays are kept inside them. Returns 0, -1 if a shot or receiver is outside of the grid or the
 * first output trace is, -3 for parameters that are not allowed, -4 if the band nxb is too small for the offset, -5 if a
 * raypath has too large a curvature (a smoother velocity is required).
 */

/* the amplitudes for the traveltimes at the source and receiver of one midpoint */
static void inv_amplitude(int nx1, int nz, int nz0, float **as, float **ag, float **bas, float **bag, float **sgs,
	float **sgg, float vs0, float vg0, float smax, float cosa, int ls, float **amp)
{
	int ix,iz;
	float temp;
	for (ix=0; ix<nx1; ++ix) {
		amp[ix][nz0] = 0;
		for (iz=nz0+1; iz<nz; ++iz) {
			/* check operator aliasing */
			if (ABS(sin(bas[ix][iz])/vs0+sin(bag[ix][iz])/vg0)>=smax) amp[ix][iz] = 0.;
			/* dip filter for imaging */
			else if (cos((as[ix][iz]+ag[ix][iz])*0.5)<cosa) amp[ix][iz] = 0.;
			else {
				temp = sqrt(vs0*sgg[ix][iz]/(vg0*sgs[ix][iz]));
				amp[ix][iz] = (cos(bas[ix][iz])*temp/vs0+cos(bag[ix][iz])/(temp*vg0))
					*ABS(cos((as[ix][iz]-ag[ix][iz])*0.5));
				/* 2.5-D amplitude correction */
				if (ls==0) amp[ix][iz] *= sqrt(sgs[ix][iz]+sgg[ix][iz]);
			}
		}
	}
}

/* the output points of one midpoint: the times ts + tg and the amplitudes amp, on the grid fx1 + ix dx, iz dz */
static void inv_accumulate(int nxo, float fxo, float dxo, int nzo, float fzo, float dzo, float fx1, float ex1, float ez,
	float odx, float odz, float odt, int nt, int nx1, int nz, float **ts, float **tg, float **amp, const float *trace,
	float **outtrace)
{
	int ixo,izo,ix,iz,nsamp;
	float xo,zo,xi,zi,sx,sz,tsd,tgd,samp,ampd;
	for (ixo=0,xo=fxo; ixo<nxo; ++ixo,xo+=dxo)
		for (izo=1,zo=fzo+dzo; izo<nzo; ++izo,zo+=dzo) {
			/* check range of output */
			if (xo-fx1<0 || xo-ex1>=0 || zo-ez>=0) continue;

			/* determine sample indices */
			xi = (xo-fx1)*odx;
			ix = xi;
			zi = zo*odz;
			iz = zi;
			if (ix<0 || ix>nx1-2 || iz<0 || iz>nz-2) continue;

			/* bilinear interpolation */
			sx = xi-ix;
			sz = zi-iz;
			tsd = (1.0-sz)*((1.0-sx)*ts[ix][iz] + sx*ts[ix+1][iz]) +
				sz*((1.0-sx)*ts[ix][iz+1] + sx*ts[ix+1][iz+1]);
			tgd = (1.0-sz)*((1.0-sx)*tg[ix][iz] + sx*tg[ix+1][iz]) +
				sz*((1.0-sx)*tg[ix][iz+1] + sx*tg[ix+1][iz+1]);
			samp = (tsd+tgd)*odt;
			if (!(samp>=0.0)) continue;
			nsamp = samp;
			if (nsamp>nt-2) continue;
			samp = samp-nsamp;
			ampd = (1.0-sz)*((1.0-sx)*amp[ix][iz] + sx*amp[ix+1][iz]) +
				sz*((1.0-sx)*amp[ix][iz+1] + sx*amp[ix+1][iz+1]);
			outtrace[ixo][izo] += ampd*(samp*trace[nsamp+1]+(1.0-samp)*trace[nsamp]);
		}
}

int su_invvxzco(int nxm, int nt, float dt, float offs, float dxm, float fxm, int nx, int nz, float dx, float fx,
	float dz, int nxd, int nxb, int nxc, int nzc, float fxo, float dxo, int nxo, float fzo, float dzo, int nzo,
	float fmax, float ang, int ls, int pert, const float *vel, const float *vp_in, const float *data, float *out)
{
	int ixm,ixo,izo,ix,iz,i,nx1,nx0,nz0,nxd1,ixd;
	float temp,temp1,xm,xs,xg,ex,ez,fx1,ex1,vs0,vg0,odx,odz,odt,smax,cosa,scal;
	float **outtrace,**vold,**vpold=NULL,**ts,**as,**sgs,**bas,**tg,**ag,**sgg,**bag,**amp,**v,**vp=NULL,
		**ts1=NULL,**tg1=NULL,**amp1=NULL;

	eikpex_error = 0;
	if (nxm<1 || nt<3 || dt<=0.0 || dxm<=0.0 || nx<2 || nz<2 || dx<=0.0 || dz<=0.0 || nxd<1 || nxb<1 || nxc<0 ||
		nzc<0 || nxo<1 || nzo<1 || fmax<=0.0) return -3;
	for (ix=0; ix<nx*nz; ix++) if (vel[ix]<=0.0) return -3;
	nxc = MIN(nxc,nxb);
	nzc = MIN(nzc,nz);
	if (fxo<fx) return -1;
	cosa = cos(ang*PI/180);

	/* check the ranges of shots and receivers */
	ex = fx+(nx-1)*dx;
	ez = (nz-1)*dz;
	for (ixm=0; ixm<nxm; ++ixm) {
		xs = fxm+ixm*dxm-0.5*offs;
		xg = xs+offs;
		if (fx>xs || ex<xs || fx>xg || ex<xg) return -1;
	}
	if (ABS(offs)>2*nxb*dx) return -4;

	odx = 1.0/dx;
	odz = 1.0/dz;
	odt = 1.0/dt;
	smax = 0.5/(dxm*fmax);

	/* allocate space */
	vold = alloc2float(nz,nx);
	if (pert) vpold = alloc2float(nz,nx);
	outtrace = alloc2float(nzo,nxo);
	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz; iz++) vold[ix][iz] = vel[(size_t)ix*nz+iz];
	if (pert) {
		/* calculate 1/2 perturbation of slowness squares */
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				vpold[ix][iz] = vp_in[(size_t)ix*nz+iz]/vold[ix][iz];
	}
	memset((void *) outtrace[0], 0, (size_t) nxo*nzo*sizeof(float));

	nx1 = 1+2*nxb;
	ts = alloc2float(nz,nx1); as = alloc2float(nz,nx1); sgs = alloc2float(nz,nx1); bas = alloc2float(nz,nx1);
	tg = alloc2float(nz,nx1); ag = alloc2float(nz,nx1); bag = alloc2float(nz,nx1); sgg = alloc2float(nz,nx1);
	amp = alloc2float(nz,nx1); v = alloc2float(nz,nx1);
	if (pert) vp = alloc2float(nz,nx1);
	if (nxd>1) {
		ts1 = alloc2float(nz,nx1); tg1 = alloc2float(nz,nx1); amp1 = alloc2float(nz,nx1);
	}
	memset((void *) ts[0],0,(size_t)nx1*nz*sizeof(float));
	memset((void *) tg[0],0,(size_t)nx1*nz*sizeof(float));
	memset((void *) amp[0],0,(size_t)nx1*nz*sizeof(float));
	if (nxd>1) {
		memset((void *) ts1[0],0,(size_t)nx1*nz*sizeof(float));
		memset((void *) tg1[0],0,(size_t)nx1*nz*sizeof(float));
		memset((void *) amp1[0],0,(size_t)nx1*nz*sizeof(float));
	}

	/* save the skipping number */
	nxd1 = nxd;

	/* the first z-sample of the velocity model which corresponds to fzo */
	nz0 = MAX(0,fzo/dz);
	if (nz0>nz-2) nz0 = nz-2;

	/* loop over midpoints */
	for (ixm=0; ixm<nxm; ixm+=nxd1) {
		xm = fxm+ixm*dxm;
		xs = xm-0.5*offs;
		xg = xs+offs;

		/* set range for traveltimes' calculation */
		fx1 = xm-nxb*dx;
		ex1 = MIN(ex+(nxd1-1)*dxm,xm+nxb*dx);
		nx1 = 1+NINT((ex1-fx1)/dx);
		if (nx1>1+2*nxb) nx1 = 1+2*nxb;
		nx0 = (fx1-fx)/dx;
		temp = (fx1-fx)/dx-nx0;

		/* transpose velocity such that the first row is at fx1 */
		for (ix=0; ix<nx1; ++ix) {
			if (ix<-nx0)
				for (iz=0; iz<nz; ++iz) v[ix][iz] = vold[0][iz];
			else if (ix+nx0>nx-2)
				for (iz=0; iz<nz; ++iz) v[ix][iz] = vold[nx-1][iz];
			else
				for (iz=0; iz<nz; ++iz) v[ix][iz] = vold[ix+nx0][iz]*(1.0-temp)+temp*vold[ix+nx0+1][iz];
		}
		if (pert) {
			for (ix=0; ix<nx1; ++ix) {
				if (ix<-nx0)
					for (iz=0; iz<nz; ++iz) vp[ix][iz] = vpold[0][iz];
				else if (ix+nx0>nx-2)
					for (iz=0; iz<nz; ++iz) vp[ix][iz] = vpold[nx-1][iz];
				else
					for (iz=0; iz<nz; ++iz) vp[ix][iz] = vpold[ix+nx0][iz]*(1.0-temp)+temp*vpold[ix+nx0+1][iz];
			}
		}

		/* change velocity values within two trianglar upper corners to avoid possible polar-turned rays */
		for (iz=0; iz<nzc; ++iz)
			for (ix=iz*nxc/nzc; ix<nxc; ++ix) {
				int ia = nxc-ix-1, ib = nxc-iz*nxc/nzc, ic = nx1-nxc+ix, id = nx1-nxc+iz*nxc/nzc-1;
				if (ia>=0 && ia<nx1 && ib>=0 && ib<nx1) v[ia][iz] = v[ib][iz];
				if (ic>=0 && ic<nx1 && id>=0 && id<nx1) v[ic][iz] = v[id][iz];
			}

		/* the velocities at the shot and receiver */
		temp = (xs-fx)/dx;
		ixd = temp;
		if (ixd>nx-2) ixd = nx-2;
		temp = temp-ixd;
		vs0 = (1-temp)*vold[ixd][0]+temp*vold[ixd+1][0];
		temp = (xg-fx)/dx;
		ixd = temp;
		if (ixd>nx-2) ixd = nx-2;
		temp = temp-ixd;
		vg0 = (1-temp)*vold[ixd][0]+temp*vold[ixd+1][0];

		if (ixm==0 || nxd1==1) {
			/* No interpolation for midpoints */
			eiktam_pert(xs,0.,nz,dz,0.,nx1,dx,fx1,v,vp,pert,ts,as,sgs,bas);
			eiktam_pert(xg,0.,nz,dz,0.,nx1,dx,fx1,v,vp,pert,tg,ag,sgg,bag);
			inv_amplitude(nx1,nz,nz0,as,ag,bas,bag,sgs,sgg,vs0,vg0,smax,cosa,ls,amp);
			inv_accumulate(nxo,fxo,dxo,nzo,fzo,dzo,fx1,ex1,ez,odx,odz,odt,nt,nx1,nz,ts,tg,amp,
				data+(size_t)ixm*nt,outtrace);
		} else {
			/* Linear interpolation for midpoints */
			eiktam_pert(xs,0,nz,dz,0,nx1,dx,fx1,v,vp,pert,ts1,as,sgs,bas);
			eiktam_pert(xg,0,nz,dz,0,nx1,dx,fx1,v,vp,pert,tg1,ag,sgg,bag);
			inv_amplitude(nx1,nz,nz0,as,ag,bas,bag,sgs,sgg,vs0,vg0,smax,cosa,ls,amp1);

			/* interpolate quantities between two midpoints */
			xm -= nxd1*dxm;
			for (i=1; i<=nxd1; ++i) {
				xm += dxm;
				xs = xm-0.5*offs;
				xg = xs+offs;
				fx1 = xm-nxb*dx;
				ex1 = MIN(ex+(nxd1-1)*dxm,xm+nxb*dx);
				nx1 = 1+NINT((ex1-fx1)/dx);
				if (nx1>1+2*nxb) nx1 = 1+2*nxb;
				temp = nxd1-i;
				temp1 = 1.0/(nxd1-i+1);
				for (ix=0; ix<nx1; ++ix)
					for (iz=nz0; iz<nz; ++iz) {
						ts[ix][iz] = (temp*ts[ix][iz]+ts1[ix][iz])*temp1;
						tg[ix][iz] = (temp*tg[ix][iz]+tg1[ix][iz])*temp1;
						amp[ix][iz] = (temp*amp[ix][iz]+amp1[ix][iz])*temp1;
					}

				/* set the right boundary for output traces */
				ex1 = MIN(ex,xm+nxb*dx);
				if (ixm-nxd1+i>=0 && ixm-nxd1+i<nxm)
					inv_accumulate(nxo,fxo,dxo,nzo,fzo,dzo,fx1,ex1,ez,odx,odz,odt,nt,nx1,nz,ts,tg,amp,
						data+(size_t)(ixm-nxd1+i)*nt,outtrace);
			}
		}

		/* recount the skipping number for the last midpoint */
		if (ixm<nxm-1 && ixm>nxm-1-nxd1) nxd1 = nxm-1-ixm;
	}

	/* scale the output traces */
	scal = 4/sqrt(PI)*dxm;
	for (ixo=0; ixo<nxo; ++ixo)
		for (izo=0; izo<nzo; ++izo)
			out[(size_t)ixo*nzo+izo] = outtrace[ixo][izo]*scal;

	free2float(vold); free2float(ts); free2float(bas); free2float(sgs); free2float(as); free2float(v);
	free2float(tg); free2float(bag); free2float(sgg); free2float(ag); free2float(amp); free2float(outtrace);
	if (pert) { free2float(vp); free2float(vpold); }
	if (nxd>1) { free2float(ts1); free2float(tg1); free2float(amp1); }
	if (eikpex_error) return -5;
	return 0;
}

static void delta_t (int na, float da, float r, float dr, 
	float uc[], float wc[], float sc[], float bc[],
	float un[], float wn[], float sn[], float bn[])
/*****************************************************************************
difference equation extrapolation of "delta t" in polar coordinates
******************************************************************************
Input:
na		number of a samples
da		a sampling interval
r		current radial distance r
dr		radial distance to extrapolate
uc		array[na] of dt/dr at current r
un		array[na] of dt/dr at next r
wc		array[na] of dt/da at current r
wn		array[na] of dt/da at next r
sc		array[na] of dt/dv at current r

Output:
sn		array[na] of dt/dv at next r 
******************************************************************************

This function implements the Crank-Nicolson finite-difference method with
boundary conditions sn[1]=sn[0] and sn[na-1]=sn[na-2].
******************************************************************************
Author:  Zhenyue Liu, Colorado School of Mines, 07/8/92
******************************************************************************/
{
	int i;
	float r1,*d,*b,*c,*e;
	
	/* allocate workspace */
	d = alloc1float(na-2);
	b = alloc1float(na-2);
	c = alloc1float(na-2);
	e = alloc1float(na-2);
	
	r1 = r+dr;
 	
	/* Crank-Nicolson */
 	for (i=0; i<na-2; ++i) {
		d[i] = (uc[i+1]+un[i+1])/(2.0*dr);
		e[i] = (wn[i+1]/(r1*r1)+wc[i+1]/(r*r))/(8.0*da);
		b[i] = (bc[i+1]+bn[i+1])*0.5-(sc[i+2]-sc[i])*e[i]
			+d[i]*sc[i+1];
		c[i] = -e[i];
	} 
	d[0] += c[0];
	d[na-3] += e[na-3]; 
	
	tripp(na-2,d,e,c,b);
	for(i=0;i<na-2; ++i) sn[i+1]=b[i];
	sn[0] = sn[1];
	sn[na-1] = sn[na-2];
	
	
	/* free workspace */
	free1float(d);
	free1float(c);
	free1float(e);
	free1float(b);
}

/* functions defined and used internally */
static void eiktam_pert (float xs, float zs, int nz, float dz, float fz, int nx, 
	float dx, float fx, float **vel, float **vel_p, int pert,
	float **time, float **angle, float **sig, float **bet)
/*****************************************************************************
Compute traveltimes t(x,z) and  propagation angle a(x,z) via eikonal equation,
 and sigma sig(x,z), incident angle bet(x,z) via Crank-Nicolson Method
******************************************************************************
Input:
xs		x coordinate of source (must be within x samples)
zs		z coordinate of source (must be within z samples)
nz		number of z samples
dz		z sampling interval
fz		first z sample
nx		number of x samples
dx		x sampling interval
fx		first x sample
vel		array[nx][nz] containing velocities
vel_p		array[nx][nz] containing velocity perturbations

Output:
time		array[nx][nz] containing first-arrival times
angle		array[nx][nz] containing propagation angles
sig  		array[nx][nz] containing sigmas
bet		array[nx][nz] containing betas
******************************************************************************
Notes:
The actual computation of times and sigmas is done in polar coordinates,
with bilinear interpolation used to map to/from rectangular coordinates.

This version takes into account a perturbation in velocity.
******************************************************************************
Revisor:  Zhenyue Liu, Colorado School of Mines, 7/8/92
******************************************************************************/
{
	int ix,iz,ia,ir,na,nr;
	float ss,a,r,da,dr,fa,fr,ex,ez,ea,rmax,rmaxs,
		**s,**sp,**tp,**up,**wp,**ap,**time_p;

	/* shift coordinates so source is at (x=0,z=0) */
	fx -= xs;
	fz -= zs;
	ex = fx+(nx-1)*dx;
	ez = fz+(nz-1)*dz;
	
	/* determine polar coordinate sampling */
	rmaxs = fx*fx+fz*fz;
	rmaxs = MAX(rmaxs,fx*fx+ez*ez);
	rmaxs = MAX(rmaxs,ex*ex+ez*ez);
	rmaxs = MAX(rmaxs,ex*ex+fz*fz);
	rmax = sqrt(rmaxs);
	dr = MIN(ABS(dx),ABS(dz));
	nr = 1+NINT(rmax/dr);
	dr = rmax/(nr-1);
	fr = 0.0;
	if (fx==0.0 && fz==0.0) {
		fa = 0.0;  ea = PI/2.0;
	} else if (fx<0.0 && fz==0.0) {
		fa = -PI/2.0;  ea = PI/2.0;
	} else if (fx==0.0 && fz<0.0) {
		fa = 0.0;  ea = PI;
	} else {
		fa = -PI;  ea = PI;
	}
	da = dr/rmax;
	na = 1+NINT((ea-fa)/da);
	da = (ea-fa)/(na-1);
	if (fa==-PI && ea==PI)
		na = na-1;
	
	/* allocate space */
	time_p = alloc2float(nz,nx);
	s = alloc2float(nz,nx);
	sp = alloc2float(na,nr);
	tp = alloc2float(na,nr);
	up = alloc2float(na,nr);
	wp = alloc2float(na,nr);
	ap = alloc2float(na,nr);
	
	/* compute slownesses */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			s[ix][iz] = 1.0/vel[ix][iz];
	
	/* convert from rectangular to polar coordinates */
	recttopolar(nz,dz,fz,nx,dx,fx,s,na,da,fa,nr,dr,fr,sp);
	
	/* average the slownesses in source region */
	for (ir=0,ss=0.0; ir<2; ++ir)
		for (ia=0; ia<na; ++ia)
			ss += sp[ir][ia];
	ss /= 2*na;

	/* compute traveltimes and derivatives in source region */
	for (ir=0,r=0; ir<2; ++ir,r+=dr) {
		for (ia=0; ia<na; ++ia) {
			up[ir][ia] = ss;
			wp[ir][ia] = 0.0;
			tp[ir][ia] = r*ss;
		}
	}

/* 	tt=cpusec();   */
	/* solve eikonal equation for remaining times and derivatives */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) {
		eikpex(na,da,r,dr,
			sp[ir],up[ir],wp[ir],tp[ir],
			sp[ir+1],up[ir+1],wp[ir+1],tp[ir+1]);
	}
	
	/* convert times from polar to rectangular coordinates */
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,time);

/*  	fprintf(stderr,"\t CPU time for traveltimes= %f \n",cpusec()-tt); 
 	tt=cpusec();   */
	
	/* compute propagation angles in polar and convert */
	for (ia=0,a=fa; ia<na; ++ia,a+=da)
		ap[0][ia] = a;
	for (ir=1,r=fr+dr; ir<nr; ++ir,r+=dr)
		for (ia=0,a=fa; ia<na; ++ia,a+=da){
		    ap[ir][ia] = a+asin(wp[ir][ia]/(sp[ir][ia]*r));
		}
	polartorect(na,da,fa,nr,dr,fr,ap,nz,dz,fz,nx,dx,fx,angle);
/*  	fprintf(stderr,"\t CPU time for propagation angles= %f\n", 	
		cpusec()-tt); 
	tt=cpusec();   */
	
	/* compute sigmas  for initial values */
	for (ir=0,r=0; ir<2; ++ir,r+=dr) 
		for (ia=0; ia<na; ++ia) tp[ir][ia] = r/ss;

	/* solve diffrence equation for remaining sigmas */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) 
 		ray_theoretic_sigma(na,da,r,dr,up[ir],wp[ir],tp[ir],
			up[ir+1],wp[ir+1],tp[ir+1]);  
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,sig);

/* 	fprintf(stderr,"\t CPU time for sigmas= %f \n",cpusec()-tt); 
	tt=cpusec(); */
	
	/* compute betas for initial values */
	for (ir=0; ir<2; ++ir) 
		for (ia=0,a=fa; ia<na; ++ia,a+=da) tp[ir][ia] = a;

	/* solve diffrence equation for remaining betas */
	for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) 
 		ray_theoretic_beta(na,da,r,dr,up[ir],wp[ir],tp[ir],
			up[ir+1],wp[ir+1],tp[ir+1]);  
	polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,bet);
	
	if(pert){
	    /* compute time perturbation by Crank-Nicolson scheme */

	    recttopolar(nz,dz,fz,nx,dx,fx,vel_p,na,da,fa,nr,dr,fr,sp);
	    /* initial values */
	    for (ir=0,r=0; ir<2; ++ir,r+=dr) 
		for (ia=0; ia<na; ++ia) 
			tp[ir][ia] = r/ss*sp[ir][ia];
	    /* solve difference equation for remaining perturbed time */
	    for (ir=1,r=dr; ir<nr-1; ++ir,r+=dr) 
 		delta_t(na,da,r,dr,up[ir],wp[ir],tp[ir],sp[ir],
			up[ir+1],wp[ir+1],tp[ir+1],sp[ir+1]);  

	    polartorect(na,da,fa,nr,dr,fr,tp,nz,dz,fz,nx,dx,fx,time_p);
	    /* sum the total time */
	    for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			time[ix][iz] += time_p[ix][iz];
	}
/* 	fprintf(stderr,"\t CPU time for incident angles= %f \n",
		cpusec()-tt); */
	
	/* free space */
	free2float(s);
	free2float(time_p);
	free2float(sp);
	free2float(tp);
	free2float(up);
	free2float(wp);
	free2float(ap);
}
