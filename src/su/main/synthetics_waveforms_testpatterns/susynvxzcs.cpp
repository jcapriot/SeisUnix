/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSYNVXZCS: $Revision: 1.17 $ ; $Date: 2015/06/02 20:15:23 $	*/

#include "su.h"
/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUSYNVXZCS - SYNthetic seismograms of common shot in V(X,Z) media via	",
" 		Kirchhoff-style modeling				",
" 									",
" susynvxzcs<vfile >outfile  nx= nz= [optional parameters]		",
" 									",
" Required Parameters:							",
" <vfile        file containing velocities v[nx][nz]			",
" >outfile      file containing seismograms of common ofset		",
" nx=           number of x samples (2nd dimension) in velocity ",
" nz=           number of z samples (1st dimension) in velocity ",
" 									",
" Optional Parameters:							",
" nt=501        	number of time samples				",
" dt=0.004      	time sampling interval (sec)			",
" ft=0.0        	first time (sec)				",
" fpeak=0.2/dt		peak frequency of symmetric Ricker wavelet (Hz)	",
" nxg=			number of receivers of input traces		",
" dxg=15		receiver sampling interval (m)			",
" fxg=0.0		first receiver (m)				",
" nxd=5         	skipped number of receivers			",
" nxs=1			number of offsets				",
" dxs=50		shot sampling interval (m)			",
" fxs=0.0		first shot (m)				",
" dx=50         	x sampling interval (m)				",
" fx=0.         	first x sample (m)				",
" dz=50         	z sampling interval (m)				",
" nxb=nx/2    	band width centered at midpoint (see note)	",
" nxc=0         hozizontal range in which velocity is changed	",
" nzc=0         vertical range in which velocity is changed	",
" pert=0        =1 calculate time correction from v_p[nx][nz]	",
" vpfile        file containing slowness perturbation array v_p[nx][nz]	",
" ref=\"1:1,2;4,2\"	reflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
" smooth=0		=1 for smooth (piecewise cubic spline) reflectors",
" ls=0			=1 for line source; =0 for point source		",
" tmin=10.0*dt		minimum time of interest (sec)			",
" ndpfz=5		number of diffractors per Fresnel zone		",
" cable=1		roll reciever spread with shot			",
" 			=0 static reciever spread			",
" verbose=0		=1 to print some useful information		",
" 									",
" Notes:								",
" This algorithm is based on formula (58) in Geo. Pros. 34, 686-703,	",
" by N. Bleistein.							",
"									",
" Traveltime and amplitude are calculated by finite difference which	",
" is done only in one of every NXD receivers; in skipped receivers, 	",
" interpolation is used to calculate traveltime and amplitude.		", 
" For each receiver, traveltime and amplitude are calculated in the 	",
" horizontal range of (xg-nxb*dx, xg+nxb*dx). Velocity is changed by 	",
" constant extropolation in two upper trianglar corners whose width is 	",
" nxc*dx and height is nzc*dz.						",
"									",
" Eikonal equation will fail to solve if there is a polar turned ray.	",
" In this case, the program shows the related geometric information. 	",
" There are three ways to remove the turned rays: smoothing velocity, 	",
" reducing nxb, and increaing nxc and nzc (if the turned ray occurs  	",
" in shallow areas). To prevent traveltime distortion from an over-	",
" smoothed velocity, traveltime is corrected based on the slowness 	",
" perturbation.								",
"									",
" More than one ref (reflector) may be specified.			",
" Note that reflectors are encoded as quoted strings, with an optional	",
" reflector amplitude: preceding the x,z coordinates of each reflector.	",
" Default amplitude is 1.0 if amplitude: part of the string is omitted.	",
"									",
NULL};
#endif


/*
 *	Author: Zhenyue Liu, 07/20/92, Center for Wave Phenomena
 *		Many subroutines borrowed from Dave Hale's program: SUSYNLV
 *
 *		Trino Salinas, 07/30/96, fixed a bug in the geometry
 *		setting to allow the spread move with the shots.
 *
 *		Chris Liner 12/10/08  added cable option, set fldr header word
 *
 * Trace header fields set: trid, counit, ns, dt, delrt,
 *				tracl. tracr, fldr, tracf,
 *				sx, gx
 */
/**************** end self doc ***********************************/

/* Library version of SUSYNVXZCS
 *
 * The main program (the parameters, the velocities and the slowness perturbations read from files, the headers of the traces) is
 * left to the caller, with the reflectors (decodeReflectors, breakReflectors, makeref), the Ricker wavelet (makericker) and the
 * half-derivative filter (mkhdiff) made by it. What is here is the loop over the receivers of one shot, which makes the traces of
 * the shot, and it holds no static state (the half-derivative filter, which the program made on its first trace, is passed in).
 *
 *	su_synvxzcs_shot:	the nxg traces of the common-shot gather of one shot
 *
 * Differences from the program:
 *  - The maximum segment length of the reflectors (dsmax), which the program works out from the velocity at the first sample, is for
 *    the caller to do.
 *  - The program has global variables (ierr, x_err, ...) with which it would report that the eikonal equation failed with a turned
 *    ray, but nothing in the library sets them, so that it never did; they are gone.
 *  - The velocities in the two upper corners (nxc, nzc) are replaced in the program with the width of the traveltime band
 *    (1+2*nxb) where the width of the model (nx) is meant, which can go past the last trace of the model; it is nx here. The model
 *    of the caller is not changed.
 */

/* parameters for half-derivative filter */
#define LHD 20
#define NHD 1+2*LHD

static void makeone (float **ts, float **as, float **sgs, float **tg,
	float **ag, float **sgg, float ex, float ez, float dx,
	float dz, float fx, float vs0, float vg0, int ls, Wavelet *w,
	int nr, Reflector *r, int nt, float dt, float ft, float *trace,
	int lhd, int nhd, float *hd);

static void delta_t (int na, float da, float r, float dr,
	float uc[], float wc[], float sc[], float bc[],
	float un[], float wn[], float sn[], float bn[]);
static void eiktam_pert (float xs, float zs, int nz, float dz, float fz, int nx,
	float dx, float fx, float **vel, float **vel_p, int pert,
	float **time, float **angle, float **sig, float **bet);

/* The traces of one shot at xs, for the receivers gx0+ixg*dxg; the velocity is vel[nx][nz] (z the fast axis) on the grid fx+ix*dx,
 * iz*dz, and pert (or NULL) the slowness perturbation, array[nx][nz]. The traces are out[nxg][nt].
 *
 * nxb: half the width, in grid samples, of the band of the traveltimes about the receiver; nxd: the number of receivers between
 * those at which the traveltimes are calculated (1 for every one); nxc, nzc: the width and height, in grid samples, of the two
 * upper corners where the velocity is replaced by constant extrapolation.
 *
 * Returns 0, -1 if the shot or a receiver is outside the grid, -3 for parameters that are not allowed.
 */
int su_synvxzcs_shot(const float *vel, const float *pert, int nx, int nz, float dx, float dz, float fx,
	int nxb, int nxd, int nxc, int nzc, float xs, int nxg, float dxg, float gx0,
	int ls, Wavelet *w, int nr, Reflector *r, int nt, float dt, float ft,
	int lhd, int nhd, float *hd, float *out)
{
	int ixg,ixd,i,ix,iz,nx1,nx0,nxd1,doperturb = (pert!=NULL);
	float temp,temp1,fx1,ex1,ex,ez,xg,vs0,vg0,
		**vold,**vpold=NULL,**told,**aold,**sgold,
		**ts,**as,**sgs,**tg,**ag,**sgg,**bas,**bag,
		**v,**vp=NULL,**tg1=NULL,**ag1=NULL,**sgg1=NULL;

	if (nx<2 || nz<2 || nxb<1 || nxd<1 || nxg<1 || nt<1 || nxc<0 || nzc<0) return -3;
	if (nxc>nxb) nxc = nxb;
	if (nzc>nz) nzc = nz;

	/* check the ranges of the shot and receivers */
	ex = fx+(nx-1)*dx;
	ez = (nz-1)*dz;
	if (fx>xs || ex<xs) return -1;
	for (ixg=0; ixg<nxg; ++ixg) {
		xg = gx0+ixg*dxg;
		if (fx>xg || ex<xg) return -1;
	}

	/* allocate space, and copy the velocities, which are changed below */
	vold = ealloc2float(nz,nx);
	aold = ealloc2float(nz,nx);
	told = ealloc2float(nz,nx);
	bas = ealloc2float(nz,nx);
	sgold = ealloc2float(nz,nx);
	memcpy(vold[0],vel,sizeof(float)*(size_t)nx*nz);
	if (doperturb) {
		vpold = ealloc2float(nz,nx);
		memcpy(vpold[0],pert,sizeof(float)*(size_t)nx*nz);
		/* calculate 1/2 perturbation of slowness squares*/
		for (ix=0; ix<nx; ++ix)
			for (iz=0; iz<nz; ++iz)
				vpold[ix][iz] = vpold[ix][iz]/vold[ix][iz];
	}

	nx1 = 1+2*nxb;
	ts = ealloc2float(nz,nx1);
	as = ealloc2float(nz,nx1);
	sgs = ealloc2float(nz,nx1);
	tg = ealloc2float(nz,nx1);
	ag = ealloc2float(nz,nx1);
	sgg = ealloc2float(nz,nx1);
	v = ealloc2float(nz,nx1);
	bag = ealloc2float(nz,nx1);
	if (doperturb) vp = ealloc2float(nz,nx1);
	if (nxd>1) {
		/* allocate space for interpolation */
		tg1 = ealloc2float(nz,nx1);
		ag1 = ealloc2float(nz,nx1);
		sgg1 = ealloc2float(nz,nx1);
  	}

	/* change velocity values within two trianglar upper corners to
 		avoid possible polar-turned rays	*/
	for(iz=0; iz<nzc; ++iz)
		for(ix=iz*nxc/nzc; ix<nxc; ++ix){
			vold[nxc-ix-1][iz] = vold[nxc-iz*nxc/nzc][iz];
			vold[nx-nxc+ix][iz] = vold[nx-nxc+iz*nxc/nzc-1][iz];
		}

	/* calculate traveltimes from shot */
	eiktam_pert(xs,0,nz,dz,0,nx,dx,fx,vold,vpold,doperturb,told,aold,sgold,bas);

	/* calculate velocity at shot 	*/
	temp = (xs-fx)/dx;
	ixd = temp;
	temp = temp-ixd;
	vs0 = (ixd+1<nx) ? (1-temp)*vold[ixd][0]+temp*vold[ixd+1][0] : vold[ixd][0];

	/* save the skipping number */
	nxd1 = nxd;

	/* loop over receivers */
	for (ixg=0; ixg<nxg; ixg +=nxd1){
		xg = gx0+ixg*dxg;

		/* set range for traveltimes' calculation */
		fx1 = xg-nxb*dx;
		nx0 = (fx1-fx)/dx;
		temp = (fx1-fx)/dx-nx0;
	       	for(ix=0;ix<nx1;++ix){
		    if(ix<-nx0)
			for(iz=0;iz<nz;++iz)
				v[ix][iz] = vold[0][iz];
		    else if(ix+nx0>nx-2)
			for(iz=0;iz<nz;++iz)
				v[ix][iz]=vold[nx-1][iz];
		    else
			for(iz=0;iz<nz;++iz){
				v[ix][iz] = vold[ix+nx0][iz]*(1.0-temp)
					+temp*vold[ix+nx0+1][iz];
				ts[ix][iz] = told[ix+nx0][iz]*(1.0-temp)
					+temp*told[ix+nx0+1][iz];
				as[ix][iz] = aold[ix+nx0][iz]*(1.0-temp)
					+temp*aold[ix+nx0+1][iz];
				sgs[ix][iz] = sgold[ix+nx0][iz]*(1.0-temp)
					+temp*sgold[ix+nx0+1][iz];
				}
		}

		if(doperturb){
		    for(ix=0;ix<nx1;++ix)
			for(iz=0;iz<nz;++iz){
			    if(ix<-nx0)
				vp[ix][iz] = vpold[0][iz];
			    else if(ix+nx0>nx-2)
				vp[ix][iz]=vpold[nx-1][iz];
			    else
				vp[ix][iz] = vpold[ix+nx0][iz]*(1.0-temp)
					+temp*vpold[ix+nx0+1][iz];
			    }
		}

		if(ixg==0 || nxd1==1){
		/* No interpolation */

			/* compute traveltimes, prop angles, sigmas from receiver */
			eiktam_pert(xg,0.,nz,dz,0.,nx1,dx,fx1,v,vp,doperturb,tg,ag,sgg,bag);

			/* calculate velocity at receiver 	*/
			vg0 = (1-temp)*vold[nx0+nxb][0]
				+temp*vold[MIN(nx0+nxb+1,nx-1)][0];

			/* make one trace */
			ex1 = MIN(ex,xg+nxb*dx);
			makeone(ts,as,sgs,tg,ag,sgg,ex1,ez,dx,dz,fx1,vs0,vg0,
				ls,w,nr,r,nt,dt,ft,out+(size_t)ixg*nt,lhd,nhd,hd);
		}
		else {
			/* Linear interpolation */

			eiktam_pert(xg,0,nz,dz,0,nx1,dx,fx1,v,vp,doperturb,tg1,ag1,sgg1,bag);

			/* calculate velocity at receiver 	*/
			vg0 = (1-temp)*vold[nx0+nxb][0]+temp*vold[MIN(nx0+nxb+1,nx-1)][0];

			/* interpolate quantities between two midpoints	*/
		    	xg -= nxd1*dxg;
		    for(i=1; i<=nxd1; ++i) {
		    	xg += dxg;
			fx1 = xg-nxb*dx;
			ex1 = xg+nxb*dx;
			temp = nxd1-i;
			temp1 = 1.0/(nxd1-i+1);
			for(ix=0;ix<nx1;++ix)
			    for(iz=0;iz<nz;++iz){
			    if(i==nxd1){
			   	tg[ix][iz] = tg1[ix][iz];
			   	ag[ix][iz] = ag1[ix][iz];
			   	sgg[ix][iz] = sgg1[ix][iz];
				}
			    else{
			   	tg[ix][iz] = (temp*tg[ix][iz]
					+tg1[ix][iz])*temp1;
			   	ag[ix][iz] = (temp*ag[ix][iz]
					+ag1[ix][iz])*temp1;
			   	sgg[ix][iz] = (temp*sgg[ix][iz]
					+sgg1[ix][iz])*temp1;
				}
			}

			/* calculate quantities at shot 	*/
	    		nx0 = (fx1-fx)/dx;
	    		temp = (fx1-fx)/dx-nx0;
	   		for(ix=MAX(0,-nx0); ix<MIN(nx1,nx-1-nx0); ++ix){
			    for(iz=0;iz<nz;++iz){
			   	ts[ix][iz] = told[ix+nx0][iz]*(1.0-temp)
					+temp*told[ix+nx0+1][iz];
			    	sgs[ix][iz] = sgold[ix+nx0][iz]*(1.0-temp)
					+temp*sgold[ix+nx0+1][iz];
			        as[ix][iz] = aold[ix+nx0][iz]*(1.0-temp)
					+temp*aold[ix+nx0+1][iz];
		   		 }
			}

			/* make one trace */
			ex1 = MIN(ex,xg+nxb*dx);
			makeone(ts,as,sgs,tg,ag,sgg,ex1,ez,dx,dz,fx1,vs0,vg0,
				ls,w,nr,r,nt,dt,ft,out+(size_t)(ixg-nxd1+i)*nt,lhd,nhd,hd);
		    }
		}
		/* recount the skipping number for the last receiver */
		if(ixg<nxg-1 && ixg>nxg-1-nxd1) nxd1 = nxg-1-ixg;
	}

	free2float(vold);
	free2float(aold);
	free2float(told);
	free2float(sgold);
	free2float(ts);
	free2float(as);
	free2float(sgs);
	free2float(bas);
	free2float(v);
	free2float(tg);
	free2float(bag);
	free2float(sgg);
	free2float(ag);
	if(doperturb) {
		free2float(vp);
		free2float(vpold);
	}
	if(nxd>1) {
		free2float(sgg1);
		free2float(tg1);
		free2float(ag1);
  	}
	return 0;
}

static void makeone (float **ts, float **as, float **sgs, 
	float **tg, float **ag, float **sgg, float ex, float ez, float dx, 
	float dz, float fx, float vs0, float vg0, int ls, Wavelet *w,
	int nr, Reflector *r, int nt, float dt, float ft, float *trace,
	int lhd, int nhd, float *hd)
/*****************************************************************************
Make one synthetic seismogram 
******************************************************************************
Input:
**ts		array[nx][nz] containing traveltimes from shot
**as		array[nx][nz] containing propagation angles from shot 
**sgs		array[nx][nz] containing sigmas from shot
**tg		array[nx][nz] containing traveltimes from receiver
**ag		array[nx][nz] containing propagation angles from receiver 
**sgg		array[nx][nz] containing sigmas from receiver
nz		number of z samples
dz		z sampling interval
nx		number of x samples
dx		x sampling interval
fx		first x sample
ex		last x sample
ez		last z sample
ls		=1 for line source amplitudes; =0 for point source
w		wavelet to convolve with trace
xs		x coordinate of source
xg		x coordinate of receiver group
nr		number of reflectors
r		array[nr] of reflectors
nt		number of time samples
dt		time sampling interval
ft		first time sample

Output:
trace		array[nt] containing synthetic seismogram
*****************************************************************************/
{
	int it,ir,is,ns,ix,iz;
	float ar,ds,xd,zd,cd,sd,xi,zi,ci,cr,time,amp,sx,sz,
		tsd,asd,sgsd,tgd,agd,sggd,
		*temp;
	ReflectorSegment *rs;
	/* zero trace */
	for (it=0; it<nt; ++it)
		trace[it] = 0.0;
	
	/* loop over reflectors */
	for (ir=0; ir<nr; ++ir) {

		/* amplitude, number of segments, segment length */
		ar = r[ir].a;
		ns = r[ir].ns;
		ds = r[ir].ds;
		rs = r[ir].rs;
	
		/* loop over diffracting segments */
		for (is=0; is<ns; ++is) {
		
			/* diffractor midpoint, unit-normal, and length */
			xd = rs[is].x;
			zd = rs[is].z;
			cd = rs[is].c;
			sd = rs[is].s;
			
			/* check range of reflector */
			if(xd<fx || xd>=ex || zd>=ez)
				continue;
			/* determine sample indices */
			xi = (xd-fx)/dx;
			ix = xi;
			zi = zd/dz;
			iz = zi;
			/* bilinear interpolation */
			sx = xi-ix;
			sz = zi-iz;
			tsd = (1.0-sz)*((1.0-sx)*ts[ix][iz] + 
						sx*ts[ix+1][iz]) +
					sz*((1.0-sx)*ts[ix][iz+1] +
						sx*ts[ix+1][iz+1]);
			asd = (1.0-sz)*((1.0-sx)*as[ix][iz] + 
						sx*as[ix+1][iz]) +
					sz*((1.0-sx)*as[ix][iz+1] +
						sx*as[ix+1][iz+1]);
			sgsd = (1.0-sz)*((1.0-sx)*sgs[ix][iz] + 
						sx*sgs[ix+1][iz]) +
					sz*((1.0-sx)*sgs[ix][iz+1] +
						sx*sgs[ix+1][iz+1]);
			tgd = (1.0-sz)*((1.0-sx)*tg[ix][iz] + 
						sx*tg[ix+1][iz]) +
					sz*((1.0-sx)*tg[ix][iz+1] +
						sx*tg[ix+1][iz+1]);
			agd = (1.0-sz)*((1.0-sx)*ag[ix][iz] + 
						sx*ag[ix+1][iz]) +
					sz*((1.0-sx)*ag[ix][iz+1] +
						sx*ag[ix+1][iz+1]);
			sggd = (1.0-sz)*((1.0-sx)*sgg[ix][iz] + 
						sx*sgg[ix+1][iz]) +
					sz*((1.0-sx)*sgg[ix][iz+1] +
						sx*sgg[ix+1][iz+1]);
			
			/* cosines of incidence and reflection angles */
			ci = cd*cos(asd)+sd*sin(asd);
			cr = cd*cos(agd)+sd*sin(agd);

			/* two-way time and amplitude */
			time = tsd+tgd;

			if (ls)
			     amp = sqrt(vs0*vg0/(sgsd*sggd));
			else
			     amp = sqrt(vs0*vg0/(sgsd*sggd*(sgsd+sggd)));
					   
			amp *= ABS(ci+cr)*ar*ds;
		
			/* add sinc wavelet to trace */
			addsinc(time,amp,nt,dt,ft,trace);
		}
	}
	
	/* allocate workspace */
	temp = ealloc1float(nt);
	
	/* apply half-derivative filter to trace */
	convolve_cwp(nhd,-lhd,hd,nt,0,trace,nt,0,temp);

	/* convolve wavelet with trace */
	convolve_cwp(w->lw,w->iw,w->wv,nt,0,temp,nt,0,trace);
	
	/* free workspace */
	free1float(temp);
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
