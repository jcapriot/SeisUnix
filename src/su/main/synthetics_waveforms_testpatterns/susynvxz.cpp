/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSYNVXZ: $Revision: 1.22 $ ; $Date: 2015/06/02 20:15:23 $	*/

#include "su.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUSYNVXZ - SYNthetic seismograms of common offset V(X,Z) media via	",
" 		Kirchhoff-style modeling				",
" 									",
" susynvxz >outfile [optional parameters]				",
" 									",
" Required Parameters:							",
" <vfile		file containing velocities v[nx][nz]		",
" nx=			number of x samples (2nd dimension)		",
" nz=			number of z samples (1st dimension)		",
" Optional Parameters:							",
" nxb=nx		band centered at midpoint			",
" nxd=1			skipped number of midponits			",
" dx=100		x sampling interval (m)				",
" fx=0.0		first x sample					",
" dz=100		z sampling interval (m)				",
" nt=101		number of time samples				",
" dt=0.04		time sampling interval (sec)			",
" ft=0.0		first time (sec)				",
" nxo=1		 	number of offsets				",
" dxo=50		offset sampling interval (m)			",
" fxo=0.0		first offset (m)				",
" nxm=101		number of midpoints				",
" dxm=50		midpoint sampling interval (m)			",
" fxm=0.0		first midpoint (m)				",
" fpeak=0.2/dt		peak frequency of symmetric Ricker wavelet (Hz)	",
" ref=\"1:1,2;4,2\"	reflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
" smooth=0		=1 for smooth (piecewise cubic spline) reflectors",
" ls=0			=1 for line source; default is point source	",
" tmin=10.0*dt		minimum time of interest (sec)			",
" ndpfz=5		number of diffractors per Fresnel zone		",
" verbose=0		=1 to print some useful information		",
" 									",
" Notes:								",
" This algorithm is based on formula (58) in Geo. Pros. 34, 686-703,	",
" by N. Bleistein.							",
" 									",
" Offsets are signed - may be positive or negative.			", 
" Traveltime and amplitude are calculated by finite differences which	",
" is done only in part of midpoints; in the skiped midpoint, interpolation",
" is used to calculate traveltime and amplitude.			", 
"									",
" More than one ref (reflector) may be specified.			",
" Note that reflectors are encoded as quoted strings, with an optional	",
" reflector amplitude: preceding the x,z coordinates of each reflector.	",
" Default amplitude is 1.0 if amplitude: part of the string is omitted.	",
"									",
NULL};
#endif


/*
 *   CWP:  Zhenyue Liu, 07/20/92
 *	Many subroutines borrowed from Dave Hale's program: SUSYNLV
 *
 * Trace header fields set: trid, counit, ns, dt, delrt,
 *				tracl. tracr,
 *				cdp, cdpt, d2, f2, offset, sx, gx
 */
/**************** end self doc ***********************************/

/* Library version of SUSYNVXZ
 *
 * The main program (the parameters, the velocities read from a file, the headers of the traces) is left to the caller, with the
 * reflectors (decodeReflectors, breakReflectors, makeref), the Ricker wavelet (makericker) and the half-derivative filter (mkhdiff)
 * made by it. What is here is the loop over the midpoints of one offset, which makes the traces of the offset, and it holds no
 * static state (the half-derivative filter, which the program made on its first trace, is passed in).
 *
 *	su_synvxz_offset:	the nxm traces of the common-offset section for one offset
 *
 * Differences from the program: the checks of the parameters return a code, and the maximum segment length of the reflectors (dsmax),
 * which the program works out from the velocity at the first sample, is for the caller to do.
 */

/* parameters for half-derivative filter */
#define LHD 20
#define NHD 1+2*LHD

static void makeone (float **ts, float **as, float **sgs, float **tg,
	float **ag, float **sgg, float ex, float ez, float dx,
	float dz, float fx, float vs0, float vg0, int ls, Wavelet *w,
	int nr, Reflector *r, int nt, float dt, float ft, float *trace,
	int lhd, int nhd, float *hd);

/* The traces of one offset xo, for the nxm midpoints fxm+ixm*dxm; the velocity is vel[nx][nz] (z the fast axis) on the grid
 * fx+ix*dx, iz*dz. The traces are out[nxm][nt].
 *
 * nxb: half the width, in grid samples, of the band of the traveltimes about the midpoint; nxd: the number of midpoints between
 * those at which the traveltimes are calculated (1 for every one).
 *
 * Returns 0, -1 if a shot or receiver is outside the grid, -2 if nxb is too small for the offset, -3 for parameters that are
 * not allowed.
 */
int su_synvxz_offset(const float *vel, int nx, int nz, float dx, float dz, float fx,
	int nxb, int nxd, float xo, int nxm, float dxm, float fxm,
	int ls, Wavelet *w, int nr, Reflector *r, int nt, float dt, float ft,
	int lhd, int nhd, float *hd, float *out)
{
	int ixm,ixd,i,ix,iz,nx1,nx0,nxd1;
	float temp,temp1,fx1,ex1,ex,ez,xm,vs0,vg0,xs,xg,
		**vold,**ts,**as,**sgs,**tg,**ag,**sgg,**bas,**bag,
		**v,**ts1=NULL,**as1=NULL,**sgs1=NULL,
		**tg1=NULL,**ag1=NULL,**sgg1=NULL;

	if (nx<2 || nz<2 || nxb<1 || nxd<1 || nxm<1 || nt<1) return -3;

	/* check the ranges of shots and receivers */
	ex = fx+(nx-1)*dx;
	ez = (nz-1)*dz;
	for (ixm=0; ixm<nxm; ++ixm) {
		xs = fxm+ixm*dxm-0.5*xo;
		xg = xs+xo;
		if (fx>xs || ex<xs || fx>xg || ex<xg) return -1;
	}
	if (ABS(xo)>nxb*dx) return -2;

	/* the velocities, with the pointers to the rows */
	vold = (float**)ealloc1(nx,sizeof(float*));
	for (ix=0; ix<nx; ++ix) vold[ix] = (float*)vel+(size_t)ix*nz;

	/* allocate space */
	nx1 = 1+2*nxb;
	ts = ealloc2float(nz,nx1);
	as = ealloc2float(nz,nx1);
	sgs = ealloc2float(nz,nx1);
	tg = ealloc2float(nz,nx1);
	ag = ealloc2float(nz,nx1);
	sgg = ealloc2float(nz,nx1);
	v = ealloc2float(nz,nx1);
	bas = ealloc2float(nz,nx1);
	bag = ealloc2float(nz,nx1);
	if(nxd>1) {
		/* allocate space for interpolation */
		ts1 = ealloc2float(nz,nx1);
		as1 = ealloc2float(nz,nx1);
		sgs1 = ealloc2float(nz,nx1);
		tg1 = ealloc2float(nz,nx1);
		ag1 = ealloc2float(nz,nx1);
		sgg1 = ealloc2float(nz,nx1);
  	}

	nxd1 = nxd;
	for (ixm=0; ixm<nxm; ixm +=nxd1){
		xm = fxm+ixm*dxm;
   		xs = xm-0.5*xo;
		xg = xs+xo;
		/* set range for traveltimes' calculation */
		fx1 = xm-nxb*dx;
		ex1 = MIN(ex+(nxd1-1)*dxm,xm+nxb*dx);
		nx1 = 1+(ex1-fx1)/dx;
		nx0 = (fx1-fx)/dx;
		temp = (fx1-fx)/dx-nx0;
		/* transpose velocity such that the first row is at fx1 */
		for(ix=0;ix<nx1;++ix)
		    for(iz=0;iz<nz;++iz){
		    	if(ix<-nx0)
			   	v[ix][iz] = vold[0][iz];
		    	else if(ix+nx0>nx-2)
				v[ix][iz]=vold[nx-1][iz];
			else
				v[ix][iz] = vold[ix+nx0][iz]*(1.0-temp)
					+temp*vold[ix+nx0+1][iz];
		    }

		if(ixm==0 || nxd1==1){
		/* No interpolation */

			/* compute traveltimes, propagation angles, sigmas
	  		 from shot and receiver respectively	*/
			eiktam(xs,0.,nz,dz,0.,nx1,dx,fx1,v,ts,as,sgs,bas);
			eiktam(xg,0.,nz,dz,0.,nx1,dx,fx1,v,tg,ag,sgg,bag);
			ixd = NINT((xs-fx)/dx);
			vs0 = vold[ixd][0];
			ixd = NINT((xg-fx)/dx);
			vg0 = vold[ixd][0];

			/* make one trace */
			ex1 = MIN(ex,xm+nxb*dx);
			makeone(ts,as,sgs,tg,ag,sgg,ex1,ez,dx,dz,fx1,vs0,vg0,
				ls,w,nr,r,nt,dt,ft,out+(size_t)ixm*nt,lhd,nhd,hd);
		}
		else {
			/* Linear interpolation */

			eiktam(xs,0,nz,dz,0,nx1,dx,fx1,v,ts1,as1,sgs1,bas);
			eiktam(xg,0,nz,dz,0,nx1,dx,fx1,v,tg1,ag1,sgg1,bag);
			ixd = NINT((xs-fx)/dx);
			vs0 = vold[ixd][0];
			ixd = NINT((xg-fx)/dx);
			vg0 = vold[ixd][0];

		    	xm -= nxd1*dxm;
		    for(i=1; i<=nxd1; ++i) {
		    	xm += dxm;
			xs = xm-0.5*xo;
			xg = xs+xo;
			fx1 = xm-nxb*dx;
			ex1 = MIN(ex+(nxd1-1)*dxm,xm+nxb*dx);
			nx1 = 1+(ex1-fx1)/dx;
			temp = nxd1-i;
			temp1 = 1.0/(nxd1-i+1);
			for(ix=0;ix<nx1;++ix)
			    for(iz=0;iz<nz;++iz){
			    if(i==nxd1){
			   	ts[ix][iz] = ts1[ix][iz];
			   	tg[ix][iz] = tg1[ix][iz];
			   	sgs[ix][iz] = sgs1[ix][iz];
			   	ag[ix][iz] = ag1[ix][iz];
			   	sgg[ix][iz] = sgg1[ix][iz];
			   	as[ix][iz] = as1[ix][iz];
				}
			    else{
			   	ts[ix][iz] = (temp*ts[ix][iz]
					+ts1[ix][iz])*temp1;
			   	tg[ix][iz] = (temp*tg[ix][iz]
					+tg1[ix][iz])*temp1;
			    	as[ix][iz] = (temp*as[ix][iz]
					+as1[ix][iz])*temp1;
			   	sgs[ix][iz] = (temp*sgs[ix][iz]
					+sgs1[ix][iz])*temp1;
			   	ag[ix][iz] = (temp*ag[ix][iz]
					+ag1[ix][iz])*temp1;
			   	sgg[ix][iz] = (temp*sgg[ix][iz]
					+sgg1[ix][iz])*temp1;
				}
			}

			/* make one trace */
			ex1 = MIN(ex,xm+nxb*dx);
			makeone(ts,as,sgs,tg,ag,sgg,ex1,ez,dx,dz,fx1,vs0,vg0,
				ls,w,nr,r,nt,dt,ft,out+(size_t)(ixm-nxd1+i)*nt,lhd,nhd,hd);
		    }
		}
		/* set skip parameter */
		if(ixm<nxm-1 && ixm>nxm-1-nxd1) nxd1 = nxm-1-ixm;
	}

	free1(vold);
	free2float(ts);
	free2float(bas);
	free2float(sgs);
	free2float(as);
	free2float(v);
	free2float(tg);
	free2float(bag);
	free2float(sgg);
	free2float(ag);
	if(nxd>1) {
		free2float(ts1);
		free2float(as1);
		free2float(sgs1);
		free2float(tg1);
		free2float(ag1);
		free2float(sgg1);
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
**v		array[nx][nz] containing velocities 
nz		number of z samples
dz		z sampling interval
nx		number of x samples
dx		x sampling interval
fx		first x sample
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
