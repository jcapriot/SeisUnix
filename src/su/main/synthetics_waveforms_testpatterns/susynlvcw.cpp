/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSYNLVCW: $Revision: 1.18 $ ; $Date: 2015/06/02 20:15:23 $	*/

#include "su.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUSYNLVCW - SYNthetic seismograms for Linear Velocity function	",
" 		for mode Converted Waves				",
"									",
" susynlvcw >outfile [optional parameters]				",
" 									",
" Optional Parameters:							",
" nt=101		number of time samples				",
" dt=0.04		time sampling interval (sec)			",
" ft=0.0		first time (sec)				",
" nxo=1			number of source-receiver offsets		",
" dxo=0.05		offset sampling interval (km)			",
" fxo=0.0		first offset (km, see notes below)		",
" xo=fxo,fxo+dxo,...	array of offsets (use only for non-uniform offsets)",
" nxm=101		number of midpoints (see notes below)		",
" dxm=0.05		midpoint sampling interval (km)		",
" fxm=0.0		first midpoint (km)				",
" nxs=101		number of shotpoints (see notes below)		",
" dxs=0.05		shotpoint sampling interval (km)		",
" fxs=0.0		first shotpoint (km)				",
" x0=0.0		distance x at which v00 is specified		",
" z0=0.0		depth z at which v00 is specified		",
" v00=2.0		velocity at x0,z0 (km/sec)			",
" gamma=1.0		velocity ratio, upgoing/downgoing		",
" dvdx=0.0		derivative of velocity with distance x (dv/dx)	",
" dvdz=0.0		derivative of velocity with depth z (dv/dz)	",
" fpeak=0.2/dt		peak frequency of symmetric Ricker wavelet (Hz)	",
" ref=\"1:1,2;4,2\"	reflector(s):  \"amplitude:x1,z1;x2,z2;x3,z3;...\"",
" smooth=0		=1 for smooth (piecewise cubic spline) reflectors",
" er=0			=1 for exploding reflector amplitudes		",
" ls=0			=1 for line source; default is point source	",
" ob=1			=1 to include obliquity factors		",
" sp=1			=1 to account for amplitude spreading		",
" 			=0 for constant amplitudes throught out		",
" tmin=10.0*dt		minimum time of interest (sec)			",
" ndpfz=5		number of diffractors per Fresnel zone		",
" verbose=0		=1 to print some useful information		",
" 									",
" Notes:								",
" 									",
" Offsets are signed - may be positive or negative.  Receiver locations	",
" are computed by adding the signed offset to the source location.	",
" 									",
" Specify either midpoint sampling or shotpoint sampling, but not both.	",
" If neither is specified, the default is the midpoint sampling above.	",
" 									",
" More than one ref (reflector) may be specified.  When obliquity factors",
" are included, then only the left side of each reflector (as the x,z	",
" reflector coordinates are traversed) is reflecting.  For example, if x",
" coordinates increase, then the top side of a reflector is reflecting.	",
" Note that reflectors are encoded as quoted strings, with an optional	",
" reflector amplitude: preceding the x,z coordinates of each reflector.	",
" Default amplitude is 1.0 if amplitude: part of the string is omitted.	",
" 									",
" Note that gamma<1 implies P-SV mode conversion, gamma>1 implies SV-P,	",
" and gamma=1 implies no mode conversion.				",
" 									",
NULL};
#endif


/*
 * based on Dave Hale's code susynlv, but modified
 * by Mohammed Alfaraj to handle mode conversion
 * Date of modification: 01/07/92
 *
 * Trace header fields set: trid, counit, ns, dt, delrt,
 *				tracl. tracr, fldr, tracf,
 *				cdp, cdpt, d2, f2, offset, sx, gx
 */
/**************** end self doc ***********************************/

/* Library version of SUSYNLVCW
 *
 * The main program (the parameters, the loops over shots or midpoints and offsets, the trace headers) is left to the caller, with
 * the reflectors (decodeReflectors, breakReflectors, makeref), the Ricker wavelet (makericker) and the half-derivative filter
 * (mkhdiff) made by it. What is here is the part that makes one seismogram, and it holds no static state: the half-derivative
 * filter, which the program made on its first trace, is passed in. The sinc table of addsinc is made by su_addsinc_table.
 *
 *	su_synlvcw:	one seismogram
 *
 * Differences from the program: the sp=0 amplitude is that of the program (1/sqrt(time)), and it is not multiplied by the
 * reflector amplitude or the segment length, so that it is only good for looking at the traveltimes.
 */

void su_synlvcw(float *data,
	float xs, float zs, float xg, float zg,
	size_t nt, float dt, float ft,
	float v00, float dvdx, float dvdz, float gamma,
	int ls, int er, int ob, int sp,
	Wavelet *w, int nr, Reflector *r, int lhd, int nhd, float *hd
)
/*****************************************************************************
Make one synthetic seismogram for linear velocity v(x,z) = v00+dvdx*x+dvdz*z, with the mode conversion
******************************************************************************
Input:
v00		velocity v at (x=0,z=0)
dvdx		derivative dv/dx of velocity v with respect to x
dvdz		derivative dv/dz of velocity v with respect to z
gamma		velocity ratio, upgoing/downgoing
ls		=1 for line source amplitudes; =0 for point source
er		=1 for exploding, =0 for normal reflector amplitudes
ob		=1 to include cos obliquity factors; =0 to omit
sp		=1 to account for amplitude spreading
w		wavelet to convolve with trace
xs, zs		source coordinates
xg, zg		receiver group coordinates
nr		number of reflectors
r		array[nr] of reflectors
nt, dt, ft	number of time samples, time sampling interval, first time sample
lhd, nhd, hd	the half-derivative filter, array[nhd], nhd = 1+2*lhd (from mkhdiff)

Output:
data		array[nt] containing synthetic seismogram
*****************************************************************************/
{
	int ir,is,ns;
	size_t it;
	float ar,ds,xd,zd,cd,sd,vs,vg,vd,cs,ss,ts,qs,cg,sg,tg,qg,
		ci,cr,time,amp,*temp,tos;
	ReflectorSegment *rs;

	/* constant depending on gamma and v00*/
	tos = 2.0/(1.0+1.0/gamma);

	/* zero trace */
	for (it=0; it<nt; ++it)
		data[it] = 0.0;
	
	/* velocities at source and receiver */
	vs = v00+dvdx*xs+dvdz*zs;
	vg = (v00+dvdx*xg+dvdz*zg)*gamma;

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
			
			/* velocity at diffractor */
			vd = (v00+dvdx*xd+dvdz*zd)*tos;

			/* ray from shot to diffractor */
			raylv2(v00,dvdx,dvdz,xs,zs,xd,zd,&cs,&ss,&ts,&qs);

			/* ray from receiver to diffractor */
			raylv2(gamma*v00,gamma*dvdx,gamma*dvdz,xg,zg,xd,zd,
			       &cg,&sg,&tg,&qg);

			/* cosines of incidence and reflection angles */
			if (ob) {
				ci = cd*cs+sd*ss;
				cr = cd*cg+sd*sg;
			} else {
				ci = 1.0;
				cr = 1.0;
			}

			/* if either cosine is negative, skip diffractor */
			if (ci<0.0 || cr<0.0) continue;

			/* two-way time and amplitude */
			time = ts+tg;
			if (sp) {
				if (er) {
					amp = sqrt(vg*vd/qg);
				} else {
					if (ls)
						amp = sqrt((vs*vd*vd*vg)/
							(qs*qg));
					else
						amp = sqrt((vs*vd*vd*vg)/
							(qs*qg*(qs+qg)));
				}
				amp *= (ci+cr)*ar*ds;
			} else 
				amp = 1.0/sqrt(time);
				
			/* add sinc wavelet to trace */
			addsinc(time,amp,nt,dt,ft,data);
		}
	}
	
	/* allocate workspace */
	temp = ealloc1float(nt);
	
	/* apply half-derivative filter to trace */
	convolve_cwp(nhd,-lhd,hd,nt,0,data,nt,0,temp);

	/* convolve wavelet with trace */
	convolve_cwp(w->lw,w->iw,w->wv,nt,0,temp,nt,0,data);
	
	/* free workspace */
	free1float(temp);
}
