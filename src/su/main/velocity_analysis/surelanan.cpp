/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SURELANAN - REsiduaL-moveout semblance ANalysis for ANisotropic media	",
"									",
" surelan refl= npicks=    [optional parameters]			",
" 									",
" Required parameters:							",
" reflector file: reflec =						",
" number of points in the reflector file =				",
" 									",
" Optional Parameters:							",
" nr1=51		number of r1-parameter samples   		",
" dr1=0.01              r1-parameter sampling interval			",
" fr1=-0.25             first value of r1-parameter			",
" nr2=51		number of r2-parameter samples   		",
" dr2=0.01              r2-parameter sampling interval			",
" fr2=-0.25             first value of r2-parameter			",
" dzratio=5             ratio of output to input depth sampling intervals",
" nsmooth=dzratio*2+1   length of semblance num and den smoothing window",
" verbose=0             =1 for diagnostic print on stderr		",
" method=linear		for linear interpolation of the interface       ",
" 			=mono for monotonic cubic interpolation of interface",
" 			=akima for Akima's cubic interpolation of interface ",
" 			=spline for cubic spline interpolation of interface ",
"									",
" Note: 								",
" 1. This program is part of Debashish Sarkar's anisotropic model building",
" technique. 								",
" 2. Input migrated traces should be sorted by cdp - surelan outputs a 	",
"    group of semblance traces every time cdp changes.  Therefore, the  ",
"    output will be useful only if cdp gathers are input.  		",
" 3. The residual-moveout semblance for cdp gathers is based		",
"	on z(h)*z(h) = z(0)*z(0) + r1*h^2 + r2*h^4/[h^2+z(0)^2] where z ",
"	depth and h is the half-offset.   				",
"									",
NULL};
#endif
/**************** end self doc *******************************************/

/* 
 * Trace header fields accessed:  ns, d1, f1, offset, cdp
 * Trace header fields modified:  ns, d1, offset, cdp
 */

/* What the program does with a trace: residual moveout of the trace for each pair of parameters (r1, r2), over the
 * nsmooth samples around the depth Z of the reflector, added to the sums.
 *
 *	z(h)^2 = z(0)^2 + r1 h^2 + r2 h^4 / (h^2 + z(0)^2),	h the half offset
 *
 * nz, dz, fz	samples, sample interval and first depth of the trace; data[nz]
 * half_offset	half the offset (the program has the integer offset divided by 2: half of an odd one is lost)
 * Z		the depth of the reflector at the CDP of the gather (the program works it out once, for the CDP of
 *		the first trace, and uses it for all of them)
 * num, den, nnz	[nr1][nr2][nsmooth] sums, to which this adds
 */
void su_relanan_accumulate(int nz, float dz, float fz, float half_offset, float Z, int nr1, float dr1, float fr1,
	int nr2, float dr2, float fr2, int nsmooth, const float *data, float *num, float *den, float *nnz)
{
	int ir1, ir2, iz, izi;
	float r1, r2, roffs2, Zmin, zn, zi, frac, temp;

	for (ir1=0,r1=fr1; ir1<nr1; ++ir1,r1+=dr1)
	for (ir2=0,r2=fr2; ir2<nr2; ++ir2,r2+=dr2) {
		size_t base = ((size_t) ir1*nr2 + ir2)*nsmooth;

		/* compute hyperbolic and nonhyperbolic terms */
		roffs2 = r1*half_offset*half_offset+r2*pow(half_offset,4)/(pow(half_offset,2)+pow(Z,2));

		/* the first depth of the window */
		Zmin=Z-((nsmooth-1)*dz/2);

		/* do rmo via linear interpolation and accumulate semblance numerator and denominator */
		for (iz=0,zn=Zmin; iz<nsmooth; ++iz,zn+=dz) {
			temp = zn*zn+roffs2;
			zi = (temp>fz*fz)?(sqrt(temp)-fz)/dz:0;
			izi = zi;
			if (izi>=0 && izi<nz-1) {
				frac = zi-izi;
				temp = (1.0-frac)*data[izi]+frac*data[izi+1];
				if (temp!=0.0) {
					num[base+iz] += temp;
					den[base+iz] += temp*temp;
					nnz[base+iz] += 1.0;
				}
			}
		}
	}
}
