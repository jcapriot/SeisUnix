/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVEL2DF: $Revision: 1.4 $ ; $Date: 2011/11/16 23:40:27 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUVEL2DF - compute stacking VELocity semblance for a single time in   ",
"			    over Vnmo and eta in 2-D			",
"									",
"    suvel2df <stdin >stdout [optional parameters]			",
"									",
" Required Parameters:							",
" tn			zero-offset time of reflection			",
" offsetm		Maximum offset considerd			",
"									",
" Optional Parameters:							",
" nv=50			number of velocities				",
" dv=50.0		velocity sampling interval			",
" fv=1500.0		first velocity					",
" nvh=50		number of horizotal velocities			",
" dvh=50.0		horizontal velocity sampling interval		",
" fvh=1500.0		first horizontal velocity			",
" xod=1.5		maximum offset-to-depth ratio to resolve	",
" dtratio=5		ratio of output to input time sampling intervals",
" nsmooth=dtratio*2+1	length of semblance num and den smoothing window",
" verbose=0		=1 for diagnostic print on stderr		",
" vavg=fv+0.5*(nv-1)*dv   average velocity used in the search		",
"									",
" Notes:								",
" Semblance is defined by the following quotient:			",
"									",
"		 n-1		 					",
"		[ sum q(t,j) ]^2					",
"		 j=0		 					",
"	s(t) = ------------------					",
"		 n-1		 					",
"		n sum [q(t,j)]^2					",
"		 j=0		 					",
"									",
" where n is the number of non-zero samples after muting.		",
" Smoothing (nsmooth) is applied separately to the numerator and denominator",
" before computing this semblance quotient.				",
"									",
" Input traces should be sorted by cdp - suvel2df outputs a group of	",
" semblance traces every time cdp changes.  Therefore, the output will	",
" be useful only if cdp gathers are input.				",
"									",
NULL};
#endif

/* Credits:
 *	CWP: Tariq Alkhalifah,  February 1997
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp.
 * Trace header fields modified:  ns, dt, offset.
 */
/**************** end self doc *******************************************/

/* What the program does with a trace: NMO by linear interpolation for each pair of velocities (vnmo, vh) over the
 * window of nsmooth samples around the time tn, added to the sums.
 *
 * nt, dt, ft	samples, sample interval and first time of the trace; data[nt]
 * offset	offset of the trace
 * tn_zero	the zero offset time of the reflection
 * nv, dv, fv	the number, increment and first of the NMO velocities
 * nvh, dvh, fvh	the same for the horizontal velocities
 * nsmooth	the length of the window
 * num, den, nnz	[nv][nvh][nsmooth] sums, to which this adds
 */
void su_vel2df_accumulate(int nt, float dt, float ft, float offset, float tn_zero, int nv, float dv, float fv, int nvh,
	float dvh, float fvh, int nsmooth, const float *data, float *num, float *den, float *nnz)
{
	int iv, ivh, it, is, ismin, ismax, iti;
	float v, vh, v2, v4, vh2, offset2, offovs, fac1, fac2, tn, ti, frac, temp;

	offset2 = offset*offset;
	for (iv=0,v=fv; iv<nv; ++iv,v+=dv) {
		for (ivh=0,vh=fvh; ivh<nvh; ++ivh,vh+=dvh) {
			size_t base = ((size_t) iv*nvh + ivh)*nsmooth;

			/* compute offset/velocity squared */
			v2=v*v;
			v4=v2*v2;
			vh2=vh*vh;
			offovs = offset2/v2;
			fac1   = (vh2-v2)*offset2*offset2/v2;
			fac2   = vh2*offset2;
			ismin = tn_zero/dt-nsmooth/2;
			ismax = ismin+nsmooth;
			if (ismin<0) ismin = 0;
			if (ismax>nt-1) ismax = nt-1;

			/* do nmo via quick and dirty linear interpolation */
			/* (accurate enough for velocity analysis) and */
			/* accumulate semblance numerator and denominator */
			for (it=ismin, is=0,tn=ft+ismin*dt; it<ismax; ++it, ++is,tn+=dt) {
				if(offset==0)
					ti = (sqrt(tn*tn+offovs)-ft)/dt;
				else
					ti = (sqrt(tn*tn+offovs-fac1/(tn*tn*v4+fac2))-ft)/dt;
				iti = ti;
				if (iti>=0 && iti<nt-1) {
					frac = ti-iti;
					temp = (1.0-frac)*data[iti]+frac*data[iti+1];
					if (temp!=0.0) {
						num[base+is] += temp;
						den[base+is] += temp*temp;
						nnz[base+is] += 1.0;
					}
				}
			}
		}
	}
}
