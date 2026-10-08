/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVELAN: $Revision: 1.20 $ ; $Date: 2011/11/16 23:40:27 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									     ",
" SUVELAN - compute stacking velocity semblance for cdp gathers		     ",
"									     ",
" suvelan <stdin >stdout [optional parameters]				     ",
"									     ",
" Optional Parameters:							     ",
" nv=50                   number of velocities				     ",
" dv=50.0                 velocity sampling interval			     ",
" fv=1500.0               first velocity				     ",
" anis1=0.0               quartic term, numerator of an extended quartic term",
" anis2=0.0               in denominator of an extended quartic term         ",
" smute=1.5               samples with NMO stretch exceeding smute are zeroed",
" dtratio=5               ratio of output to input time sampling intervals   ",
" nsmooth=dtratio*2+1     length of semblance num and den smoothing window   ",
" verbose=0               =1 for diagnostic print on stderr		     ",
" pwr=1.0                 semblance value to the power      		     ",
"									     ",
" Notes:								     ",
" Velocity analysis is usually a two-dimensional screen for optimal values of",
" the vertical two-way traveltime and stacking velocity. But if the travel-  ",
" time curve is no longer close to a hyperbola, the quartic term of the      ",
" traveltime series should be considered. In its easiest form (with anis2=0) ",
" the optimizion of all parameters requires a three-dimensional screen. This ",
" is done by a repetition of the conventional two-dimensional screen with a  ",
" variation of the quartic term. The extended quartic term is more accurate, ",
" though the function is no more a polynomial. When screening for optimal    ",
" values the theoretical dependencies between these paramters can be taken   ",
" into account. The traveltime function is defined by                        ",
"                                                                            ",
"                1            anis1                                          ",
" t^2 = t_0^2 + --- x^2 + ------------- x^4                                  ",
"               v^2       1 + anis2 x^2                                      ",
"                                                                            ",
" The coefficients anis1, anis2 are assumed to be small, that means the non- ",
" hyperbolicity is assumed to be small. Triplications cannot be handled.     ",
"                                                                            ",
" Semblance is defined by the following quotient:			     ",
"									     ",
"                 n-1                 					     ",
"               [ sum q(t,j) ]^2      					     ",
"                 j=0                 					     ",
"       s(t) = ------------------     					     ",
"                 n-1                 					     ",
"               n sum [q(t,j)]^2      					     ",
"                 j=0                 					     ",
"									     ",
" where n is the number of non-zero samples after muting.		     ",
" Smoothing (nsmooth) is applied separately to the numerator and denominator ",
" before computing this semblance quotient.				     ",
"									     ",
" Then, the semblance is set to the power of the parameter pwr. With pwr > 1 ",
" the difference between semblance values is stretched in the upper half of  ",
" the range of semblance values [0,1], but compressed in the lower half of   ",
" it; thus, the few large semblance values are enhanced. With pwr < 1 the    ",
" many small values are enhanced, thus more discernible against background   ",
" noise. Of course, always at the expanse of the respective other feature.   ",
"									     ",
" Input traces should be sorted by cdp - suvelan outputs a group of	     ",
" semblance traces every time cdp changes.  Therefore, the output will	     ",
" be useful only if cdp gathers are input.				     ",
NULL};
#endif


/* Credits:
 *	CWP, Colorado School of Mines:
 *           Dave Hale (everything except ...)
 *           Bjoern Rommel (... the quartic term)
 *      SINTEF, IKU Petroleumsforskning
 *           Bjoern Rommel (... the power-of-semblance function)
 *
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp
 * Trace header fields modified:  ns, dt, offset, cdp
 */
/**************** end self doc *******************************************/

/* Library version of SUVELAN
 *
 * The main program is not part of the library. It reads the gathers (the traces that have the same cdp), and writes nv
 * semblance traces for each of them. Doing that is left to the caller, and what is here is the work of the program on
 * the traces and the sums:
 *
 *	su_velan_accumulate:	NMO of a trace at each of the velocities, added to the sums for the semblance
 *	su_velan_semblance:	the semblance of one velocity from its sums
 *
 * All the arrays are the caller's, none are allocated here. The sums are num[nv*nt], den[nv*nt] and nnz[nv*nt] (the row of
 * each velocity after the other), that start at 0 for each gather.
 *
 * Differences from the program:
 *  - The smoothing window of the semblance has nsmooth samples, where the program leaves out the last one (SURELAN,
 *    which does the same job, has them all).
 *  - A sample whose NMO time is not real (a velocity that has no moveout, see the return value) is left out, where the
 *    program takes the square root of a negative number and uses what comes out as an index.
 */

/* NMO of a trace of nt samples that start at time ft and are dt apart, at offset `offset` and at the velocities
 * fv + iv * dv, with linear interpolation (accurate enough for velocity analysis) and the sums for the semblance.
 *
 * anis1, anis2: the quartic term, and the extension of it in the denominator. smute: samples with NMO stretch exceeding
 * smute are zeroed.
 *
 * Returns -1 if anis2 is too small for this offset (the program stops with an error), 1 if there is a velocity that has no
 * moveout (the program warns: check anis1 and anis2), otherwise 0.
 */
int su_velan_accumulate(int nt, float dt, float ft, float offset, int nv, float dv, float fv, float anis1,
	float anis2, float smute, const float *data, float *num, float *den, float *nnz)
{
	int iv,it,itmute,iti;
	float v,offovs,offan,tnmute,tn,ti,frac,temp,tsq;
	int nomoveout = 0;

	if ((1.0 + offset*offset*anis2) <= 0.0) return -1;
	offan = (offset*offset*offset*offset*anis1) / (1.0 + offset*offset*anis2);

	for (iv=0,v=fv; iv<nv; ++iv,v+=dv) {
		float *pnum = num + (size_t) iv*nt;
		float *pden = den + (size_t) iv*nt;
		float *pnnz = nnz + (size_t) iv*nt;

		/* compute offset/velocity squared */
		offovs = (offset*offset)/(v*v) + offan;
		if (offovs < 0.0) nomoveout = 1;

		/* determine mute time after nmo */
		tnmute = sqrt(offovs/(smute*smute-1.0));
		if (tnmute > ft) {
			itmute = (tnmute-ft)/dt;
		} else {
			itmute = 0;
		}

		/* do nmo via quick and dirty linear interpolation
		   (accurate enough for velocity analysis) and
		   accumulate semblance numerator and denominator
		*/
		for (it=itmute,tn=ft+itmute*dt; it<nt; ++it,tn+=dt) {
			tsq = tn*tn+offovs;
			if (tsq < 0.0) continue;
			ti = (sqrt(tsq)-ft)/dt;
			iti = ti;
			if (iti>=0 && iti<nt-1) {
				frac = ti-iti;
				temp = (1.0-frac)*data[iti]+frac*data[iti+1];
				if (temp!=0.0) {
					pnum[it] += temp;
					pden[it] += temp*temp;
					pnnz[it] += 1.0;
				}
			}
		}
	}
	return nomoveout;
}

/* The semblance of one velocity (or residual moveout), from its sums num[nt], den[nt] and nnz[nt]: ntout = 1+(nt-1)/dtratio
 * samples, each from the sums of the nsmooth samples about it, to the power pwr.
 */
void su_velan_semblance(int nt, int ntout, int dtratio, int nsmooth, float pwr, const float *num, const float *den,
	const float *nnz, float *sem)
{
	int itout,it,is,ismin,ismax;
	float nsum,dsum;

	for (itout=0; itout<ntout; ++itout) {
		it = itout*dtratio;
		ismin = it-nsmooth/2;
		ismax = it+nsmooth/2;
		if (ismin<0) ismin = 0;
		if (ismax>nt-1) ismax = nt-1;
		nsum = dsum = 0.0;
		for (is=ismin; is<=ismax; ++is) {
			nsum += num[is]*num[is];
			dsum += nnz[is]*den[is];
		}
		sem[itout] = (dsum!=0.0?nsum/dsum:0.0);
	}

	/* powering the semblance */
	if (pwr != 1.0)
		for (itout=0; itout<ntout; ++itout)
			sem[itout] = pow(sem[itout], pwr);
}
