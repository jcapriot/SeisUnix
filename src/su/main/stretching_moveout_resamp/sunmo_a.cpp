/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNMO: $Revision: 1.5 $ ; $Date: 2012/01/03 19:45:15 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									     ",
" SUNMO_a - NMO for an arbitrary velocity function of time and CDP with	     ",
"		experimental Anisotropy options				     ",
"  sunmo <stdin >stdout [optional parameters]				     ",
"									     ",
" Optional Parameters:							     ",
" tnmo=0,...		NMO times corresponding to velocities in vnmo	     ",
" vnmo=1500,..		NMO velocities corresponding to times in tnmo	     ",
" anis1=0		two anisotropy coefficients making up quartic term   ",
" anis2=0		in traveltime curve, corresponding to times in tnmo  ",
" cdp=			CDPs for which vnmo & tnmo are specified (see Notes) ",
" smute=1.5		samples with NMO stretch exceeding smute are zeroed  ",
" lmute=25		length (in samples) of linear ramp for stretch mute  ",
" sscale=1		=1 to divide output samples by NMO stretch factor    ",
" invert=0		=1 to perform (approximate) inverse NMO		     ",
" upward=0		=1 to scan upward to find first sample to kill	",
"									     ",
" Notes:								     ",
" For constant-velocity NMO, specify only one vnmo=constant and omit tnmo.   ",
"									     ",
" NMO interpolation error is less than 1% for frequencies less than 60% of   ",
" the Nyquist frequency.						     ",
"									     ",
" Exact inverse NMO is impossible, particularly for early times at large     ",
" offsets and for frequencies near Nyquist with large interpolation errors.  ",
" 								     	     ",
" The \"offset\" header field must be set.				     ",
" Use suazimuth to set offset header field when sx,sy,gx,gy are all	     ",
" nonzero. 							   	     ",
"									     ",
" Anisotropy option:							     ",
" Caveat, this is an experimental option,				     ",
"									     ",
" The anisotropy coefficients anis1, anis2 permit non-hyperbolicity due	     ",
" to layering, mode conversion, or anisotropy. Default is isotropic NMO.     ",
"									     ",
" Moveout is defined by							     ",
"									     ",
"   1		 anis1							     ",
"  --- x^2 + ------------- x^4.						     ",
"  v^2	     1 + anis2 x^2						     ",
"									     ",
" Caveat:								     ",
" Nmo cannot handle negative moveout as in triplication caused by	     ",
" anisotropy. But negative moveout happens necessarily for negative anis1 at ",
" sufficiently large offsets. Then the error-negative moveout- is printed.   ",
" Check anis1. An error (anis2 too small) is also printed if the	     ",
" denominator of the quartic term becomes negative. Check anis2. These errors",
" are prompted even if they occur in traces which would not survive the	     ",
" NMO-stretch threshold. Chop off enough far-offset traces (e.g. with suwind)",
" if anis1, anis2 are fine for near-offset traces.			     ",
"									     ",
NULL};
#endif

/* Credits:
 *	SEP: Shuki, Chuck Sword
 *	CWP: Shuki, Jack, Dave Hale, Bjoern Rommel
 *      Modified: 08/08/98 - Carlos E. Theodoro - option for lateral offset
 *      Modified: 07/11/02 - Sang-yong Suh -
 *	  added "upward" option to handle decreasing velocity function.
 *      CWP: Sept 2010: John Stockwell
 *	  replaced Carlos Theodoro's fix
 *	  and added the instruction in the selfdoc to use suazimuth to set
 *	    offset so that it accounts for lateral offset
 *
 * Technical Reference:
 *	The Common Depth Point Stack
 *	William A. Schneider
 *	Proc. IEEE, v. 72, n. 10, p. 1238-1254
 *	1984
 *
 * Trace header fields accessed: ns, dt, delrt, offset, cdp, scalel
 */
/**************** end self doc *******************************************/

/* The tables of the program, as su_nmo_tables of sunmo.cpp makes them for sunmo, with the anisotropy coefficients
 * (the quartic term of the traveltime). The NMO of a trace with them is su_nmo (sunmo.cpp).
 *
 * oa1t[nt], oa2t[nt]: anis1 and anis2 at the sample times.
 *
 * Returns 0, 1 if there is negative moveout (check anis1: the program stops with an error), or 2 if 1 + anis2 x^2
 * is not positive (check anis2).
 *
 * As in su_nmo_tables, the first sample that is not muted of the inverse NMO is kept from being negative (the
 * program uses it as an index).
 */
int su_nmo_a_tables(int nt, float dt, float ft, float offset, const float *ovvt, const float *oa1t, const float *oa2t,
	float smute, int upward, int invert, int sscale, float *ttn, float *atn, float *tnt, float *at, int *itmute_out)
{
	int it, itmute;
	float tn, temp, tsq, osmute;

	/* compute time t(tn) (normalized) */
	temp = ((float) offset*offset)/(dt*dt);
	for (it=0,tn=ft/dt; it<nt; ++it,tn+=1.0) {
		tsq = temp*ovvt[it] + oa1t[it]*temp*temp / (1.0+oa2t[it]*temp);
		if (tsq<0.0) return 1;
		if ((1.0+oa2t[it]*temp)<=0.0) return 2;
		ttn[it] = sqrt (tn*tn + tsq);
	}

	/* compute inverse of stretch factor a(tn) */
	atn[0] = ttn[1]-ttn[0];
	for (it=1; it<nt; ++it)
		atn[it] = ttn[it]-ttn[it-1];

	/* determine index of first sample to survive mute */
	osmute = 1.0/smute;
	if(!upward) {
		for (it=0; it<nt-1 && atn[it]<osmute; ++it);
	} else {
		/* scan samples from bottom to top */
		for (it=nt-1; it>0 && atn[it]>=osmute; --it);
	}
	itmute = it;

	/* if inverse NMO will be performed */
	if (invert) {
		/* compute tn(t) from t(tn) */
		yxtoxy(nt-itmute,1.0,ft/dt+itmute,&ttn[itmute],
			nt-itmute,1.0,ft/dt+itmute,
			ft/dt-nt,ft/dt+nt,&tnt[itmute]);

		/* adjust mute time */
		itmute = 1.0+ttn[itmute]-ft/dt;
		itmute = MIN(nt-2,itmute);
		itmute = MAX(0,itmute);	/* (the program does not have this) */

		/* compute a(t) */
		if (sscale) {
			for (it=itmute+1; it<nt; ++it)
				at[it] = tnt[it]-tnt[it-1];
			at[itmute] = at[itmute+1];
		}
	}
	*itmute_out = itmute;
	return 0;
}
