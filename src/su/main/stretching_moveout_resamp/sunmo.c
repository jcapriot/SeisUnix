/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUNMO: $Revision: 1.31 $ ; $Date: 2013/03/06 20:35:27 $		*/
 
#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									     ",
" SUNMO - NMO for an arbitrary velocity function of time and CDP	     ",
"									     ",
"  sunmo <stdin >stdout [optional parameters]				     ",
"									     ",
" Optional Parameters:							     ",
" tnmo=0,...		NMO times corresponding to velocities in vnmo	     ",
" vnmo=1500,...		NMO velocities corresponding to times in tnmo	     ",
" cdp=			CDPs for which vnmo & tnmo are specified (see Notes) ",
" smute=1.5		samples with NMO stretch exceeding smute are zeroed  ",
" lmute=25		length (in samples) of linear ramp for stretch mute  ",
" sscale=1		=1 to divide output samples by NMO stretch factor    ",
" invert=0		=1 to perform (approximate) inverse NMO		     ",
" upward=0		=1 to scan upward to find first sample to kill	     ",
" voutfile=		if set, interplolated velocity function v[cdp][t] is ",
"			output to named file.			     	     ",
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
" For NMO with a velocity function of time only, specify the arrays	     ",
"	   vnmo=v1,v2,... tnmo=t1,t2,...				     ",
" where v1 is the velocity at time t1, v2 is the velocity at time t2, ...    ",
" The times specified in the tnmo array must be monotonically increasing.    ",
" Linear interpolation and constant extrapolation of the specified velocities",
" is used to compute the velocities at times not specified.		     ",
" The same holds for the anisotropy coefficients as a function of time only. ",
"									     ",
" For NMO with a velocity function of time and CDP, specify the array	     ",
"	   cdp=cdp1,cdp2,...						     ",
" and, for each CDP specified, specify the vnmo and tnmo arrays as described ",
" above. The first (vnmo,tnmo) pair corresponds to the first cdp, and so on. ",
" Linear interpolation and constant extrapolation of 1/velocity^2 is used    ",
" to compute velocities at CDPs not specified.				     ",
"									     ",
" The format of the output interpolated velocity file is unformatted C floats",
" with vout[cdp][t], with time as the fast dimension and may be used as an   ",
" input velocity file for further processing.				     ",
"									     ",
" Note that this version of sunmo does not attempt to deal with	anisotropy.  ",
" The version of sunmo with experimental anisotropy support is \"sunmo_a\"   ",
"									     ",
NULL};
#endif


/* Credits:
 *	SEP: Shuki Ronen, Chuck Sword
 *	CWP: Shuki Ronen, Jack, Dave Hale, Bjoern Rommel
 *      Modified: 08/08/98 - Carlos E. Theodoro - option for lateral offset
 *      Modified: 07/11/02 - Sang-yong Suh -
 *	  added "upward" option to handle decreasing velocity function.
 *      CWP: Sept 2010: John Stockwell
 *	  1. replaced Carlos Theodoro's fix 
 *	  2. added  the instruction in the selfdoc to use suazimuth to set 
 *	      offset so that it accounts for lateral offset. 
 *        3. removed  Bjoren Rommel's anisotropy stuff. sunmo_a is the 
 *           version with the anisotropy parameters left in.
 *        4. note that scalel does not scale the offset field in
 *           the segy standard.
 * Technical Reference:
 *	The Common Depth Point Stack
 *	William A. Schneider
 *	Proc. IEEE, v. 72, n. 10, p. 1238-1254
 *	1984
 *
 * Trace header fields accessed: ns, dt, delrt, offset, cdp, scalel
 */
/**************** end self doc *******************************************/


/* Library version of SUNMO
 *
 * The main program is not part of the library. It reads the velocity functions and works out the sloth (1/v^2)
 * at the sample times of each trace, for the CDP of the trace. Doing that is left to the caller, and what is here
 * is what is done to a trace once the sloth, and the offset, are known:
 *
 *	su_nmo_tables:	times t(tn), the stretch, and the mute (these only change if the sloth or the offset do)
 *	su_nmo:		the NMO (or the inverse NMO) of a trace, with those tables
 *
 * All the arrays are the caller's, none are allocated here (the tables and the work array are made by
 * the main program).
 *
 * Differences from the program:
 *  - In the inverse NMO the first sample that is not muted can come out negative, and is used as an index. It is
 *    kept from being negative.
 *
 * nt must be at least 2 (4 for the inverse NMO).
 */

/* The tables for NMO, or the inverse NMO, of traces of nt samples that start at time ft and are dt apart, at offset
 * `offset`, and where the sloth (1/v^2) is ovvt[nt].
 *
 * smute: samples with NMO stretch exceeding smute are zeroed. upward: scan upward to find the first sample to kill.
 * sscale: the stretch is going to be applied to the amplitudes (so the table for it is made for the inverse NMO).
 *
 * ttn[nt]: time t(tn) for NMO, in samples. atn[nt]: amplitude a(tn) for NMO (the inverse of the stretch factor).
 * tnt[nt]: time tn(t) for the inverse NMO, in samples, and at[nt] amplitude a(t); (only used for the inverse NMO)
 * itmute: samples with indices less than this are zeroed
 */
void su_nmo_tables(int nt, float dt, float ft, float offset, const float *ovvt, float smute, int upward,
	int invert, int sscale, float *ttn, float *atn, float *tnt, float *at, int *itmute_out)
{
	int it, itmute;
	float tn, temp, tsq, osmute;

	/* compute time t(tn) (normalized) */
	temp = ((float) offset*offset)/(dt*dt);
	for (it=0,tn=ft/dt; it<nt; ++it,tn+=1.0) {
		tsq = temp*ovvt[it];
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
}

/* NMO (or inverse NMO) of data[nt], in place, with the tables from su_nmo_tables.
 *
 * lmute: length (in samples) of the linear ramp for the stretch mute (forward NMO only)
 * sscale: divide the output samples by the NMO stretch factor
 * q[nt]: scratch space
 */
void su_nmo(float *data, int nt, float dt, float ft, int itmute, int lmute, int sscale, int invert,
	const float *ttn, const float *atn, const float *tnt, const float *at, float *q)
{
	int it;

	/* if forward (not inverse) nmo */
	if (!invert) {
		/* do nmo via 8-point sinc interpolation */
		ints8r(nt,1.0,ft/dt,data,0.0,0.0,
			nt-itmute,(float *) &ttn[itmute],&q[itmute]);

		/* apply mute */
		for (it=0; it<itmute; ++it)
			q[it] = 0.0;

		/* apply linear ramp */
		for (it=itmute; it<itmute+lmute && it<nt; ++it)
			q[it] *= (float)(it-itmute+1)/(float)lmute;

		/* if specified, scale by the NMO stretch factor */
		if (sscale)
			for (it=itmute; it<nt; ++it)
				q[it] *= atn[it];

	/* else inverse nmo */
	} else {
		/* do inverse nmo via 8-point sinc interpolation */
		ints8r(nt,1.0,ft/dt,data,0.0,0.0,
			nt-itmute,(float *) &tnt[itmute],&q[itmute]);

		/* apply mute */
		for (it=0; it<itmute; ++it)
			q[it] = 0.0;

		/* if specified, undo NMO stretch factor scaling */
		if (sscale)
			for (it=itmute; it<nt; ++it)
				q[it] *= at[it];
	}

	/* copy the result to the trace */
	memcpy( (void *) data, (const void *) q, nt*sizeof(float));
}
