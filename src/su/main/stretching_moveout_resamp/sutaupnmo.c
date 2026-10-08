/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTAUPNMO: $Revision: 1.3 $ ; $Date: 2011/11/16 23:21:55 $		*/
 
#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUTAUPNMO - NMO for an arbitrary velocity function of tau and CDP	",
"									",
"  sutaupnmo <stdin >stdout [optional parameters]			",
"									",
" Optional Parameters:							",
" tnmo=0,...		NMO times corresponding to velocities in vnmo	",
" vnmo=1500,...		NMO velocities corresponding to times in tnmo	",
" cdp=			CDPs for which vnmo & tnmo are specified (see Notes) ",
" smute=1.5		samples with NMO stretch exceeding smute are zeroed  ",
" lmute=25		length (in samples) of linear ramp for stretch mute  ",
" sscale=1		=1 to divide output samples by NMO stretch factor    ",
"									",
" Notes:								",
"									",
" For constant-velocity NMO, specify only one vnmo=constant and omit tnmo.",
"									",
" For NMO with a velocity function of tau only, specify the arrays	",
"	   vnmo=v1,v2,... tnmo=t1,t2,...				",
" where v1 is the velocity at tau t1, v2 is the velocity at tau t2, ...    ",
" The taus specified in the tnmo array must be monotonically increasing.    ",
" Linear interpolation and constant extrapolation of the specified velocities",
" is used to compute the velocities at taus not specified.		",
"									",
" For NMO with a velocity function of tau and CDP, specify the array	",
"	   cdp=cdp1,cdp2,...						",
" and, for each CDP specified, specify the vnmo and tnmo arrays as described ",
" above. The first (vnmo,tnmo) pair corresponds to the first cdp, and so on. ",
" Linear interpolation and constant extrapolation of velocity^2 is used	 ",
" to compute velocities at CDPs not specified.				",
"									",
" Moveout is defined by							",
"									",
"  tau^2 + tau^2.p^2.vel^2						",
"									",
" Note: In general, the user should set the cdp parameter.  The default is   ",
"	to use tr.cdp from the first trace and assume only one cdp.	 ",
" Caveat:								",
" Taunmo should handle triplication					",
"									",
" NMO interpolation error is less than 1% for frequencies less than 60% of   ",
" the Nyquist frequency.						",
"									",
" Exact inverse NMO is not implemented, nor has anisotropy		",
" Example implementation:						",
"   sutaup dx=25 option=2 pmin=0 pmax=0.0007025 < cmpgather.su |	",
"   supef minlag=0.2 maxlag=0.8 |					",
"   sutaupnmo tnmo=0.5,2,4 vnmo=1500,2000,3200 smute=1.5 |		",
"   sumute key=tracr mode=1 ntaper=20 xmute=1,30,40,50,85,15  		",
"				 tmute=7.8,7.8,4.5,3.5,2.0,0.35 |	",
"   sustack key=cdp | ... [...]						",
"									",
NULL};
#endif

/* Credits:
 *	 Durham, Richard Hobbs modified from SUNMO credited below
 *	SEP: Shuki Ronen, Chuck Sword
 *	CWP: Shuki Ronen, Jack K. Cohen , Dave Hale
 *
 * Technical Reference:
 *	van der Baan papers in geophysics (2002 & 2004)
 *
 * Trace header fields accessed: ns, dt, delrt, offset, cdp, sy
 */
/**************** end self doc *******************************************/

/* Library version of SUTAUPNMO
 *
 * The main program is not part of the library. It reads the velocity functions, makes v^2 at the sample times of each
 * trace for the CDP of the trace, and takes the ray parameter of the trace from its header (f2 + (tracr - 1) * d2).
 * Doing that is left to the caller, and what is here is what is done to a trace once v^2 and p are known:
 *
 *	su_taupnmo_tables:	times t(tn), the stretch, and the mute (these only change if v^2 or p do)
 *	su_taupnmo:		the NMO of a trace, with those tables
 *
 * All the arrays are the caller's, none are allocated here.
 *
 * Differences from the program:
 *  - The program scales the samples by the stretch factor from the first sample that is muted, which does nothing
 *    (they are zero). It scales the samples that are not muted.
 *  - The linear ramp of the mute starts lmute samples before the first muted sample, and the program writes
 *    before the start of the trace when that is before sample 0. It is kept at sample 0.
 *
 * nt must be at least 2.
 */

/* The tables for the NMO of traces of nt samples that start at time ft and are dt apart, for the ray parameter p
 * (in time units per unit of x, the units of dt) where v^2 is vvt[nt].
 *
 * smute: samples with NMO stretch exceeding smute are zeroed.
 *
 * ttn[nt]: time t(tn) for NMO, in samples. atn[nt]: amplitude a(tn) for NMO (the inverse of the stretch factor).
 * itmute: samples with indices of this or more are zeroed
 */
void su_taupnmo_tables(int nt, float dt, float ft, float p, const float *vvt, float smute, float *ttn, float *atn,
	int *itmute_out)
{
	int it, itmute;
	float tn, tsq, osmute;

	/* compute time t(tn) (normalized) */
	for (it=0,tn=ft/dt; it<nt; ++it,tn+=1.0) {
		tsq = tn*tn - p*p*tn*tn*vvt[it];
		if (tsq < 0.0) {
			ttn[it] = 0;
		} else {
			ttn[it] = sqrt(tsq);
		}
	}

	/* compute inverse of stretch factor a(tn) */
	atn[0] = ttn[1]-ttn[0];
	for (it=1; it<nt; ++it)
		atn[it] = ttn[it]-ttn[it-1];

	/* determine index of first sample to be muted */
	osmute = 1.0/smute;
	for (it=0; it<nt-1 && atn[it]>osmute; ++it)
		;
	itmute = it;
	*itmute_out = itmute;
}

/* NMO of data[nt], in place, with the tables from su_taupnmo_tables.
 *
 * lmute: length (in samples) of the linear ramp before the muted samples
 * sscale: divide the output samples by the NMO stretch factor
 * q[nt]: scratch space
 */
void su_taupnmo(float *data, int nt, float dt, float ft, int itmute, int lmute, int sscale,
	const float *ttn, const float *atn, float *q)
{
	int it;

	/* do nmo via 8-point sinc interpolation */
	ints8r(nt,1.0,ft/dt,data,0.0,0.0,itmute,(float *) &ttn[0],&q[0]);

	/* apply mute */
	for (it=itmute; it<nt; ++it)
		q[it] = 0.0;

	/* apply linear ramp */
	for (it=(itmute-lmute > 0 ? itmute-lmute : 0); it<itmute && it<nt; ++it)
		q[it] *= (float)(itmute-it)/(float)lmute;

	/* if specified, scale by the NMO stretch factor */
	if (sscale)
		for (it=0; it<itmute; ++it)
			q[it] *= atn[it];

	/* copy the result to the trace */
	memcpy( (void *) data, (const void *) q, nt*sizeof(float));
}
