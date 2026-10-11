/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUDMOTX: $Revision: 1.17 $ ; $Date: 2015/08/07 22:01:42 $		*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUDMOTX - DMO via T-X domain (Kirchhoff) method for common-offset gathers",
"									",
" sudmotx <stdin >stdout cdpmin= cdpmax= dxcdp= noffmix= [optional parms]",
"									",
" Required Parameters:							",
" cdpmin                  minimum cdp (integer number) for which to apply DMO",
" cdpmax                  maximum cdp (integer number) for which to apply DMO",
" dxcdp                   distance between successive cdps		",
" noffmix                 number of offsets to mix (see notes)		",
"									",
" Optional Parameters:							",
" offmax=3000.0           maximum offset				",
" tmute=2.0               mute time at maximum offset offmax		",
" vrms=1500.0             RMS velocity at mute time tmute		",
" verbose=0               =1 for diagnostic print			",
" tmpdir=	if non-empty, use the value as a directory path	prefix	",
"		for storing temporary files; else if the CWP_TMPDIR	",
"		environment variable is set use	its value for the path;	",
"		else use tmpfile()					",
"									",
"									",
" Notes:								",
" Input traces should be sorted into common-offset gathers.  One common-",
" offset gather ends and another begins when the offset field of the trace",
" headers changes.							",
" 									",
" The cdp field of the input trace headers must be the cdp bin NUMBER, NOT",
" the cdp location expressed in units of meters or feet.		",
"									",
" The number of offsets to mix (noffmix) should typically equal the ratio of",
" the shotpoint spacing to the cdp spacing.  This choice ensures that every",
" cdp will be represented in each offset mix.  Traces in each mix will	",
" contribute through DMO to other traces in adjacent cdps within that mix.",
"									",
" The defaults for offmax and vrms are appropriate only for metric units.",
" If distances are measured in feet, then these parameters should be	",
" specified explicitly.							",
"									",
" offmax, tmute, and vrms need not be specified precisely.		",
" If these values are unknown, then one should overestimate offmax	",
" and underestimate tmute and vrms.					",
"									",
" No muting is actually performed.  The tmute parameter is used only to	",
" determine parameters required to perform DMO.				",
NULL};
#endif


/* Credits:
 *	CWP: Dave Hale
 *
 * Technical Reference:
 *      A non-aliased integral method for dip-moveout
 *      Dave Hale
 *      submitted to Geophysics, June, 1990
 *
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp.
 */
/**************** end self doc *******************************************/

/* Library version of SUDMOTX
 *
 * The main program (the parameters, the traces and their headers, the offsets and the mixes of them) is left to the caller.
 *
 *	su_dmotx_tables:	the time shifts and amplitudes of the operator (maketa of the program), and the shaping filter
 *				that corrects its horizontal reflection response (makeds)
 *	su_dmotx_trace:		the contribution of one input trace to the cdps that it can change (the loop of the program
 *				over the output traces, with dmotx)
 *	su_dmotx_shape:		the shaping filter applied to the accumulated traces of all of the cdps
 *
 * Differences from the program:
 *  - makeds had arrays of 400 samples on the stack, and wrote before the start of the one of the response when a time
 *    shift was larger than that (a trace with a small dt); it has the length that it needs.
 *  - The first sample to accumulate was a float made an int without a bound, which is not defined past the range of the
 *    int (an offset that is large against the cdp interval); it is limited to nt.
 *  - The program counted a change of offset only for a trace that is in the range of cdps, so that a trace out of the
 *    range with another offset made a new mix of every trace that followed it; a change of offset is counted for every
 *    trace (in the caller).
 */

#define DMOTX_LDSMAX 300

static void maketa (float dx, float dt, float offmax, float tmute, float vrms,
	int nsmax, int *ns, float *ts, float *as)
/*****************************************************************************
make time shifts and amplitudes for non-aliased (t,x) domain DMO
******************************************************************************
Input:
dx		midpoint sampling interval (see notes)
dt		time sampling interval (see notes)
offmax		maximum offset (see notes)
tmute		mute time at maximum offset (see notes)
vrms		RMS velocity at mute time (see notes)
nsmax		maximum number of time shifts

Output:
ns		number of shifts
ts		array[ns] of time shifts (normalized by sampling interval dt)
as		array[ns] of amplitudes corresponding to time shifts
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/04/90
*****************************************************************************/
{
	int nsi,is,itaper,ltaper;
	float alpha,dalpha,dalphan,tsmax;

	/* ensure positive arguments */
	dx = ABS(dx);
	dt = ABS(dt);
	offmax = ABS(offmax);
	tmute = ABS(tmute);
	vrms = ABS(vrms);

	/* determine nominal increment in alpha = sqrt(ts) */
	dalphan = MIN(sqrt(2.0*(tmute/dt)*pow((dx/offmax),2.0)),1.0);

	/* determine maximum time shift for steepest reflection slope */
	tsmax = (tmute/dt)*(1.0-1.0/sqrt(1.0+pow(offmax/(vrms*tmute),2.0)));

	/* recursively compute time shifts and amplitudes */
	alpha = 0.0;
	dalpha = dalphan;
	for (nsi=0; alpha*alpha<=tsmax && nsi<nsmax; ++nsi) {
		ts[nsi] = alpha*alpha;
		as[nsi] = dalpha;
		alpha = alpha+dalpha;
		dalpha = sqrt(alpha*alpha+1.0)-alpha;
		if (dalpha>dalphan) dalpha = dalphan;
	}

	/* taper the amplitudes with a raised cosine */
	ltaper = nsi/3;
	for (is=nsi-ltaper,itaper=1; is<nsi; ++is,++itaper)
		as[is] *= 0.54+0.46*cos(PI*itaper/ltaper);
	*ns = nsi;
}

static void makeds (int ns, float *ts, float *as, int lds, int ifds, float *ds)
/*****************************************************************************
make shaping filter to correct DMO horizontal reflection response
******************************************************************************
Input:
ns		number of shifts
ts		array[ns] of time shifts (normalized by sampling interval)
as		array[ns] of amplitudes corresponding to time shifts
lds		length of shaping filter
ifds		index of first sample of shaping filter (see notes)

Output:
ds		array[lds] containing shaping filter
******************************************************************************
Notes:
Reasonable values for lds and ifds are lds=125 and ifds=-100.
This function must be kept consistent with that used to perform DMO.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 11/04/90
*****************************************************************************/
{
	int ld=400,ifd,i,is;
	float one=1.0,*d,*dd,*di,*work;

	/* the response is as long as the largest time shift needs */
	if (ns>0 && (int)ts[ns-1]+2>ld) ld = (int)ts[ns-1]+2;
	ifd = 1-ld;
	d = ealloc1float(ld);
	dd = ealloc1float(lds);
	di = ealloc1float(lds);
	work = ealloc1float(lds);

	/* compute dmo horizontal reflection response d(t) */
	for (i=0; i<ld; ++i)
		d[i] = 0.0;
	d[ld-1] = as[0];
	for (is=1; is<ns; ++is) {
		i = (int)ts[is];
		d[ld-i-1] += 2.0*as[is];
	}

	/* compute autocorrelation of d(t) */
	xcor(ld,ifd,d,ld,ifd,d,lds,0,dd);

	/* compute crosscorrelation of d(t) and desired impulse i(t) */
	xcor(ld,ifd,d,1,0,&one,lds,ifds,di);

	/* solve symmetric toeplitz system of equations for filter */
	stoepf(lds,dd,di,ds,work);

	free1float(d);
	free1float(dd);
	free1float(di);
	free1float(work);
}

/* the tables of the operator; ts and as have nsmax values, ds has lds values; returns the number of shifts, or a negative
   code if the parameters are not allowed */
int su_dmotx_tables(float dxcdp, float dt, float offmax, float tmute, float vrms, int nsmax, int lds, int ifds,
	float *ts, float *as, float *ds)
{
	int ns;

	if (dxcdp==0.0 || dt==0.0 || offmax==0.0 || tmute==0.0 || vrms==0.0 || nsmax<2 || lds<1 || lds>DMOTX_LDSMAX) return -3;
	maketa(dxcdp,dt,offmax,tmute,vrms,nsmax,&ns,ts,as);
	if (ns<1) return -3;
	makeds(ns,ts,as,lds,ifds,ds);
	return ns;
}

/* the contribution of one trace p[nt] at cdp index icdp of ncdp to the traces q[ncdp][nt] */
void su_dmotx_trace(int ns, const float *ts, const float *as, float offset, int icdp, int ncdp, float dxcdp, int nt, float dt,
	float ft, const float *p, float *q)
{
	int it1,it2,its,it,is,jcdp,jcdplo,jcdphi,itmute=0;
	float h,xb,ftodt,x1,x2,scale1,scale2,asi,tsi,x;

	jcdplo = MAX(0,icdp-0.5*ABS(offset/dxcdp));
	jcdphi = MIN(ncdp-1,icdp+0.5*ABS(offset/dxcdp));

	/* compute half-offset */
	h = offset/2.0;

	/* normalize the time of first sample */
	ftodt = ft/dt;

	/* loop over potentially modified output traces */
	for (jcdp=jcdplo; jcdp<=jcdphi; ++jcdp) {
		float *qj = q+(size_t)jcdp*nt;

		/* midpoint distance from input trace to output trace */
		x = (jcdp-icdp)*dxcdp;

		/* round midpoint distance to nearest multiple of sampling interval */
		xb = ABS(NINT(x/dxcdp)*dxcdp);

		/* if rounded midpoint distance is zero */
		if (xb==0.0) {

			/* compute midpoint distance at right edge of midpoint bin */
			x1 = xb+0.5*ABS(dxcdp);

			/* compute scale factor used to determine index it1 below */
			scale1 = 1.0/(1.0-sqrt(MAX(0.0,1.0-x1*x1/MAX(h*h,x1*x1))));

			/* accumulate output for first (zero) time shift */
			asi = as[0];
			for (it=itmute; it<nt; ++it)
				qj[it] += asi*p[it];

			/* loop over non-zero time shifts */
			for (is=1; is<ns; ++is) {
				float first;

				tsi = ts[is];
				asi = as[is];

				/* round the normalized time shift to lower int */
				its = tsi;

				/* compute index of first time sample for this shift (bounded) */
				first = scale1*tsi-ftodt;
				it1 = first>=nt ? nt : (int)first;
				if (it1<itmute+its) it1 = itmute+its;

				/* if index is beyond last time sample, break */
				if (it1>=nt) break;

				/* accumulate shifted and weighted p(t) in q(t) */
				for (it=it1; it<nt; ++it)
					qj[it-its] += 2.0*asi*p[it];
			}

		/* else if non-zero time shifts exist */
		} else if (ns>1) {

			/* compute midpoint distances at right and left edges of bin */
			x1 = xb+0.5*ABS(dxcdp);
			x2 = xb-0.5*ABS(dxcdp);

			/* if midpoint distance exceeds half-offset, return */
			if (x2*x2>h*h) continue;

			/* compute scale factors for it1 and it2 below */
			scale1 = 1.0/(1.0-sqrt(MAX(0.0,1.0-x1*x1/MAX(h*h,x1*x1))));
			scale2 = 1.0/(1.0-sqrt(MAX(0.0,1.0-x2*x2/MAX(h*h,x2*x2))));

			/* loop over time shifts */
			for (is=1; is<ns; ++is) {
				float first,last;

				tsi = ts[is];
				asi = as[is];

				/* round normalized time shift to lower int */
				its = tsi;

				/* compute indices of first and last samples */
				first = scale1*tsi-ftodt;
				it1 = first>=nt ? nt : (int)first;
				if (it1<itmute+its) it1 = itmute+its;
				last = scale2*tsi-ftodt;
				it2 = last>=nt ? nt : (int)last;
				if (it2>nt) it2 = nt;

				/* if first sample index is out of bounds, break */
				if (it1>=nt) break;

				/* compute contribution to q(t) for this shift */
				for (it=it1; it<it2; ++it)
					qj[it-its] += asi*p[it];
			}
		}
	}
}

/* the shaping filter ds[lds], first sample at ifds, applied to each of the ncdp traces of q[ncdp][nt], in place */
void su_dmotx_shape(int ncdp, int nt, int lds, int ifds, const float *ds_in, float *q)
{
	int icdp;
	float *temp = ealloc1float(nt), *ds = ealloc1float(lds);

	memcpy(ds,ds_in,sizeof(float)*lds);
	for (icdp=0; icdp<ncdp; ++icdp) {
		float *qi = q+(size_t)icdp*nt;

		convolve_cwp(lds,ifds,ds,nt,0,qi,nt,0,temp);
		memcpy(qi,temp,nt*sizeof(float));
	}
	free1float(temp);
	free1float(ds);
}
