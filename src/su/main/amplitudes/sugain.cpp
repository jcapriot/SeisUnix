/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUGAIN: $Revision: 1.63 $ ; $Date: 2019/07/22 16:31:00 $		*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <signal.h>
#include <float.h>

/*********************** self documentation *****************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUGAIN - apply various types of gain				  	",
"									",
" sugain <stdin >stdout [optional parameters]			   	",
"									",
" Required parameters:						  	",
"	none (no-op)						    	",
"									",
" Optional parameters:						  	",
"	panel=0	        =1  gain whole data set (vs. trace by trace)	",
"	tpow=0.0	multiply data by t^tpow			 	",
"	epow=0.0	multiply data by exp(epow*t)		    	",
"	etpow=1.0	multiply data by exp(epow*t^etpow)	    	",
"	gpow=1.0	take signed gpowth power of scaled data	 	",
"	agc=0	   flag; 1 = do automatic gain control	     		",
"	gagc=0	  flag; 1 = ... with gaussian taper			",
"	wagc=0.5	agc window in seconds (use if agc=1 or gagc=1)  ",
"	trap=none	zero any value whose magnitude exceeds trapval  ",
"	clip=none	clip any value whose magnitude exceeds clipval  ",
"	pclip=none	clip any value greater than clipval  		",
"	nclip=none	clip any value less than  clipval 		",
"	qclip=1.0	clip by quantile on absolute values on trace    ",
"	qbal=0	  flag; 1 = balance traces by qclip and scale     	",
"	pbal=0	  flag; 1 = bal traces by dividing by rms value   	",
"	mbal=0	  flag; 1 = bal traces by subtracting the mean    	",
"	maxbal=0	flag; 1 = balance traces by subtracting the max ",
"	scale=1.0	multiply data by overall scale factor	   	",
"	norm=0.0	divide data by overall scale factor	     	",
"	bias=0.0	bias data by adding an overall bias value	",
"	jon=0	   	flag; 1 means tpow=2, gpow=.5, qclip=.95	",
"	verbose=0	verbose = 1 echoes info				",
"	mark=0		apply gain only to traces with tr.mark=0	",
"			=1 apply gain only to traces with tr.mark!=0    ",
"	vred=0	  reducing velocity of data to use with tpow		",
"									",
" 	tmpdir=		if non-empty, use the value as a directory path	",
"			prefix for storing temporary files; else if the ",
"			the CWP_TMPDIR environment variable is set use  ",
"			its value for the path; else use tmpfile()	",
"									",
" Operation order:							",
" if (norm) scale/norm						  	",
"									",
" out(t) = scale * BAL{CLIP[AGC{[t^tpow * exp(epow * t^tpow) * ( in(t)-bias )]^gpow}]}",
"									",
" Notes:								",
"	The jon flag selects the parameter choices discussed in		",
"	Claerbout's Imaging the Earth, pp 233-236.			",
"									",
"	Extremely large/small values may be lost during agc. Windowing  ",
"	these off and applying a scale in a preliminary pass through	",
"	sugain may help.						",
"									",
"	Sugain only applies gain to traces with tr.mark=0. Use sushw,	",
"	suchw, suedit, or suxedit to mark traces you do not want gained.",
"	See the selfdocs of sushw, suchw, suedit, and suxedit for more	",
"	information about setting header fields. Use \"sukeyword mark\" ",
"	for more information about the mark header field.		",
"									",
"      debias data by using mbal=1					",
"									",
"      option etpow only becomes active if epow is nonzero		",
NULL};
#endif


/* Credits:
 *	SEP: Jon Claerbout
 *	CWP: Jack K. Cohen, Brian Sumner, Dave Hale
 *
 * Note: Have assumed tr.deltr >= 0 in tpow routine.
 *
 * Technical Reference:
 *	Jon's second book, pages 233-236.
 *
 * Trace header fields accessed: ns, dt, delrt, mark, offset
 */
/**************** end self doc *******************************************/


/* Library version of SUGAIN
 *
 * The main program is not part of the library, the functions below are what it does to a trace.
 *
 * Compared to sugain.c the lookup tables and scratch arrays are arguments, instead of function level
 * statics that are allocated and filled by the first trace to come through (those could not be shared
 * between traces of different lengths, or used from more than one thread). The caller owns all the memory.
 * The reducing velocity is handled by the caller, by choosing the reduced time `tred` of su_gain_tpow_table().
 */

#define EPS     3.8090232	/* exp(-EPS*EPS) = 5e-7, "noise" level  */

static float quant(float *a, int k, int n);
static void do_trap(float *data, float trap, int nt);
static void do_clip(float *data, float clip, int nt);
static void do_nclip(float *data, float nclip, int nt);
static void do_pclip(float *data, float pclip, int nt);
static void do_qclip(float *data, float qclip, int nt, float *absdata);
static void do_qbal(float *data, float qclip, int nt, float *absdata);
static void do_agc(float *data, int iwagc, int nt, float *agcdata, float *d2);
static void do_gagc(float *data, int iwagc, int nt, float *agcdata, float *w, float *d2, float *s);

/* Table of t^tpow for nt samples starting at time tmin, with sample interval dt.
 * tred is the reduced time (|offset|/vred), or 0 for no reduction. */
void su_gain_tpow_table(float *tpowfac, int nt, float tmin, float dt, float tpow, float tred)
{
	int i;

	/* protect against negative tpow */
	tpowfac[0] = (tmin == 0.0) ? 0.0 : pow(tmin, tpow);
	for (i = 1; i < nt; ++i)
		tpowfac[i] = pow(tmin + tred + i*dt, tpow);
}

/* Table of exp(epow * t^etpow) for nt samples starting at time tmin, with sample interval dt */
void su_gain_epow_table(float *epowfac, int nt, float tmin, float dt, float epow, float etpow)
{
	int i;
	float etpowfac;

	/* protect against negative etpow */
	etpowfac = (tmin == 0.0) ? 0.0 : pow(tmin, etpow);
	epowfac[0] = exp(epow * etpowfac);
	for (i = 1; i < nt; ++i) {
		etpowfac = pow(tmin + i*dt, etpow);
		epowfac[i] = exp(epow * etpowfac);
	}
}

/* Apply all the various gains to data[nt]. This is gain() from sugain.c.
 *
 * tpowfac, epowfac: tables from the functions above (only used if tpow, or epow, are not 0)
 *
 * The scratch arrays, which the program has as statics, are:
 * absdata[nt]	used by qclip and qbal
 * agcdata[nt]	used by agc and gagc
 * d2[nt]	used by agc and gagc
 * w[iwagc]	used by gagc
 * s[nt]	used by gagc
 */
void su_gain(float *data, float tpow, float epow, float gpow,
	  int agc, int gagc, int qbal, int pbal, int mbal, float scale, float bias,
	  float trap, float clip, float qclip, int iwagc,
	  int nt, int maxbal, float pclip, float nclip,
	  const float *tpowfac, const float *epowfac,
	  float *absdata, float *agcdata, float *d2, float *w, float *s)
{
	float f_two  = 2.0;
	float f_one  = 1.0;
	float f_half = 0.5;
	int i;

	if (bias) {
		for (i = 0; i < nt; ++i)  data[i]+=bias ;
	}

	if (tpow) {
		for (i = 0; i < nt; ++i)  data[i] *= tpowfac[i];
	}

	if (epow) {
		for (i = 0; i < nt; ++i)  data[i] *= epowfac[i];
	}

	if (!CLOSETO(gpow, f_one)) {
		float val;
		if (CLOSETO(gpow, f_half)) {
			for (i = 0; i < nt; ++i) {
				val = data[i];
				data[i] = (val >= 0.0) ?
					sqrt(val) :
				-sqrt(-val);
			}
		} else if (CLOSETO(gpow, f_two)) {
			for (i = 0; i < nt; ++i) {
				val = data[i];
				data[i] = val * ABS(val);
			}
		} else {
			for (i = 0; i < nt; ++i) {
				val = data[i];
				data[i] = (val >= 0.0) ?
					pow(val, gpow) :
				-pow(-val, gpow);
			}
		}
	}

	if (agc)		   do_agc(data, iwagc, nt, agcdata, d2);
	if (gagc)		   do_gagc(data, iwagc, nt, agcdata, w, d2, s);
	if (trap > 0.0)	    do_trap(data, trap, nt);
	if (clip > 0.0)	    do_clip(data, clip, nt);
	if (pclip < FLT_MAX )	do_pclip(data, pclip, nt);
	if (nclip > -FLT_MAX )     do_nclip(data, nclip, nt);
	if (qclip < 1.0 && !qbal)  do_qclip(data, qclip, nt, absdata);
	if (qbal)		  do_qbal(data, qclip, nt, absdata);

	if (pbal) {
		int i;
		float val;
		float rmsq = 0.0;

		/* rmsq = sqrt (SUM( a()*a() ) / nt) */
		for (i = 0; i < nt; ++i) {
			val = data[i];
			rmsq += val * val;
		}
		rmsq = sqrt(rmsq / nt);
		if (rmsq) {
			for (i = 0; i < nt; ++i)
				data[i] /= rmsq;
		}
	}

	if (mbal) {
		int i;
		float mean = 0.0;

		/* mean = SUM (data[i] / nt) */
		for (i = 0; i < nt; ++i) {
			mean+=data[i];
		}

		/* compute the mean */
		mean/=nt;

		/* subtract the mean from each sample */
		if (mean) {
			for (i = 0; i < nt; ++i)
				data[i]-=mean;
		}
	}

	if (maxbal) {
		int i;
		float max = data[0];

		/* max */
		for (i = 0; i < nt; ++i) {
			if( data[i] > max ) max = data[i];
		}

		/* subtract max */
		for (i = 0; i < nt; ++i) data[i]-=max;
	}

	if (!CLOSETO(scale, f_one)) {
		int i;
		for (i = 0; i < nt; ++i)  data[i] *= scale;
	}
}

/* Zero out outliers */
static void do_trap(
	float *data,		/* the data			*/
	float trap,    /* zero if magnitude > trap     */
	int nt	 /* number of samples	    */
)
{
	float *dataptr = data;

	while (nt--) {
		if (ABS(*dataptr) > trap) *dataptr = 0.0;
		dataptr++;
	}
}

/* Hard clip outliers */
static void do_clip(
	float *data,		/* the data				*/
	float clip,    /* hard clip if magnitude > clip	*/
	int nt	 /* number of samples		    */
)
{
	float *dataptr = data;
	float mclip = -clip;

	while (nt--) {
		if (*dataptr > clip) {
			*dataptr = clip;
		} else if (*dataptr < mclip) {
			*dataptr = mclip;
		}
		dataptr++;
	}
}

/* Hard clip maxima */
static void do_pclip(
	float *data,		/* the data				*/
	float pclip,    /* hard clip if magnitude > clip	*/
	int nt	 /* number of samples		    */
)
{
	float *dataptr = data;

	while (nt--) {
		if (*dataptr > pclip) {
			*dataptr = pclip;
		}
		dataptr++;
	}
}

/* Hard clip minima */
static void do_nclip(
	float *data,		/* the data				*/
	float nclip,    /* hard clip if magnitude > clip	*/
	int nt	 /* number of samples		    */
)
{
	float *dataptr = data;

	while (nt--) {
		if (*dataptr < nclip) {
			*dataptr = nclip;
		}
		dataptr++;
	}
}

/* Quantile clip on magnitudes of trace values (absdata is scratch space for nt floats) */
static void do_qclip(
	float *data,	/* the data			*/
	float qclip,    /* quantile at which to clip    */
	int nt,	  /* number of sample points	*/
	float *absdata
)
{
	int i;
	int iq;			/* index of qclipth quantile    */
	float clip;		/* ... value of rank[iq]	*/

	iq = (int) (qclip * nt - 0.5); /* round, don't truncate */

	/* Clip on value corresponding to qth quantile */
	for (i = 0; i < nt; ++i)  absdata[i] = ABS(data[i]);
	clip = quant(absdata, iq, nt);
	do_clip(data, clip, nt);
}

/* Quantile balance (absdata is scratch space for nt floats) */
static void do_qbal(
	float *data,	/* the data			*/
	float qclip,    /* quantile at which to clip    */
	int nt,	  /* number of sample points	*/
	float *absdata
)
{
	int i;
	int iq;			/* index of qclipth quantile    */
	float bal;		/* value used to balance trace  */

	if (qclip == 1.0) { /* balance by max magnitude on trace */
		bal = ABS(data[0]);
		for (i = 1; i < nt; ++i)  bal = MAX(bal, ABS(data[i]));
		if (bal == 0.0) {
			return;
		} else {
			for (i = 0; i < nt; ++i)  data[i] /= bal;
			return;
		}
	}

	iq = (int) (qclip * nt - 0.5); /* round, don't truncate */

	/* Balance by quantile value (qclip < 1.0) */
	for (i = 0; i < nt; ++i)  absdata[i] = ABS(data[i]);
	bal = quant(absdata, iq, nt);
	if (bal == 0.0) {
		return;
	} else {
		for (i = 0; i < nt; ++i)  data[i] /= bal;
		do_clip(data, 1.0, nt);
		return;
	}
}

/* Automatic Gain Control--standard box (agcdata and d2 are scratch space for nt floats) */
static void do_agc(float *data, int iwagc, int nt, float *agcdata, float *d2)
{
	int i,j;
	float val;
	float sum;
	int nwin;
	float rms;

	/* Compute square of data */
	for (i = 0; i < nt; ++i) {
		val = data[i];
		d2[i] = val * val;
	}

	/* intialize first half window and gain first sample */
	sum = 0.0;
	for (i = 0; i < iwagc; ++i) {
		sum += d2[i];
	}
	nwin = iwagc;
	rms = sum/nwin;

	/* rms = 0 implies data[0]=0 */
	if (rms == 0) {
		agcdata[0]=0;
	} else {
		/* (The original has data[i] here, where i == iwagc is left over from the loop above.
		 * That is not the first sample, so it was a typo for data[0]) */
		agcdata[0]=data[0]/sqrt(rms);
	}

	/* ramping on : increase sum and nwin & gain data until reaching 2*iwagc-1 window */
	/* processing samples from 1 to iwagc-1 */
	for (i = 1; i < iwagc; ++i) {
		sum += d2[i+iwagc-1];
		++nwin;
		rms = sum/nwin;

		/* rms = 0 implies data[i]=0 */
		if (rms == 0) {
			agcdata[i]=0;
		} else {
			agcdata[i]=data[i]/sqrt(rms);
		}
	}

	/*  full 2*iagc rms window -- gain data */
	/* compute sum from 0 at each sample -- decreasing sum give inaccurate results, even negative RMS */
	/* processing samples from iwagc to nt-iwagc-1 */
	++nwin;
	for (i = iwagc; i < nt-iwagc; ++i) {
		sum=0;
		for (j = i-iwagc; j < i+iwagc; ++j) {
			sum +=d2[j];
		}
		rms = sum/nwin;
		if (rms == 0) {
			agcdata[i]=0;
		} else {
			agcdata[i]=data[i]/sqrt(rms);
		}
	}

	/* ramping off -- decrease nwin -- gain data */
	/* compute sum from 0 at each sample -- decreasing sum give inaccurate results, even negative RMS */
	/* processing samples from nt-iwagc to nt-1 */
	for (i = nt-iwagc; i < nt; ++i) {
		sum=0;
		for (j = i-iwagc; j < nt; ++j) {
			sum +=d2[j];
		}
		--nwin;
		rms = sum/nwin;
		if (data[i] == 0) {
			agcdata[i]=0;
		} else {
			agcdata[i]=data[i]/sqrt(rms);
		}
	}

	/* copy data back into trace */
	memcpy( (void *) data, (const void *) agcdata, nt*FSIZE);
	return;
}

/* Automatic Gain Control--gaussian taper
 * (agcdata, d2 and s are scratch space for nt floats, w for iwagc floats) */
static void do_gagc(float *data, int iwagc, int nt, float *agcdata, float *w, float *d2, float *s)
{
	float u;		/* related to reciprocal of std dev     */
	float usq;		/* u*u				  */

	/* Compute Gaussian window weights */
	u = EPS / ((float) iwagc);
	usq = u*u;
	{
		int i;
		float floati;
		for (i = 1; i < iwagc; ++i) {
			floati = (float) i;
			w[i] = exp(-(usq*floati*floati));
		}
	}

	/* Agc the trace */
	{
		int i, j, k;
		float val;
		float wtmp;
		float stmp;

		/* Put sum of squares of data in d2 and */
		/* initialize s to d2 to get center point set */
		for (i = 0; i < nt; ++i) {
			val = data[i];
			s[i] = d2[i] = val * val;
		}

		/* Compute weighted sum s; use symmetry of Gaussian */
		for (j = 1; j < iwagc; ++j) {
			wtmp = w[j];
			for (i = j; i < nt; ++i)  s[i] += wtmp*d2[i-j];
			k = nt - j;
			for (i = 0; i < k; ++i)   s[i] += wtmp*d2[i+j];
		}

		for (i = 0; i < nt; ++i) {
			stmp = s[i];
			agcdata[i] = (!stmp) ? 0.0 : data[i]/sqrt(stmp);
		}

		/* Copy data back into trace */
		memcpy( (void *) data, (const void *) agcdata, nt*FSIZE);
	}
	return;
}

/*
 * QUANT - find k/n th quantile of a[]
 *
 * Works by reordering a so a[j] < a[k] if j < k.
 *
 * Parameters:
 *    a	 - data
 *    k	 - indicates quantile
 *    n	 - number of points in data
 *
 * This is Hoare's algorithm worked over by SEP (#10, p100) and Brian.
 */
static float quant(float *a, int k, int n)
{
	int i, j;
	int low, hi;
	float ak, aa;

	low = 0; hi = n-1;
	while (low < hi) {
		ak = a[k];
		i = low;
		j = hi;
		do {
			while (a[i] < ak) i++;
			while (a[j] > ak) j--;
			if (i <= j) {
				aa = a[i]; a[i] = a[j]; a[j] = aa;
				i++;
				j--;
			}
		} while (i <= j);
		if (j < k) low = i;
		if (k < i) hi = j;
	}
	return(a[k]);
}
