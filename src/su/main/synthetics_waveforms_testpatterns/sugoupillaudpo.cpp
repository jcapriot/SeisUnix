#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									  ",
" SUGOUPILLAUDPO - calculate Primaries-Only impulse response of a lossless",
"	      GOUPILLAUD medium for plane waves at normal incidence	",
"									  ",
" sugoupillaudpo < stdin > stdout [optional parameters]		      ",
"									    ",
" Required parameters:							     ",
"	none								",
"      									     ",
" Optional parameters:						       ",
"	l=1	   source layer number; 1 <= l <= tr.ns		  ",
"		      Source is located at the top of layer l.		     ",
"	k=1	   receiver layer number; 1 <= k			 ",
"		      Receiver is located at the top of layer k.	    ",
"	tmax	  number of output time-samples;			",
"		      default: long enough to capture all primaries	 ",
"	pV=1	  flag for vector field seismogram		      ",
"		      (displacement, velocity, acceleration);	       ",
"		      =-1 for pressure seismogram.			  ",
"	verbose=0     silent operation, =1 list warnings		    ",
"									     ",
" Input: Reflection coefficient series:				      ",
"									    ",
"			       impedance[i]-impedance[i+1]		   ",
"		       r[i] = -----------------------------		  ",
"			       impedance[i]+impedance[i+1]		   ",
"									    ",
"	r[0]= surface refl. coef. (as seen from above)		      ",
"	r[n]= refl. coef. of the deepest interface			  ",
"									     ",
" Input file is to be in SU format, i.e., binary floats with a SU header.    ",
"									    ",
" Remarks:								   ",
" 1. For vector fields, a buried source produces a spike of amplitude 1      ",
" propagating downwards and a spike of amplitude -1 propagating upwards.     ",
" A buried pressure source produces spikes of amplitude 1 both in the up-    ",
" and downward directions.						   ",
"    A surface source induces only a downgoing spike of amplitude 1 at the   ",
" top of the first layer (both for vector and pressure fields).	      ",
" 2. The sampling interval dt in the header of the input reflectivity file   ",
" is interpreted as a two-way traveltime thicknes of the layers. The sampling",
" interval of the output seismogram is the same as that of the input file.   ",
NULL};
#endif


/* 
 * Credits:
 *	CWP: Albena Mateeva, April 2001.
 *
 */

/**************** end self doc ***********************************/

/* Library version of SUGOUPILLAUDPO
 *
 * The main program is not part of the library. It reads the reflection coefficients of the first trace of its input, and writes
 * one trace, the seismogram. Doing that is left to the caller, and what is here is what is worked out from the reflectivity:
 *
 *	su_goupillaudpo_tmax:	the default length of the seismogram, long enough for all of the primaries
 *	su_goupillaudpo:	the seismogram
 *
 * All the arrays are the caller's, none are allocated here.
 *
 * Differences from the program:
 *  - The reflectivity is not changed (the program multiplies it by pV in place).
 *  - The receiver layer must be at most n+1 (the program reads beyond the reflectivity for a deeper one).
 */

static int gpo_imin(int c1, int c2) { return c1 <= c2 ? c1 : c2; }
static int gpo_imax(int c1, int c2) { return c1 >= c2 ? c1 : c2; }

/* The default number of output samples, for n interfaces (n+1 reflection coefficients) and the source and receiver layers l, k */
int su_goupillaudpo_tmax(int n, int l, int k)
{
	int n1 = gpo_imin(k,l);
	int n2 = gpo_imax(k,l);

	return (n2-n1+2+2*gpo_imax(gpo_imax(n1-1,n+1-n2),n+1-n1))/2;
}

/* The primaries-only impulse response of a lossless Goupillaud medium for plane waves at normal incidence.
 *
 * n: number of subsurface interfaces; r[n+1]: the reflection coefficients, r[0] that of the surface
 * l, k: the source and the receiver layers (1 to n+1; the layer 1 is at the surface)
 * tmax: number of output samples
 * pV: 1 for a vector field (displacement, velocity, acceleration), -1 for pressure
 * x[2*tmax]: scratch. out[tmax]: the seismogram.
 * *odd: set to 1 if the seismogram is shifted by half a sample (k-l is odd), otherwise 0
 *
 * Returns 0, or -1 if a reflection coefficient is not from -1 to 1, -2 if the seismogram is too short to see any signal,
 * -3 for parameters that are not allowed.
 */
int su_goupillaudpo(int n, const float *r, int l, int k, int tmax, int pV, float *x, float *out, int *odd)
{
	int i, rmax, skl, n1, n2;
	float effsd = 1;	/* effective source strength in downward direction */
	float transm1 = 1;	/* one-way transmission coefficient from source to receiver */
	float transm2 = 1;	/* two-way transmission coefficient between reflector and receiver */

	if (n < 0 || k < 1 || l < 1 || l > n+1 || k > n+1 || tmax < 0 || !(pV == 1 || pV == -1)) return -3;
	for (i=0; i<=n; ++i)
		if (r[i] > 1. || r[i] < -1.) return -1;

	if (k>l) skl=1;
	else if (k<l) skl=-1;
	else skl=0;
	n1 = gpo_imin(k,l);
	n2 = gpo_imax(k,l);

	/* Trivial case */
	if (n2-n1 >= 2*tmax) return -2;

	/* Initial zeroing of the seismogram */
	memset((void *) x, 0, 2*tmax*FSIZE);

	/* one-way transmission coef. between source and receiver */
	for (i=n1; i<n2; ++i) transm1 *= 1+skl*(pV*r[i]);

	/* direct arrival */
	if (pV == -1)    x[n2-n1] = transm1;	/* r[i]=pV*r[i] already */
	else if (l == 1) x[n2-n1] = transm1;	/* downgoing 1 for pV=1 */
	else		 x[n2-n1] = skl*transm1;	/* upgoing -1  for pV=1 */

	/* effective source in downward direction */
	if (!(l==1)) effsd = 1+pV*(pV*r[l-1]);

	/* index of the deepest observable reflector */
	i = (n2+n1-1)/2;
	rmax = gpo_imin(n, tmax-1+i);

	/* primary reflections from below the receiver (and below the source if l>k) */
	for (i=n2; i<=rmax; ++i) {
		x[2*(i+1)-n2-n1] = effsd*transm1*(pV*r[i])*transm2;
		transm2 *= 1-(pV*r[i])*(pV*r[i]);
	}

	/* primary reflections from above the receiver (and above the source if l<k) */
	if (!(l==1))
		for (i=n1-1, transm2=1; i>=0 && n2-n1+2*(n1-1-i)<2*tmax; --i) {
			x[n2-n1+2*(n1-1-i)] += pV*transm1*(pV*r[i])*transm2;
			transm2 *= 1-(pV*r[i])*(pV*r[i]);
		}

	/* the seismogram */
	i = (k-l)/2;
	if (2*i == k-l) {
		for (i=0; i<tmax; ++i) out[i] = x[2*i];
		*odd = 0;
	} else {
		for (i=0; i<tmax; ++i) out[i] = x[2*i+1];
		*odd = 1;
	}
	return 0;
}
