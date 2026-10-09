/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSYNCZ: $Revision: 1.16 $ ; $Date: 2011/11/12 00:40:42 $	*/

#include "su.h"
#include "segy.h" 

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUSYNCZ - SYNthetic seismograms for piecewise constant V(Z) function	",
"	   True amplitude (primaries only) modeling for 2.5D		",
" 									",
"  susyncz > outfile [parameters]					",
" 									",
" Required parameters:							",
" none									",
" 									",
" Optional Parameters:							",
" ninf=4        number of interfaces (not including upper surface)	",
" dip=5*i       dips of interfaces in degrees (i=1,2,3,4)		",
" zint=100*i    z-intercepts of interfaces at x=0 (i=1,2,3,4)		",
" v=1500+ 500*i velocities below surface & interfaces (i=0,1,2,3,4)	",
" rho=1,1,1,1,1 densities below surface & interfaces (i=0,1,2,3,4)	",
" nline=1	number of (identical) lines				",
" ntr=32        number of traces					",
" dx=10         trace interval						",
" tdelay=0      delay in recording time after source initiation		",
" dt=0.004      time interval						",
" nt=128        number of time samples					",
" 			 						",
" Notes:								",
" The original purpose of this code was to create some nontrivial	",
" data for Brian Sumner's CZ suite.					",
" 			 						",
" The program produces zero-offset data over dipping reflectors.	",
" 			 						",
" In the original fortran code, some arrays had the index		",
" interval 1:ninf, as a natural way to index over the subsurface	",
" reflectors.  This indexing was preserved in this C translation.	",
" Consequently, some arrays in the code do not use the 0 \"slot\".	",
" 			 						",
" Example:								",
"	susyncz | sufilter | sugain tpow=1 | display_program		",
" 			 						",
" Trace header fields set: tracl, ns, dt, delrt, ntr, sx, gx		",
NULL};
#endif


/*
 * Credits:
 * 	CWP: Brian Sumner, 1983, 1985, Fortran design and code 
 *      CWP: Stockwell & Cohen, 1995, translation to C 
 *
 */

/**************** end self doc *******************************************/

/* Library version of SUSYNCZ
 *
 * The main program is not part of the library. It reads the parameters, and writes the traces with their headers. Doing that is
 * left to the caller, and what is here is the work on the model:
 *
 *	su_syncz_tables:	the tables for a model (once)
 *	su_syncz_trace:		the zero-offset trace at one position
 *
 * All the arrays are the caller's, none are allocated here. The tables with two indices (theta[j][i], m, b, d, k, w) are
 * arrays of (ninf+1)*(ninf+1) numbers, the row j after the row j-1; those with one index have ninf+1 numbers. The index 0 is
 * the surface, and the interfaces are 1 to ninf.
 *
 * Differences from the program:
 *  - The table theta is made for all of the interfaces. The program leaves the deepest one out (it makes it for j < ninf, and
 *    uses it for j <= ninf), so it works with numbers that it never set.
 *  - The velocities and densities of all of the layers are checked (the program leaves out the deepest).
 *  - The check that the interfaces do not meet in the region that is observed goes through all of the pairs of interfaces. The
 *    program has the wrong counter in the loop that follows the checks of the first interface (j instead of j2), which
 *    stops the loop after it, and works with the table m beyond what is set.
 */

#define SYNCZ_TORADS (PI/180.0)
#define SYNCZ_FOURPI (4.0*PI)
#define SYNCZ_AT(a, n, j, i) ((a)[(size_t) (j) * ((n) + 1) + (i)])

/* The tables of a model of ninf interfaces under the surface, for ntr traces dx apart.
 *
 * zint[ninf+1], dip[ninf+1]: the intercepts of the interfaces with x = 0, and their dips in degrees (the first of each, the surface,
 *	is 0). v[ninf+1], rho[ninf+1]: the velocities and densities of the layers (the first is that above the first interface)
 * dipr[ninf+1], xl[ninf+1], xr[ninf+1]: scratch
 * theta, m, b, d, k, w: the tables [ (ninf+1)*(ninf+1) ]. meet[ninf+1], trcoefs[ninf+1]: the tables. dorow[ninf+1]: 1 where the
 *	interfaces i and i+1 meet, 0 where they are parallel
 *
 * Returns 0, or the number of the error of the program: 1 a dip that is not from -90 to 90 degrees, 2 a change of dip of more
 * than 90 degrees from an interface to the next, 3 intercepts that do not go down, 4 a velocity that is not positive, 5 a density
 * that is not positive, 6 no specular ray to an interface, 7 interfaces that meet in the region that is observed. *info is the
 * interface that it is about.
 */
int su_syncz_tables(int ninf, int ntr, float dx, const float *zint, const float *dip, const float *v, const float *rho,
	float *dipr, float *xl, float *xr, float *theta, float *m, float *b, float *d, float *k, float *w, float *meet,
	float *trcoefs, int *dorow, int *info)
{
	int i, j, j2;
	float xmin, xmax;

	/* check dips and intercepts */
	dipr[0] = 0.0;
	for (i = 1; i <= ninf; ++i) {
		*info = i;
		if (dip[i] < -90.0 || dip[i] > 90.0) return 1;
		if (ABS(dip[i] - dip[i-1]) > 90.0) return 2;
		dipr[i] = dip[i] * SYNCZ_TORADS;
		if (zint[i] <= zint[i-1]) return 3;
	}
	for (i = 0; i <= ninf; ++i) {
		*info = i;
		if (v[i] <= 0.0) return 4;
		if (rho[i] <= 0.0) return 5;
	}

	/* Table theta[j][i], the angle at which the specular to layer j leaves interface i. */
	for (j = 1; j <= ninf; ++j) {
		SYNCZ_AT(theta, ninf, j, j-1) = -dipr[j];
		for (i = j-1; i > 0; --i) {
			float temp;
			temp = v[i-1]*sin(SYNCZ_AT(theta, ninf, j, i)+dipr[i])/v[i];
			if (ABS(temp) > 1.0) { *info = j; return 6; }
			SYNCZ_AT(theta, ninf, j, i-1) = asin(temp) - dipr[i];
		}
	}

	/* Table m and b, which are used to find the x coordinate of the specular at the next interface from its x coordinate at
	 * the previous interface. */
	for (j = 1; j <= ninf; ++j) {
		for (i = 0; i < j; ++i) {
			float s, temp;
			s = sin(SYNCZ_AT(theta, ninf, j, i));
			temp = 1.0 - tan(dipr[i+1])*s;
			SYNCZ_AT(m, ninf, j, i) = (1.0 - tan(dipr[i])*s)/temp;
			SYNCZ_AT(b, ninf, j, i) = (zint[i+1]-zint[i])*s/temp;
		}
	}

	/* Make the final checks on the substructure specification. The strategy is to check the x coordinates of the interface
	 * intercepts against x coordinates of specular-interface intercepts. Only the rays from the left and right endpoints need
	 * to be checked, since they define the area being "observed". */
	xmin = 0.0;
	xmax = dx * (ntr - 1);
	for (j = 1; j <= ninf; ++j) {
		xl[j] = xmin;
		xr[j] = xmax;
	}
	for (j=0; j<ninf; ++j) {
		for (j2=j+1; j2<=ninf; ++j2) {
			if (dipr[j2] != dipr[j]) {
				float intercept;
				intercept = (zint[j2]-zint[j])/(tan(dipr[j])-tan(dipr[j2]));
				if (intercept > xmin && intercept < xmax) { *info = j2; return 7; }
			}
		}
		for (j2=j+1; j2<=ninf; ++j2) {
			xl[j2] = SYNCZ_AT(m, ninf, j2, j)*xl[j2] + SYNCZ_AT(b, ninf, j2, j);
			xr[j2] = SYNCZ_AT(m, ninf, j2, j)*xr[j2] + SYNCZ_AT(b, ninf, j2, j);
			if (xl[j2] < xmin) xmin = xl[j2];
			if (xr[j2] > xmax) xmax = xr[j2];
		}
	}

	/* Table the arrays meet and k if two adjacent interfaces intersect, otherwise d if they are parallel. Meet is the x
	 * coordinate of intersection if there is an intersection, k is a value which will be used to find the distance travelled
	 * in that layer. If there is no intersection d will be a constant for all j and can be stored now. */
	for (i = 0; i < ninf; ++i) {
		if (dipr[i+1] == dipr[i]) {
			dorow[i] = 0;
			for (j = i + 1; j <= ninf; ++j)
				SYNCZ_AT(d, ninf, j, i) = (zint[i+1]-zint[i])*cos(dipr[i]) / cos(SYNCZ_AT(theta, ninf, j, i) + dipr[i]);
		} else {
			dorow[i] = 1;
			meet[i] = (zint[i+1]-zint[i])/(tan(dipr[i])-tan(dipr[i+1]));
			for (j = i + 1; j <= ninf; ++j)
				SYNCZ_AT(k, ninf, j, i) = sin(dipr[i]-dipr[i+1]) /
					(cos(SYNCZ_AT(theta, ninf, j, i)+dipr[i+1])*cos(dipr[i]));
		}
	}

	/* Table trcoefs, the transmission coefficients, as the product of the reflection coefficient of the interface and the
	 * transmission coefficients (1 - r^2) of the interfaces above it. */
	for (j = 1; j <= ninf; ++j) {
		float t1, t2, p;
		t1 = rho[j]*v[j];
		t2 = rho[j-1]*v[j-1];
		p = (t1 - t2)/(t1 + t2);
		for (i = 1; i < j; ++i) {
			float r;
			t1 = rho[i]*v[i]*cos(SYNCZ_AT(theta, ninf, j, i-1)+dipr[i]);
			t2 = rho[i-1]*v[i-1]*cos(SYNCZ_AT(theta, ninf, j, i)+dipr[i]);
			r = (t1 - t2)/(t1 + t2);
			p *= 1.0 - r*r;
		}
		trcoefs[j] = p;
	}

	/* Table w, the coefficients for the spreading factor calculation. */
	for (j = 1; j <= ninf; ++j) {
		SYNCZ_AT(w, ninf, j, 0) = 1.0;
		for (i = 1; i < j; ++i) {
			float c1, c2;
			c1 = cos(SYNCZ_AT(theta, ninf, j, i-1) + dipr[i]);
			c2 = cos(SYNCZ_AT(theta, ninf, j, i) + dipr[i]);
			SYNCZ_AT(w, ninf, j, i) = SYNCZ_AT(w, ninf, j, i-1)*c1*c1/(c2*c2);
		}
	}
	return 0;
}

/* The trace at the position x, which is the sum of the spike responses of the interfaces, from the tables of su_syncz_tables.
 *
 * nt, dt, tdelay: the number of samples, the sample interval, and the delay of the recording time
 * tout[nt]: the times of the samples, it*dt
 * data[nt]: scratch. trace[nt]: the trace
 */
void su_syncz_trace(int ninf, float x, int nt, float dt, float tdelay, const float *v, const float *m, const float *b,
	const float *d, const float *k, const float *w, const float *meet, const float *trcoefs, const int *dorow,
	const float *tout, float *data, float *trace)
{
	int i, j, it;

	memset((void *) trace, 0, nt * FSIZE);
	for (j = 1; j <= ninf; ++j) {
		float ex, t0, time, tmp, p1, p2, dist, amp[1];

		ex = x;
		if (dorow[0]) dist = (meet[0]-ex)*SYNCZ_AT(k, ninf, j, 0);
		else          dist = SYNCZ_AT(d, ninf, j, 0);
		t0 = dist/v[0];
		tmp = dist*v[0];
		p1 = tmp;
		p2 = tmp;
		for (i = 1; i < j; ++i) {
			ex = SYNCZ_AT(m, ninf, j, i-1)*ex + SYNCZ_AT(b, ninf, j, i-1);
			if (dorow[i]) dist = (meet[i] - ex)*SYNCZ_AT(k, ninf, j, i);
			else          dist = SYNCZ_AT(d, ninf, j, i);
			t0 += dist/v[i];
			tmp = dist*v[i];
			p1 += tmp;
			p2 += tmp*SYNCZ_AT(w, ninf, j, i);
		}

		/* set strength and time of spike response */
		amp[0] = v[0]*trcoefs[j]/(dt*SYNCZ_FOURPI*2.0*sqrt(p1*p2));
		time = 2.0*t0 - tdelay;

		/* distribute spike response over full trace, and add it to the trace */
		ints8r(1,dt,time,amp,0,0,nt,(float *) tout,data);
		for (it = 0; it < nt; ++it)
			trace[it] += data[it];
	}
}
