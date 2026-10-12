/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUEIPOFI: $Revision: 1.6 $ ; $Date: 2011/11/16 22:58:31 $	*/

#include "su.h"
#include <vector>
#include "segy.h"
#include "header.h"

/*********************** self documentation *****************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                       ",
" SUEIPOFI - EIgenimage (SVD) based POlarization FIlter for             ",
"            three-component data                                       ",
"                                                                       ",
" sueipofi <stdin >stdout [optional parameters]                         ",
"                                                                       ",
" Required parameters:                                                  ",
"    none                                                               ",
"                                                                       ",
" Optional parameters:                                                  ",
"    dt=(from header)  time sampling intervall in seconds               ",
"    wl=0.1            SVD time window length in seconds                ",
"    pwr=1.0           exponent of filter weights                       ",
"    interp=cubic      interpolation between initially calculated       ",
"                      weights, choose \"cubic\" or \"linear\"          ",
"    verbose=0         1 = echo additional information                  ",
"                                                                       ",
"    file=polar        base name for additional output file(s) of       ",
"                      filter weights (see flags below)                 ",
"    rl1=0             1 = rectilinearity along first principal axis    ",
"    rl2=0             1 = rectilinearity along second principal axis   ",
"    pln=0             1 = planarity                                    ",
"                                                                       ",
"                                                                       ",
" Notes:                                                                ",
"    Three adjacent traces are considered as one three-component        ",
"    dataset.                                                           ",
"                                                                       ",
"    The filter is the sum of the first two eigenimages of the singular ",
"    value decomposition (SVD) of the signal matrix (time window).      ",
"    Weighting functions depending on linearity and planarity of the    ",
"    signal are applied, additionally. To avoid edge effects, these are ",
"    interpolated linearily or via cubic splines between initially      ",
"    calculated values of non-overlapping time windows.                 ",
"    The algorithm is based on the assumption that the particle motion  ",
"    trajectory is essentially 2D (elliptical polarization).            ",
"                                                                       ",
" Caveat:                                                               ",
"    Cubic spline interpolation may result in filter weights exceeding  ",
"    the set of values of initial weights. Weights outside the valid    ",
"    interval [0.0, 1.0] are clipped.                                   ",
"                                                                       ",
NULL};
#endif

/* 
 * Author: Nils Maercklin, 
 *         GeoForschungsZentrum (GFZ) Potsdam, Germany, 2001.
 *         E-mail: nils@gfz-potsdam.de
 *
 *
 * References:
 *    Franco, R. de, and Musacchio, G., 2000: Polarization Filter with
 *       Singular Value Decomposition, submitted to Geophysics and
 *       published electronically in Geophysics online (www.geo-online.org).
 *    Jurkevics, A., 1988: Polarization analysis of three-comomponent
 *       array data, Bulletin of the Seismological Society of America, 
 *       vol. 78, no. 5.
 *    Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P.
 *       1996: Numerical Recipes in C - The Art of Scientific Computing,
 *       Cambridge University Press, Cambridge.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 */
/**************** end self doc *******************************************/



/* rectilinearity, calculated from singular values (after Jurkevics, 1988, Franco and Musacchio, 2000)
   w[3] = array containing singular values in descending order */
static float calc_svdrl(float *w, int axis, float pwr)
{
	float rl;
	if (w[axis]*w[axis]) {
		rl = 1.0 - (w[2]*w[2]) / (w[axis]*w[axis]);
	}
	else {
		rl = 0.0;
	}
	return pow(rl,pwr);
}

/* planarity, calculated from singular values (Jurkevics, 1988) */
static float calc_planarity(float *w, float pwr)
{
	float p;
	if (w[0] + w[1]) {
		p = 1.0 - (2.0 * w[2]*w[2]) / (w[0]*w[0] + w[1]*w[1]);
	}
	else {
		p = 0.0;
	}
	return pow(p,pwr);
}

/* the values of the attribute iny at the centers of the windows (inx0 + ix din), at all of the samples of the trace, by
   cubic splines (the weights are clipped to 0 to 1 by the caller) or by linear interpolation (the program reads
   iny[nin] for the value above the last center: iny[nin-1] here) */
static void interpolate_weights(const float *iny, float *outy, int nin, int din, int inx0, int nout, int cubic)
{
	int ix;
	std::vector<float> inx(nin), outx(nout), y(iny, iny+nin);
	std::vector<float> yd((size_t) 4*nin);
	for (ix=0; ix<nin; ix++) inx[ix] = (float) (inx0 + ix * din);
	for (ix=0; ix<nout; ix++) outx[ix] = (float) ix;
	if (cubic) {
		csplin(nin,inx.data(),y.data(),(float (*)[4]) yd.data());
		intcub(0,nin,inx.data(),(float (*)[4]) yd.data(),nout,outx.data(),outy);
	} else {
		intlin(nin,inx.data(),y.data(),y[0],y[nin-1],nout,outx.data(),outy);
	}
}

/* The filter of one three-component station of the program.
 *
 * nt, iwl		samples of a trace, and the length of the (non overlapping) windows of the SVD: iwl >= 3
 * pwr			the exponent of the polarization attributes (the filter weights)
 * cubic		cubic spline (1) or linear (0) interpolation of the weights between the centers of the windows
 * data[3][nt]		the three components, which are replaced by the filtered ones in the first (nt / iwl) * iwl samples
 *			(the rest are left as they are)
 * polar[3][nt]		the weights: rectilinearity along the first and second axes, and planarity
 *
 * Needs nt / iwl >= 1 windows, and 2 for the cubic spline.
 */
void su_eipofi_station(int nt, int iwl, float pwr, int cubic, float *data, float *polar)
{
	int i, j, it, jt, iw;
	const int nw = nt / iwl;
	std::vector<float*> u(iwl), v(3);
	std::vector<float> ubuf((size_t) iwl*3), vbuf(9), sval(3);
	std::vector<float> eigen((size_t) 3*3*nt, 0.0f), initpa((size_t) 3*nw);

	for (i=0; i<iwl; i++) u[i] = &ubuf[(size_t) i*3];
	for (i=0; i<3; i++) v[i] = &vbuf[(size_t) i*3];

	for (iw=0,jt=0; iw<nw; iw++,jt+=iwl) {
		/* copy data into matrix for SVD */
		for (i=0; i<3; i++)
			for (it=0; it<iwl; it++)
				u[it][i] = data[(size_t) i*nt + it+jt];

		/* a window of zeros has no singular values (compute_svd gives it undefined ones) */
		bool nonzero = false;
		for (i=0; i<3 && !nonzero; i++)
			for (it=0; it<iwl; it++)
				if (u[it][i] != 0.0f) { nonzero = true; break; }
		if (nonzero) {
			/* perform singular value decomposition (SVD), and sort singular values in descending order */
			compute_svd(u.data(), iwl, 3, sval.data(), v.data());
			svd_sort(u.data(), sval.data(), v.data(), 3, iwl);

			/* construct and save eigenimages */
			for (j=0; j<3; j++)
				for (i=0; i<3; i++)
					for (it=jt; it<jt+iwl; it++)
						eigen[((size_t) j*3 + i)*nt + it] = u[it-jt][j] * v[i][j] * sval[j];
		} else {
			sval[0] = sval[1] = sval[2] = 0.0f;
		}

		/* calculate polarization attributes */
		for (j=0; j<2; j++) initpa[(size_t) j*nw + iw] = calc_svdrl(sval.data(),j,pwr);
		initpa[(size_t) 2*nw + iw] = calc_planarity(sval.data(),pwr);
	}

	/* interpolation of all polarization attributes: the initial values are at the centers of the windows */
	for (i=0; i<3; i++) {
		float *pol = polar + (size_t) i*nt;
		interpolate_weights(&initpa[(size_t) i*nw], pol, nw, iwl, iwl/2, nt, cubic);
		if (cubic) {
			for (it=0; it<nt; it++) {
				if (pol[it]>1.0) pol[it] = 1.0;
				if (pol[it]<0.0) pol[it] = 0.0;
			}
		}
	}

	/* filter: polarization attribute weighted eigenimage stack */
	for (i=0; i<3; i++) {
		for (it=0; it<nw*iwl; it++) {
			float x = eigen[((size_t) 0*3 + i)*nt + it]*polar[(size_t) 0*nt + it] +
				eigen[((size_t) 1*3 + i)*nt + it]*polar[(size_t) 1*nt + it];
			data[(size_t) i*nt + it] = x * polar[(size_t) 2*nt + it];
		}
	}
}
