/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUVELAN_UCCS : $Revision: 1.7 $ ; $Date: 2011/11/16 23:40:27 $	*/

#include "su.h"
#include <vector>
#include <algorithm>
#include "segy.h"
#include "string.h"

/*********************** self documentation **********************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									     ",
" SUVELAN_UCCS - compute stacking VELocity panel for cdp gathers	     ",
"		using UnNormalized CrossCorrelation Sum 	             ",
"									     ",
" suvelan_uccs <stdin >stdout [optional parameters]			     ",
"									     ",
" Optional Parameters:							     ",
" nx=tr.cdpt              number of traces in cdp 			     ",
" nv=50                   number of velocities				     ",
" dv=50.0                velocity sampling interval			     ",
" fv=1500.0               first velocity				     ",
" smute=1.5               samples with NMO stretch exceeding smute are zeroed",
" dtratio=5               ratio of output to input time sampling intervals   ",
" nsmooth=dtratio*2+1     length of smoothing window                         ",
" verbose=0               =1 for diagnostic print on stderr		     ",
" pwr=1.0                 semblance value to the power      		     ",
"									     ",
"Notes:									     ",
" Unnormalized crosscorrelation sum: sum all possible crosscorrelation trace ",
" pairs in a CMP gather for each trial velocity and zero-offset two-way      ",
" travel time inside a time window. This unnormalized coherency measure      ",
" produces large spectral amplitudes for strong reflections and small        ",
" spectral amplitudes for weaker ones. If M is the number of traces in the   ",
" CMP gather M(M-1)/2 is the total number of crosscorrelations for each trial",
" velocity and zero-offset two-way traveltime.			 	     ",
NULL};
#endif

/* 
 * Credits: CWP: Valmore Celis, Sept 2002	
 * 
 * Based on the original code: suvelan.c 
 *    Colorado School of Mines:  Dave Hale c. 1989
 *
 *
 * Reference: Neidell, N.S., and Taner, M.T., 1971, Semblance and other 
 *            coherency measures for multichannel data: Geophysics, 36, 498-509.
 *
 * Trace header fields accessed:  ns, dt, delrt, offset, cdp
 * Trace header fields modified:  ns, dt, offset, cdp
 */
/**************** end self doc *******************************************/

/* The four programs suvelan_uccs, suvelan_nccs, suvelan_usel and suvelan_nsel: the stacking velocity panel of a CDP
 * gather from the sum of the cross-correlations of the pairs of its traces, after NMO, in windows of nsmooth samples:
 *
 *	uccs	unnormalized sum of all of the pairs (selective = 0, normalize = 0)
 *	nccs	each pair divided by the geometric mean of the energies of its traces (+ 1) in the window, and the sum
 *		by the number of pairs: 2 / (M (M - 1)) (selective = 0, normalize = 1)
 *	usel	unnormalized sum of the pairs of a large enough differential moveout (selective = 1, normalize = 0)
 *	nsel	the same, normalized (selective = 1, normalize = 1)
 *
 * nt, dt, ft	samples, sample interval and first time of the traces
 * nx		number of traces; data[nx][nt]; offsets[nx]
 * nv, dv, fv	number, increment and first of the velocities
 * smute	samples with NMO stretch exceeding smute are zeroed
 * dtratio, nsmooth	ratio of output to input sampling intervals, and the length of the window
 * tau		the threshold of the selection, from 0 (all pairs) to 1 (only the pair of the shortest and the
 *		longest offset): a pair is used if the difference of the squares of its offsets (the parabolic
 *		differential moveout) is at least tau times that of the whole gather
 * pwr		each value is taken to this power (with its sign)
 * sem[nv][ntout]	the panel; ntout = 1 + (nt - 1) / dtratio
 *
 * The programs output one trace (the one of the last velocity) for each gather, from a gather that is overwritten
 * by the next one, with the velocities before it left in the NMO of the muted samples, and a selection from index
 * formulas that need regularly spaced offsets sorted by the (unsigned) offset: it is made here from the offsets. The
 * programs also add the squares (the pairs of a trace with itself) to the sum, which their documentation does not
 * have, and do not use pwr.
 */
void su_velan_cc(int nt, float dt, float ft, int nx, const float *offsets, const float *data, int nv, float dv, float fv,
	float smute, int dtratio, int nsmooth, int normalize, int selective, float tau, float pwr, float *sem)
{
	const int ntout = 1 + (nt-1)/dtratio;
	int iv, ix, jx, it, itout, is, ismin, ismax, itmute, iti;
	float v, offovs, tnmute, tn, ti, frac, ff, sum, range = 0.0, a, b;
	std::vector<int> order(nx);
	std::vector<float> aoff(nx), datn((size_t) nx*nt), norm(nx);
	std::vector<char> pair((size_t) nx*nx, 0);
	double beta = 0.0;

	/* the traces by increasing offset */
	for (ix=0; ix<nx; ++ix) order[ix] = ix;
	std::stable_sort(order.begin(), order.end(), [&](int p, int q) { return fabs(offsets[p]) < fabs(offsets[q]); });
	for (ix=0; ix<nx; ++ix) aoff[ix] = offsets[order[ix]];

	/* the pairs that are summed */
	if (nx > 1) {
		range = aoff[nx-1]*aoff[nx-1] - aoff[0]*aoff[0];
		for (jx=1; jx<nx; ++jx) {
			for (ix=0; ix<jx; ++ix) {
				bool use = !selective || (aoff[jx]*aoff[jx] - aoff[ix]*aoff[ix] >= tau*range);
				pair[(size_t) jx*nx + ix] = use;
				if (use) beta += 1.0;
			}
		}
	}

	for (iv=0,v=fv; iv<nv; ++iv,v+=dv) {
		float *out = sem + (size_t) iv*ntout;
		std::fill(datn.begin(), datn.end(), 0.0f);
		for (ix=0; ix<nx; ++ix) {
			const float *x = data + (size_t) order[ix]*nt;
			float *y = &datn[(size_t) ix*nt];

			/* compute offset/velocity squared */
			offovs = (aoff[ix]*aoff[ix])/(v*v);

			/* determine mute time after nmo */
			tnmute = sqrt(offovs/(smute*smute-1.0));
			itmute = (tnmute > ft) ? (int) ((tnmute-ft)/dt) : 0;

			/* do nmo via quick and dirty linear interpolation (accurate enough for velocity analysis) */
			for (it=itmute,tn=ft+itmute*dt; it<nt; ++it,tn+=dt) {
				ti = (sqrt(tn*tn+offovs)-ft)/dt;
				iti = ti;
				if (iti>=0 && iti<nt-1) {
					frac = ti-iti;
					y[it] = (1.0-frac)*x[iti]+frac*x[iti+1];
				}
			}
		}
		for (itout=0; itout<ntout; ++itout) {
			it = itout*dtratio;
			ismin = it-nsmooth/2;
			ismax = it+nsmooth/2;
			if (ismin<0) ismin = 0;
			if (ismax>nt-1) ismax = nt-1;
			for (ix=0; ix<nx; ++ix) {
				sum = 0.0;
				for (is=ismin; is<ismax; ++is) sum += datn[(size_t) ix*nt+is]*datn[(size_t) ix*nt+is];
				norm[ix] = normalize ? sqrt(sum)+1 : 1.0;
			}
			ff = 0.0;
			for (is=ismin; is<ismax; ++is) {
				for (jx=1; jx<nx; ++jx) {
					b = datn[(size_t) jx*nt+is]/norm[jx];
					for (ix=0; ix<jx; ++ix) {
						if (pair[(size_t) jx*nx + ix]) {
							a = datn[(size_t) ix*nt+is]/norm[ix];
							ff += a*b;
						}
					}
				}
			}
			if (normalize) {
				ff = beta > 0.0 ? ff/beta : 0.0;
				if (ff < 0.0) ff = 0.0;
			}
			if (pwr != 1.0) ff = (ff < 0.0 ? -1.0 : 1.0) * pow(fabs(ff), pwr);
			out[itout] = ff;
		}
	}
}
