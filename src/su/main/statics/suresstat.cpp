/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURESSTAT: $Revision: 1.18 $ ; $Date: 2011/11/16 23:16:23 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>
#include <algorithm>

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SURESSTAT - Surface consistent source and receiver statics calculation",
" 									",
"   suresstat <stdin [optional parameters]				",
" 									",
" Required parameters: 							",
" ssol=		output file source statics				",
" rsol=		output file receiver statics				",
" ntraces=	number of traces in input data set (must be correct!)	",
" 									",
" Optional parameters:							",
" ntpick=50 	maximum static shift (samples)         			",
" niter=5 	number of iterations					",
" nshot=240 	largest shot number (fldr=1 to nshot)			",
" nr=335 	largest receiver number (tracf=1 to nr)			",
" nc=574 	maximum number of cmp's (for array allocation)		",
" sfold=96 	maximum shot gather fold				",
" rfold=96 	maximum receiver gather fold				",
" cfold=48 	maximum cmp gather fold					",
" sub=0 	subtract super trace 1 from super trace 2 (=1)		",
" 		sub=0 strongly biases static to a value of 0		",
" mode=0 	use global maximum in cross-correllation window		",
"		=1 choose the peak perc=percent smaller than the global max.",
" perc=10. 	percent of global max (used only for mode=1)		",
" verbose=0 	print diagnostic output (verbose=1)                     ",
" 									",
" Notes:								",
" Estimates surface-consistent source and receiver statics, meaning that",
" there is one static correction value estimated for each shot and receiver",
" position.								",
" 									",
" The method employed here is based on the method of Ronen and Claerbout:",
" Geophysics 50, 2759-2767 (1985).					",
"  									",
" The code builds a supertrace1 and supertrace2, which are subsequently	",
" cross-correllated. The program then picks the time lag associated with",
" the largest peak in the cross-correllation according to two possible	",
" criteria set by the parameter \"mode\". If mode=0, the maximum of the	",
" cross-correllation window is chosen. If mode=1, the program will pick ",
" a peak which is up to perc=percent smaller than the global maximum, but",
" closer to zero lag than the global maximum.	(Choosing mode=0 is	",
" recommended.)								",
"									",
" To apply the static corrections, use sustatic with hdrs=3		",
NULL};
#endif

/* Reference:
 *
 *  Ronen, J. and Claerbout, J., 1985, Surface-consistent residual statics
 *      estimation  by stack-power maximization: Geophysics, vol. 50,
 *      2759-2767.
 *
 * Credits:
 *	CWP: Timo Tjan, 4 October 1994
 *
 *      rewritten by Thomas Pratt, USGS, Feb. 2000.
 *
 * Trace header fields accessed: ns, dt, tracf, fldr, cdp
 */
/**************** end self doc *******************************************/

/* Pick the time lag of the peak of the cross-correlation of 2 * ntpick + 1 lags (corr[ntpick] is the zero lag).
 *
 * mode 0 is the global maximum (the first one, and zero lag if nothing is positive); mode 1 is the lag closest to zero
 * of those that are above perc percent below the global maximum (the program gives -(ntpick + 1) if there are none:
 * it is zero lag here). */
static int pick_lag(const float *corr, int mode, float perc, int n)
{
	float maxamp = 0.;
	int zero = (n - 1) / 2;
	int globmax = zero, locmax = zero;
	int j;
	float frac = 1.0 - (perc/100);
	bool found = false;

	for (j = 0; j < n; j++) {
		if (corr[j] > maxamp) {
			maxamp = corr[j];
			globmax = j;
		}
	}
	if (mode == 1) {
		/* the max within perc% of global max, but closer to zero lag */
		for (j = 0; j < n; j++) {
			if (corr[j] > frac*maxamp && (!found || ABS(j - zero) < ABS(locmax - zero))) {
				locmax = j;
				found = true;
			}
		}
		return locmax - zero;
	}
	return globmax - zero;
}

/* The algorithm of the program on traces that are in memory.
 *
 * nt, ntraces		samples per trace, number of traces
 * data[ntraces][nt]	the traces
 * icmp, irec, ishot	[ntraces] the index (from 0) of the CMP, receiver and shot of each trace
 * ncmp, nrec, nshot	the numbers of CMPs, receivers and shots
 * ntpick		the greatest static shift (samples)
 * niter		the number of iterations
 * mode, perc, sub	the options of the program
 * tsstat[nshot]	the source statics (samples, to apply to the traces of the shot: positive is earlier)
 * trstat[nrec]		the receiver statics (samples)
 *
 * Differences from the program, which are its bugs: the samples that it reads from a trace are 1 to nt (the first is
 * dropped, and the one after the last is used); the pick of mode 1 can be a lag of -(ntpick + 1); and its supertraces
 * are as long as for the greatest number of the shots and the receivers, whatever the fold is (they are only as long as
 * the fold needs here: the rest is zeros that do not change the sums).
 */
void su_resstat(int nt, int ntraces, const float *data_in, const int *icmp, const int *irec, const int *ishot,
	int ncmp, int nrec, int nshot, int ntpick, int niter, int mode, float perc, int sub, float *tsstat,
	float *trstat)
{
	int i, j, it, iter, lag, s;
	const int nt_super = nt + 2*ntpick;	/* trace length plus a buffer on each end */
	const int ntout = 2*ntpick + 1;
	std::vector<float> data(data_in, data_in + (size_t) ntraces*nt);
	std::vector<float> model((size_t) ncmp*nt), corr(ntout), shifted(nt);
	std::vector<int> sstat(nshot, 0), rstat(nrec, 0);
	std::vector<std::vector<int> > cmp_loc(ncmp), rec_loc(nrec), shot_loc(nshot);

	for (i = 0; i < ntraces; i++) {
		cmp_loc[icmp[i]].push_back(i);
		rec_loc[irec[i]].push_back(i);
		shot_loc[ishot[i]].push_back(i);
	}
	std::fill(tsstat, tsstat + nshot, 0.0f);
	std::fill(trstat, trstat + nrec, 0.0f);

	/* the sum of the traces of each CMP */
	auto stack = [&]() {
		std::fill(model.begin(), model.end(), 0.0f);
		for (int c = 0; c < ncmp; c++)
			for (size_t k = 0; k < cmp_loc[c].size(); k++)
				for (int t = 0; t < nt; t++)
					model[(size_t) c*nt + t] += data[(size_t) cmp_loc[c][k]*nt + t];
	};

	/* the static of one gather (shot or receiver): the lag that the supertrace of its traces correlates best with the
	   supertrace of the stacks of their CMPs */
	auto estimate = [&](const std::vector<int> &traces) -> int {
		const int fold = (int) traces.size();
		const size_t total = (size_t) fold*nt_super;
		std::vector<float> g(total + ntpick + 1, 0.0f), mt(total + ntpick + 1, 0.0f);
		for (int ichan = 0; ichan < fold; ichan++) {
			size_t jj = (size_t) ichan*nt_super + ntpick;
			int c = icmp[traces[ichan]];
			for (int t = 0; t < nt; t++) {
				g[jj + t] = data[(size_t) traces[ichan]*nt + t];
				mt[jj + t] = model[(size_t) c*nt + t];
			}
		}
		/* subtract the gather supertrace from the CMP supertrace if desired */
		if (sub == 1)
			for (size_t k = 0; k < total; k++) mt[k] -= g[k];
		/* cross-correlate */
		for (int k = 0; k < ntout; k++) {
			int tlag = k - ntpick;
			float sum = 0.0f;
			for (size_t jj = ntpick; jj <= total - ntpick; jj++)
				sum += g[jj + tlag]*mt[jj];
			corr[k] = sum;
		}
		return pick_lag(corr.data(), mode, perc, ntout);
	};

	/* the traces of a gather moved by its static: samples that come from outside are 0 */
	auto apply = [&](const std::vector<int> &gather_of_trace, const std::vector<int> &stat) {
		for (int tr = 0; tr < ntraces; tr++) {
			s = stat[gather_of_trace[tr]];
			if (s == 0) continue;
			float *x = &data[(size_t) tr*nt];
			std::fill(shifted.begin(), shifted.end(), 0.0f);
			for (it = std::max(0, -s); it < std::min(nt, nt - s); it++) shifted[it] = x[it + s];
			std::copy(shifted.begin(), shifted.end(), x);
		}
	};
	std::vector<int> shot_of(ishot, ishot + ntraces), rec_of(irec, irec + ntraces);

	for (iter = 1; iter <= niter; iter++) {
		stack();
		for (int sh = 0; sh < nshot; sh++) {
			if (shot_loc[sh].size() > 1) {
				lag = estimate(shot_loc[sh]);
				sstat[sh] = lag;
				tsstat[sh] += lag;
			}
		}
		apply(shot_of, sstat);

		/* the stacks of the corrected traces */
		stack();
		for (int r = 0; r < nrec; r++) {
			if (rec_loc[r].size() > 1) {
				lag = estimate(rec_loc[r]);
				rstat[r] = lag;
				trstat[r] += lag;
			}
		}
		apply(rec_of, rstat);
	}
	(void) j;
}
