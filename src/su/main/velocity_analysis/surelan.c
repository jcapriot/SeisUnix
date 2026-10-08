/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURELAN: $Revision: 1.8 $ ; $Date: 2011/11/16 23:40:27 $		*/

#include "su.h"
#include "segy.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SURELAN - compute residual-moveout semblance for cdp gathers based	",
"	on z(h)*z(h) = z(0)*z(0) + r*h*h where z depth and h offset.	",
"									",
" surelan <stdin >stdout   [optional parameters]			",
"									",
" Optional Parameters:							",
" nr=51			number of r-parameter samples   		",
" dr=0.01               r-parameter sampling interval			",
" fr=-0.25               first value of b-parameter			",
" smute=1.5             samples with RMO stretch exceeding smute are zeroed",
" dzratio=5             ratio of output to input depth sampling intervals",
" nsmooth=dzratio*2+1   length of semblance num and den smoothing window",
" verbose=0             =1 for diagnostic print on stderr		",
"									",
" Note: 								",
" 1. This program is part of Zhenyue Liu's velocity analysis technique.	",
" 2. Input migrated traces should be sorted by cdp - surelan outputs a 	",
"    group of semblance traces every time cdp changes.  Therefore, the  ",
"    output will be useful only if cdp gathers are input.  		",
" 3. The parameter r may take negative values. The range of r can be 	",
"     controlled by maximum of (z(h)*z(h)-z(0)*z(0))/(h*h)   		",
NULL};
#endif

/**************** end self doc *******************************************/

/* Library version of SURELAN
 *
 * The main program is not part of the library. It reads the gathers (the traces that have the same cdp), and writes nr
 * semblance traces for each of them. Doing that is left to the caller, and what is here is the work of the program on the
 * traces: su_relan_accumulate, the residual moveout of a trace at each of the r parameters, added to the sums for the
 * semblance. The semblance of each r parameter is made from the sums by su_velan_semblance (the smoothing is the same).
 *
 * All the arrays are the caller's, none are allocated here. The sums are num[nr*nz], den[nr*nz] and nnz[nr*nz] (the row of
 * each r parameter after the other), that start at 0 for each gather.
 */

/* Residual moveout z(h)^2 = z(0)^2 + r h^2 of a trace of nz samples that start at depth fz and are dz apart, at offset h =
 * `offset`, and at the parameters fr + ir * dr, with linear interpolation and the sums for the semblance.
 *
 * smute: samples with residual moveout stretch exceeding smute are zeroed.
 */
void su_relan_accumulate(int nz, float dz, float fz, float offset, int nr, float dr, float fr, float smute,
	const float *data, float *num, float *den, float *nnz)
{
	int ir,iz,izmute,izi;
	float r,roffs2,znmute,zn,zi,frac,temp;

	for (ir=0,r=fr; ir<nr; ++ir,r+=dr) {
		float *pnum = num + (size_t) ir*nz;
		float *pden = den + (size_t) ir*nz;
		float *pnnz = nnz + (size_t) ir*nz;

		/* compute r multiplied by offset squared  */
		roffs2 = r*offset*offset;

		/* determine mute depth after rmo */
		znmute = roffs2/(smute*smute-1.0);
		znmute = (znmute>0)?sqrt(znmute):sqrt(-znmute);
		izmute = (znmute-fz)/dz;
		if(izmute<0) izmute = 0;

		/* do rmo via linear interpolation and  */
		/* accumulate semblance numerator and denominator */
		for (iz=izmute,zn=fz+izmute*dz; iz<nz; ++iz,zn+=dz) {
			temp = zn*zn+roffs2;
			zi = (temp>fz*fz)?(sqrt(temp)-fz)/dz:0;
			izi = zi;
			if (izi>=0 && izi<nz-1) {
				frac = zi-izi;
				temp = (1.0-frac)*data[izi]+frac*data[izi+1];
				if (temp!=0.0) {
					pnum[iz] += temp;
					pden[iz] += temp*temp;
					pnnz[iz] += 1.0;
				}
			}
		}
	}
}
