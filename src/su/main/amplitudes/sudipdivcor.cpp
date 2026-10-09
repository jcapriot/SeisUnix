/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */


/* SUDIPDIVCOR: $Revision: 1.4 $ ; 					*/

#include "su.h"
#include "segy.h"
#include "header.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 								",
" SUDIPDIVCOR - Dip-dependent Divergence (spreading) correction	",
" 								",
"	sudipdivcor <stdin >stdout  [optional parms]		",
" 								",
" Required Parameters:						",
"	dxcdp	distance between sucessive cdps	in meters	",
" 								",
" Optional Parameters:						",
"	np=50		number of slopes			",
"	tmig=0.0	times corresponding to rms velocities in vmig",
"	vmig=1500.0	rms velocities corresponding to times in tmig",
"	vfile=binary	(non-ascii) file containing velocities vmig(t) ",
"	conv=0		=1 to apply the conventional divergence correction",
"	trans=0		=1 to include transmission factors 	",
"	verbose=0	=1 for diagnostic print			",
" 								",
" Notes:								",
" The tmig, vmig arrays specify an rms velocity function of time.",
" Linear interpolation and constant extrapolation is used to determine",
" rms velocities at times not specified.  Values specified in tmig",
" must increase monotonically.					",
" 								",
" Alternatively, rms velocities may be stored in a binary file",
" containing one velocity for every time sample.  If vfile is	",
" specified, then the tmig and vmig arrays are ignored.		",
" The time of the first sample is assumed to be constant, and is",
" taken as the value of the first trace header field delrt. 	",
" 								",
" Whereas the conventional divergence correction (sudivcor) is	",
" valid only for horizontal reflectors, which have zero reflection",
" slope, the dip-dependent divergence correction is valid for any",
" reflector dip or reflection slope.  Only the conventional	",
" correction will be applied to the data if conv=1 is specified. ",
" Note that the conventional correction over-amplifies		",
" reflections from dipping beds					",
" 								",
" The transmission factor should be applied when the divergence ",
" corrected data is to be migrated with a reverse time migration ",
" based on the constant density wave equation.			",
"								",
" Trace header fields accessed:  ns, dt, delrt			",
" 								",
NULL};
#endif
/**************** end self doc**************************************/

/*  The hidden variable "norm" turns on and off the normalization of the
divergence correction by the correction at zero reflection slope and the
first time sampling interval. 	The dip-dependent divergence correction is
equivalent to the conventional correction (sudivcor) only when the
normalization is turned on (the default).  The dip-dependent divergence
correction is designed to correct for geometrical spreading along the
downgoing zero-offset raypath to the reflection point.   Conventional
correction, however, compensates for spreading along both the incident and
the reflected zero-offset raypath.  The conventional correction is equal to
twice the dip-dependent divergence correction for zero reflection slope.
When the dip-dependent divergence correction is normalized, this constant
scale factor of two cancels. 
	If the normalization is turned off and transmission is turned on;
zero-offset `susynlv' amplitudes after dip-dependent divergence correction
will equal exploding reflector `susynlv' amplitudes.
*/

/* Library version of SUDIPDIVCOR
 *
 * The main program is not part of the library. It reads the velocity function, Fourier transforms the traces in x
 * (the traces are the panel), loops over the wavenumbers, and transforms back. Doing that is left to the caller, and
 * what is here is the work of the program that is not the transforms in x:
 *
 *	su_dipdivcor_table:	the table of divergence corrections, for np reflection slopes
 *	su_dipdivcor_filter:	the dip filter of the divergence correction, for one wavenumber (of p(t,k))
 *
 * All the arrays are the caller's, none are allocated here.
 *
 * Differences from the program:
 *  - The table does not scale the velocity in place (the program leaves vt[] multiplied by 0.0005, and then uses
 *    vt[0] to work out the slope sampling dpx and the highest wavenumber): the caller gets those from
 *    su_dipdivcor_scale.
 *  - In the dip filter, the upper frequency for a slope of 0 is not worked out by dividing by 0 (it was replaced just
 *    after).
 */

/* The factor that the velocities are multiplied by: v/2 to compensate for v(two-way t), and from m/s to km/ms */
#define DIPDIVCOR_VSCALE 0.0005

float su_dipdivcor_scale(void)
{
	return DIPDIVCOR_VSCALE;
}

/* The table of divergence corrections for the np slopes, dpx = 1/(v0 scaled)/(np-1) apart.
 *
 * nt, dt: number of time samples and the sample interval (s). tt[nt]: the times, it*dt. vt[nt]: the velocity (m/s)
 * at those times. vs[nt] and vind[nt][4]: scratch (vs is the scaled velocity, vind the spline coefficients)
 * divcor[np*nt]: the table (np slopes of nt samples each)
 * trans: include the transmission factors. norm: normalize by the correction at zero slope.
 */
void su_dipdivcor_table(int nt, int np, float dt, const float *tt, const float *vt, float *vs, float (*vind)[4],
	float *divcor, int trans, int norm)
{
	int it,ip,flag;
	float dpx,px,px2,vel,velold,gamma,gammaold,pz,tau,alpha,v0,
		denom,vel1,vel1old,vel2,q,qold,p,pold,
		sigma,pz1;

	/* use v/2 to compensate for v(two-way t)*/
	for (it=0; it<nt; it++)
		vs[it] = vt[it]*DIPDIVCOR_VSCALE;

	/* establish spline coefficients*/
	cmonot(nt,(float *) tt,vs,vind);

	/* clear divergence correction array */
	for (ip=0; ip<np; ip++)
		for (it=0; it<nt; it++)
			divcor[ip*nt+it] = 0.0;

	v0 = vind[0][0];
	dpx = 1.0/v0/(np-1);

	/* evaluate divergence correction */
	for (ip=0; ip<np; ip++){

		/* calculate px */
		px = ip*dpx;

		/* initialize variables */
		vel = v0;
		vel1 = vind[0][1];
		vel2 = vind[0][2];
		tau=0.0;
		px2= px*px;
		gamma = px2*(pow(vel1/vel,2.0) - vel2/vel);
		p=1.0/v0;
		q=0.0;
		sigma=0.0;
		if (ip==np-1) {pz=0.0; flag=0;}
		else {pz = sqrt(1/vel/vel - px2); flag=1;}

		/* ray tracing */
		for (it=1; it<nt; it++){

			tau += vel*pz*dt;

			if (tau>=0.0){
				velold = vel;
				vel1old = vel1;
				intcub(0,nt,(float *) tt,vind,1,&tau,&vel);
				intcub(1,nt,(float *) tt,vind,1,&tau,&vel1);
				intcub(2,nt,(float *) tt,vind,1,&tau,&vel2);
				if (flag==1) {
					/* pz1(it=1) from oldvel oldpz*/
					pz1 = pz-vel1old/velold/velold*dt;
					/* pz(it=1) based on vel */
					pz = sqrt(1/vel/vel - px2);
					if (pz1<0.0 || (px*vel>1.0)) {
						flag=0; pz = pz1;}
				}
				else pz += -vel1old/velold/velold*dt;

				sigma += velold*velold*dt;

				/* Crank-Nicolson on p&q */
				gammaold = gamma;
				gamma = px2*(pow(vel1/vel,2.0) - vel2/vel);
				alpha = 1.0 - gamma*pow((dt/2*vel),2.0);
				qold = q;
				pold = p;
				q = pold*(dt/2/alpha)*(vel*vel + velold*velold) + qold/alpha*(1.0 + pow((dt*vel/2),2.0)*gammaold);
				p = qold*(dt/2/alpha)*(gamma + gammaold) + pold/alpha*(1.0 + pow((dt*velold/2),2.0)*gamma);

/*factor of 4 for vrms agreement */
/*correct zero-offset to exp. reflector; one-way w/transmission */

				if (trans==1){
					divcor[ip*nt+it] =sqrt(2*sigma*ABS(q)/vel);
				} else {
					divcor[ip*nt+it]=sqrt(2*sigma*ABS(q)/v0);
				}
			} else{
				vel = v0;
				sigma += vel*vel*dt;
				gammaold = gamma;
				gamma = 0.0;
				alpha = 1.0;
				qold = q;
				pold = p;
				q = pold*(dt*vel*vel)
					+ qold*(1.0
						+ pow((dt*vel/2),2.0)*gammaold);
				p = qold*dt/2*gammaold + pold;
				if (trans==1){
					divcor[ip*nt+it]=sqrt(2*sigma*ABS(q)/vel);
				} else {
					divcor[ip*nt+it]=sqrt(2*sigma*ABS(q)/v0);
				}
			}

		}/* end t*/
	} /* end p */

	/* normalize dip-dependent correction */
	if (norm==1){
		denom=1.0/divcor[1];

		for (ip=np-1; ip>=0; ip--)
			for (it=1; it<nt; it++)
				divcor[ip*nt+it] *= denom;
	}
}

/*********************************************************************
Jacubowicz filter to apply dip-dependent divergence correction
**********************************************************************
Input:
k		wavenumber
dpx		dip sampling interval
dt		time sampling interval
np		number of reflection slopes
nt		number of time samples
nw		number of frequency samples (npfa(nt))
div		amplitude table [np*nt]
p		array[nt] containing input p(t,k)
kq, qq		scratch arrays of nw complex numbers

Output:
q		array[nt] containing divergence corrected output q(t,k)
		(p and q may be the same array)
*********************************************************************/
void su_dipdivcor_filter(float k, float dpx, float dt, int np, int nw, int nt, const float *div,
	const complex *p, complex *q, complex *kq, complex *qq)
{
	int ip,iw,it,iwl,iwh;
	float dw,wny,pscl,fftscl,pmin,ph,pm,pl;

	dw = 2.0*PI/(nw*dt);
	wny = PI/dt;
	pscl = 0.5;
	fftscl = 1.0/nw;
	pmin = k/wny;

	/* initialize qq */
	for (iw=0; iw<nw; iw++){
		qq[iw] = 0.0f;
	}

	for (ip=np-1; ip>=0; ip--){

		ph=dpx*(ip+1);
		pm=dpx*ip;
		pl=dpx*(ip-1);

		/* define frequency range */
		iwl = k/ph/dw + 1;
		iwh = pl>0.0 ? k/pl/dw : 0;

		if (pl<1.01*pmin) iwh = (nw%2 ? (nw-1)/2 : nw/2-1);
		if (pm<1.01*pmin) iwh = (nw%2 ? nw/2 : nw/2-1);

		/* sum over frequency */
		if (iwh>=iwl){
			for (it=0; it<nt; it++){
				kq[it] = p[it]*div[ip*nt+it];
			}
			for (it=nt; it<nw; it++){
				kq[it] = 0.0f;
			}
			pfacc(1,nw,kq);

			/* dip filter positive frequencies */
			for (iw=iwl; iw<=iwh; iw++){
				qq[iw] += kq[iw]*pscl;
			}

			/* dip filter negative frequencies */
			iwl=nw-iwl;
			iwh=nw-iwh;

			for (iw=iwh; iw<=iwl; iw++){
				qq[iw] += kq[iw]*pscl;
			}
		}
		if (pm<1.01*pmin) break;
	}

	/* Fourier transform w to t */
	pfacc(-1,nw,qq);
	for (it=0; it<nt; it++){
		q[it] = qq[it]*fftscl;
	}
}
