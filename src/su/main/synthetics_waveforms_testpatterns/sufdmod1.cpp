/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFDMOD1: $Revision: 1.4 $ ; $Date: 2011/11/12 00:40:42 $	*/

#include "par.h"
#include "su.h"
#include "segy.h"
/*********************** self documentation ********************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUFDMOD1 - Finite difference modelling (1-D 1rst order) for the	",
" acoustic wave equation						"
"									",
" sufdmod1 <vfile >sfile nz= tmax= sz= [optional parameters]		",
"									",
" Required parameters :							",
" <vfile or vfile=	binary file containing velocities[nz]		",
" >sfile or sfile=	SU file containing seimogram[nt]		",
" nz=		 number of z samples				   	",
" tmax=		maximum propagation time				",
" sz=		 z coordinate of source					",
"									",
" Optional parameters :							",
" dz=1	   z sampling interval						",
" fz=0.0	 first depth sample					",
" rz=1	   coordinate of receiver					",
" sz=1	   coordinate of source						",
" dfile=	 binary input file containing density[nz]		",
" wfile=	 output file for wave field (snapshots in a SU trace panel)",
" abs=0,1	absorbing conditions on top and bottom			",
" styp=0	 source type (0: gauss, 1: ricker 1, 2: ricker 2)	",
" freq=15.0	approximate source center frequency (Hz)		",
" nt=1+tmax/dt   number od time samples (dt determined for numerical	",
" stability)								",
" zt=1	   trace undersampling factor for trace and snapshots	 	",
" zd=1	   depth undersampling factor for snapshots		   	",
" press=1	to record the pressure field; 0 records the particle	",
"		velocity						",
" verbose=0	=1 for diagnostic messages				",
"									",
" Notes :								",
"  This program uses a first order explicit velocity/pressure  finite	",
"  difference equation.							",
"  The source function is applied on the pressure component.		",
"  If no density file is given, constant density is assumed	 	",
"  Wavefield  can be easily viewed with suximage, user must provide f2=0",
"  to the ximage program in order to  get correct time labelling	",
"  Seismic trace is shifted in order to get a zero phase source		",
"  Source begins and stop when it's amplitude is 10^-4 its maximum	",
"  Time and depth undersampling only modify the output trace and snapshots.",
"  These parameters are useful for keeping snapshot file small and	",
"  the number of samples under SU_NFLTS.				",
"									",
NULL  };
#endif


/* Library version of SUFDMOD1
 *
 * The main program (the parameters, the files, the headers of the traces) is left to the caller. What is here are the two steps of
 * the computation, with no static state, and all of the arrays are owned by the caller:
 *
 *	su_fdmod1_plan:	the time step for stability, the number of time steps, and the time shift of the source
 *	su_fdmod1_run:	the propagation, giving the seismogram and optionally the snapshots
 *
 * Differences from the program: only the first of the two values of abs was read (the second was whatever was in memory); both are
 * used here. The checks of the parameters return a code.
 */

static float source (float t, int styp, float dt, float dz, float t0, float alpha);

/* The time step dt (for numerical stability), the number of time steps nt (nt_in if it is not 0, otherwise 1+tmax/dt), the time shift
 * t0 of the source to make it near causal, and the number of time steps ies for which the source is added.
 * Returns 0, or -3 for parameters that are not allowed.
 */
int su_fdmod1_plan(const float *rv, int nz, float dz, float tmax, int nt_in, int styp, float freq,
	float *dt, int *nt, float *t0, int *ies)
{
	int iz;
	float vmax=0.0,alpha,epst0,t,dtl;

	if (nz<2 || dz<=0.0 || freq<=0.0 || tmax<0.0 || nt_in<0 || styp<0 || styp>2) return -3;
	for (iz=0; iz<nz; iz++) if (rv[iz]>vmax) vmax=rv[iz];
	if (vmax<=0.0) return -3;
	dtl = dz/1.414/vmax/2;
	*dt = dtl;
	*nt = (nt_in==0) ? (int)(1+tmax/dtl) : nt_in;

	/* source parameter computation */
	alpha=2*9.8696*freq*freq;

	/* time shift to get a t0 centered source */
	if ((styp==0) || (styp == 2)) epst0=fabs(source (0, styp, dtl, dz, 0, alpha) / 1e4);
	else epst0=fabs(source (1/sqrt(2*alpha), styp, dtl, dz, 0, alpha)) / 1e4;

	t=tmax+dtl;
	do t=t-dtl; while (t>0.0 && fabs(source(t, styp, dtl, dz, 0, alpha))<epst0);
	*t0=t;
	*ies=(int)(2*t/dtl);
	return 0;
}

/* The propagation of the pressure p and the particle velocity v (first order, explicit) in the velocities rv[nz] and densities
 * rd[nz], for nt+1 time steps of dt, with the source at the sample isz and the receiver at irz, abs[2] the absorbing conditions at
 * the top and bottom, td and zd the undersampling in time and depth, press=1 to record the pressure and 0 the particle velocity.
 *
 * sismo: array[nt/td+1], the seismogram; snaps: NULL, or array[nt/td+1][nz/zd], the snapshots.
 * Returns the number of samples of the seismogram, or -3 for parameters that are not allowed.
 */
int su_fdmod1_run(const float *rv, const float *rd, int nz, float dz, float dt, int nt, float t0, int ies,
	int isz, int irz, const int *abs, int styp, float freq, int td, int zd, int press,
	float *sismo, float *snaps)
{
	float *p,*v,*rd1_5;
	float alpha,t;
	int iz,it,itsis=0,nzs;

	if (nz<2 || nt<1 || td<1 || zd<1 || isz<0 || isz>=nz || irz<0 || irz>=nz || dz<=0.0 || dt<=0.0) return -3;
	nzs = nz/zd;
	p=ealloc1float(nz);
	v=ealloc1float(nz);
	rd1_5=ealloc1float(nz);
	for (iz=0; iz<nz-1; iz++) rd1_5[iz]=(rd[iz]+rd[iz+1])/2;
	rd1_5[nz-1]=rd[nz-1];
	alpha=2*9.8696*freq*freq;

	/* array initialization */
	for (iz=0; iz<nz; iz++) {  v[iz]=0; p[iz]=0;  }

	/* propagation computation */
	for (it=0; it<=nt; it++) {
		t=it*dt;
		if (abs[0]==1) p[0]=(p[0]*(1-rv[0]*dt/dz)+2*rd[0]*rv[0]*rv[0]*dt/dz*v[0])/(1+rv[0]*dt/dz);
		else p[0]=0;
		for (iz=1; iz<nz; iz++) p[iz]=p[iz]+rd[iz]*rv[iz]*rv[iz]*dt/dz*(v[iz]-v[iz-1]);
		if (abs[1]!=1) p[nz-1]=0;
		if (it<ies) {
		p[isz]=p[isz]+source(t, styp, dt, dz, t0, alpha);
		}

		for (iz=0; iz<nz-1; iz++) v[iz]=v[iz]+dt/rd1_5[iz]/dz*(p[iz+1]-p[iz]);

		if (abs[1] != 1) v[nz-1]=0;
		else
		v[nz-1]=((rd1_5[nz-1]*dz-dt*rd[nz-1]*rv[nz-1])*v[nz-1]-2*dt*p[nz-1])/(rd1_5[nz-1]*dz+dt*rd[nz-1]*rv[nz-1]);

		if (it % td == 0) {
			sismo[itsis] = (press==1) ? p[irz] : v[irz];
			if (snaps!=NULL) {
				float *snap = snaps+(size_t)itsis*nzs;
				for (iz=0; iz<nzs; ++iz) snap[iz] = (press==1) ? p[iz*zd] : v[iz*zd];
			}
			itsis++;
		}
	}
	free1float(p);
	free1float(v);
	free1float(rd1_5);
	return itsis;
}


static float source (float t, int styp, float dt, float dz, float t0, float alpha)
{
	float x=0.0, sou=0.0;
	x=-alpha*(t-t0)*(t-t0);
	if (x>-40) {
	 if (styp==0) sou=exp(x);
	 	else if (styp==1) sou=-2*alpha*(t-t0)*exp(x);
	 	else if (styp==2) sou=2*alpha*(1+2*x)*exp(x);
		}
	else sou=0;
	return sou/dz*dt*1e8;
}
