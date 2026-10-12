/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTIVEL: $Revision: 1.5 $ ; $Date: 2011/11/16 23:40:27 $        */

#include "su.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = { 
"									",
"  SUTIVEL -  SU Transversely Isotropic velocity table builder		",
"	computes vnmo or vphase as a function of Thomsen's parameters and",
"	theta and optionally interpolate to constant increments in slowness",
"									",
" Optional Parameters:							",
" a=2500.		alpha (vertical p velocity)			",
" b=1250.		beta (vertical sv velocity)			",
" e=.20			epsilon (horiz p-wave anisotropy)		",
" d=.10			delta (strange parameter)			",
" maxangle=90.0		max angle in degrees				",
" nangle=9001		number of angles to compute			",
" verbose=0		set to 1 to see full listing			",
" np=8001		number of slowness values to output		",
" option=1		1=output vnmo(p) (result used for TI DMO)	",
"			2=output vnmo(theta) in degrees			",
"			3=output vnmo(theta) in radians			",
"			4=output vphase(p)				",
"			5=output vphase(theta) in degrees		",
"			6=output vphase(theta) in radians		",
"			7=output first derivative vphase(p)		",
"			8=output first derivative vphase(theta) in degrees",
"			9=output first derivative vphase(theta) in radians",
"			10=output second derivative vphase(p)		",
"			11=output second derivative vphase(theta) in degrees",
"			12=output second derivative vphase(theta) in radians",
"			13=( 1/vnmo(0)^2 -1/vnmo(theta)^2 )/p^2 test vs theta",
"			   (result should be zero for all theta for d=e)",
"			14=return vnmo(p) for weak anisotropy		",
" normalize=0		=1 means scale vnmo by cosine and scale vphase by",
" 			    1/sqrt(1+2*e*sin(theta)*sin(theta)		",
"	 		   (only useful for vphase when d=e for constant",
"				result)					",
"			=0 means output vnmo or vphase unnormalized	",
"									",
" Output on standard output is ascii text with:				",
" line   1: number of values						",
" line   2: abscissa increment (p or theta increment, always starts at zero)",
" line 3-n: one value per line						",
"									",
NULL};
#endif

/*
 * Author: (visitor to CSM form Mobil) John E. Anderson, Spring 1994
 */
/**************** end self doc ********************************/

static void vget( float a, float b, float e, float d, float theta, float *vel);
static float vnmoweakTI( float a, float b, float e, float d, float p);
static float DweakTI( float a, float b, float e, float d, float p);

static void vget( float a, float b, float e, float d, float theta, float *vel)
/*

This routine returns phase and NMO velocity information given
an input angle and Thomsen's parameters for a TI medium.

input parameters:
----------------
a      = alpha = vertical p-wave phase velocity = vrms/sqrt(1+2*d)    
b      = beta = vertical shear wave phase velocity
e      = epsilon = anisotropy factor for horizontally propagating p waves  
d      = delta = 0.5*(e + ds/(1-(b*b)/(a*a)))	
theta  = angle in radians 

returned parameters:
-------------------
vel[0] = p-wave phase velocity
vel[1] = first derivative of p-wave phase velocity with respect to theta
vel[2] = second derivative of p-wave phase velocity with respect to theta
vel[3] = NMO velocity
vel[4] = ray parameter p = sin(theta)/vphase
vel[5] = NMO velocity using weak anisotropy assumption

*/
{
	float p, psi, sint, cost, sin2t, cos2t, tant, eps, f, g;
	float vnmo,v,dv,d2v,gamma,dgamma,d2gamma,sqgamma;
	eps=1.0e-7;
	sint=sin(theta);
	cost=cos(theta);
	sin2t=sint*sint;
	cos2t=cost*cost;
	psi=1.-(b*b)/(a*a);
	f = psi*(2*d-e);
	g = (psi+e)*e;
	gamma=0.25*psi*psi+f*sin2t*cos2t+g*sin2t*sin2t;
	dgamma=f*2*sint*cost*cos2t +
		( 4*g - 2*f )*sint*sin2t*cost;
	d2gamma=f*2*cos2t*cos2t +
		( g - f )*12*cos2t*sin2t +
		( 2*f - 4*g )*sin2t*sin2t;	
	sqgamma=sqrt(gamma);
	v=a*sqrt(1.+e*sin2t-0.5*psi+sqgamma);
	dv=0.5*a*a*(2*e*sint*cost + 0.5*dgamma/sqgamma)/v;
	d2v=0.5*a*a*(2*e*(cos2t-sin2t)+0.5*d2gamma/sqgamma
		-0.25*dgamma*dgamma/(sqgamma*gamma) )/v - dv*dv/v;
	if(cost<eps) cost=eps;
	tant=sint/cost;
	vnmo = ( v*sqrt( 1 + d2v/v) )/( cost * (1-tant*dv/v) );
	p=sint/v;
	vel[0]=v;
	vel[1]=dv;
	vel[2]=d2v;
	vel[3]=vnmo;
	vel[4]=p;
	return;
}
static float DweakTI( float a, float b, float e, float d, float p)
/*
	returns D(p) for the weak anisotropy approximation

input parameters:
----------------
a      = alpha = vertical p-wave phase velocity = vrms/sqrt(1+2*d)    
b      = beta = vertical shear wave phase velocity
e      = epsilon = anisotropy factor for horizontally propagating p waves  
d      = delta = 0.5*(e + ds/(1-(b*b)/(a*a)))	
p      = slowness = sin(theta)/v

*/
{
	float y;
	float vnmo0;
	float D;

	vnmo0=a*sqrt(1+2*d);
	y=p*vnmo0;
	y=y*y;
	D=p*p*( 1 + 2*(e-d)*(4*y*y-9*y+6) );
	return (D);
}

static float vnmoweakTI( float a, float b, float e, float d, float p)
/*
	returns vnmo(p) for the weak anisotropy approximation
	or a value of -1 for the evanescent region

input parameters:
----------------
a      = alpha = vertical p-wave phase velocity = vrms/sqrt(1+2*d)    
b      = beta = vertical shear wave phase velocity
e      = epsilon = anisotropy factor for horizontally propagating p waves  
d      = delta = 0.5*(e + ds/(1-(b*b)/(a*a)))	
p      = slowness = sin(theta)/v

*/
{
	float y;
	float D;
	float vnmo;

	D=DweakTI(a,b,e,d,p);
	y=a*a*(1+2*d);
	y=1/y;
	if(y>D) {
		vnmo=1./sqrt(y-D);
	}else{
		vnmo=-1;
	}
	return (vnmo);
}

/* The values of the program for the angles from 0 to maxangle (degrees), nangle of them:
 *
 * p[nangle]	the ray parameter sin(theta) / vphase
 * v, dv, d2v	the phase velocity and its first and second derivatives with respect to theta (radians)
 * vnmo		the NMO velocity
 *
 * normalize: scale vnmo by the cosine and vphase by 1 / sqrt(1 + 2 e sin^2(theta)) (the program does that to the
 * velocity but not to its derivatives).
 */
void su_tivel_table(float a, float b, float e, float d, float maxangle, int nangle, int normalize, float *p, float *v,
	float *dv, float *d2v, float *vnmo)
{
	int iangle;
	float vel[6], dangle, theta;

	if (nangle < 2) nangle = 2;
	dangle = maxangle/(nangle-1);
	for (iangle=0; iangle<nangle; iangle++) {
		theta = iangle*dangle*PI/180.;
		vget(a,b,e,d,theta,vel);
		if (normalize) {
			vel[0] /= sqrt(1 + 2*e*sin(theta)*sin(theta));
			vel[3] *= cos(theta);
		}
		v[iangle] = vel[0];
		dv[iangle] = vel[1];
		d2v[iangle] = vel[2];
		vnmo[iangle] = vel[3];
		p[iangle] = vel[4];
	}
}

/* The NMO velocity of the weak anisotropy approximation at the ray parameter p, or -1 in the evanescent region */
float su_tivel_vnmoweak(float a, float b, float e, float d, float p)
{
	return vnmoweakTI(a, b, e, d, p);
}
