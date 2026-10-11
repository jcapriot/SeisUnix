/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUTIHALEDMO: $Revision: 1.5 $ ; $Date: 2011/11/16 17:51:02 $        */

#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = { 
" SUTIHALEDMO - TI Hale Dip MoveOut (based on Hale's PhD thesis)	",
"									",
"  sutihaledmo <infile >outfile [optional parameters]			",
"									",
"									",
" Required Parameters:							",
" nxmax		  maximum number of midpoints in common offset gather",
"									",
" Optional Parameters:							",
" option=1		1 = traditional Hale DMO (from PhD thesis)	",
"			2 = Bleistein's true amplitude DMO		",
"			3 = Bleistein's cos*cos weighted DMO		",
"			4 = Zhang's DMO					",
"			5 = Tsvankin's anisotropic DMO			",
"			6 = Tsvankin's VTI DMO weak anisotropy approximation",
" dx=50.		 midpoint sampling interval between traces	",
"			in a common offset gather.  (usually shot	",
"			interval in meters)				",
" v=1500.0		velocity (in meters/sec)			",
"			(must enter a positive value for option=3)	",
"			(for excluding evanescent energy)		",
" h=200.0		source-receiver half-offset (in meters)		",
" ntpad=0		number of time samples to pad			",
" nxpad=h/dx		number of midpoints to pad			",
" file=vnmo		name of file with vnmo as a function of p	",
"			used for option=5--otherwise not used		",
"			(Generate this file by running program		",
"			sutivel with appropriate list of Thomsen's	",
"			parameters.)					",
" e=0.			Thompsen's epsilon				",
" d=0.			Thompsen's delta				",
"									",
"Note:									",
"									",
" This module assumes a single common offset gather after NMO is	",
" to be input, DMO corrected, and output.  It is useful for computing	",
" theoretical DMO impulse responses.  The Hale algorithm is		",
" computationally intensive and not commonly used for bulk processing	",
" of all of the offsets on a 2-D line as there are cheaper alternative	",
" algorithms.  The Hale algorithm is commonly used in theoretical studies.",
" Bulk processing for multiple common offset gathers is typically done	",
" using other modules.							",
"									",
" Test run:   suspike | sutihaledmo nxmax=32 option=1 v=1500 | suxwigb & ",
NULL};
#endif
 

/*
 * Author:  (Visitor to CSM from Mobil) John E. Anderson Spring 1994
 * References: Anderson, J.E., and Tsvankin, I., 1994, Dip-moveout by
 *	Fourier transform in anisotropic media, CWP-146
 */
/**************** end self doc ********************************/

static void haledmo(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,float v)
/*  hale fk dmo

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
v	velocity in m/sec (used to limit aperature of DMO response)

*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	float hkpw;		/*	h*k/w				*/
	float hkpwt;		/*	h*k/(w*tn)			*/
	float A;		/*	Hale's A			*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float eps=1.0e-20;	/*	zero offset*wavenumber thresh	*/
	float scale;		/*	forward/inverse fft scale factor*/
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int iwtest;		/*	test for pos or neg freq	*/
	

	hk = h*k;
	if(hk<eps) return;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	for(iw=1;iw<ntfft;iw++) {
		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;
		csum=complex(0.0f,0.0f);
		if( (k*v) < (2*fabs(w)) ) {
			hkpw=hk/(w);
			for(it=1;it<nt;it++) {
				tn=it*dt;
				hkpwt=hkpw/tn;
				A = sqrt(1. + hkpwt*hkpwt);
				phase = w*tn*A;
				cfac = complex( cos(phase)/A, sin(phase)/A ); 		
				csum += p[it]*cfac;

			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}
static void bleisteindmo(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,float v)
/*  bleistein true amplitude fk dmo

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
v	velocity in m/sec (used to limit aperature of DMO response)

*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	float hkpw;		/*	h*k/w				*/
	float hkpwt;		/*	h*k/(w*tn)			*/
	float A;		/*	Hale's A			*/
	float B;		/*	Bleistein's scale factor	*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float eps=1.0e-20;	/*	zero offset*wavenumber thresh	*/
	float scale;		/*	forward/inverse fft scale factor*/
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int iwtest;		/*	test for pos or neg freq	*/
	

	hk = h*k;
	if(hk<eps) return;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	for(iw=1;iw<ntfft;iw++) {
		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;
		csum=complex(0.0f,0.0f);
		if( (k*v) < (2*fabs(w)) ) {
			hkpw=hk/(w);
			for(it=1;it<nt;it++) {
				tn=it*dt;
				hkpwt=hkpw/tn;
				A = sqrt(1. + hkpwt*hkpwt);
				B = sqrt(1. + 2*hkpwt*hkpwt);
				phase = w*tn*A;
				cfac = complex( B*cos(phase)/A, B*sin(phase)/A ); 		
				csum += p[it]*cfac;
			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}
static void bleisteindmo2(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,float v)
/*  bleistein cos*cos weighted amplitude fk dmo

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
v	velocity in m/sec (used to limit aperature of DMO response)

*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	float hkpw;		/*	h*k/w				*/
	float hkpwt;		/*	h*k/(w*tn)			*/
	float hkpwt2;		/*	square of above			*/
	float hpv;		/* 	2*h/v				*/
	float hpvt2;		/*	[(2*h)/(v*tn)]^2		*/
	float A;		/*	Hale's A			*/
	float B;		/*	Bleistein's scale factor	*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float scale;		/*	forward/inverse fft scale factor*/
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int iwtest;		/*	test for pos or neg freq	*/
	

	hk = h*k;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	hpv=2.*h/v;
	for(iw=1;iw<ntfft;iw++) {
		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;
		csum=complex(0.0f,0.0f);
		if( (k*v) < (2*fabs(w)) ) {
			hkpw=hk/(w);
			for(it=1;it<nt;it++) {
				tn=it*dt;
				hkpwt=hkpw/tn;
				hkpwt2=hkpwt*hkpwt;
				hpvt2=hpv/tn;
				hpvt2=hpvt2*hpvt2;
				A = sqrt(1. + hkpwt2);
				B = (1.+2*hkpwt2)*(1+hkpwt2)/(1+hpvt2);
				phase = w*tn*A;
				cfac = complex( B*cos(phase)/A, B*sin(phase)/A ); 		
				csum += p[it]*cfac;
			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}
static void zhangdmo(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,float v)
/*  zhang fk dmo

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
v	velocity in m/sec (used to limit aperature of DMO response)

*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	float hkpw;		/*	h*k/w				*/
	float hkpwt;		/*	h*k/(w*tn)			*/
	float hkpwt2;		/* 	square of above			*/
	float A;		/*	Hale's A			*/
	float B;		/*	Zhang's scale factor		*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float eps=1.0e-20;	/*	zero offset*wavenumber thresh	*/
	float scale;		/*	forward/inverse fft scale factor*/
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int iwtest;		/*	test for pos or neg freq	*/
	
	hk = h*k;
	if(hk<eps) return;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	for(iw=1;iw<ntfft;iw++) {
		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;
		csum=complex(0.0f,0.0f);
		if( (k*v) < (2*fabs(w)) ) {
			hkpw=hk/(w);
			for(it=1;it<nt;it++) {
				tn=it*dt;
				hkpwt=hkpw/tn;
				hkpwt2=hkpwt*hkpwt;
				A = sqrt(1. + hkpwt2);
				B = (1.+2*hkpwt2)/(1+hkpwt2);
				phase = w*tn*A;
				cfac = complex( B*cos(phase)/A, B*sin(phase)/A ); 		
				csum += p[it]*cfac;
			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}
static void tsvankindmo(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,int np, float dp, float *vnmo)
/*  tsvankin fk dmo

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
np	number of slowness values for which vnmo has been tabled
dp	slowness increment for vnmo table
vnmo	nmo velocity in m/sec as a function of p 

*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	double v;		/*	vnmo(p)				*/
	double v0;		/*	vnmo(0)				*/
	float A;		/*	Hale's A			*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float tn2;
	float eps=1.0e-20;	/*	zero offset*wavenumber thresh	*/
	float scale;		/*	forward/inverse fft scale factor*/
	float pp;		/*	ray parameter			*/
	double fac;
	float rp;
	float wgt;
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int it1;		/* 	first stable time index		*/
	int iwtest;		/*	test for pos or neg freq	*/
	int ip;
	int npm;

	hk = h*k;
	if(hk<eps) return;
	npm=np-2;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	for(iw=1;iw<ntfft;iw++) {

		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;

		csum=complex(0.0f,0.0f);
		pp=fabs(k/(2.*w));
		rp=pp/dp;
		ip=rp;
		if(ip<npm) {
			wgt=rp-ip;
			v=(1.-wgt)*vnmo[ip] + wgt*vnmo[ip+1];
			v0=vnmo[0];
			fac=4*h*h*(1./(v0*v0) - 1./(v*v) );
			if(fac<0.) {
				it1=2+sqrt(fabs(fac))/dt;
			}else{
				it1=1;
			}
			for(it=it1;it<nt;it++) {
				tn=it*dt;
				tn2=tn*tn;
				A=sqrt(1. + fac/tn2);
				phase = w*tn*A;
				cfac = complex( cos(phase)/A, sin(phase)/A ); 		
				csum += p[it]*cfac;	
			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}
static void weakdmo(float h,float k,complex *p,complex *pwork,
		int nt,int ntfft,float dt,float e,float d, float v)
/*  tsvankin fk dmo weak approximation

h	half offset
k	wavenumber
p	input NMO corrected data in wavenumber domain, p(tn,k,h)
pwork	complex array of lenght ntfft for work space
nt	number of time samples
ntfft	lenght of temporal fft
dt	temporal sample rate in seconds
e	epsilon
d	delta
v	v0 velocity
*/
{
	complex csum;		/* complex variable for sum 		*/
	complex cfac;	
	float w;		/* 	angular freq 			*/
	float hk;		/* 	h*k				*/
	float A;		/*	Hale's A			*/
	float phase;		/*	phase				*/
	float dw;		/*	angular freq inc		*/
	float tn;		/*	input event time after NMO	*/
	float tn2;
	float eps=1.0e-20;	/*	zero offset*wavenumber thresh	*/
	float scale;		/*	forward/inverse fft scale factor*/
	float pp;		/*	ray parameter			*/
	float y;
	double fac;
	int iw;			/*	freq index			*/
	int it;			/*	tn time index			*/
	int it1;		/* 	first stable time index		*/
	int iwtest;		/*	test for pos or neg freq	*/

	hk = h*k;
	if(hk<eps) return;
	dw = 2*PI/(ntfft*dt);
	scale=1./ntfft;
	iwtest=ntfft/2+1;
	pwork[0]=complex(0.0f,0.0f);
	for(iw=1;iw<ntfft;iw++) {

		if(iw<iwtest)
			w=iw*dw;
		else
			w=(iw-ntfft)*dw;

		csum=complex(0.0f,0.0f);
		pp=fabs(k/(2.*w));
		if(pp<(1/v)) {
			y=pp*pp*v*v;
			fac=4*h*h*pp*pp*(1+2*(e-d)*(4*y*y-9*y+6));
			if(fac<0.) {
				it1=2+sqrt(fabs(fac))/dt;
			}else{
				it1=1;
			}
			for(it=it1;it<nt;it++) {
				tn=it*dt;
				tn2=tn*tn;
				A=sqrt(1. + fac/tn2);
				phase = w*tn*A;
				cfac = complex( cos(phase)/A, sin(phase)/A ); 		
				csum += p[it]*cfac;	
			}
		}
		pwork[iw]=csum;
	}
	pfacc(-1,ntfft,pwork);
	for(it=0;it<nt;it++) p[it]=pwork[it]*scale;
	return;
}

/* Library version of SUTIHALEDMO
 *
 * The main program (the parameters, the traces and their headers, the file of the VND library that held the transform of
 * the traces) is left to the caller.
 *
 *	su_tihaledmo:	DMO of one common-offset gather data[nx][nt] (NMO corrected) in the f-k domain, with the operators
 *			option 1 to 6 of the program: Hale, Bleistein's true amplitude, Bleistein's cos*cos weighted,
 *			Zhang, Tsvankin (with vnmo(p)), and Tsvankin's VTI weak anisotropy approximation
 *
 * Differences from the program:
 *  - The program kept nxmax + 1 traces: the last could be past the length of the transform (nxmax + nxpad is not enough
 *    when nxpad is small) and be lost; the number of traces is nx.
 *  - The arrays are in memory, as the program's file of the transform was a means of not having it so.
 */
int su_tihaledmo(int nx, int nt, float dt, float dx, float h, float v, int option, int ntpad, int nxpad, int np, float dp,
	const float *vnmo, float e, float d, const float *data, float *out)
{
	int ntfft, nxfft, nk, ik, it, ix;
	float k, dk, scale;
	float *rt;
	complex *crt, *ctemp, *ck;

	if (nx<1 || nt<2 || dt<=0.0 || dx==0.0 || h<0.0 || option<1 || option>6 || ntpad<0 || nxpad<0) return -3;
	if ((option==3 || option==6) && v<=0.0) return -3;
	if (option==5 && (np<3 || dp<=0.0 || vnmo==NULL)) return -3;

	ntfft = npfao(nt+ntpad,2*(nt+ntpad));
	nxfft = npfar(nx+nxpad);
	if (nxfft<nx) return -3;
	dk = 2.0*PI/(nxfft*dx);
	scale = 1.0/nxfft;
	nk = nxfft/2+1;

	rt = ealloc1float(nxfft+2);
	crt = ealloc1complex(nk > nt ? nk : nt);
	ctemp = ealloc1complex(ntfft);
	ck = ealloc1complex((size_t)nk*nt);

	/* forward FFT from x to k, for each time sample */
	for (it=0; it<nt; ++it) {
		for (ix=0; ix<nx; ++ix) rt[ix] = data[(size_t)ix*nt+it]*scale;
		for (ix=nx; ix<nxfft; ++ix) rt[ix] = 0.0;
		pfarc(1,nxfft,rt,crt);
		for (ik=0; ik<nk; ++ik) ck[(size_t)ik*nt+it] = crt[ik];
	}

	/* DMO correction of each wavenumber component */
	for (ik=0; ik<nk; ++ik) {
		complex *p = ck+(size_t)ik*nt;

		k = ik*dk;
		if (option==1) haledmo(h,k,p,ctemp,nt,ntfft,dt,v);
		else if (option==2) bleisteindmo(h,k,p,ctemp,nt,ntfft,dt,v);
		else if (option==3) bleisteindmo2(h,k,p,ctemp,nt,ntfft,dt,v);
		else if (option==4) zhangdmo(h,k,p,ctemp,nt,ntfft,dt,v);
		else if (option==5) tsvankindmo(h,k,p,ctemp,nt,ntfft,dt,np,dp,(float*)vnmo);
		else weakdmo(h,k,p,ctemp,nt,ntfft,dt,e,d,v);
	}

	/* inverse FFT from k to x */
	for (it=0; it<nt; ++it) {
		for (ik=0; ik<nk; ++ik) crt[ik] = ck[(size_t)ik*nt+it];
		pfacr(-1,nxfft,crt,rt);
		for (ix=0; ix<nx; ++ix) out[(size_t)ix*nt+it] = rt[ix];
	}

	free1float(rt);
	free1complex(crt);
	free1complex(ctemp);
	free1complex(ck);
	return 0;
}
