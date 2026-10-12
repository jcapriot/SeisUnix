/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUINTERP: $Revision: 1.15 $ ; $Date: 2011/11/16 22:12:22 $                */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>
#include <algorithm>

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = { 
"									",
" SUINTERP - interpolate traces using automatic event picking		",
"									",
"           suinterp < stdin > stdout					",
"									",
" ninterp=1    number of traces to output between each pair of input traces",
" nxmax=500    maximum number of input traces				",
" freq1=4.     starting corner frequency of unaliased range		",
" freq2=20.    ending corner frequency of unaliased range		",
" deriv=0      =1 means take vertical derivative on pick section        ",
"              (useful if interpolating velocities instead of seismic)  ",
" linear=0     =0 means use 8 point sinc temporal interpolation         ",
"              =1 means use linear temporal interpolation               ",
"              (useful if interpolating velocities instead of seismic)  ",
" lent=5       number of time samples to smooth for dip estimate	",
" lenx=1       number of traces to smooth for dip estimate		",
" lagc=400     number of ms agc for dip estimate			",
" xopt=0       0 compute spatial derivative via FFT			",
"                 (assumes input traces regularly spaced and relatively	",
"                  noise-free)						",
"              1 compute spatial derivative via differences		",
"                 (will work on irregulary spaced data)			",
" iopt=0     0 = interpolate",
"            1 = output low-pass model: useful for QC if interpolator failing",
"            2 = output dip picks in units of samples/trace		",
"									",
" verbose=0	verbose = 1 echoes information				",
"									",
" tmpdir= 	 if non-empty, use the value as a directory path	",
"		 prefix for storing temporary files; else if the	",
"	         the CWP_TMPDIR environment variable is set use		",
"	         its value for the path; else use tmpfile()		",
" 									",
" Notes:								",
" This program outputs 'ninterp' interpolated traces between each pair of",
" input traces.  The values for lagc, freq1, and freq2 are only used for",
" event tracking. The output data will be full bandwidth with no agc.  The",
" default parameters typically will do a satisfactory job of interpolation",
" for dips up to about 12 ms/trace.  Using a larger value for freq2 causes",
" the algorithm to do a better job on the shallow dips, but to fail on the",
" steep dips.  Only one dip is assumed at each time sample between each pair",
" of input traces.							",
" 									",
" The key assumption used here is that the low frequency data are unaliased",
" and can be used for event tracking. Those dip picks are used to interpolate",
" the original full-bandwidth data, giving some measure of interpolation",
" at higher frequencies which otherwise would be aliased.  Using iopt equal",
" to 1 allows you to visually check whether the low-pass picking model is",
" aliased.								",
" 									",
" Trace headers for interpolated traces are not updated correctly.	",
" The output header for an interpolated traces equals that for the preceding",
" trace in the original input data.  The original input traces are passed",
" through this module without modification.				",
"									",
" The place this code is most likely to fail is on the first breaks.	",
"									",
" Example run:    suplane | suinterp | suxwigb &			",
"									",
NULL};
#endif 

/*
 * Credit: John Anderson (visiting scholar from Mobil) July 1994
 *
 * Trace header fields accessed: ns, dt
 */

/**************** end self doc ********************************/

static void runav(int n, int len, const float *a, float *b);

/* The program on traces in memory: interpolate input data in space, placing ninterp synthetic traces between each
 * pair of the original traces.
 *
 * nt, nx		number of time samples, and of traces; data[nx][nt]
 * ninterp		number of traces to interpolate between each input trace
 * freq1, freq2		the low and high end frequencies (Hz) for picking
 * lagc			length of the AGC operator for picking (ms)
 * lent, lenx		length of the time and space smoothers (samples) for the picker
 * xopt			1 = differences for the spatial derivative, 0 = FFT (differences for less than 16 traces)
 * dt			sample interval (s)
 * iopt			0 = interpolate: 1 + (nx - 1) (1 + ninterp) traces (nx if there is nothing to interpolate)
 *			1 = the low-pass model: nx traces
 *			2 = the dip picks in samples per trace: nx traces
 * deriv		1 = take the vertical derivative for the picking section
 * linear		1 = linear (not 8 point sinc) interpolation in time
 * out[nout][nt]	the traces
 *
 * Returns nout. The original traces are passed through unchanged.
 */
int su_interp_traces(int nt, int nx, int ninterp, const float *data, float freq1, float freq2, int lagc, int lent,
	int lenx, int xopt, float dt, int iopt, int deriv, int linear, float *out)
{
	int ntfft,ntfftny,nxfft,nxfftny,j,k,ix,it,ixm;
	float df,dff,wa,wb,dxx,eps=1.0e-30,f,fcl,fch;

	lent=1+2*(lent/2);
	lenx=1+2*(lenx/2);
	lagc=1 + lagc*0.001/dt;

	ntfft=npfar(nt);
	ntfftny=1+ntfft/2;
	nxfft=npfar(nx);
	nxfftny=1+nxfft/2;
	df=1./(ntfft*dt);

	if(nx<2 || (iopt==0 && ninterp==0)) {
		memcpy((void *) out, (const void *) data, (size_t) nx*nt*sizeof(float));
		return nx;
	}

	std::vector<complex> crt(std::max(ntfftny, nxfftny)), ccrt(ntfftny);
	std::vector<float> rt(std::max(ntfft, nxfft)+2), rrt(ntfft+2);
	std::vector<float> a(std::max(nx,nt)), b(std::max(nx,nt)), p(nt), time(nt), aa(std::max(nx,nt)),
		bb(std::max(nx,nt)), tin(nt);
	std::vector<float> vnda((size_t) nx*nt), vndb((size_t) nx*nt);
	for (j=0; j<nt; j++) tin[j]=j;

	/* loop computing filtered data for picking purposes in vnda */
	/* compute time derivative of filtered data in vndb */
	dff=2.*PI/ntfft;
	for(ix=0;ix<nx;ix++) {
		for(j=0;j<nt;j++) rt[j]=data[(size_t) ix*nt+j];
		for(j=0;j<nt;j++) a[j]=fabs(rt[j]);
		runav(nt,lagc,a.data(),b.data());
		runav(nt,lagc,b.data(),a.data());
		for(j=0;j<nt;j++) rt[j]=rt[j]/(a[j]+eps);
		for(j=nt;j<ntfft;j++) rt[j]=0.;
		if(deriv) {
			for(j=nt-1;j>0;j--){
				rt[j]=rt[j]-rt[j-1];
			}
			rt[0]=rt[1];
		}
		pfarc(1,ntfft,rt.data(),crt.data());
		if(freq1>0.){
			for(j=0;j<ntfftny;j++){
				f=j*df;
				fcl=(f/freq1);
				fcl=fcl*fcl*fcl*fcl;
				fch=(f/freq2);
				fch=fch*fch*fch*fch;
				f=fcl/( (1.+fcl)*(1.+fch) );
				crt[j]=crt[j]*f;
				ccrt[j]=complex(crt[j].imag()*(j*dff), -crt[j].real()*(j*dff));
			}
		}else{
			for(j=0;j<ntfftny;j++){
				f=j*df;
				fch=(f/freq2);
				f=1./(1.+fch*fch*fch*fch);
				crt[j]=crt[j]*f;
				ccrt[j]=complex(crt[j].imag()*(j*dff), -crt[j].real()*(j*dff));
			}
		}
		pfacr(-1,ntfft,crt.data(),rt.data());
		memcpy(&vnda[(size_t) ix*nt], rt.data(), nt*sizeof(float));
		pfacr(-1,ntfft,ccrt.data(),rrt.data());
		memcpy(&vndb[(size_t) ix*nt], rrt.data(), nt*sizeof(float));
	}

	if(iopt==1){
		memcpy((void *) out, (const void *) vnda.data(), (size_t) nx*nt*sizeof(float));
		return nx;
	}

	/* loop computing spatial derivative of data for picking purposes*/
	dxx=2.*PI/(nxfft*nxfft);
	if(nx<16) xopt=1;
	for(it=0;it<nt;it++) {
		for(j=0;j<nx;j++) rt[j]=vnda[(size_t) j*nt+it];
		if(xopt) {
			for(j=0;j<nx-1;j++) rt[j]=rt[j+1]-rt[j];
			rt[nx-1]=rt[nx-2];
		}else{
			for(j=nx;j<nxfft;j++) rt[j]=0.;
			pfarc(1,nxfft,rt.data(),crt.data());
			for(j=0;j<nxfftny;j++){
				crt[j]=complex(crt[j].imag()*(j*dxx), -crt[j].real()*(j*dxx));
			}
			pfacr(-1,nxfft,crt.data(),rt.data());
		}
		for(j=0;j<nx;j++) vnda[(size_t) j*nt+it]=rt[j];
	}

	/* compute dot products and smooth over time */
	for(ix=0;ix<nx;ix++) {
		for(it=0;it<nt;it++) {
			a[it]=vnda[(size_t) ix*nt+it];
			b[it]=vndb[(size_t) ix*nt+it];
		}
		for(it=0;it<nt;it++) {
			aa[it]=a[it]*b[it];
			bb[it]=b[it]*b[it];
		}
		runav(nt,lent,aa.data(),a.data());
		runav(nt,lent,a.data(),aa.data());
		runav(nt,lent,bb.data(),b.data());
		runav(nt,lent,b.data(),bb.data());
		memcpy(&vnda[(size_t) ix*nt], aa.data(), nt*sizeof(float));
		memcpy(&vndb[(size_t) ix*nt], bb.data(), nt*sizeof(float));
	}

	/* smooth dot products in x */
	if(lenx>1){
		for(it=0;it<nt;it++) {
			for(ix=0;ix<nx;ix++) {
				a[ix]=vnda[(size_t) ix*nt+it];
				b[ix]=vndb[(size_t) ix*nt+it];
			}
			runav(nx,lenx,a.data(),aa.data());
			runav(nx,lenx,aa.data(),a.data());
			runav(nx,lenx,b.data(),bb.data());
			runav(nx,lenx,bb.data(),b.data());
			for(ix=0;ix<nx;ix++) {
				vnda[(size_t) ix*nt+it]=a[ix];
				vndb[(size_t) ix*nt+it]=b[ix];
			}
		}
	}

	/* loop computing p, interpolating, and outputting results */
	for(it=0;it<nt;it++) a[it]=data[it];
	for(ix=1;ix<nx;ix++) {
		ixm=ix-1;
		for(it=0;it<nt;it++) {
			p[it] = - vnda[(size_t) ixm*nt+it]/( vndb[(size_t) ixm*nt+it] + eps );
		}
		for(it=0;it<nt;it++) b[it]=data[(size_t) ix*nt+it];
		if(iopt==2) {
			/* don't output dip picks except on original traces */
			memcpy(&out[(size_t) ixm*nt], p.data(), nt*sizeof(float));
		}else{
			memcpy(&out[(size_t) ixm*(ninterp+1)*nt], a.data(), nt*sizeof(float));
			for(k=0;k<ninterp;k++){
				wa=(1.+k)/(1+ninterp);
				wb=1.-wa;
				for(it=0;it<nt;it++) time[it] = ((float ) it) - p[it]*wa;
				if(linear){
					intlin(nt,tin.data(),a.data(),a[0],a[nt-1],nt,time.data(),aa.data());
				}else{
					ints8r(nt,1.0,0.,a.data(),0.0,0.0,nt,time.data(),aa.data());
				}
				for(it=0;it<nt;it++) time[it] = ((float ) it) + p[it]*wb;
				if(linear){
					intlin(nt,tin.data(),b.data(),b[0],b[nt-1],nt,time.data(),bb.data());
				}else{
					ints8r(nt,1.0,0.,b.data(),0.0,0.0,nt,time.data(),bb.data());
				}
				for(it=0;it<nt;it++)
					aa[it]=wb*aa[it]+wa*bb[it];
				memcpy(&out[(size_t) (k+1+ixm*(ninterp+1))*nt], aa.data(), nt*sizeof(float));
			}
		}
		std::swap(a, b);
	}
	if(iopt==2) {
		memcpy(&out[(size_t) (nx-1)*nt], p.data(), nt*sizeof(float));
		return nx;
	}
	memcpy(&out[(size_t) (nx-1)*(ninterp+1)*nt], a.data(), nt*sizeof(float));
	return 1+(nx-1)*(ninterp+1);
}

/* compute a boxcar running average filter: n samples in a[] and b[], len the length of the average */
static void runav(int n, int len, const float *a, float *b)
{
	float sum=0.;
	int j,lenh=len/2;

	if(len<=1) {
		for(j=0;j<n;j++) b[j]=a[j];
		return;
	}
	for(j=0;j<MIN(len,n);j++) sum+=a[j];
	for(j=0;j<MIN(lenh+1,n);j++) b[j]=sum;
	for(j=lenh+1;j<n-lenh;j++) {
		sum=sum+a[j+lenh]-a[j-lenh-1];
		b[j]=sum;
	}
	for(j=MAX(0,n-lenh);j<n;j++) b[j]=sum;
	sum=1./len;
	for(j=0;j<n;j++) b[j]*=sum;
	return;
}
