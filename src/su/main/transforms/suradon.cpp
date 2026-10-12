/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SURADON: $Revision: 1.21 $ ; $Date: 2015/08/07 22:26:35 $	*/

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>
#include <algorithm>

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SURADON - compute forward or reverse Radon transform or remove multiples",
"           by using the parabolic Radon transform to estimate multiples",
"           and subtract.						",
"									",
"     suradon <stdin >stdout [Optional Parameters]			",
"									",
" Optional Parameters:							",
" choose=0    0  Forward Radon transform				",
"             1  Compute data minus multiples				",
"             2  Compute estimate of multiples				",
"             3  Compute forward and reverse transform			",
"             4  Compute inverse Radon transform			",
" igopt=1     1  parabolic transform: g(x) = offset**2			",
"             2  Foster/Mosher psuedo hyperbolic transform		",
"                   g(x) = sqrt(depth**2 + offset**2)			",
"             3  Linear tau-p: g(x) = offset				",
"             4  abs linear tau-p: g(x) = abs(offset)			",
" offref=2000.    reference maximum offset to which maximum and minimum	",
"                 moveout times are associated				",
" interoff=0.     intercept offset to which tau-p times are associated	",
" pmin=-200       minimum moveout in ms on reference offset		",
" pmax=400        maximum moveout in ms on reference offset		",
" anderson=1	  p values in ms, =0 pvalues in horizontal slowness	",
" dp=16           moveout increment in ms on reference offset		",
" pmula=80        moveout in ms on reference offset where multiples begin",
"                     at maximum time					",
" pmulb=200       moveout in ms on reference offset where multiples begin",
"                     at zero time					",
" depthref=500.   Reference depth for Foster/Mosher hyperbolic transform",
" nwin=1          number of windows to use through the mute zone	",
" f1=60.          High-end frequency before taper off			",
" f2=80.          High-end frequency					",
" prewhite=0.1    Prewhitening factor in percent.			",
" cdpkey=cdp      name of header word for defining ensemble		",
" offkey=offset   name of header word with spatial information		",
" nxmax=240       maximum number of input traces per ensemble		",
" ltaper=7	  taper (integer) for mute tapering function		",
"									",
" Optimizing Parameters:						",
" The following parameters are occasionally used to avoid spatial aliasing",
" problems on the linear tau-p transform.  Not recommended for other	",
" transforms...								",
" ninterp=0      number of traces to interpolate between each input trace",
"                   prior to computing transform			",
" freq1=4.0      low-end frequency in Hz for picking (good default: 3 Hz)",
"                (Known bug: freq1 cannot be zero) ",
" freq2=20.0     high-end frequency in Hz for picking (good default: 20 Hz)",
" lagc=400       length of AGC operator for picking (good default: 400 ms)",
" lent=5         length of time smoother in samples for picker		",
"                     (good default: 5 samples)				",
" lenx=7         length of space smoother in samples for picker		",
"                     (good default: 1 sample)				",
" xopt=1         1 = use differences for spatial derivative		",
"                        (works with irregular spacing)			",
"                0 = use FFT derivative for spatial derivatives		",
"                      (more accurate but requires regular spacing and	",
"                      at least 16 input tracs--will switch to differences",
"                      automatically if have less than 16 input traces)	",
"									",
NULL};
#endif

/* Credits:
 *	CWP: John Anderson (visitor to CSM from Mobil) Spring 1993
 *
 * Multiple removal notes:
 *	Usually the input data are NMO corrected CMP gathers.  The
 *	first pass is to compute a parabolic Radon transform and
 * 	identify the multiples in the transform domain.  Then, the
 * 	module is run on all the data using "choose=1" to estimate
 * 	and subtract the multiples.  See the May, 1993 CWP Project
 *	Review for more extensive documentation.
 *
 * NWIN notes:
 *	The parabolic transform runs with higher resolution if the
 * 	mute zone is honored.  When "nwin" is specified larger than
 *   	one (say 6), then multiple windows are used through the mute
 * 	zone.  It is assumed in this case that the input data are
 * 	sorted by the offkey header item from small offset to large
 * 	offset.  This causes the code to run 6 times longer.  The
 *      mute time is taken from the "muts" header word.
 *      You may have to manually set this header field yourself, if
 *      it is not already set.
 *
 * References:
 * Anderson, J. E., 1993, Parabolic and linear 2-D, tau-p transforms
 *       using the generalized radon tranform, in May 11-14, 1993
 *       Project Review, Consortium Project on Seismic Inverse methods
 *       for Complex Structures, CWP-137, Center for Wave Phenomena
 *       internal report.
 * Other References cited in above paper:
 * Beylkin, G,.1987, The discrete Radon transform: IEEE Transactions
 *       of Acoustics, Speech, and Signal Processing, 35, 162-712.
 * Chapman, C.H.,1981, Generalized Radon transforms and slant stacks:
 *       Geophysical Journal of the Royal Astronomical Society, 66,
 *       445-453.
 * Foster, D. J. and Mosher, C. C., 1990, Multiple supression
 *       using curvilinear Radon transforms: SEG Expanded Abstracts 1990,
 *       1647-1650.
 * Foster, D. J. and Mosher, C. C., 1992, Suppression of multiples
 *       using the Radon transform: Geophysics, 57, No. 3, 386-395.
 * Gulunay, N., 1990, F-X domain least-squares Tau-P and Tau-Q: SEG
 *       Expanded Abstracts 1990, 1607-1610.
 * Hampson, D., 1986, Inverse velocity stacking for multiple elimination:
 *       J. Can. Soc. Expl. Geophs., 22, 44-55.
 * Hampson, D., 1987, The discrete Radon transform: a new tool for image
 *       enhancement and noise suppression: SEG Expanded Abstracts 1978,
 *       141-143.
 * Johnston, D.E., 1990, Which multiple suppression method should I use?
 *       SEG Expanded Abstracts 1990, 1750-1752.
 *
 * Trace header words accessed: ns, dt, cdpkey, offkey, muts
 */
/**************** end self doc ********************************/

/* The functions of the program on gathers that are in memory (the program keeps them in temporary files of its own
   virtual memory, VND, and does the choices of the main program for each ensemble: that is done by the caller).

   Layouts: a gather is [trace][sample], a tau-p panel is [p][tau]. All arrays are the caller's. */

static void compute_r(float w, int nx, const float *g, int np, float dp, complex *r);
static void compute_rhs(float w, int nx, const float *g, const complex *data, int np, float pmin, float dp,
		complex *rhs);
static int ctoep(int n, const complex *r, complex *a, const complex *b, complex *f, complex *g);
static int ctoephcg(int niter, int n, const complex *r, complex *x, const complex *y, complex *s, complex *ss,
		complex *g, complex *rr);
static float rcdot(int n, const complex *a, const complex *b);
static void htmul(int n, const complex *a, const complex *x, complex *y);
static void runav(int n, int len, const float *a, float *b);

float su_radon_gofx(int igopt, float offset, float intercept_off, float refdepth)
/*******************************************************************
return g(x) for various options
******************************************************************
igopt		1 = parabolic transform
		2 = Foster/Mosher pseudo hyperbolic option
		3 = linear tau-p
		4 = linear tau-p using absolute value of offset
offset		offset in m
intercept_off	offset corresponding to intercept time
refdepth	reference depth in m for igopt=2
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	offset=offset-intercept_off;
	if(igopt==1) {
		return(offset*offset);
	}
	if(igopt==2) {
		return( sqrt(refdepth*refdepth + offset*offset) );
	}
	if(igopt==3) {
		return(offset);
	}
	if(igopt==4) {
		return(fabs(offset));
	}
	return(offset);
}

static float freqweight(int j, float df, float f1, float f2)
/* weight for each frequency: j freq index, df freq increment, f1 taper off freq, f2 freq beyond which all are 0 */
{
	float w;
	float f=j*df;
	if(f<=f1) return (1.);
	if(f>=f2) return (0.);
	w = (f2-f)/(f2-f1);
	return (w);
}

void su_radon_taupmute(int ip, int ipa, int ipb, int nt, int ltap, float *rt)
/*******************************************************************
do simple tau-p mute to elliminate multiples
******************************************************************
ip		current ray parameter index
ipa		max ray parameter primary  pick at maximum time
ipb		max ray parmater primary pick at minimum time
nt		number of time samples
ltap		length of mute taper in samples
rt[nt]		tau-p data for all tau values

The program's taper is (k-j)/ltap with integers, which is 0 for all of the
samples that it is applied to: it is a ramp here, as the intent is.
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	int j,k;
	float w;
	if(ip>=ipb) return;
	if(ip<=ipa) {
		for(k=0;k<nt;k++) rt[k]=0;
		return;
	}
	w=MAX(ipb-ipa,1);
	w=(ipb-ip)/w;
	j=w*nt;
	for(k=0;k<j;k++) rt[k]=0.;
	for(k=j;k<MIN(nt,j+ltap);k++) rt[k]*=(float)(k-j)/ltap;
}

int su_radon_forward(int nx, int nt, const float *data, const float *g, float dt, int ntfft, int np, float pmin,
	float dp, const float *mutetime, int nk, float f1, float f2, float prewhite, float *ptau)
/*******************************************************************
do forward generalized radon transform
******************************************************************
nx		number of input traces
nt		number of input time samples
data[nx][nt]	input data in time-space domain
g[nx]		g[ix]=offset**2 for parabolic transform
dt		sample rate in seconds
ntfft		length of time fft and output tau-p data in samples
np		number of generalized ray parameters
pmin		minimum generalized ray parameter
dp		generalized ray parameter increment
mutetime[nx]	array of mute times in ms (used for nk > 1)
nk		number of offset ranges to transform separately
		through the mute zone
f1		max freq without taper
f2		max non-zero freq component
prewhite	0.01 means prewhiten 1 percent
ptau[np][ntfft]	the transform

key assumption: offsets are sorted to increase with index
Returns the number of coefficients that the least squares solutions gave
for the last frequency that was solved, or 0.
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	int ix,ip,it,iw,j,k,nxx,nxxinc,ik,ik2,ntfftny;
	float w,dw,rsum,fac,wa,wb,rk[2],rit[2],df;
	int status=0;
	const complex czero(0.0f, 0.0f);

	fac=1./ntfft;
	ntfftny=1+ntfft/2;
	df=1./(ntfft*dt);
	dw=2.*PI*df;
	nxxinc=1+(nx-1)/nk;

	std::vector<float> kindex(nk), tindex(nk);
	if(nk>1) {
		/* a set of (mute time, group index) pairs */
		for(k=0;k<nk;k++) {
			kindex[k]=k;
			nxx=MIN(nx-1,(int)((k+0.75)*nxxinc));
			tindex[k]=0.001*mutetime[nxx]/dt;
		}
	}

	std::vector<float> rt(ntfft+2);
	std::vector<complex> spectra((size_t) nx*ntfftny), crt(std::max(nx, ntfftny)), ccrt((size_t) (nk+1)*np),
		r(np), rhs(np), wrk1(np), wrk2(np), wrk3(np), wrk4(np);
	/* the frequency by (np*nk) panel, and its transform to tau */
	std::vector<complex> byfreq((size_t) ntfftny*np*nk);
	std::vector<float> bytau((size_t) np*nk*ntfft);

	/* do forward time to frequency fft */
	for(ix=0;ix<nx;ix++) {
		for(it=0;it<nt;it++) rt[it]=data[(size_t) ix*nt+it]*fac;
		for(j=nt;j<ntfft;j++) rt[j]=0.;
		pfarc(1,ntfft,rt.data(),&spectra[(size_t) ix*ntfftny]);
	}

	/* do radon transform, frequency by frequency, for multiple spatial windows */
	for(iw=0;iw<ntfftny;iw++) {
		wa=freqweight(iw,df,f1,f2);
		if(wa>0.) {
		    w=iw*dw;
		    for(ix=0;ix<nx;ix++) crt[ix]=spectra[(size_t) ix*ntfftny+iw];
		    if(wa<1.) {
		    	for(ix=0;ix<nx;ix++) crt[ix]=crt[ix]*wa;
		    }
		    for(k=0;k<nk;k++) {
			nxx=MIN(nx,(k+1)*nxxinc);
			compute_rhs(w,nxx,g,crt.data(),np,pmin,dp,rhs.data());
			compute_r(w,nxx,g,np,dp,r.data());
			r[0]=complex(r[0].real()*(1.+prewhite), r[0].imag());
			for(rsum=0.,j=1;j<np;j++)
				rsum +=sqrt(r[j].real()*r[j].real() + r[j].imag()*r[j].imag());
			rsum=rsum/r[0].real();
			if (rsum>1.+np/5) {
				j=ctoephcg(np/7,np,r.data(),&ccrt[(size_t) k*np],rhs.data(),
					wrk1.data(),wrk2.data(),wrk3.data(),wrk4.data());
			}else{
				j=ctoep(np,r.data(),&ccrt[(size_t) k*np],rhs.data(),wrk1.data(),wrk2.data());
			}
			status=j;
		    }
		}else{
		    for(ip=0;ip<np*nk;ip++) ccrt[ip]=czero;
		}
		for(ip=0;ip<np*nk;ip++) byfreq[(size_t) ip*ntfftny+iw]=ccrt[ip];
	}

	/* do fourier transform from frequency to tau */
	for(ip=0;ip<np*nk;ip++) {
		pfacr(-1,ntfft,&byfreq[(size_t) ip*ntfftny],&bytau[(size_t) ip*ntfft]);
	}

	/* merge appropriate tau-p transform for each window using mute zone information */
	if(nk>1) {
		for(it=0;it<ntfft;it++) {
			rit[0]=it;
			intlin(nk,tindex.data(),kindex.data(),kindex[0],kindex[nk-1],1,rit,rk);
			ik=rk[0];
			ik2=MIN(ik+1,nk-1);
			wb=rk[0]-ik;
			wa=1-wb;
			for(ip=0;ip<np;ip++) {
				ptau[(size_t) ip*ntfft+it]=wa*bytau[(size_t) (ik*np+ip)*ntfft+it]
					+ wb*bytau[(size_t) (ik2*np+ip)*ntfft+it];
			}
		}
	} else {
		memcpy((void *) ptau, (const void *) bytau.data(), (size_t) np*ntfft*sizeof(float));
	}
	return status;
}

void su_radon_inverse(int np, int ntfft, int nx, const float *ptau, const float *g, float dt, float pmin, float dp,
	int ip1, float f1, float f2, float *out)
/*******************************************************************
do inverse generalized radon transform
******************************************************************
np		number of generalized ray parameters
ntfft		length of time fft and input tau-p data in samples
nx		number of output traces
ptau[np][ntfft]	the data in the tau-p domain
g[nx]		g[ix]=offset**2 for parabolic transform
dt		sample rate in seconds
pmin		minimum generalized ray parameter
dp		generalized ray parameter increment
ip1		starting generalized ray parameter index to compute
		(use as 0 for full inversion, positive integer to
		just invert multiples)
f1		max freq without taper
f2		max non-zero freq component
out[nx][ntfft]	the data in the time-space domain
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	int ip,iw,ntfftny,ix,it;
	float w,dw,p,rsum,isum,dr,di,tr,ti,fac,wa,df;
	const complex czero(0.0f, 0.0f);

	ntfftny=1+ntfft/2;
	df=1./(ntfft*dt);
	dw=2.*PI*df;

	std::vector<float> rt(ntfft+2);
	std::vector<complex> pspectra((size_t) np*ntfftny), xspectra((size_t) nx*ntfftny), ctemp(np);

	fac=1./ntfft;
	for(ip=0;ip<np;ip++) {
		for(it=0;it<ntfft;it++) rt[it]=ptau[(size_t) ip*ntfft+it]*fac;
		pfarc(1,ntfft,rt.data(),&pspectra[(size_t) ip*ntfftny]);
	}

	fac=1./np;
	for(iw=0;iw<ntfftny;iw++) {
		wa=freqweight(iw,df,f1,f2);
		if(wa>0.) {
			w=iw*dw;
			for(ip=0;ip<np;ip++) {
				complex c=pspectra[(size_t) ip*ntfftny+iw];
				if(wa<1.) c=c*wa;
				ctemp[ip]=c;
			}
			for(ix=0;ix<nx;ix++) {
			    rsum = isum = 0.;
			    for(ip=ip1;ip<np;ip++) {
				p = pmin + ip*dp;
				tr = cos(w*p*g[ix]);
				ti = sin(w*p*g[ix]);
				dr = ctemp[ip].real();
				di = ctemp[ip].imag();
				rsum += tr*dr - ti*di;
				isum += tr*di + ti*dr;
			    }
			    xspectra[(size_t) ix*ntfftny+iw]=complex(fac*rsum, fac*isum);
			}
		}else{
			for(ix=0;ix<nx;ix++) xspectra[(size_t) ix*ntfftny+iw]=czero;
		}
	}
	for(ix=0;ix<nx;ix++) {
		pfacr(-1,ntfft,&xspectra[(size_t) ix*ntfftny],&out[(size_t) ix*ntfft]);
	}
}

static void compute_r(float w, int nx, const float *g, int np, float dp, complex *r)
/*******************************************************************
Compute the top row of the Hermitian Toeplitz Matrix
			+
		  R = B B

		  i w p g(x)
where B = (1/np) e	    for equal increments in p as

     +           -i w p g(x)
and B = (1/nx) e

as used for the Discrete Radon Transform computation for
linear or parabolic tau-p.
******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
******************************************************************/
{
	int j,k;
	float rsum, isum, fac;
	fac= 1./(nx*np);
	for(j=0;j<np;j++) {
		rsum=0.;
		isum=0.;
		for(k=0;k<nx;k++) {
			rsum = rsum+cos( w*j*dp*g[k] );
			isum = isum+sin( w*j*dp*g[k] );
		}
		r[j]=complex(fac*rsum, fac*isum);
	}
}

static void compute_rhs(float w, int nx, const float *g, const complex *data, int np, float pmin, float dp,
		complex *rhs)
/*********************************************************************
				     +
Compute the right-hand-side vector  B  data(x)

	+	    -i w p g(x)
where B   = (1/nx) e	        for equal increments in p as
used for the Discrete Radon Transform computation for
linear or parabolic tau-p.
*********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*********************************************************************/
{
	int ip, ix;
	float p, rsum, isum, dr, di, tr, ti, fac;
	fac=1./nx;
	for(ip=0;ip<np;ip++) {
		p = pmin + ip*dp;
		rsum = isum = 0.;
		for(ix=0;ix<nx;ix++) {
			tr = cos(w*p*g[ix]);
			ti = -sin(w*p*g[ix]);
			dr = data[ix].real();
			di = data[ix].imag();
			rsum += tr*dr - ti*di;
			isum += tr*di + ti*dr;
		}
		rhs[ip]=complex(fac*rsum, fac*isum);
	}
}

static int ctoep(int n, const complex *r, complex *a, const complex *b, complex *f, complex *g)
/***********************************************************************
Complex Hermitian Toeplitz Solver for

N-1
Sum  R	     A  = B      for i=0,1,2,...,N-1
j=0   (i-j)   j    i

where R is Hermitian Toeplitz and A and B are complex.

n     dimension of system
*r    provides the top row of the Hermitian Toeplitz matrix R
*a    returns the complex solution vector A
*b    input as complex vector B (not changed during call)
*f    returns the complex spiking filter F
*g    work space of length n complex values

The function value returns as the number of successfully
computed complex filter coefficients (up to n) if successful or
0 if no coefficients could be computed.
***********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
***********************************************************************/
{
	float er, ei, vr, vi, cr, ci, vsq;
	int j;	  	/*  for the jth iteration, j=0,n-1 	*/
	int k;			/*  for the kth component, k=0,j-1 	*/
	int jmk;		/*  j-k 				*/
	if (r[0].real()==0.) return 0;

	f[0]=complex(1.0/r[0].real(), 0.0);
	a[0]=complex(b[0].real()/r[0].real(), b[0].imag()/r[0].real());
	vr=1.;
	vi=0.;
	vsq=1.;

	for(j=1;j<n;j++) {     	/* iteration loop for iteration j	*/
	/*  	Compute spiking filter that outputs {v,0,0,...}
		for this iteration step j			*/
		f[j]=complex(0.0f, 0.0f);
		er=ei=0.;
		for(k=0;k<j;k++) {
			jmk=j-k;
			er+=r[jmk].real()*f[k].real()+r[jmk].imag()*f[k].imag();
			ei+=r[jmk].real()*f[k].imag()-r[jmk].imag()*f[k].real();
		}
		cr  = (er*vr - ei*vi)/vsq;
		ci  = (er*vi + ei*vr)/vsq;
		vr  = vr - (cr*er+ci*ei);
		vi  = vi + (cr*ei-ci*er);
		vsq =  vr*vr + vi*vi;
		if (vsq <= 0.) break;
		for(k=0;k<=j;k++) {
			jmk=j-k;
			g[k]=complex(f[k].real() - cr*f[jmk].real() - ci*f[jmk].imag(),
				f[k].imag() + cr*f[jmk].imag() - ci*f[jmk].real());
		}
		for(k=0;k<=j;k++) {
			f[k]=g[k];
		}

		/*  Compute shaping filter for this iteration */
		a[j]=complex(0.0f, 0.0f);
		er=ei=0.;
		for(k=0;k<j;k++) {
			jmk=j-k;
			er+=r[jmk].real()*a[k].real()+r[jmk].imag()*a[k].imag();
			ei+=r[jmk].real()*a[k].imag()-r[jmk].imag()*a[k].real();
		}
		er  = er-b[j].real();
		ei  = ei-b[j].imag();
		cr  = (er*vr - ei*vi)/vsq;
		ci  = (er*vi + ei*vr)/vsq;
		for(k=0;k<=j;k++) {
			jmk=j-k;
			a[k]=complex(a[k].real() - cr*f[jmk].real() - ci*f[jmk].imag(),
				a[k].imag() + cr*f[jmk].imag() - ci*f[jmk].real());
		}
	}

	/* (the program goes on to normalize the spiking filter f, which is not used by the caller) */
	return (j);
}

static int ctoephcg(int niter, int n, const complex *a, complex *x, const complex *y,
	complex *s, complex *ss, complex *g, complex *rr)
/*********************************************************************
Hestenes and Stiefel conjugate gradient algorithm
specialized for solving Hermitian Toeplitz
system.  a[] is input as a vector defining the only the
top row of A.  x[] is the solution vector returned.
y[] is input.  niter is the maximum number of conjugate
gradient steps to compute.  The function returns as
the number of steps actually computed.  The other
vectors provide workspace.
********************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*********************************************************************/
{
	int j, iter;
	const complex czero(0.0f, 0.0f);
	float alpha, beta, gamma, gammam, rsq, rp, test;
	float eps=1.0e-6;
	rp   = rcdot(n,y,y);
	test = n*eps*eps*rp;

	for(j=0;j<n;j++) {
		x[j]=czero;
		rr[j]=y[j];
	}
	htmul(n,a,rr,g);	   /*  adjoint matrix multiply */

	for(j=0;j<n;j++) s[j]=g[j];
	gammam=rcdot(n,g,g);

	for(iter=0;iter<niter;iter++) { /* forward matrix multiply  */
		htmul(n,a,s,ss);
		alpha  = gammam/rcdot(n,ss,ss);
		for(j=0;j<n;j++) {
			x[j] =x[j]+s[j]*alpha;
			rr[j]=rr[j]-ss[j]*alpha;
		}
		rsq = rcdot(n,rr,rr);
		if ( iter>0 && ( rsq==rp || rsq<test ) ) return(iter-1);
		rp = rsq;

		htmul(n,a,rr,g);   /*  adjoint matrix multiply  */
		gamma  = rcdot(n,g,g);
		if (gamma<eps) break;
		beta   = gamma/gammam;
		gammam = gamma;

		for(j=0;j<n;j++) {
			s[j] =g[j]+s[j]*beta;
		}
	}
	return(iter);
}

static float rcdot(int n, const complex *a, const complex *b)
/* the real part of a complex dot product where the first vector is the one complex conjugated */
{
	int j;
	float sum=0.;
	for(j=0;j<n;j++) sum += a[j].real() * b[j].real() + a[j].imag() * b[j].imag();
	return(sum);
}

static void htmul(int n, const complex *a, const complex *x, complex *y)
/* Hermitian Toeplitz matrix multiply y = A x, where A is defined by the vector a giving the top row of A */
{
	int j,irow;

	for(irow=0;irow<n;irow++) {
		float yr=0.0f, yi=0.0f;
		for(j=0;j<irow;j++) {
			/* conjg(a[irow-j])*x[j] + y */
			float ar=a[irow-j].real(), ai=-a[irow-j].imag(), xr=x[j].real(), xi=x[j].imag();
			yr = (ar*xr - ai*xi) + yr;
			yi = (ar*xi + ai*xr) + yi;
		}
		for(j=irow;j<n;j++) {
			float ar=a[j-irow].real(), ai=a[j-irow].imag(), xr=x[j].real(), xi=x[j].imag();
			yr = (ar*xr - ai*xi) + yr;
			yi = (ar*xi + ai*xr) + yi;
		}
		y[irow]=complex(yr, yi);
	}
}

void su_radon_xinterp(int nt, int nx, int ninterp, const float *data, float freq1, float freq2, int lagc,
		int lent, int lenx, int xopt, float dt, float *out)
/*******************************************************************
interpolate input data in space placing "ninterp" synthetic traces
between each pair of original input traces
******************************************************************
nt			number of time samples
nx			number of input traces
ninterp			number of traces to interpolate between each input trace
data[nx][nt]		the input data
freq1			low-end frequency in Hz for picking (good default: 3 Hz)
freq2			high-end frequency in Hz for picking (good default: 20 Hz)
lagc			length of AGC operator for picking (good default: 400 ms)
lent			length of time smoother in samples for picker
lenx			length of space smoother in samples for picker
xopt			1 = use differences for spatial derivative
			0 = use FFT derivative for spatial derivatives
			    (switches to differences with less than 16 traces)
dt			sample rate in sec
out[1+(nx-1)*(1+ninterp)][nt]	the output: original plus interpolated data

The place this code is most likely to fail is on the first breaks.
(The options iopt=1 and 2 of the program's routine, which output the model
for picking and the dips, are not here.)
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
{
	int	ntfft,ntfftny,nxfft,nxfftny,j,k,ix,it,ixm;
	float	df,dff,wa,wb,dxx,eps,f,fcl,fch;

	lent=1+2*(lent/2);
	lenx=1+2*(lenx/2);
	lagc=1 + lagc*0.001/dt;

	ntfft=npfar(nt);
	ntfftny=1+ntfft/2;
	nxfft=npfar(nx);
	nxfftny=1+nxfft/2;

	df=1./(ntfft*dt);

	if(nx<2 || ninterp==0) {
		memcpy((void *) out, (const void *) data, (size_t) nx*nt*sizeof(float));
		return;
	}

	std::vector<complex> crt(std::max(ntfftny, nxfftny)), ccrt(ntfftny);
	std::vector<float> rt(std::max(ntfft, nxfft)+2), rrt(ntfft+2);
	std::vector<float> a(std::max(nx,nt)), b(std::max(nx,nt)), p(nt), time(nt), aa(std::max(nx,nt)),
		bb(std::max(nx,nt));
	/* the data for picking, and its time derivative */
	std::vector<float> vnda((size_t) nx*nt), vndb((size_t) nx*nt);

	eps = 0;
	for(ix=0;ix<nx;ix++)
		for(j=0;j<nt;j++) eps = MAX(eps, fabs(data[(size_t) ix*nt+j]));
	eps = MAX(1e-2*eps, 1e-30);

	/* loop computing filtered data for picking purposes in vnda */
	/* compute time derivative of filtered data in vndb */
	dff=2.*PI/ntfft;
	for(ix=0;ix<nx;ix++) {
		for(j=0;j<nt;j++) rt[j]=data[(size_t) ix*nt+j];
		for(j=0;j<nt;j++) a[j]=fabs(rt[j]);
		runav(nt,lagc,a.data(),b.data());
		runav(nt,lagc,b.data(),a.data());
		for(j=0;j<nt;j++) rt[j]=rt[j]/(a[j]+ eps);
		for(j=nt;j<ntfft;j++) rt[j]=0.;
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

	eps = 0;
	for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++) eps = MAX(eps, vndb[(size_t) ix*nt+it]);
	eps = MAX(1e-5*eps, 1e-30);

	/* loop computing p, interpolating, and outputting results */
	for(it=0;it<nt;it++) a[it]=data[it];
	for(ix=1;ix<nx;ix++) {
		ixm=ix-1;
		for(it=0;it<nt;it++) {
			p[it] = - vnda[(size_t) ixm*nt+it]/( vndb[(size_t) ixm*nt+it] + eps );
		}
		for(it=0;it<nt;it++) b[it]=data[(size_t) ix*nt+it];
		memcpy(&out[(size_t) ixm*(ninterp+1)*nt], a.data(), nt*sizeof(float));
		for(k=0;k<ninterp;k++){
			wa=(1.+k)/(1+ninterp);
			wb=1.-wa;
			for(it=0;it<nt;it++) time[it] = it - p[it]*wa;
			ints8r(nt,1.0,0.,a.data(),0.0,0.0,nt,time.data(),aa.data());
			for(it=0;it<nt;it++) time[it] = it + p[it]*wb;
			ints8r(nt,1.0,0.,b.data(),0.0,0.0,nt,time.data(),bb.data());
			for(it=0;it<nt;it++)
				aa[it]=wb*aa[it]+wa*bb[it];
			memcpy(&out[(size_t) (k+1+ixm*(ninterp+1))*nt], aa.data(), nt*sizeof(float));
		}
		std::swap(a, b);
	}
	memcpy(&out[(size_t) (nx-1)*(ninterp+1)*nt], a.data(), nt*sizeof(float));
}

static void runav(int n, int len, const float *a, float *b)
/******************************************************************
compute a boxcar running average filter
*******************************************************************
n   	number of samples in a[] and b[]
len 	length of running average in samples
a[n]	input array
b[n]	output array
*******************************************************************
Author: John Anderson (visitor to CSM from Mobil) Spring 1993
*******************************************************************/
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
