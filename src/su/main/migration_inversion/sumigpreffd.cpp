/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGPREFFD: $Vision: 1.00 $ ; $Date: 2015/08/07 22:19:43 $	*/

#include "su.h"
#include <vector>
#include "segy.h"
#include "header.h"

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
"SUMIGPREFFD - The 2-D prestack common-shot Fourier finite-difference	",
"		depth  migration.					",
"									",
"  sumigpreffd <indata >outfile [parameters]				", 
"									",
" Required Parameters:							",  
" nxo=	   number of total horizontal output samples			",
" nxshot=	number of shot gathers to be migrated			",
" nz=		number of depth samples					",
" dx=		horizontal sampling interval				",
" dz=		depth sampling interval					",
" vfile=	 velocity profile, it must be binary format.		",
"									",  
" Optional Parameters:							",
" fmax=25	the peak frequency of Ricker wavelet used as source wavelet",
" f1=5,f2=10,f3=40,f4=50	frequencies to build a Hamming window	",
" lpad=9999,rpad=9999		number of zero traces padded on both	",
"				sides of depth section to determine the ",
"				migration aperature, the default	",
"				values are using the full aperature.	",
" verbose=0		silent, =1 additional runtime information	",
"									",  
" Notes:								",
" The input velocity file consists of C-style binary floats.		",  
" The structure of this file is vfile[iz][ix]. Note that this means that",
" the x-direction is the fastest direction instead of z-direction! Such a",
" structure is more convenient for the downward continuation type	",
" migration algorithm than using z as fastest dimension as in other SU  ", 
" programs.								",
"									",
" Because most of the tools in the SU package (such as  unif2, unisam2, ", 
" and makevel) produce output with the structure vfile[ix][iz], you will",
" need to transpose the velocity files created by these programs. You may",
" use the SU program \'transp\' in SU to transpose such files into the  ",
" required vfile[iz][ix] structure.					",
" (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  		",
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	",
" programmers may expect.)						", 
"									",
" Also, sx must be monotonically increasing throughout the dataset, and ",
" and gx must be monotonically increasing within a shot. You may resort ",
" your data with \'susort\', accordingly.				",
"									",
" The scalco header field is honored so this field must be set correctly.",
" See selfdocs of \'susort\', \'suchw\'. Also:   sukeyword scalco	",
"									",
NULL};
#endif

/*
 * Credits: CWP, Baoniu Han, bhan@dix.mines.edu, April 19th, 1998
 *
 *	  Modified: Chris Stolk, 11 Dec 2005, - changed data input
 *		    to remove erroneous time delay.
 *	  Modified: CWP, John Stockwell 26 Sept 2006 - replaced Han's
 *	  "goto-loop" with  "do { }while loops".
 *	  Fixed it so that sx, gx, and scalco are honored.
 *
 *
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */

/**************** end self doc *******************************************/

/* the Ricker wavelet of the program: 2*N+1 samples, with the peak in the middle */
static float *ricker(float Freq,float dt,int *Npoint)
{
	int i;
	float Bpar,t,u,*Amp;
	int Np1,N;

	if(Freq==0.0)Freq=30.0;
	if(dt==0.0)dt=0.004;
	Bpar=sqrt(6.0)/(PI*Freq);
	N=ceil(1.35*Bpar/dt);
	Np1=N;
	*Npoint=2*N+1;

	Amp=alloc1float(*Npoint);

	Amp[Np1]=1.0;

	for(i=1;i<=N;i++) {
		t=dt*(float)i;
		u=2.0*sqrt(6.0)*t/Bpar;
		Amp[Np1+i]=Amp[Np1-i]=0.5*(2.0-u*u)*exp(-u*u/4.0);
	}

	return Amp;
}

/* the tridiagonal solve of the finite differences */
static void retris(complex *data,complex *a,complex *c,complex *b,complex endl,complex endr,int nx,complex *d)
{
	int ix;
	complex den;
	std::vector<complex> e(nx), f(nx);

	e[0] = -a[0]/endl;
	f[0] = d[0]/endl;
	for (ix=1; ix<nx-1; ++ix) {
		den = b[ix]+c[ix]*e[ix-1];
		e[ix] = -a[ix]/den;
		f[ix] = (d[ix]-f[ix-1]*c[ix])/den;
	}
	data[nx-1] = (d[nx-1]-f[nx-2]*c[nx-2])/(endr+c[nx-2]*e[nx-2]);
	for (ix=nx-2; ix>-1; --ix)
		data[ix] = data[ix+1]*e[ix]+f[ix];
}

/* The 15 degree finite-difference downward continuation of the program: step (coefa, coefb) of the dip, for the
 * frequencies fw + iw dw, and v[nx] the velocities */
static void fdstep_ffd(complex **cp, int nx, int nw, const float *v, float fw, float dw, float dz, float dx, float dt,
	const float *p, float vc, float coefa, float coefb)
{
	int iw,ix;
	float w,v1,vn,trick=0.1;
	complex cp2,cp3,cpnm1,cpnm2;
	complex a1,a2,b1,b2;
	complex endl,endr;
	std::vector<float> s1(nx), s2(nx);
	std::vector<complex> data(nx), d(nx), a(nx), b(nx), c(nx);

	v1 = v[0];
	vn = v[nx-1];

	for (iw=0,w=fw; iw<nw; iw++,w+=dw) {
		if (fabs(w)<=1.0e-10) w = 1.0e-10/dt;

		for (ix=0; ix<nx; ix++) {
			s1[ix] = (v[ix]*v[ix])*p[ix]*coefb/(dx*dx*w*w)+trick;
			s2[ix] = -(1-vc/v[ix])*v[ix]*dz*coefa/(w*dx*dx)*0.5;
		}

		for (ix=0; ix<nx; ix++)
			data[ix] = cp[iw][ix];

		cp2 = data[1];
		cp3 = data[2];
		cpnm1 = data[nx-2];
		cpnm2 = data[nx-3];

		a1 = cp2*std::conj(cp3)*2.0f;
		b1 = cp2*std::conj(cp2)+cp3*std::conj(cp3);
		if (b1.real()==0.0 && b1.imag()==0.0)
			a1 = std::exp(complex(0.0f,-w*dx*0.5/v1));
		else
			a1 = a1/b1;
		if (a1.imag()>0.0) a1 = std::exp(complex(0.0f,-w*dx*0.5/v1));

		a2 = cpnm1*std::conj(cpnm2)*2.0f;
		b2 = cpnm1*std::conj(cpnm1)+cpnm2*std::conj(cpnm2);
		if (b2.real()==0.0 && b2.imag()==0.0)
			a2 = std::exp(complex(0.0f,-w*dx*0.5/vn));
		else
			a2 = a2/b2;
		if (a2.imag()>0.0) a2 = std::exp(complex(0.0f,-w*dx*0.5/vn));

		for (ix=0; ix<nx; ix++) {
			a[ix] = complex(s1[ix],s2[ix]);
			b[ix] = complex(1.0-2.0*s1[ix],-2.0*s2[ix]);
		}

		for (ix=1; ix<nx-1; ix++)
			d[ix] = data[ix+1]*a[ix+1]+data[ix-1]*a[ix-1]+data[ix]*b[ix];
		d[0] = (b[0]+a[0]*a1)*data[0]+data[1]*a[1];
		d[nx-1] = (b[nx-1]+a[nx-1]*a2)*data[nx-1]+data[nx-2]*a[nx-2];

		for (ix=0; ix<nx; ix++) {
			data[ix] = complex(s1[ix],-s2[ix]);
			b[ix] = complex(1.0-2.0*s1[ix],2.0*s2[ix]);
		}

		endl = b[0]+data[0]*a1;
		endr = b[nx-1]+data[nx-1]*a2;

		for (ix=1; ix<nx-1; ix++) {
			a[ix] = data[ix+1];
			c[ix] = data[ix-1];
		}
		a[0] = data[1];
		c[nx-1] = data[nx-2];

		retris(data.data(),a.data(),c.data(),b.data(),endl,endr,nx,d.data());

		for (ix=0; ix<nx; ix++)
			cp[iw][ix] = data[ix];
	}
}

/* The Fourier finite-difference correction of SUMIGPREFFD: one 45 (or 65) degree step, for the velocity v and the
 * reference velocity vc */
static void fdmig(complex **cp, int nx, int nw, const float *v, float fw, float dw, float dz, float dx, float dt,
	float vc, int dip)
{
	float coefa,coefb;
	std::vector<float> p(nx);
	int ix;
	for (ix=0; ix<nx; ix++) {
		p[ix] = vc/v[ix];
		p[ix] = (p[ix]*p[ix]+p[ix]+1.0);
	}
	if (dip!=65) {
		coefa = 0.5; coefb = 0.25;
	} else {
		coefa = 0.4784689;
		coefb = 0.37607656;
	}
	fdstep_ffd(cp,nx,nw,v,fw,dw,dz,dx,dt,p.data(),vc,coefa,coefb);
}

/* Library version of SUMIGPREFFD
 *
 * The main program (the parameters, the traces and their headers, the velocity file, the loop over shots) is left to the
 * caller. The Fourier finite-difference depth migration of a shot.
 *
 *	su_migpreffd_shot:	migrates one shot gather, p[nxo][nt] (the traces of the gather at the positions of their receivers, zero
 *		where there is no trace), whose shot is at the position ixshot, and adds its image to cresult[nxo][nz]
 *		(the sum of the shots). ix2 and ix3 are the first and last positions of the migrated aperture. The
 *		velocities are vp[nz][nxo]. Frequencies f1, f2, f3 and f4 are those of the Hamming window.
 *
 * Differences from the program, which are its bugs: the Hamming window of the program was shifted by the index of the first
 * frequency of the band (its taper was reversed); the number of frequencies and the aperture are checked.
 * Returns 0, or -3 for parameters that are not allowed.
 */

int su_migpreffd_shot(int nt, int nxo, int nz, float dt, float dx, float dz, int dip, float fmax, float f1, float f2, float f3,
	float f4, int ix2, int ix3, int ixshot, const float *pin, const float *vpin, float *cresult_out)
{
	int iz,iw,ix,it,ik;
	int ntfft,nxfft=0,nw,truenw,nk=0,ntw,nx;
	int nf1,nf2,nf3,nf4;
	float dw,dk=0.0,fw,fk=0.0,w,k,v1=0.0,vmin=0.0;
	double kz1,kz2,phase1;
	complex cshift1,cshift2;
	float **p,**v,*wl,*wtmp;
	complex *wlsp,**cp,**cp1,**cq,**cq1=NULL;


	if (nt<2 || nxo<4 || nz<1 || dt<=0.0 || dx<=0.0 || dz<=0.0 || fmax<0.0) return -3;
	if (ix2<0 || ix3>nxo-1 || ix3-ix2+1<4 || ixshot<ix2 || ixshot>ix3) return -3;
	for (iz=0; iz<nz; iz++)
		for (ix=0; ix<nxo; ix++)
			if (vpin[(size_t)iz*nxo+ix] <= 0.0) return -3;


	/* determine frequency sampling interval */
	ntfft = npfar(nt);
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);

	/* compute the index of the frequency to be migrated */
	nf1 = 2.0*PI*f1/dw+0.5;
	nf2 = 2.0*PI*f2/dw+0.5;
	nf3 = 2.0*PI*f3/dw+0.5;
	nf4 = 2.0*PI*f4/dw+0.5;
	if (nf1<0 || nf2<=nf1 || nf3<nf2 || nf4<=nf3 || nf4>nw-1) return -3;

	/* the number of frequencies to migrated */
	truenw = nf4-nf1+1;
	fw = 0.0+nf1*dw;

	/* the Ricker wavelet, centered, and its Fourier transform */
	wl = alloc1float(ntfft);
	wlsp = alloc1complex(nw);
	wtmp = ricker(fmax,dt,&ntw);
	memset((void *) wl, 0, ntfft*FSIZE);
	for (it=0; it<ntw; it++)
		wl[(it-ntw/2+ntfft) % ntfft] = wtmp[it];
	free1float(wtmp);
	pfarc(-1,ntfft,wl,wlsp);

	/* the traces of the shot, in frequency */
	p = alloc2float(ntfft,nxo);
	cq = alloc2complex(nw,nxo);
	memset((void *) p[0], 0, ntfft*nxo*FSIZE);
	for (ix=0; ix<nxo; ix++)
		for (it=0; it<nt; it++) p[ix][it] = pin[(size_t)ix*nt+it];
	pfa2rc(1,1,ntfft,nxo,p[0],cq[0]);

	/* the aperture */
	nx = ix3-ix2+1;
	nw = truenw;

	/* determine wavenumber sampling (for complex to complex FFT) */
	nxfft = npfa(nx);
	nk = nxfft;
	dk = 2.0*PI/(nxfft*dx);
	fk = -PI/dx;


	v = alloc2float(nx,nz);
	for (iz=0; iz<nz; ++iz)
		for (ix=0; ix<nx; ++ix)
			v[iz][ix] = vpin[(size_t)iz*nxo+ix+ix2];

	cp = alloc2complex(nx,nw);
	cp1 = alloc2complex(nx,nw);

	/* transpose the frequency domain data from data[ix][iw] to data[iw][ix] and apply a Hamming window at the same time */
	for (ix=0; ix<nx; ++ix) {
		for (iw=0; iw<nw; iw++) {
			float tmpp;
			if (iw>=0 && iw<=(nf2-nf1)) {
				tmpp = 0.54+0.46*cos(PI/(nf2-nf1)*iw-PI);
				cp[iw][ix] = cq[ix+ix2][iw+nf1]*tmpp;
			} else if (iw>=(nf3-nf1) && iw<=(nf4-nf1)) {
				tmpp = 0.54+0.46*cos(PI/(nf4-nf3)*(iw-(nf3-nf1)));
				cp[iw][ix] = cq[ix+ix2][iw+nf1]*tmpp;
			} else {
				cp[iw][ix] = cq[ix+ix2][iw+nf1];
			}
			cp1[iw][ix] = complex(0.0f,0.0f);
		}
	}

	/* the source: the wavelet at the shot */
	for (iw=0; iw<nw; iw++)
		cp1[iw][ixshot-ix2] = wlsp[iw+nf1];

	free2float(p);
	free2complex(cq);
	free1float(wl);
	free1complex(wlsp);

	cq = alloc2complex(nxfft,nw);
	cq1 = alloc2complex(nxfft,nw);

	/* loops over depth */
	for (iz=0; iz<nz; ++iz) {

		/* the imaging condition */
		for (ix=0; ix<nx; ++ix) {
			for (iw=0,w=fw; iw<nw; w+=dw,iw++) {
				complex tmp;
				float ratio=10.0;
				if (abs(ix+ix2-ixshot)*dx<ratio*iz*dz)
					tmp = cp[iw][ix]*cp1[iw][ix];
				else
					tmp = complex(0.0f,0.0f);
				cresult_out[(size_t)(ix+ix2)*nz+iz] += tmp.real()/ntfft;
			}
		}


		/* the reference velocity */
		/* get the minimum velocity */
		vmin = v[iz][0];
		for (ix=0; ix<nx; ix++)
			if (v[iz][ix]<vmin) vmin = v[iz][ix];

		/* shift the wavefield */
		for (ik=0; ik<nx; ++ik) {
			for (iw=0; iw<nw; ++iw) {
				cq[iw][ik] = ik%2 ? -cp[iw][ik] : cp[iw][ik];
				cq1[iw][ik] = ik%2 ? -cp1[iw][ik] : cp1[iw][ik];
			}
		}

		/* zero out parts of the cq[][] and cq1[][] arrays */
		for (ik=nx; ik<nk; ++ik) {
			for (iw=0; iw<nw; ++iw) {
				cq[iw][ik] = complex(0.0f,0.0f);
				cq1[iw][ik] = complex(0.0f,0.0f);
			}
		}

		/* FFT to W-K domain */
		pfa2cc(-1,1,nk,nw,cq[0]);
		pfa2cc(-1,1,nk,nw,cq1[0]);

		v1 = vmin;

		/* apply phase shift */
		for (ik=0,k=fk; ik<nk; ++ik,k+=dk) {
			for (iw=0,w=fw; iw<nw; ++iw,w+=dw) {
				if (w==0.0) w = 1.0e-10/dt;

				kz1 = 1.0-pow(v1*k/w,2.0);
				if (kz1>0.15) {
					phase1 = -w*sqrt(kz1)*dz/v1;
					cshift1 = complex(cos(phase1), sin(phase1));
					cq[iw][ik] = cq[iw][ik]*cshift1;
					cq1[iw][ik] = cq1[iw][ik]*cshift1;
				} else {
					cq[iw][ik] = cq1[iw][ik] = complex(0.0f,0.0f);
				}
			}
		}

		/* fourier transform */
		pfa2cc(1,1,nk,nw,cq[0]);
		pfa2cc(1,1,nk,nw,cq1[0]);

		for (ix=0; ix<nx; ++ix) {
			for (iw=0,w=fw; iw<nw; w+=dw,++iw) {
				cq[iw][ix] = cq[iw][ix]*(1.0f/nxfft);
				cq[iw][ix] = ix%2 ? -cq[iw][ix] : cq[iw][ix];
				kz2 = (1.0/v1-1.0/v[iz][ix])*w*dz;
				cshift2 = complex(cos(kz2),sin(kz2));
				cp[iw][ix] = cq[iw][ix]*cshift2;

				cq1[iw][ix] = cq1[iw][ix]*(1.0f/nxfft);
				cq1[iw][ix] = ix%2 ? -cq1[iw][ix] : cq1[iw][ix];
				cp1[iw][ix] = cq1[iw][ix]*cshift2;
			}
		}

		/* apply fdmig algorithm */
		fdmig(cp,nx,nw,v[iz],fw,dw,dz,dx,dt,v1,dip);
		fdmig(cp1,nx,nw,v[iz],fw,dw,dz,dx,dt,v1,dip);


	}

	free2complex(cp);
	free2complex(cp1);
	free2float(v);

	free2complex(cq);
	free2complex(cq1);

	return 0;
}
