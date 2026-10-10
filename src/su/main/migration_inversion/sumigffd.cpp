/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPSMIGFFD: $Revision: 1.9 $ ; $Date: 2015/08/07 22:19:43 $	*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUMIGFFD - Fourier finite difference depth migration for		",
"	    zero-offset data. This method is a hybrid migration which	",
"	    combines the advantages of phase shift and finite difference", 
"	    migrations.							",
"									",
" sumigffd <infile >outfile vfile= [optional parameters]		",
"									",
" Required Parameters:						  	",
" nz=		   number of depth samples			 	", 
" dz=		   depth sampling interval			 	",
" vfile=		name of file containing velocities	      	",
"									",
" Optional Parameters:						  	",
" dt=from header(dt) or .004    time sampling interval		  	",
" dx=from header(d2) or 1.0     midpoint sampling interval	  	",
" ft=0.0			first time sample			",
" fz=0.0			first depth sample		      	",
"									",
" tmpdir=	if non-empty, use the value as a directory path		",
"		prefix for storing temporary files; else if the		",
"		the CWP_TMPDIR environment variable is set use		",
"		its value for the path; else use tmpfile()		",
"									", 
" The input velocity file \'vfile\' consists of C-style binary floats.  ",  
" The structure of this file is vfile[iz][ix]. Note that this means that",
" the x-direction is the fastest direction instead of z-direction! Such a",
" structure is more convenient for the downward continuation type	",
" migration algorithm than using z as fastest dimension as in other SU  ", 
" programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  ",
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	",
" programmers may expect.)						", 
"									",
" Because most of the tools in the SU package (such as  unif2, unisam2, ", 
" and makevel) produce output with the structure vfile[ix][iz], you will",
" need to transpose the velocity files created by these programs. You may",
" use the SU program \'transp\' in SU to transpose such files into the  ",
" required vfile[iz][ix] structure.					",
"									",
"									",
NULL};
#endif


/*
 * Credits: CWP Baoniu Han, July 21th, 1997
 *
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */
/**************** end self doc *******************************************/


/* Library version of SUMIGFFD
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migffd:	migrates the zero-offset section p[nx][nt] (the traces in order) to cresult[nx][nz], with the
 *			velocities v[nz][nx]
 *
 * The velocities are not changed. Returns 0, or -3 for parameters that are not allowed.
 */

static void retris(complex *data,complex *a,complex *c,complex *b,
		complex endl,complex endr, int nx, complex *d)
{
	int ix;
	complex *e,den;
	complex *f;

	e=alloc1complex(nx);
	f=alloc1complex(nx);
	e[0]=-a[0]/endl;
	f[0]=d[0]/endl;

	for(ix=1;ix<nx-1;++ix){
		den=b[ix]+c[ix]*e[ix-1];
		e[ix]=-a[ix]/den;
		f[ix]=(d[ix]-f[ix-1]*c[ix])/den;
	}

	data[nx-1]=(d[nx-1]-f[nx-2]*c[nx-1])/(endr+c[nx-1]*e[nx-2]);

	for(ix=nx-2;ix>-1;--ix)
		data[ix]=data[ix+1]*e[ix]+f[ix];

	free1complex(e);
	free1complex(f);
}

static void fdmig( complex **cp, int nx, int nw, float *v,float fw,float
		dw,float dz,float dx,float dt,float vc)
{
	int iw,ix;
	float *p,*s1,*s2,w,a0=2.0,v1,vn,trick=0.14;
	complex endl,endr;
	complex cp2,cp3,cpnm1,cpnm2;
	complex a1,a2,b1,b2;
	complex *data,*d,*a,*b,*c;
	const complex czero(0.0f,0.0f);

	p=alloc1float(nx);
	s1=alloc1float(nx);
	s2=alloc1float(nx);

	data=alloc1complex(nx);
	d=alloc1complex(nx);
	a=alloc1complex(nx);
	b=alloc1complex(nx);
	c=alloc1complex(nx);

	for(ix=0;ix<nx;ix++){
		p[ix]=vc*2.0/v[ix];
		p[ix]=0.5*(p[ix]*p[ix]+p[ix]+1.0);
	}

	v1=v[0];vn=v[nx-1];

	for(iw=1,w=fw+dw;iw<nw;iw++,w+=dw){
		if(w==0)w=1.0e-10/dt;

		for(ix=0;ix<nx;ix++){
			s1[ix]=p[ix]*v[ix]*v[ix]/(4.0*a0*dx*dx*w*w)+trick;
			s2[ix]=(1-vc*2.0/v[ix])*v[ix]*dz*0.5/(2.0*w*a0*dx*dx);
		}

		for(ix=0;ix<nx;ix++){
			data[ix]=cp[ix][iw];
		}

		cp2=data[1];
		cp3=data[2];
		cpnm1=data[nx-2];
		cpnm2=data[nx-3];
		a1=2.0f*(cp2*std::conj(cp3));
		b1=cp2*std::conj(cp2)+cp3*std::conj(cp3);

		if(b1==czero)
			a1=std::exp(complex(0.0f,-w*dx*0.5/v1));
		else a1=a1/b1;

		if(a1.imag()>0.0)a1=std::exp(complex(0.0f,-w*dx*0.5/v1));

		a2=2.0f*(cpnm1*std::conj(cpnm2));
		b2=cpnm1*std::conj(cpnm1)+cpnm2*std::conj(cpnm2);

		if(b2==czero)
			a2=std::exp(complex(0.0f,-w*dx*0.5/vn));
		else a2=a2/b2;
		if(a2.imag()>0.0)a2=std::exp(complex(0.0f,-w*dx*0.5/vn));

		for(ix=0;ix<nx;ix++){
			a[ix]=complex(s1[ix],-s2[ix]);
			b[ix]=complex(1.0-2.0*s1[ix],2.0*s2[ix]);
		}

		d[0]=(b[0]+a[0]*a1)*data[0]+data[1]*a[1];
		d[nx-1]=(b[nx-1]+a[nx-1]*a2)*data[nx-1]+data[nx-2]*a[nx-2];

		for(ix=1;ix<nx-1;ix++){
			d[ix]=data[ix+1]*a[ix+1]+data[ix-1]*a[ix-1]+data[ix]*b[ix];
		}

		for(ix=0;ix<nx;ix++){
			data[ix]=complex(s1[ix],s2[ix]);
			b[ix]=complex(1.0-2.0*s1[ix],-2.0*s2[ix]);
		}

		endl=b[0]+data[0]*a1;
		endr=b[nx-1]+data[nx-1]*a2;

		for(ix=1;ix<nx-1;ix++) {
			a[ix]=data[ix+1];
			c[ix]=data[ix-1];
		}
		a[0]=data[1];
		c[nx-1]=data[nx-2];

		retris(data,a,c,b,endl,endr,nx,d);

		for(ix=0;ix<nx;ix++){
			cp[ix][iw]=data[ix];
		}
	}

	free1complex(data);
	free1complex(d);
	free1float(p);
	free1complex(b);
	free1complex(c);
	free1complex(a);
	free1float(s1);
	free1float(s2);
}

int su_migffd(int nx, int nt, int nz, float dz, float dt, float dx, const float *pin, const float *vin, float *cresult_out)
{
	int ik,iz,iw,ix,it,nxfft,ntfft,nk,nw;
	float dk,dw,fk,fw,k,w,vmin,v1;
	double kz1,kz2,phase1;
	float **p,**cresult,**v;
	complex cshift1,cshift2,**cp,**cq;

	if (nx<4 || nt<2 || nz<1 || dz<=0.0 || dt<=0.0 || dx<=0.0) return -3;

	/* determine frequency sampling (for real to complex FFT) */
	ntfft = npfar(nt);
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);
	fw = 0.0;

	/* determine wavenumber sampling (for complex to complex FFT) */
	nxfft = npfa(nx);
	nk = nxfft;
	dk = 2.0*PI/(nxfft*dx);
	fk = -PI/dx;

	/* allocate space */
	p = alloc2float(ntfft,nx);
	cp = alloc2complex(nw,nx);
	cq = alloc2complex(nw,nk);
	cresult = alloc2float(nz,nx);
	v=alloc2float(nx,nz);

	for(ix=0;ix<nx;ix++) {
		for(it=0;it<nt;it++) p[ix][it] = pin[(size_t)ix*nt+it];
		for (it=nt; it<ntfft; it++) p[ix][it] = 0.0;
	}
	for(iz=0;iz<nz;iz++)
		for(ix=0;ix<nx;ix++) {
			v[iz][ix] = vin[(size_t)iz*nx+ix];
			if (v[iz][ix] <= 0.0) { free2float(p); free2complex(cp); free2complex(cq); free2float(cresult); free2float(v); return -3; }
		}

	pfa2rc(1,1,ntfft,nx,p[0],cp[0]);

	/* loop over depth*/
	for(iz=0;iz<nz;++iz){

		for(ix=0;ix<nx;ix++){
			cresult[ix][iz] =0.0;
			for(iw=1;iw<nw;iw++)
				cresult[ix][iz]+=cp[ix][iw].real()/ntfft;
		}

		vmin=v[iz][0];
		for(ix=0;ix<nx;ix++){
			if(v[iz][ix]<=vmin) vmin=v[iz][ix];
		}

		for (ik=0; ik<nx; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = ik%2 ? -cp[ik][iw] : cp[ik][iw];

		for (ik=nx; ik<nk; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = complex(0.0f,0.0f);

		/* FFT to W-K domain */
		pfa2cc(-1,2,nw,nk,cq[0]);

		v1=vmin*0.5;

		for(ik=0,k=fk;ik<nk;++ik,k+=dk)
			for(iw=1,w=fw+dw;iw<nw;++iw,w+=dw){
				if(w==0.0)w=1.0e-10/dt;
				kz1=1.0-pow(v1*k/w,2.0);
				if(kz1>0){
					phase1 = -w*sqrt(kz1)*dz/v1;
					cshift1 = complex(cos(phase1), sin(phase1));
					cq[ik][iw] = cq[ik][iw]*cshift1;
				} else {
					cq[ik][iw] = complex(0.0f,0.0f);
				}
			}

		pfa2cc(1,2,nw,nk,cq[0]);
		for(ix=0;ix<nx;++ix)
			for(iw=1,w=fw+dw;iw<nw;w+=dw,++iw){
				cq[ix][iw] = cq[ix][iw]*(1.0f/nxfft);
				cp[ix][iw] = ix%2 ? -cq[ix][iw] : cq[ix][iw];
			}

		for(ix=0;ix<nx;++ix)
			for(iw=1,w=fw+dw;iw<nw;w+=dw,++iw){
				kz2=(1.0/v1-2.0/v[iz][ix])*w*dz;
				cshift2=complex(cos(kz2),sin(kz2));
				cp[ix][iw]=cp[ix][iw]*cshift2;
			}

		fdmig( cp, nx, nw,v[iz],fw,dw,dz,dx,dt,v1);
	}

	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++) cresult_out[(size_t)ix*nz+iz] = cresult[ix][iz];

	free2float(p);
	free2complex(cp);
	free2complex(cq);
	free2float(cresult);
	free2float(v);
	return 0;
}
