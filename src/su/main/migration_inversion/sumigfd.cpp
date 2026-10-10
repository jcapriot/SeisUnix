/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGFD: $Revision: 1.13 $ ; $Date: 2015/08/07 22:19:43 $	*/

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"									",
" SUMIGFD - 45-90 degree Finite difference depth migration for		",
"           zero-offset data.						",
"									",
"   sumigfd <infile >outfile vfile= [optional parameters]		",
"									",
" Required Parameters:							",
" nz=		number of depth samples					",
" dz=		depth sampling interval					",
" vfile=	name of file containing velocities			",
" 		(see Notes below concerning format of this file)	",
"									",
" Optional Parameters:							",
" dt=from header(dt) or .004    time sampling interval			",
" dx=from header(d2) or 1.0	midpoint sampling interval		",
" dip=45,65,79,80,87,89,90  	Maximum angle of dip reflector		",
"									",
" tmpdir=	if non-empty, use the value as a directory path		",
"		prefix for storing temporary files; else if the		",
"		the CWP_TMPDIR environment variable is set use		",
"		its value for the path; else use tmpfile()		",
"									", 
" Notes:								", 
" The computation cost by dip angle is 45=65=79<80<87<89<90		",
"									", 
" The input velocity file \'vfile\' consists of C-style binary floats.	", 
" The structure of this file is vfile[iz][ix]. Note that this means that",
" the x-direction is the fastest direction instead of z-direction! Such a",
" structure is more convenient for the downward continuation type	",
" migration algorithm than using z as fastest dimension as in other SU	",
" programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  ",
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	",
" programmers may expect.)						", 
"									", 
" Because most of the tools in the SU package (such as  unif2, unisam2,	",
" and makevel) produce output with the structure vfile[ix][iz], you will",
" need to transpose the velocity files created by these programs. You may",
" use the SU program \'transp\' in SU to transpose such files into the	",
" required vfile[iz][ix] structure.					",
"									",
NULL};
#endif


/* 
 * Credits: CWP Baoniu Han, April 20th, 1998
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */
/**************** end self doc *******************************************/


/* Library version of SUMIGFD
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migfd:	migrates the zero-offset section p[nx][nt] (the traces in order) to cresult[nx][nz], with the
 *			velocities v[nz][nx] (in the section's x direction fast)
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

	data[nx-1]=(d[nx-1]-f[nx-2]*c[nx-2])/(endr+c[nx-2]*e[nx-2]);

	for(ix=nx-2;ix>-1;--ix)
		data[ix]=data[ix+1]*e[ix]+f[ix];

	free1complex(e);
	free1complex(f);
}

static void fdmig( complex **cp, int nx, int nw, float *v,float fw,float
		dw,float dz,float dx,float dt,int dip)
{
	int iw,ix,step=1;
	float *s1,*s2,w,coefa[5],coefb[5],v1,vn,trick=0.1;
	complex cp2,cp3,cpnm1,cpnm2;
	complex a1,a2,b1,b2;
	complex endl,endr;
	complex *data,*d,*a,*b,*c;
	const complex czero(0.0f,0.0f);

	s1=alloc1float(nx);
	s2=alloc1float(nx);

	data=alloc1complex(nx);
	d=alloc1complex(nx);
	a=alloc1complex(nx);
	b=alloc1complex(nx);
	c=alloc1complex(nx);

	if(dip!=45&&dip!=65&&dip!=79&&dip!=80&&dip!=87&&dip!=89&&dip!=90)
	dip=79;

	if(dip==45){
	coefa[0]=0.5;coefb[0]=0.25;
	step=1;
	}

	if(dip==65){
	coefa[0]=0.478242060;coefb[0]=0.376369527;
	step=1;
	}

	if(dip==79){
	coefa[0]=coefb[0]=0.4575;
	step=1;
	}

	if(dip==80){
	coefa[1]=0.040315157;coefb[1]=0.873981642;
	coefa[0]=0.457289566;coefb[0]=0.222691983;
	step=2;
	}

	if(dip==87){
	coefa[2]=0.00421042;coefb[2]=0.972926132;
	coefa[1]=0.081312882;coefb[1]=0.744418059;
	coefa[0]=0.414236605;coefb[0]=0.150843924;
	step=3;
	}

	if(dip==89){
	coefa[3]=0.000523275;coefb[3]=0.994065088;
	coefa[2]=0.014853510;coefb[2]=0.919432661;
	coefa[1]=0.117592008;coefb[1]=0.614520676;
	coefa[0]=0.367013245;coefb[0]=0.105756624;
	step=4;
	}

	if(dip==90){
	coefa[4]=0.000153427;coefb[4]=0.997370236;
	coefa[3]=0.004172967;coefb[3]=0.964827992;
	coefa[2]=0.033860918;coefb[2]=0.824918565;
	coefa[1]=0.143798076;coefb[1]=0.483340757;
	coefa[0]=0.318013812;coefb[0]=0.073588213;
	step=5;
	}

	v1=v[0];vn=v[nx-1];

	while (step>0) {
	step--;

	for(iw=0,w=fw;iw<nw;iw++,w+=dw){

		if(fabs(w)<=1.0e-10)w=1.0e-10/dt;

		for(ix=0;ix<nx;ix++){
			s1[ix]=(v[ix]*v[ix])*coefb[step]/(dx*dx*w*w)+trick;
			s2[ix]=-v[ix]*dz*coefa[step]/(w*dx*dx)*0.5;
		}

		for(ix=0;ix<nx;ix++){
			data[ix]=cp[ix][iw];
		}

		cp2=data[0];
		cp3=data[1];
		cpnm1=data[nx-1];
		cpnm2=data[nx-2];
		a1=cp2*std::conj(cp3);
		b1=cp3*std::conj(cp3);
		if(b1==czero)
			a1=std::exp(complex(0.0f,-w*dx*0.5/v1));
		else
			a1=a1/b1;

		if(a1.imag()>0.0)a1=std::exp(complex(0.0f,-w*dx*0.5/v1));

		a2=cpnm1*std::conj(cpnm2);
		b2=cpnm2*std::conj(cpnm2);

		if(b2==czero)
			a2=std::exp(complex(0.0f,-w*dx*0.5/vn));
		else
			a2=a2/b2;

		if(a2.imag()>0.0)a2=std::exp(complex(0.0f,-w*dx*0.5/vn));

		for(ix=0;ix<nx;ix++){
			a[ix]=complex(s1[ix],s2[ix]);
			b[ix]=complex(1.0-2.0*s1[ix],-2.0*s2[ix]);
		}

		for(ix=1;ix<nx-1;ix++){
			d[ix]=data[ix+1]*a[ix+1]+data[ix-1]*a[ix-1]+data[ix]*b[ix];
		}

		d[0]=(b[0]+a[0]*a1)*data[0]+data[1]*a[1];
		d[nx-1]=(b[nx-1]+a[nx-1]*a2)*data[nx-1]+data[nx-2]*a[nx-2];

		for(ix=0;ix<nx;ix++){
			data[ix]=complex(s1[ix],-s2[ix]);
			b[ix]=complex(1.0-2.0*s1[ix],2.0*s2[ix]);
		}

		endl=b[0]+data[0]*a1;
		endr=b[nx-1]+data[nx-1]*a2;

		for(ix=1;ix<nx-1;ix++){
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
	}

	free1complex(data);
	free1complex(d);
	free1complex(b);
	free1complex(c);
	free1complex(a);
	free1float(s1);
	free1float(s2);
}

int su_migfd(int nx, int nt, int nz, float dz, float dt, float dx, int dip, const float *pin, const float *vin, float *cresult_out)
{
	int iz,iw,ix,it,ntfft,nw;
	float dw,fw=0.0,w;
	double kz2;
	float **p,**cresult,**v;
	complex cshift2,**cp;

	if (nx<3 || nt<2 || nz<1 || dz<=0.0 || dt<=0.0 || dx<=0.0) return -3;

	/* determine frequency sampling (for real to complex FFT) */
	ntfft = npfar(nt);
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);

	/* allocate space */
	p = alloc2float(ntfft,nx);
	cp = alloc2complex(nw,nx);
	cresult = alloc2float(nz,nx);
	v=alloc2float(nx,nz);

	/* Zero all arrays */
	memset((void *) p[0], 0, FSIZE*(size_t)ntfft*nx);
	memset((void *) cp[0], 0, sizeof(complex)*(size_t)nw*nx);
	memset((void *) cresult[0], 0, FSIZE*(size_t)nz*nx);

	/* load the traces and the velocities */
	for(ix=0;ix<nx;ix++)
		for(it=0;it<nt;it++) p[ix][it] = pin[(size_t)ix*nt+it];
	for(iz=0;iz<nz;iz++)
		for(ix=0;ix<nx;ix++) {
			v[iz][ix] = vin[(size_t)iz*nx+ix];
			if (v[iz][ix] <= 0.0) { free2float(p); free2complex(cp); free2float(cresult); free2float(v); return -3; }
		}

	/* Fourier transform */
	pfa2rc(1,1,ntfft,nx,p[0],cp[0]);

	/* loop over depth*/
	for(iz=0;iz<nz;++iz){
		for(ix=0;ix<nx;ix++){
			cresult[ix][iz] =0.0;
			v[iz][ix]=v[iz][ix]/2.0;
			for(iw=1;iw<nw;iw++)
				cresult[ix][iz]+=cp[ix][iw].real()/ntfft;
		}

		for(ix=0;ix<nx;++ix)
		for(iw=1,w=fw+dw;iw<nw;w+=dw,++iw) {
			kz2=-1.0/v[iz][ix]*w*dz;
			cshift2=complex(cos(kz2),sin(kz2));
			cp[ix][iw]=cp[ix][iw]*cshift2;
		}

		fdmig( cp, nx, nw,v[iz],fw,dw,dz,dx,dt,dip);
	}

	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++) cresult_out[(size_t)ix*nz+iz] = cresult[ix][iz];

	free2float(p);
	free2complex(cp);
	free2float(cresult);
	free2float(v);
	return 0;
}
