/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGPSPI: $Revision: 1.12 $ ; $Date: 2015/08/07 22:19:43 $       */

#include "su.h"
#include "segy.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                       ",
" SUMIGPSPI - Gazdag's phase-shift plus interpolation depth migration   ",
"            for zero-offset data, which can handle the lateral         ",
"            velocity variation.                                        ",
"                                                                       ",
" sumigpspi <infile >outfile vfile= [optional parameters]               ",
"                                                                       ", 
" Required Parameters:							",
" nz=		number of depth samples					",
" dz=		depth sampling interval					",
" vfile=	name of file containing velocities			",
"		(Please see Notes below concerning the format of vfile)	",
"									",
" Optional Parameters:                                                  ",
" dt=from header(dt) or .004    time sampling interval                  ",
" dx=from header(d2) or 1.0     midpoint sampling interval              ",
"                                                                       ",
" tmpdir=        if non-empty, use the value as a directory path        ",
"                prefix for storing temporary files; else if the        ",
"                the CWP_TMPDIR environment variable is set use         ",
"                its value for the path; else use tmpfile()             ",
"                                                                       ",
" Notes:								",
" The input velocity file 'vfile' consists of C-style binary floats.	",
" The structure of this file is vfile[iz][ix]. Note that this means that",
" the x-direction is the fastest direction instead of z-direction! Such a",
" structure is more convenient for the downward continuation type	",
" migration algorithm than using z as fastest dimension as in other SU	",
" programs. (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]	",
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	",
" programmers may expect.)						",
"									",
" Because most of the tools in the SU package (such as  unif2, unisam2,	",
" and makevel) produce output with the structure vfile[ix][iz], you will",
" need to transpose the velocity files created by these programs. You may",
" use the SU program 'transp' in SU to transpose such files into the	",
" required vfile[iz][ix] structure.					",
"									",
"									",
NULL};
#endif


/*
 * Credits: CWP, Baoniu Han, April 20th, 1998
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */
/**************** end self doc *******************************************/


/* Library version of SUMIGPSPI
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migpspi:	migrates the zero-offset section p[nx][nt] (the traces in order) to cresult[nx][nz], with the
 *			velocities v[nz][nx]
 *
 * The velocities are not changed. Returns 0, or -3 for parameters that are not allowed, -4 if there is not memory for the
 * reference velocities (the program stopped).
 */

int su_migpspi(int nx, int nt, int nz, float dz, float dt, float dx, const float *pin, const float *vin, float *cresult_out)
{
	const int L=10;
	int Bz=0;
	float c[41], *V=NULL, P[40], Y[41];
	float Sz=0;
	int ik,iz,iw,ix,it;
	int nxfft,ntfft;
	int nk,nw;
	float dk,dw;
	float fk,fw;
	float k,w;
	float **p,**cresult;
	float vmin=0.0,vmax=0.0,v1,v2,lvmin,lvmax;
	double kz1;
	float **v;
	double phase1;
	complex cshift1;
	complex **cp,**cq,**cq1,***cq2;
	double a,a1,a2,theta,theta1,theta2;

	if (nx<4 || nt<2 || nz<1 || dz<=0.0 || dt<=0.0 || dx<=0.0) return -3;

	/* determine frequency sampling (for real to complex FFT) */
	ntfft = npfar(nt);
	nw = ntfft/2+1;
	dw = 2.0*PI/(ntfft*dt);
	fw = 0;

	/* determine wavenumber sampling (for complex to complex FFT) */
	nxfft = npfa(nx);
	nk = nxfft;
	dk = 2.0*PI/(nxfft*dx);
	fk = -PI/dx;

	/* allocate space */
	p = alloc2float(ntfft,nx);
	cp = alloc2complex(nw,nx);
	cq = alloc2complex(nw,nk);
	cq1 = alloc2complex(nw,nk);
	cresult = alloc2float(nz,nx);
	v=alloc2float(nx,nz);

	for(ix=0;ix<nx;ix++) {
		for(it=0;it<nt;it++) p[ix][it] = pin[(size_t)ix*nt+it];
		for(it=nt;it<ntfft;it++) p[ix][it] = 0.0;
	}
	for(iz=0;iz<nz;iz++)
		for(ix=0;ix<nx;ix++) {
			v[iz][ix] = vin[(size_t)iz*nx+ix];
			if (v[iz][ix] <= 0.0) {
				free2float(p); free2complex(cp); free2complex(cq); free2complex(cq1); free2float(cresult); free2float(v);
				return -3;
			}
		}

	vmax=v[0][0];vmin=v[0][0];
	for(iz=0;iz<nz;++iz)
		for(ix=0;ix<nx;ix++) {
			if(v[iz][ix]>=vmax) vmax=v[iz][ix];
			if(v[iz][ix]<=vmin) vmin=v[iz][ix];
		}

	pfa2rc(1,1,ntfft,nx,p[0],cp[0]);

	/*loops over depth*/
	for(iz=0;iz<nz;++iz){

		for(ix=0;ix<nx;ix++){
			cresult[ix][iz] =0.0;
			for(iw=0;iw<nw;iw++)
				cresult[ix][iz]+=cp[ix][iw].real()/ntfft;
		}

		lvmax=v[iz][0];
		lvmin=v[iz][0];

		for(ix=0;ix<nx;ix++){
			if(v[iz][ix]>=lvmax) lvmax=v[iz][ix];
			if(v[iz][ix]<=lvmin) lvmin=v[iz][ix];
		}

		for (ik=0; ik<nx; ++ik)
			for (iw=0,w=fw; iw<nw;w+=dw, ++iw){
				cp[ik][iw]=cp[ik][iw]*std::exp(complex(0.0f,-w*dz*2.0/v[iz][ik]));
				cq[ik][iw] = ik%2 ? -cp[ik][iw] : cp[ik][iw];
			}

		for (ik=nx; ik<nk; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = complex(0.0f,0.0f);

		/* FFT to W-K domain */
		pfa2cc(-1,2,nw,nk,cq[0]);

		/* The second time phase shift */
		v1=lvmin*0.5;
		v2=lvmax*0.5;

		if((v2-v1)/v1<0.01){

			for(ik=0,k=fk;ik<nk;++ik,k+=dk)
				for(iw=0,w=fw;iw<nw;++iw,w+=dw){

					if(w==0.0)w=1.0e-10/dt;
					kz1=1.0-pow(v1*k/w,2.0);

					if(kz1>0){
						phase1 = -w*sqrt(kz1)*dz/v1+w*dz/v1;
						cshift1 = complex(cos(phase1), sin(phase1));
						cq1[ik][iw] = cq[ik][iw]*cshift1;
					} else {
						phase1 = -w*sqrt(-kz1)*dz/v1;
						cshift1=std::exp(complex(phase1,w*dz/v1));
						cq1[ik][iw] = cq[ik][iw]*cshift1;
					}
				}

			pfa2cc(1,2,nw,nk,cq1[0]);

			for(ix=0;ix<nx;++ix)
				for(iw=0;iw<nw;++iw){
					cq1[ix][iw] = cq1[ix][iw]*(1.0f/nxfft);
					cp[ix][iw] = ix%2 ? -cq1[ix][iw] : cq1[ix][iw];
				}
		}
		else{

			for(ik=0;ik<=L;ik++)
				c[ik]=vmin+ik*1.0*(vmax-vmin)/(L*1.0);

			for(ik=0;ik<L;ik++)
				P[ik]=0.0;

			for(ix=0;ix<nx;ix++){
				for(ik=0;ik<L;ik++){
					if(((v[iz][ix]>=c[ik])&&(v[iz][ix]<c[ik+1]))||((ik==L-1)&&(v[iz][ix]==vmax))){
						P[ik]+=1.0/nx; break;
					}
				}
			}

			Sz=0.0;
			for(ik=0;ik<L;ik++)
				if(P[ik]!=0.00) Sz=Sz-P[ik]*log(P[ik]);

			Bz=exp(Sz)+0.5;
			Y[0]=0.0; Y[L]=1.0;

			for(ik=1;ik<L;ik++){
				Y[ik]=0.0;
				for(ix=0;ix<ik;ix++)
					Y[ik]=Y[ik]+P[ix];
			}

			V=alloc1float(Bz+1);
			for(ix=0;ix<=Bz;ix++) V[ix]=vmax;

			V[0]=vmin;

			for(ix=1;ix<=Bz;ix++){
				for(ik=0;ik<L;ik++){
					if((ix*1.0/Bz>Y[ik])&&(ix*1.0/Bz<=Y[ik+1])){
						V[ix]=c[ik]+(ix*1.0/Bz-Y[ik])*(c[ik+1]-c[ik])/(Y[ik+1]-Y[ik]);
						break;
					}
				}
			}
			V[Bz]=V[Bz]*1.005;

			cq2=ealloc3complex(nw,nk,Bz+1);

			for(ix=0;ix<Bz+1;ix++){
				for(iw=0,w=fw;iw<nw;++iw,w+=dw)
					for(ik=0,k=fk;ik<nk;++ik,k+=dk){

						if(w==0.0)w=1.0e-10/dt;

						kz1=1.0-pow(V[ix]/2.0*k/w,2.0);
						if(kz1>=0.00){
							phase1 =-w*sqrt(kz1)*dz*2.0/V[ix]+w*dz*2.0/V[ix];
							cshift1 = std::exp(complex(0.0f,phase1));
							cq2[ix][ik][iw] = cq[ik][iw]*cshift1;
						} else {
							phase1 =-w*sqrt(-kz1)*dz*2.0/V[ix];
							cshift1 = std::exp(complex(phase1,w*dz*2.0/V[ix]));
							cq2[ix][ik][iw] = cq[ik][iw]*cshift1;
						}
					}

				pfa2cc(1,2,nw,nk,cq2[ix][0]);

				for(ik=0;ik<nx;++ik)
					for(iw=0,w=fw;iw<nw;w+=dw,++iw){
						float ga=0.015,g=1.0;
						int I=20;

						if(ik<=I)g=exp(-ga*(I-ik)*(I-ik));
						if(ik>=nx-I)g=exp(-ga*(-nx+I+ik)*(-nx+I+ik));

						cq2[ix][ik][iw] = cq2[ix][ik][iw]*(g*1.0f/nxfft);
						cq2[ix][ik][iw] = ik%2 ? -cq2[ix][ik][iw] : cq2[ix][ik][iw];
					}
			}

			for(ix=0;ix<nx;++ix)
				for(ik=0;ik<Bz;++ik){

					if(((v[iz][ix]>=V[ik])&&(v[iz][ix]<V[ik+1]))) {

						v1=V[ik];v2=V[ik+1];

						for(iw=0,w=fw;iw<nw;w+=dw,++iw){

							a1=cq2[ik][ix][iw].real();a2=cq2[ik+1][ix][iw].real();
							theta1=cq2[ik][ix][iw].imag();theta2=cq2[ik+1][ix][iw].imag();

							a= a1*(v2-v[iz][ix])/(v2-v1)+a2*(v[iz][ix]-v1)/(v2-v1);
							theta=theta1*(v2-v[iz][ix])/(v2-v1)+theta2*(v[iz][ix]-v1)/(v2-v1);

							cp[ix][iw] =complex(a,theta);
						}

						break;
					}
				}
			free3complex(cq2);
			free1float(V);
		}
	}

	for(ix=0; ix<nx; ix++)
		for(iz=0; iz<nz; iz++) cresult_out[(size_t)ix*nz+iz] = cresult[ix][iz];

	free2float(p);
	free2complex(cp);
	free2complex(cq);
	free2complex(cq1);
	free2float(cresult);
	free2float(v);
	return 0;
}
