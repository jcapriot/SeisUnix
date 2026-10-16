/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUMIGSPLIT: $Revision: 1.9 $ ; $Date: 2011/11/16 22:14:43 $       */

#include "su.h"
#include "segy.h"
#include "header.h"
#include <vector>

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                       ",
" SUMIGSPLIT - Split-step depth migration for zero-offset data.         ",
"                                                                       ",
" sumigsplit <infile >outfile vfile= [optional parameters]              ",
"                                                                       ",
" Required Parameters:                                                  ",
" nz=                   number of depth samples                         ",
" dz=                   depth sampling interval                         ",
" vfile=                name of file containing velocities              ",
"                                                                       ",
" Optional Parameters:                                                  ",
" dt=from header(dt) or .004    time sampling interval                  ",
" dx=from header(d2) or 1.0     midpoint sampling interval              ",
" ft=0.0                        first time sample                       ",
" fz=                           first depth sample                      ",
"                                                                       ",
" tmpdir=        if non-empty, use the value as a directory path        ",
"                prefix for storing temporary files; else if the        ",
"                the CWP_TMPDIR environment variable is set use         ",
"                its value for the path; else use tmpfile()             ",
"                                                                       ",
"                                                                       ",
" Notes:                                                                ",
" The input velocity file \'vfile\' consists of C-style binary floats.  ",
" The structure of this file is vfile[iz][ix]. Note that this means that",
" the x-direction is the fastest direction instead of z-direction! Such a",
" structure is more convenient for the downward continuation type       ",
" migration algorithm than using z as fastest dimension as in other SU  ",
" programs.                                                     	",
"                                                                       ",
" Because most of the tools in the SU package (such as  unif2, unisam2, ",
" and makevel) produce output with the structure vfile[ix][iz], you will",
" need to transpose the velocity files created by these programs. You may",
" use the SU program \'transp\' in SU to transpose such files into the  ",
" required vfile[iz][ix] structure.                                     ",
" (In C  v[iz][ix] denotes a v(x,z) array, whereas v[ix][iz]  		",
" denotes a v(z,x) array, the opposite of what Matlab and Fortran	",
" programmers may expect.)						", 
"									",
NULL};
#endif
  
/*
 * Credits: CWP Baoniu Han, July 21th, 1997
 *
 * Reference:
 *  Stoffa, P. L. and Fokkema, J. T. and Freire, R. M. and Kessinger, W. P.,
 *  1990, Split-step Fourier migration, Geophysics, 55, 410-421.
 *
 *
 * Trace header fields accessed: ns, dt, delrt, d2
 * Trace header fields modified: ns, dt, delrt
 */

/**************** end self doc *******************************************/


/* Library version of SUMIGSPLIT
 *
 * The main program (the parameters, the traces and their headers, the temporary files, the velocity file) is left to the caller.
 *
 *	su_migsplit:	migrates the zero-offset section pin[nx][nt] (the traces in order) to cresult_out[nx][nz], with the
 *			velocities vin[nz][nx], by split-step Fourier depth migration
 *
 * Returns 0, or -3 for parameters that are not allowed.
 */
int su_migsplit(int nx, int nt, int nz, float dz, float dt, float dx, const float *pin, const float *vin,
	float *cresult_out)
{
	int ik,iz,iw,ix,it;
	int nxfft,ntfft;
	int nk,nw;
	float dk,dw;
	float fk,fw;
	float k,w;
	float vmin,v1;
	double kz1,kz2;
	double phase1;
	std::complex<float> cshift1,cshift2;
	float **p,**cresult,**v;
	complex **cp,**cq,**cq1;

	if (nx<4 || nt<2 || nz<1 || dz<=0.0 || dt<=0.0 || dx<=0.0) return -3;
	for (iz=0; iz<nz; iz++)
		for (ix=0; ix<nx; ix++)
			if (vin[(size_t)iz*nx+ix] <= 0.0) return -3;

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
	v = alloc2float(nx,nz);

	for (ix=0; ix<nx; ix++) {
		for (it=0; it<nt; it++) p[ix][it] = pin[(size_t)ix*nt+it];
		for (it=nt; it<ntfft; it++) p[ix][it] = 0.0;
	}
	for (iz=0; iz<nz; iz++)
		for (ix=0; ix<nx; ix++) v[iz][ix] = vin[(size_t)iz*nx+ix];

	/* Fourier transform t to w */
	pfa2rc(1,1,ntfft,nx,p[0],cp[0]);

	/* loops over depth */
	for (iz=0; iz<nz; ++iz) {
		for (ix=0; ix<nx; ix++) {
			cresult[ix][iz] = 0.0;
			for (iw=0; iw<nw; iw++)
				cresult[ix][iz] += cp[ix][iw].real()/ntfft;

			for (iw=0; iw<nw; iw++)
				cp[ix][iw] = complex(cp[ix][iw].real()-cresult[ix][iz]/ntfft, cp[ix][iw].imag());
		}

		vmin = 0;
		for (ix=0; ix<nx; ix++)
			vmin += 1.0/v[iz][ix]/nx;
		vmin = 1.0/vmin*0.5;

		for (ik=0; ik<nx; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = ik%2 ? -cp[ik][iw] : cp[ik][iw];

		for (ik=nx; ik<nk; ++ik)
			for (iw=0; iw<nw; ++iw)
				cq[ik][iw] = complex(0.0f,0.0f);

		/* FFT to W-K domain */
		pfa2cc(-1,2,nw,nk,cq[0]);

		v1 = vmin;

		for (ik=0,k=fk; ik<nk; ++ik,k+=dk)
			for (iw=0,w=fw; iw<nw; ++iw,w+=dw) {
				if (w==0.0) w = 1.0e-10/dt;
				kz1 = 1.0-pow(v1*k/w,2.0);
				if (kz1>0) {
					phase1 = -w*sqrt(kz1)*dz/v1;
					cshift1 = complex(cos(phase1), sin(phase1));
					cq1[ik][iw] = cq[ik][iw]*cshift1;
				} else {
					cq1[ik][iw] = complex(0.0f,0.0f);
				}
			}

		pfa2cc(1,2,nw,nk,cq1[0]);

		for (ix=0; ix<nx; ++ix)
			for (iw=0,w=fw; iw<nw; w+=dw,++iw) {
				cq1[ix][iw] = cq1[ix][iw]*(1.0f/nxfft);
				cq1[ix][iw] = ix%2 ? -cq1[ix][iw] : cq1[ix][iw];
				kz2 = (1.0/v1-2.0/v[iz][ix])*w*dz;
				cshift2 = complex(cos(kz2),sin(kz2));
				cp[ix][iw] = cq1[ix][iw]*cshift2;
			}
	}

	for (ix=0; ix<nx; ix++)
		for (iz=0; iz<nz; iz++) cresult_out[(size_t)ix*nz+iz] = cresult[ix][iz];

	free2float(p);
	free2complex(cp);
	free2complex(cq);
	free2complex(cq1);
	free2float(cresult);
	free2float(v);
	return 0;
}
