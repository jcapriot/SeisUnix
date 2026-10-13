/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

#include "su.h"
#include <vector>
#include "segy.h"
#include "header.h"
#include <signal.h>

/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[]={
"									",
" SUDATUMFD - 2D zero-offset Finite Difference acoustic wave-equation	",
"		 DATUMing    						",
"									",
" sudatumfd <stdin > stdout [optional parameters]			",
"									",
" Required parameters:						   	",
"									",
" nt=	   number of time samples on each trace	       			",
" nx=	   number of receivers per shot gather				",
" nsx=	  number of shot gathers				    	",
" nz=	   number of downward continuation depth steps			",
" dz=	   depth sampling interval (in meters)				",
" mx=	   number of horizontal samples in the velocity model		",
" mz=	   number of vertical samples in the velocity model		",
" vfile1=       velocity file used for thin-lens term	    		",
" vfile2=       velocity file used for diffraction term			",
" dx=           horizontal sampling interval (in meters)                ",
"									",
" Optional parameters:						   	",
"									",
" dt=.004       time sampling interval (in seconds)			",
" buff=5	number of zero traces added to each side of each   	",
"	     shot gather as a pad			       		",
" tap_len=5     taper length (in number of traces)			",
" x_0=0.0       x coordinate of leftmost position in velocity model     ",
"									",
" Notes:								",
" The algorithm is a 45-degree implicit-finite-difference scheme based  ",
" on the one-way wave equation.  It works on poststack (zero-offset)    ",
" data only.  The two velocity files, vfile1 and vfile2, are binary     ",
" files containing floats with the format v[ix][iz].  There are two     ",
" potentially different velocity files for the thin-lens and            ",
" diffraction terms to allow for the use of a zero-velocity layer       ",
" which allows for datuming from an irregular surface.                  ",
"                                                                       ",
" Source and receiver locations must be set in the header values in     ",
" order for the datuming to work properly.  The leftmost position of    ",
" of the velocity models given in vfile1 and vfile2 must also be given. ",
" 									",
" 									",
NULL};
#endif

/* 
 * Author:  Chris Robinson, 10/16/00, CWP, Colorado School of Mines
 *
 *
 * References:
 *  Beasley, C., and Lynn, W., 1992, The zero-velocity layer: migration
 *    from irregular surfaces: Geophysics, 57, 1435-1443.
 *
 *  Claerbout, J. F., 1985, Imaging the earth's interior:  Blackwell
 *    Scientific Publications.
 *
 */

/**************** end self doc *******************************************/

/**************** end self doc *******************************************/

/* The datuming of the program for the shot gathers in memory.
 *
 * nt, dt		samples and sample interval of the traces
 * nx, nsx		receivers per shot gather, and number of shot gathers
 * nz, dz		number of downward continuation steps, and their depth (m)
 * dx			horizontal sampling interval (m) of the velocity models
 * x_0			x coordinate (m) of the leftmost position in the velocity models
 * buff			number of zero traces added to each side of each shot gather as a pad
 * tap_len		taper length (in traces)
 * mx, mz		number of horizontal and vertical samples of the velocity models
 * vel1, vel2		[mx][mz] the velocities for the thin lens term and for the diffraction term
 * gx[nsx][nx]		the x coordinate of each receiver (that selects the column of the velocity)
 * in[nsx][nx][nt]	the traces
 * out[nsx][nx][nt]	the datumed traces
 *
 * Returns 0, or 1 if a receiver is outside of the velocity models (its column is the nearest one then).
 *
 * Differences from the program, which are its bugs: the frequency sampling is 2 pi / ((ntfft - 1) dt) where the
 * transform has 2 pi / (ntfft dt); and a receiver outside of the velocity model reads outside of it (the nearest column
 * here).
 */
int su_datumfd(int nt, float dt, int nx, int nsx, int nz, float dz, float dx, float x_0, int buff, int tap_len, int mx,
	int mz, const float *vel1, const float *vel2, const float *gx, const float *in, float *out)
{
	typedef std::complex<float> cf;
	int status = 0;
	const int ntfft = npfar(2*nt);
	const int nw = ntfft/2 + 1;
	const float dw = 2.0*PI/(ntfft*dt);
	const int nx_pad = nx + 2*buff;
	const float sigma = 0.5, tau = 0.25;
	int isx, ix, it, iw, n, j, k, m, iz;

	std::vector<float> p(ntfft);
	std::vector<cf> cp(nw);
	std::vector<cf> ifr1((size_t) nx_pad*nw), ifr2((size_t) nx_pad*nw);
	std::vector<cf> e(nx_pad), f(nx_pad), d(nx_pad), ohm(nx_pad), alpha(nx_pad), beta(nx_pad), gamma(nx_pad);
	std::vector<float> tap_func(nx_pad);
	std::vector<int> col(nx_pad);

	/* the tapering function */
	if (tap_len == 0) {
		for (m=0; m<nx_pad; m++) tap_func[m] = 1.0;
	} else {
		for (m=0; m<nx_pad; m++) {
			float v = 0.0;
			if (m<buff) v = 0.0;
			if (m>=buff && m<(buff+tap_len-1)) v = sin(2.0*PI*(m+1-buff)/(4.0*tap_len));
			if (m>=(buff+tap_len-1) && m<=(buff+nx-1-(tap_len-1))) v = 1.0;
			if (m>(buff+nx-1-(tap_len-1)) && m<=(buff+nx-1))
				v = cos(2.0*PI*(m-(buff+nx-1-(tap_len-1)))/(4.0*tap_len));
			if (m>(buff+nx-1)) v = 0.0;
			tap_func[m] = v;
		}
	}

	for (isx=0; isx<nsx; isx++) {
		const float *gxs = gx + (size_t) isx*nx;

		/* the column of the velocity of each trace position (the padding takes those of the ends) */
		for (n=0; n<nx_pad; n++) {
			int a;
			if (n < buff) a = (int) gxs[0];
			else if (n > nx_pad-1-buff) a = (int) gxs[nx-1];
			else a = (int) gxs[n-buff];
			j = (int) ((a-x_0)/dx);
			if (j < 0) { j = 0; status = 1; }
			if (j > mx-1) { j = mx-1; status = 1; }
			col[n] = j;
		}

		/* temporal Fourier transform of the traces, padded in x */
		std::fill(ifr1.begin(), ifr1.end(), cf(0.0f, 0.0f));
		for (ix=0; ix<nx; ix++) {
			const float *x = in + ((size_t) isx*nx + ix)*nt;
			for (it=0; it<nt; it++) p[it] = x[it];
			for (it=nt; it<ntfft; it++) p[it] = 0.0;
			pfarc(1, ntfft, p.data(), cp.data());
			for (iw=0; iw<nw; iw++) ifr1[(size_t) (ix+buff)*nw + iw] = cp[iw];
		}

		for (iz=0; iz<nz; iz++) {
			/* taper */
			for (n=0; n<nx_pad; n++)
				for (iw=0; iw<nw; iw++) ifr2[(size_t) n*nw + iw] = tap_func[n]*ifr1[(size_t) n*nw + iw];

			/* thin lens term */
			for (n=0; n<nx_pad; n++) {
				float v = 0.5*vel1[(size_t) col[n]*mz + iz];
				float ig = (-1.0/v)*dz;
				for (k=0; k<nw; k++) {
					cf mult(cos(k*dw*ig), sin(k*dw*ig));
					ifr1[(size_t) n*nw + k] = ifr2[(size_t) n*nw + k]*mult;
				}
			}

			/* diffraction term */
			for (k=0; k<nw; k++) {
				float w = k*dw;
				if (w < 5.0) w = 5.0;
				for (n=0; n<nx_pad; n++) {
					float v = 0.5*vel2[(size_t) col[n]*mz + iz];
					alpha[n] = cf(0.0, -1.0*v*sigma*dz/(dx*dx*2.0*w));
					beta[n] = alpha[n]*(cf(1.0, 0.0) + cf(0.0, -2.0*v*tau/(w*sigma*dz)));
					gamma[n] = alpha[n]*(cf(1.0, 0.0) + cf(0.0, 2.0*v*tau/(w*sigma*dz)));
				}
				auto u = [&](int row) -> cf & { return ifr1[(size_t) row*nw + k]; };

				/* the d[n]'s */
				d[0] = (cf(1.0, 0.0) + cf(-2.0, 0.0)*gamma[0])*u(0) + gamma[0]*u(1);
				ohm[0] = cf(-1.0, 0.0)*(d[0]/beta[0]);
				for (n=1; n<nx_pad-1; n++) {
					d[n] = gamma[n]*u(n-1) + (cf(1.0, 0.0) + cf(-2.0, 0.0)*gamma[n])*u(n) + gamma[n]*u(n+1);
					ohm[n] = cf(-1.0, 0.0)*(d[n]/beta[n]);
				}
				d[nx_pad-1] = gamma[nx_pad-1]*u(nx_pad-2) +
					(cf(1.0, 0.0) + cf(-2.0, 0.0)*gamma[nx_pad-1])*u(nx_pad-1);
				ohm[nx_pad-1] = cf(-1.0, 0.0)*(d[nx_pad-1]/beta[nx_pad-1]);

				/* the e's and f's */
				e[1] = cf(1.0, 0.0)/beta[1] + cf(2.0, 0.0);
				f[1] = ohm[1];
				for (m=2; m<nx_pad-1; m++) {
					e[m] = (cf(1.0, 0.0)/beta[m] + cf(2.0, 0.0)) + cf(-1.0, 0.0)/e[m-1];
					f[m] = ohm[m] + f[m-1]/e[m-1];
				}
				e[0] = cf(0.0, 0.0);
				e[nx_pad-1] = cf(0.0, 0.0);

				/* downward continue one frequency of one shot gather */
				ifr2[(size_t) (nx_pad-1)*nw + k] = cf(0.0, 0.0);
				ifr2[(size_t) (nx_pad-2)*nw + k] = cf(-1.0, 0.0)*(f[nx_pad-2]/e[nx_pad-2]);
				for (j=nx_pad-3; j>0; j--)
					ifr2[(size_t) j*nw + k] = (ifr2[(size_t) (j+1)*nw + k] + cf(-1.0, 0.0)*f[j])/e[j];
				ifr2[(size_t) 0*nw + k] = cf(0.0, 0.0);
			}
			ifr1 = ifr2;
		}

		/* temporal inverse Fourier transform of the unpadded traces */
		for (ix=0; ix<nx; ix++) {
			float *y = out + ((size_t) isx*nx + ix)*nt;
			for (iw=0; iw<nw; iw++) cp[iw] = ifr1[(size_t) (ix+buff)*nw + iw];
			pfacr(-1, ntfft, cp.data(), p.data());
			for (it=0; it<nt; it++) y[it] = p[it]/(1.0*ntfft);
		}
	}
	return status;
}
