 /* Copyright (c) Colorado School of Mines, 2011.*/
 /* All rights reserved.                       */

 /* SUMIGPREPSPI: $Vision: 1.00 $ ; $Date: 2015/08/07 22:19:43 $*/

 #include "su.h"
#include <vector>
 #include "segy.h"
 #include "header.h"
 

 /*********************** self documentation ******************************/
 #if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
 "									",
 " SUMIGPREPSPI --- The 2-D PREstack commom-shot Phase-Shift-Plus 	",
 "			interpolation depth MIGration.			",
 "									",
 "   sumigprepspi <indata >outfile [parameters] 			", 
 "									",
 " Required Parameters:						   	",
 "									",
 " nxo=     number of total horizontal output samples			",
 " nxshot=  number of shot gathers to be migrated		   	",
 " nz=      number of depth samples				 	",
 " dx=      horizontal sampling interval			  	",   
 " dz=      depth sampling interval				 	",
 " vfile=   velocity profile, it must be binary format.                 ",
 " 									",
 " Optional Parameters:						   	",
 " fmax=25    the peak frequency of Ricker wavelet used as source wavelet",
 " f1=5,f2=10,f3=40,f4=50     frequencies to build a Hamming window     ",
 " lpad=9999,rpad=9999        number of zero traces padded on both	",
 "                            sides of depth section to determine the	",
 "                            migration aperture, the default values    ",
 "                            are using the full aperture.              ",
 " nflag=0    normalization of cross-correlation:                       ",
 "            0: none, 1: by source wave field                          ",
 " verbose=0  silent, =1 additional runtime information	                ",
 "									",  
 " Notes:								",
 " The input velocity file \'vfile\' consists of C-style binary floats.	",  
 " The structure of this file is vfile[iz][ix]. Note that this means that",
 " the x-direction is the fastest direction instead of z-direction! Such a",
 " structure is more convenient for the downward continuation type	",
 " migration algorithm than using z as fastest dimension as in other SU  ", 
 " programs.						   		",
 "									",
 " Because most of the tools in the SU package (such as  unif2, unisam2, ", 
 " and makevel) produce output with the structure vfile[ix][iz], you will",
 " need to transpose the velocity files created by these programs. You may",
 " use the SU program \'transp\' in SU to transpose such files into the  ",
 " required vfile[iz][ix] structure.					",
 "									",
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
  *	  Modified: Chris Stolk, 11 Dec 2005, - changed data input
  *		    to remove erroneous time delay.
  *	  Modified: CWP, John Stockwell 26 Sept 2006 - replaced Han's
  *	  "goto-loop" in two places with "do { }while loops".
  *	  Fixed it so that sx, gx, and scalco are honored.
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

/* Library version of SUMIGPREPSPI
 *
 * The main program (the parameters, the traces and their headers, the velocity file, the loop over shots) is left to the
 * caller. The phase shift plus interpolation depth migration of a shot, with the normalization of the image by the source
 * wavefield if nflag is not zero.
 *
 *	su_migprepspi_shot:	migrates one shot gather, p[nxo][nt] (the traces of the gather at the positions of their receivers, zero
 *		where there is no trace), whose shot is at the position ixshot, and adds its image to cresult[nxo][nz]
 *		(the sum of the shots). ix2 and ix3 are the first and last positions of the migrated aperture. The
 *		velocities are vp[nz][nxo]. Frequencies f1, f2, f3 and f4 are those of the Hamming window.
 *
 * Differences from the program, which are its bugs: the Hamming window of the program was shifted by the index of the first
 * frequency of the band (its taper was reversed); the number of frequencies and the aperture are checked.
 * Returns 0, or -3 for parameters that are not allowed.
 */

int su_migprepspi_shot(int nt, int nxo, int nz, float dt, float dx, float dz, int nflag, float fmax, float f1, float f2, float f3,
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
	const int L=10;
	int Bz=0;
	float c[101],P[100],Y[101];
	float Sz=0,vmax,v2;
	std::vector<float> V;
	double a,a1,a2,theta,theta1,theta2;
	complex ***cq2=NULL,***cq3=NULL;


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

		float eps=0.0f, rtmp;
		if (nflag) {
			/* stabilizing parameter for normalization by source wave field */
			eps = 0.0f;
			for (ix=0; ix<nx; ix++) {
				for (iw=0; iw<nw; iw++) {
					rtmp = std::abs(cp1[iw][ix]);
					eps = eps > rtmp ? eps : rtmp;
				}
			}
			eps *= 10.; /* ad-hoc scaling of stabilizing parameter */
			eps = eps * eps;
		}

		/* the imaging condition */
		for (ix=0; ix<nx; ix++) {
			for (iw=0,w=fw; iw<nw; w+=dw,iw++) {
				complex tmp, tmpd(0.0f,0.0f);
				float ratio=10.0;
				if (abs(ix+ix2-ixshot)*dx<ratio*iz*dz) {
					tmp = cp[iw][ix]*cp1[iw][ix];
					if (nflag) tmpd = cp1[iw][ix]*std::conj(cp1[iw][ix]);
				} else {
					tmp = complex(0.0f,0.0f);
					if (nflag) tmpd = complex(1.0f,0.0f);
				}
				if (nflag) /* normalization by source wave field */
					cresult_out[(size_t)(ix+ix2)*nz+iz] += tmp.real()/(tmpd.real()+eps)/ntfft;
				else /* no normalization  */
					cresult_out[(size_t)(ix+ix2)*nz+iz] += tmp.real()/ntfft;
			}
		}


		/* shift the wavefield with the velocity of each trace */
		for (ik=0; ik<nx; ++ik) {
			for (iw=0,w=fw; iw<nw; w+=dw,++iw) {
				cp[iw][ik] = cp[iw][ik]*std::exp(complex(0.0f,-w*dz/v[iz][ik]));
				cp1[iw][ik] = cp1[iw][ik]*std::exp(complex(0.0f,-w*dz/v[iz][ik]));
				cq[iw][ik] = ik%2 ? -cp[iw][ik] : cp[iw][ik];
				cq1[iw][ik] = ik%2 ? -cp1[iw][ik] : cp1[iw][ik];
			}
		}
		for (ik=nx; ik<nk; ++ik) {
			for (iw=0; iw<nw; ++iw) {
				cq[iw][ik] = complex(0.0f,0.0f);
				cq1[iw][ik] = complex(0.0f,0.0f);
			}
		}

		/* FFT to W-K domain */
		pfa2cc(-1,1,nk,nw,cq[0]);
		pfa2cc(-1,1,nk,nw,cq1[0]);

		vmax = v[iz][0]; vmin = v[iz][0];
		for (ix=0; ix<nx; ++ix) {
			if (v[iz][ix]>=vmax) vmax = v[iz][ix];
			if (v[iz][ix]<=vmin) vmin = v[iz][ix];
		}

		/* The second time phase shift */
		v1 = vmin;
		v2 = vmax;

		if ((v2-v1)/v1<0.01) {
			for (ik=0,k=fk; ik<nk; ++ik,k+=dk) {
				for (iw=0,w=fw; iw<nw; ++iw,w+=dw) {
					if (w==0.0) w = 1.0e-10/dt;
					kz1 = 1.0-pow(v1*k/w,2.0);
					if (kz1>=0.0) {
						phase1 = -w*sqrt(kz1)*dz/v1+w*dz/v1;
						cshift1 = complex(cos(phase1), sin(phase1));
						cq[iw][ik] = cq[iw][ik]*cshift1;
						cq1[iw][ik] = cq1[iw][ik]*cshift1;
					} else {
						phase1 = -w*sqrt(-kz1)*dz/v1;
						cshift1 = std::exp(complex(phase1,w*dz/v1));
						cq[iw][ik] = cq[iw][ik]*cshift1;
						cq1[iw][ik] = cq1[iw][ik]*cshift1;
					}
				}
			}

			pfa2cc(1,1,nk,nw,cq[0]);
			pfa2cc(1,1,nk,nw,cq1[0]);

			for (ix=0; ix<nx; ++ix) {
				for (iw=0; iw<nw; ++iw) {
					cq[iw][ix] = cq[iw][ix]*(1.0f/nxfft);
					cp[iw][ix] = ix%2 ? -cq[iw][ix] : cq[iw][ix];
					cq1[iw][ix] = cq1[iw][ix]*(1.0f/nxfft);
					cp1[iw][ix] = ix%2 ? -cq1[iw][ix] : cq1[iw][ix];
				}
			}
		} else {
			for (ik=0; ik<=L; ik++)
				c[ik] = vmin+ik*1.0*(vmax-vmin)/((float) L);

			memset((void *) P, 0, sizeof(P));
			for (ix=0; ix<nx; ix++) {
				for (ik=0; ik<L; ik++) {
					if (((v[iz][ix]>=c[ik]) && (v[iz][ix]<c[ik+1])) || ((ik==L-1) && (v[iz][ix]==vmax))) {
						P[ik] += 1.0/nx;
						break;
					}
				}
			}

			Sz = 0.0;
			for (ik=0; ik<L; ik++) {
				if (P[ik]!=0.00) Sz = Sz-P[ik]*log(P[ik]);
			}

			Bz = exp(Sz)+0.5;
			Y[0] = 0.0;
			Y[L] = 1.0;

			for (ik=1; ik<L; ik++) {
				Y[ik] = 0.0;
				for (ix=0; ix<ik; ix++)
					Y[ik] = Y[ik]+P[ix];
			}

			V.assign(Bz+1, vmax);
			V[0] = vmin;
			for (ix=1; ix<=Bz; ix++) {
				for (ik=0; ik<L; ik++) {
					if ((ix*1.0/Bz>Y[ik]) && (ix*1.0/Bz<=Y[ik+1])) {
						V[ix] = c[ik]+(ix*1.0/Bz-Y[ik])*(c[ik+1]-c[ik])/(Y[ik+1]-Y[ik]);
						break;
					}
				}
			}
			V[Bz] = V[Bz]*1.005;

			cq2 = ealloc3complex(nk,nw,Bz+1);
			cq3 = ealloc3complex(nk,nw,Bz+1);

			for (ix=0; ix<Bz+1; ix++) {
				for (iw=0,w=fw; iw<nw; ++iw,w+=dw)
					for (ik=0,k=fk; ik<nk; ++ik,k+=dk) {
						if (w==0.0) w = 1.0e-10/dt;
						kz1 = 1.0-pow(V[ix]*k/w,2.0);
						if (kz1>=0.00) {
							phase1 = -w*sqrt(kz1)*dz/V[ix]+w*dz/V[ix];
							cshift1 = std::exp(complex(0.0f,phase1));
						} else {
							phase1 = -w*sqrt(-kz1)*dz/V[ix];
							cshift1 = std::exp(complex(phase1,w*dz/V[ix]));
						}
						cq2[ix][iw][ik] = cq[iw][ik]*cshift1;
						cq3[ix][iw][ik] = cq1[iw][ik]*cshift1;
					}
				pfa2cc(1,1,nk,nw,cq2[ix][0]);
				pfa2cc(1,1,nk,nw,cq3[ix][0]);

				for (ik=0; ik<nx; ++ik)
					for (iw=0,w=fw; iw<nw; w+=dw,++iw) {
						cq2[ix][iw][ik] = cq2[ix][iw][ik]*(1.0f/nxfft);
						cq2[ix][iw][ik] = ik%2 ? -cq2[ix][iw][ik] : cq2[ix][iw][ik];
						cq3[ix][iw][ik] = cq3[ix][iw][ik]*(1.0f/nxfft);
						cq3[ix][iw][ik] = ik%2 ? -cq3[ix][iw][ik] : cq3[ix][iw][ik];
					}
			}

			for (ix=0; ix<nx; ++ix) {
				for (ik=0; ik<Bz; ++ik) {
					if (((v[iz][ix]>=V[ik]) && (v[iz][ix]<V[ik+1]))) {
						v1 = V[ik]; v2 = V[ik+1];
						for (iw=0,w=fw; iw<nw; w+=dw,++iw) {
							a1 = cq2[ik][iw][ix].real();
							a2 = cq2[ik+1][iw][ix].real();
							theta1 = cq2[ik][iw][ix].imag(); theta2 = cq2[ik+1][iw][ix].imag();
							a = a1*(v2-v[iz][ix])/(v2-v1)+a2*(v[iz][ix]-v1)/(v2-v1);
							theta = theta1*(v2-v[iz][ix])/(v2-v1)+theta2*(v[iz][ix]-v1)/(v2-v1);
							cp[iw][ix] = complex(a,theta);

							a1 = cq3[ik][iw][ix].real(); a2 = cq3[ik+1][iw][ix].real();
							theta1 = cq3[ik][iw][ix].imag();
							theta2 = cq3[ik+1][iw][ix].imag();
							a = a1*(v2-v[iz][ix])/(v2-v1)+a2*(v[iz][ix]-v1)/(v2-v1);
							theta = theta1*(v2-v[iz][ix])/(v2-v1)+theta2*(v[iz][ix]-v1)/(v2-v1);
							cp1[iw][ix] = complex(a,theta);
						}
						break;
					}
				}
			}
			free3complex(cq2);
			free3complex(cq3);
		}

	}

	free2complex(cp);
	free2complex(cp1);
	free2float(v);

	free2complex(cq);
	free2complex(cq1);

	return 0;
}
