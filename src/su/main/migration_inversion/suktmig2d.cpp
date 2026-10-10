/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.		       */

/* SUKTMIG2D: $Revision: 1.7 $ ; $Date: 2013/08/14 18:34:35 $*/

#include "su.h"
#include "segy.h"
/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 								",
" SUKTMIG2D - prestack time migration of a common-offset	",
"	section with the double-square root (DSR) operator	",
"								",
" 								",
"   suktmig2d < infile vfile= [parameters]  > outfile		",
" 								",
" Required Parameters:						",
" vfile=	rms velocity file (units/s) v(t,x) as a function",
"		of time						",
" dx=		distance (units) between consecutive traces	",
"								",
" Optional parameters:						",
" fcdpdata=tr.cdp	first cdp in data			",
" firstcdp=fcdpdata	first cdp number in velocity file	",
" lastcdp=from header	last cdp number in velocity file	",
" dcdp=from header	number of cdps between consecutive traces",
" angmax=40	maximum aperture angle for migration (degrees)	",
" hoffset=.5*tr.offset		half offset (m)			",
" nfc=16	number of Fourier-coefficients to approximate	",
"		low-pass					",
"		filters. The larger nfc the narrower the filter	",
" fwidth=5 	high-end frequency increment for the low-pass	",
" 		filters						",
" 		in Hz. The lower this number the more the number",
"		of lowpass filters to be calculated for each 	",
"		input trace.					",
"								",
" Caveat: this code may need some work				",
" Notes:							",
" Data must be preprocessed with sufrac to correct for the	",
" wave-shaping factor using phasefac=.25 for 2D migration.	",
"								",
" Input traces must be sorted into offset and cdp number.	",
" The velocity file consists of rms velocities for all CMPs as a",
" function of vertical time and horizontal position v(t,x)	",
" in C-style binary floating point numbers.  It's easiest to 	",
" supply v(t,x) that has the same dimensions as the input data to",
" be migrated. Note that time t is the fast dimension in these  ",
" the input velocity file.					",
" 								",
" The units may be feet or meters, as long as these are		",
" consistent.							",
" Antialias filter is performed using (Gray,1992, Geoph. Prosp), ",
" using nc low- pass filtered copies of the data. The cutoff	",
" frequencies are calculated  as fractions of the Nyquist	",
" frequency.							",
"								",
" The maximum allowed angle is 80 degrees(a 10 degree taper is ",
" applied to the end of the aperture)				",
NULL};
#endif


/* Library version of SUKTMIG2D
 *
 * The main program (the parameters, the traces and their headers, the velocity file, the cdp bookkeeping) is left to the caller.
 *
 *	su_ktmig2d:	migrates the common-offset section data[ntr][nt] with the rms velocities vel[ntr][nt] (those at the cdp of
 *			each trace: the caller selects them) to mig[ntr][nt]
 *
 * Differences from the program:
 *  - The low-pass filters are made once for the cut frequencies (the program made them again for each trace).
 *  - On the left of the midpoint the program used 'fphi=++fplo', which moved the low cut filter along with the high one;
 *    it is fplo and fplo+1 on both sides, as it is on the right.
 *  - The checks of the parameters return a code.
 */

static void lpfilt(int nfc, int nfft, float dt, float fhi, float *filter)
{
	int i,j;  /* counters */
	int nf;   /* Number of frequencies (including Nyquist) */
	float onfft;  /* reciprocal of nfft */
	float fn; /* Nyquist frequency */
	float df; /* frequency interval */
	float dw; /* frequency interval in radians */
	float whi;/* cut-frequency in radians */
	float w;  /* radian frequency */

	nf= nfft/2 + 1;
	onfft=1.0/nfft;
	fn=1.0/(2*dt);
	df=onfft/dt;
	whi=fhi*PI/fn;
	dw=df*PI/fn;

	for(i=0; i<nf; ++i){
		filter[i]= whi/PI;
		w=i*dw;

		for(j=1; j<nfc; ++j){
			float c= sin(whi*j)*sin(PI*j/nfc)*2*nfc/(PI*PI*j*j);
			filter[i] +=c*cos(j*w);
		}
	}
}

/* one image point of the operator: the amplitude of the input at the time t, from the low-pass filtered versions */
static void kt_add(float **mig, int iip, int it, float geoms, float obliq, float angtaper, float t, int fplo, int nc, float wlo,
	float whi, float dt, int nt, float **lowpass)
{
	float datalo[8], datahi[8], amplo, amphi, firstt;
	int itb, ite, k, fphi = fplo+1;

	itb = MAX((int)ceil(t/dt)-3,0);
	ite = MIN(itb+8,nt);
	firstt = (itb-1)*dt;

	/* Move energy from CMP to CIP */
	if (fplo>=nc) {
		for (k=0; k<8; ++k) datalo[k] = 0.0;
		for (k=itb; k<ite; ++k)
			datalo[k-itb]=lowpass[nc][k];
		ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
		mig[iip][it] += geoms*obliq*angtaper*amplo;
	} else {
		for (k=0; k<8; ++k) datalo[k] = datahi[k] = 0.0;
		for (k=itb; k<ite; ++k) {
			datalo[k-itb]=lowpass[fplo][k];
			datahi[k-itb]=lowpass[fphi][k];
		}
		ints8r(8,dt,firstt,datalo,0.0,0.0,1,&t,&amplo);
		ints8r(8,dt,firstt,datahi,0.0,0.0,1,&t,&amphi);
		mig[iip][it] += geoms*obliq*angtaper*(wlo*amplo + whi*amphi);
	}
}

int su_ktmig2d(int ntr, int nt, float dt, float dx, float h, float angmax, int nfc, int fwidth,
	const float *data_in, const float *vel_in, float *mig_out)
{
	int imp,iip,it,ifc,i;
	int nc,nfft,nf;
	float x,mp,ip,t,t0,tmax,fnyq,v,ang,angtaper=0.0,geoms,obliq,pmin,p;
	float **data,**vel,**mig,**lowpass,**filters,*fc,*rtin,*rtout;
	complex *ct;

	if (ntr<2 || nt<2 || dt<=0.0 || dx<=0.0 || nfc<1 || fwidth<1 || angmax<0.0 || angmax>80.0) return -3;
	tmax = (nt-1)*dt;

	/* Set up FFT parameters */
	nfft = npfaro(nt, 2*nt);
	if (nfft >= 720720) return -3;
	nf = nfft/2 + 1;

	/* Determine number of filters for antialiasing */
	fnyq = 1.0/(2*dt);
	nc = (int)ceil(fnyq/fwidth);

	data = alloc2float(nt,ntr);
	vel = alloc2float(nt,ntr);
	mig = alloc2float(nt,ntr);
	lowpass = alloc2float(nt,nc+1);
	filters = alloc2float(nf,nc+1);
	fc = alloc1float(nc+1);
	rtin = ealloc1float(nfft);
	rtout = ealloc1float(nfft);
	ct = ealloc1complex(nf);
	memcpy(data[0],data_in,sizeof(float)*(size_t)ntr*nt);
	memcpy(vel[0],vel_in,sizeof(float)*(size_t)ntr*nt);
	memset(mig[0],0,sizeof(float)*(size_t)ntr*nt);
	memset(lowpass[0],0,sizeof(float)*(size_t)nt*(nc+1));
	memset(rtin,0,sizeof(float)*nfft);
	fc[0] = 0.0;

	/* the cut frequencies of the low-pass filters */
	for (i=1; i<nc+1; ++i) {
		fc[i] = fnyq*i/nc;
		lpfilt(nfc,nfft,dt,fc[i],filters[i]);
	}
	for (i=0; i<nf; ++i) filters[0][i] = 0.0;

	/* Loop over input mid-points */
	for (imp=0; imp<ntr; ++imp) {
		mp=imp*dx;

		/* low-pass filtered versions of the data, for antialiasing */
		for (it=0; it<nt; ++it) rtin[it] = data[imp][it];
		for (it=nt; it<nfft; ++it) rtin[it] = 0.0;
		for (ifc=1; ifc<nc+1; ++ifc) {
			memset((void *) rtout, 0, nfft*FSIZE);
			pfarc(1,nfft,rtin,ct);
			for (it=0; it<nf; ++it)
				ct[it] *= filters[ifc][it];
			pfacr(-1,nfft,ct,rtout);
			for (it=0; it<nt; ++it)
				lowpass[ifc][it] = rtout[it];
		}

		/* Loop over vertical traveltimes */
		for (it=0; it<nt; ++it) {
			int lx,ux,side;

			t0=it*dt;
			v=vel[imp][it];
			if (v<=0.0) continue;
			{
				float xmax=tan((angmax+10.0)*PI/180.0)*v*t0;
				lx=MAX(0,imp - (int)ceil(xmax/dx));
				ux=MIN(ntr,imp + (int)ceil(xmax/dx));
			}

			/* the image points to the left of the midpoint, then to the right */
			for (side=0; side<2; ++side) {
				int first = side==0 ? imp : imp+1;
				int step = side==0 ? -1 : 1;
				for (iip=first; side==0 ? iip>lx : iip<ux; iip+=step) {
					float ts,tr,ref,wlo,whi;
					int fplo=0;

					ip=iip*dx;
					x=ip-mp;
					ts=sqrt( pow(t0/2,2) + pow((x+h)/v,2) );
					tr=sqrt( pow(t0/2,2) + pow((h-x)/v,2) );
					t= ts + tr;
					if(t>=tmax) break;
					geoms=sqrt(1/(t*v));
					obliq=sqrt(.5*(1 + (t0*t0/(4*ts*tr))
						- (1/(ts*tr))*sqrt(ts*ts - t0*t0/4)*sqrt(tr*tr - t0*t0/4)));
					ang=180.0*fabs(acos(t0/t))/PI;
					if(ang<=angmax) angtaper=1.0;
					if(ang>angmax) angtaper=cos((ang-angmax)*PI/20);

					/* Evaluate migration operator slowness p to determine */
					/* the low-pass filtered trace for antialiasing */
					pmin=1/(2*dx*fnyq);
					p=fabs((x+h)/(pow(v,2)*ts) + (x-h)/(pow(v,2)*tr));
					if(p>0){fplo=(int)floor(nc*pmin/p);}
					if(p==0){fplo=nc;}
					ref=fmod(nc*pmin,p);
					wlo=1-ref;
					whi=ref;
					kt_add(mig,iip,it,geoms,obliq,angtaper,t,fplo,nc,wlo,whi,dt,nt,lowpass);
				}
			}
		}
	}

	memcpy(mig_out,mig[0],sizeof(float)*(size_t)ntr*nt);

	free2float(data); free2float(vel); free2float(mig); free2float(lowpass); free2float(filters);
	free1float(fc); free1float(rtin); free1float(rtout); free1complex(ct);
	return 0;
}

/* The tables of the sinc interpolators (ints8r and ints8c keep them in static variables, made at the first call, which is not
 * safe from several threads at once): call it once while nothing else is running. */
void su_ints8_tables(void)
{
	float yr = 0.0, xr = 0.0, outr;
	complex yc(0.0f,0.0f), outc;
	ints8r(1,1.0,0.0,&yr,0.0,0.0,1,&xr,&outr);
	ints8c(1,1.0,0.0,&yc,yc,yc,1,&xr,&outc);
}
