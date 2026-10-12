/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUALFORD: $Revision: 1.6 $ ; $Date: 2011/11/16 22:58:31 $	*/

#include "su.h"
#include <vector>
#include "segy.h" 
#include "header.h"
/*********************** self documentation ******************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUALFORD - trace by trace Alford Rotation of shear wave data volumes  ",
" 									",
" sualford inS11=file1 inS22=file2 inS12=file3 inS21=file4		",
" outS11=file5 outS22=file6 outS12=file7 outS21=file8 [optional         ",
" parameters]                                                           ",
"                                                                       ",
" Required Parameters:                                                  ",
" inS11=	input data volume for the 11 component			",
" inS12=	input data volume for the 12 component			",
" inS21=	input data volume for the 21 component			",
" inS22=	input data volume for the 22 component			",
" outS11=	output data volume for the 11 component			",
" outS12=	output data volume for the 11 component			",
" outS21=	output data volume for the 11 component			",
" outS22=	output data volume for the 11 component			",
"                                                                       ",
" Optional parameters:                                                  ",
" angle_inc=               sets the increment to the angle by which	",
"                         the data sets are rotated. The minimum is     ",
"                         set to be 1 degree and default is 5.          ",
" Az_key=                  to set the header field storing the azimuths	",
"                         for the fast shear wave on the output volumes ",
" Q_key=                   to set the header field storing the quality	",
"                         factors of performed optimum rotations        ",
" lag_key=                 to set the header field storing the lag in	",
"                         miliseconds the fast and slow shear components",
" xcorr_key=               to set the header field containing the maxi-	",
"                         mum normalized cross-correlation between the	",
"                         and slow shear waves.                         ",
" taper=		  2*taper+1 is the length of the sample overlap ",
"			  between the unrotated data with the rotated   ",
"			  data on the traces. The boundary between them ",
" 			  is defined by time windowning.                ",
"				taper = -1, for no-overlap		",
"				taper = 0, for overlap of one sample	",
"				taper =>1, for use of cosine scale to   ",
"					   to interpolate between the 	",
"					   unrotated and rotated parts	",
"					   of the traces		",
"									",
" taperwin=               another taper used to taper the data within   ",
"			  the window of analysis to diminish the effect ",
"                         of data near the window edges.In this way one ",
"                         can focus on a given reflector. Also given in ",
"                         number of samples                             ",
"									",
" maxlag=		  maximum limit in ms for the lag between fast  ",
" 			  and slow shear waves. If this threshold is 	",
"			  attained or surpassed, the quality factor for	",
"			  the rotation is zeroed as well as all the     ",
"			  parameters found for that certain rotation 	",
"			 						",
"									",
" ntraces=		  number of traces to be used per computation   ",
"			  ntraces=3 will use three adjacent traces to   ",
"		          compute the angle of rotation                 "
"                                                                       ",
" Notes:                                                                ",
"                                                                       ",
" The Alford Rotation is a method to rotate the four components         ",
" of a shear wave survey into its natural coordinate system, where      ",
" the fast and slow shear correspond to the inline to inline shear (S11)",
" and xline to xline (S22) volumes, respectively.                       ",
"                                                                       ",
" This Alford Rotation code tries to maximize the energy in the         ",
" diagonal volumes, i.e., S11 and S22, while minimizing the energy      ",
" in the off-diagonals, i.e., in volumes S12 and S21, in a trace by     ",
" trace manner. It then returns the new rotated volumes, saving the     ",
" the quality factor for the rotation and azimuth angle of the fast     ",
" shear wave direction for each trace headers of the new rotated S11    ",
" volume.                                                               ",
"                                                                       ",
" The fields in the header containing the Azimuth and Quality factor    ",
" and the sample lag between fast and slow shear are otrav, grnolf and  ",
" grnofr, respectively, by default. The values are multiplied by ten in ",
" the case of the angles and by a thousand for quality factors. To      ",
" change this defaults use the optional parameters Az_key, Q_key and    ",
" lag_key                                                             	",
"                                                                       ",
"			   						", 
" modified header fields:                                               ",
" the ones specified by Az_key, Q_key, lag_key and xcorr_key. By default",
" these are otrav, grnlof, tstat and grnors, respectively.            	",
NULL};
#endif

/* Credits:
 *	CWP: Rodrigo Felicio Fuck
 *      Code translated and adapted from original version in Fortran
 *      by Ted Schuck (1993)
 *
 *
 * Schuck, E. L. , 1993, Multicomponent, three dimensional seismic 
 * characterization of a fractured coalbed methane reservoir, 
 * Cedar Hill Field, San Juan County, New Mexico, Ph.D. Thesis,
 * Colorado School of Mines
 *
 */
/**************** end self doc *******************************************/




/* RMS amplitude of a trace */
static float RMS_amp(const float *_tr, int nt)
{
	float _amp = 0;
	int i;
	for (i=0; i<=nt-1; ++i) _amp += _tr[i]*_tr[i];
	return sqrt(_amp/nt);
}

/* The quality factor of a rotation: the energy of the diagonal traces against that of the off-diagonal ones, as defined by
   Schuck (1993). mode 0: 1 - (S12 + S21) / (S11 + S22) (amplitudes) */
static float Quality_Factor(const float *_S11, const float *_S12, const float *_S21, const float *_S22, int nt, int mode)
{
	float amp_S11, amp_S12, amp_S21, amp_S22;
	float Q, min, max, num, den;
	float noise = 1.0E-30;

	amp_S11 = RMS_amp(_S11,nt) + noise;
	amp_S12 = RMS_amp(_S12,nt) + noise;
	amp_S21 = RMS_amp(_S21,nt) + noise;
	amp_S22 = RMS_amp(_S22,nt) + noise;
	if (amp_S11 >= amp_S22) { max = amp_S11; min = amp_S22; }
	else { max = amp_S22; min = amp_S11; }
	num = amp_S12 + amp_S21;
	den = amp_S11 + amp_S22;
	if (mode==0) Q = (1.0-(num/den));
	else Q = (min/max)*(1.0 - (num/den));
	return Q;
}

/* rotate the samples itmin to itmax of the four traces by _theta (degrees) into _S1, _D1, _D2, _S2 (nt = itmax - itmin + 1
   samples), with a cosine taper of taper samples at the ends (when it is more than 1).
   (The program's taper at the lower end is at nt - i, one sample after the trace: it is nt - 1 - i here.) */
static void ROT4TR(const float *_S11, const float *_S12, const float *_S21, const float *_S22,
	float *_S1, float *_D1, float *_D2, float *_S2, float _theta, int itmin, int itmax, int taper)
{
	float rad, cos2, sin2, sincos, temp, env=0.0, f;
	int i, nt;

	rad    = _theta*(PI/180);
	cos2   = pow(cos(rad),2);
	sin2   = pow(sin(rad),2);
	sincos = cos(rad)*sin(rad);
	nt     = itmax - itmin +1;

	for (i=itmin;i <= itmax;i++) {
		temp = sincos*(_S12[i]+ _S21[i]);
		_S1[i-itmin] = _S11[i]*cos2 + _S22[i]*sin2 + temp;
		_S2[i-itmin] = _S11[i]*sin2 + _S22[i]*cos2 - temp;
		temp         = sincos*(_S22[i]- _S11[i]);
		_D1[i-itmin] =  _S12[i]*cos2 - _S21[i]*sin2 + temp;
		_D2[i-itmin] = -_S12[i]*sin2 + _S21[i]*cos2 + temp;
	}

	/* tapering the data */
	if (taper > 1) {
		for (i=0; i<taper && i<nt; i++) {
			f = (float)i/(float)taper;
			env = 0.5*(1.0-cos(PI*f));
			_S1[i] *= env; _S2[i] *= env; _D1[i] *= env; _D2[i] *= env;
			_S1[nt-1-i] *= env; _S2[nt-1-i] *= env; _D1[nt-1-i] *= env; _D2[nt-1-i] *= env;
		}
	}
}

/* quadratic interpolation of three values, the second of which is the biggest: the fractional distance of the maximum
   from the middle one, and its value in max (Schuck, 1993) */
static float QuadFit(const float *y, float *max)
{
	float a, b, c, x;
	c = y[1];
	b = (y[2]-y[0])*0.5;
	a = y[0]+b-c;
	x = -b/(2.0*a);
	*max = c + x*(b+x*a);
	return x;
}

/* the index of the maximum of y[nt], and the value */
static int FindPeak(const float *y, int nt, float *val)
{
	float max = -1E30;
	int i, index=0;
	for (i=0; i <= nt-1; i++) {
		if (y[i]> max) { max = y[i]; index = i; }
	}
	*val= max;
	return index;
}

/* the normalized cross-correlation of x and y (xcor divided by the product of their norms) */
static void nxcor(int nx, int x1, float *x, int ny, int y1, float *y, int nz, int z1, float *z)
{
	float sumX = 0, sumY = 0, N;
	int i;

	xcor(nx, x1, x, ny, y1, y, nz, z1, z);
	for (i=0; i<=nx-1; i++) sumX += x[i]*x[i];
	for (i=0; i<=ny-1; i++) sumY += y[i]*y[i];
	N = sqrt(sumX) * sqrt(sumY);
	for (i=0; i<=nz-1; i++) z[i] = z[i]/N;
}

/* The optimum angle for the Alford rotation of the samples itmin to itmax of four (mixed) traces: the rotation of the
 * window that makes the energy of the diagonal traces most different from the off-diagonal ones (the quality factor), for
 * angles from 0 to 90 degrees in steps of ang_inc (at least 1, 5 if less), refined by a quadratic fit of the quality factors
 * around the maximum; and the time lag between the fast and the slow shear waves (from the cross-correlation of the
 * rotated traces, which are resampled to 4 times as many samples), which decides which of the two is the fast one.
 *
 * S11, S12, S21, S22		the traces
 * dt				sample interval
 * itmin, itmax			the window of analysis
 * ang_inc			angle increment (degrees)
 * taper			taper (samples) of the window for the rotations that are searched
 * quality, inputquality	the quality factor of the best rotation and of no rotation
 * lag				the time lag between the fast and the slow shear waves (the same unit as dt, not negative)
 * xcorr			the maximum of the normalized cross-correlation of the fast and the slow shear waves
 *
 * Returns the angle (degrees, from 0 to 180) of the fast shear wave (the program can give an angle a little under 0).
 *
 * Differences from Find_Angle of the program, which are its bugs: it does not rotate the traces for the window that is the
 * whole trace, in a swapped way that is not related to the angle that it returns (the rotation of the traces is made with
 * su_alford_rotate, by the angle that it returns).
 */
float su_alford_angle(const float *S11, const float *S12, const float *S21, const float *S22, float dt, int itmin,
	int itmax, float ang_inc, int taper, float *quality, float *inputquality, float *lag, float *xcorr)
{
	const int nt = itmax - itmin + 1;
	int i, imax=0, ns_angle, peak_index, ntres, ntemp, zero_lag;
	float qmax, y[3], angle_frac, AMPS1, AMPS2, dtres, angle;

	if (ang_inc < 1.0) ang_inc = 5.0;

	std::vector<float> S1(nt), S2(nt), D1(nt), D2(nt);
	ns_angle = 1 + (int) floor(90/ang_inc);
	std::vector<float> angles(ns_angle), Q(ns_angle);

	/* the table of angles */
	angles[0] = 0;
	for (i=1; i<=ns_angle-1; i++) angles[i] = angles[i-1] + ang_inc;

	/* rotate for each of the angles, and the quality factors */
	for (i=0; i<=ns_angle-1; i++) {
		ROT4TR(S11,S12,S21,S22,S1.data(),D1.data(),D2.data(),S2.data(),angles[i],itmin,itmax,taper);
		Q[i] = Quality_Factor(S1.data(),D1.data(),D2.data(),S2.data(),nt,0);
	}
	*inputquality = Q[0];

	/* the maximum */
	qmax = -1.0E30;
	for (i=0; i<=ns_angle-1; i++) {
		if (Q[i] > qmax) { qmax = Q[i]; imax = i; }
	}

	/* the three values for the quadratic interpolation */
	if (imax!=0 && imax != ns_angle-1) {
		for (i=0; i<=2; i++) y[i] = Q[imax-1+i];
	} else {
		/* a maximum on the first or last of them: use the symmetry of the problem (95 degrees is the same as 5, and
		   -5 the same as 85) */
		y[1] = Q[imax];
		if (imax == 0) {
			y[0] = Q[ns_angle > 1 ? ns_angle-2 : 0];
			y[2] = Q[ns_angle > 1 ? imax+1 : 0];
		} else {
			y[0] = Q[imax-1];
			y[2] = Q[1];
		}
	}
	angle_frac = QuadFit(y,&qmax);
	if (!(angle_frac == angle_frac)) angle_frac = 0.0;   /* (a flat maximum) */
	angle = angles[imax] + ang_inc*angle_frac;
	if (angle < 0.0) angle += 180;   /* (a maximum on the first angle can be a little before it) */

	/* rotate the window by it, and cross-correlate the fast and slow traces (at 4 times as many samples) */
	ROT4TR(S11,S12,S21,S22,S1.data(),D1.data(),D2.data(),S2.data(),angle,itmin,itmax,taper);

	ntres = 4*nt;
	dtres = dt/4;
	std::vector<float> S1res(ntres), S2res(ntres), t(ntres);
	for (i=0; i<ntres; i++) t[i] = 0.25f*i;
	ints8r(nt,1.0,0.0,S1.data(),0.0,0.0,ntres,t.data(),S1res.data());
	ints8r(nt,1.0,0.0,S2.data(),0.0,0.0,ntres,t.data(),S2res.data());
	ntemp = 2*ntres-1;
	std::vector<float> temp(ntemp);
	nxcor(ntres,0,S1res.data(),ntres,0,S2res.data(),ntemp,-ntres+1,temp.data());
	peak_index = FindPeak(temp.data(), ntemp, xcorr);
	zero_lag = ntres-1;
	*lag = dtres*(peak_index - zero_lag);

	/* a peak before the zero lag: S2 is delayed in relation to S1 the other way: S2 is the fast shear wave and the angle
	   is 90 degrees more */
	if (peak_index < zero_lag) {
		angle += 90;
		if (angle >= 180) angle -= 180;
		*lag *= -1;
	}

	/* no lag, but bigger amplitudes in S2: also 90 degrees from the angle that was found */
	AMPS1 = RMS_amp(S1.data(),nt);
	AMPS2 = RMS_amp(S2.data(),nt);
	if (peak_index == zero_lag && AMPS2 > AMPS1) {
		angle += 90;
		if (angle >= 180) angle -= 180;
	}

	*quality = Quality_Factor(S1.data(), D1.data(), D2.data(), S2.data(), nt, 0);
	return angle;
}

/* Rotate four whole traces by an angle (degrees), from the sample itmin on, in place. The part above is not rotated, and
 * is joined to the rotated part by an overlap of 2 * taper + 1 samples, around itmin, with a cosine scale (taper = -1 has no
 * overlap, taper = 0 an average of the one sample).
 *
 * (The program's scale has the integer division i / (overlap - 1), which makes it 1 for all but the last sample of the
 * overlap; it also needs itmin > taper and does not rotate a window that is the whole trace.)
 */
void su_alford_rotate(float *S11, float *S12, float *S21, float *S22, float theta, int nt, int itmin, int taper)
{
	float rad, cos2, sin2, sincos, temp, c1, c2;
	int i, first, overlap;
	std::vector<float> S1(S11, S11+nt), S2(S22, S22+nt), D1(S12, S12+nt), D2(S21, S21+nt);

	rad = theta*(PI/180);
	cos2 = pow(cos(rad),2);
	sin2 = pow(sin(rad),2);
	sincos = cos(rad)*sin(rad);

	/* the overlap and the part below are rotated */
	overlap = taper < 0 ? 0 : 2*taper+1;
	first = itmin - (taper < 0 ? 0 : taper);
	if (first < 0) first = 0;
	for (i=first; i<=nt-1; i++) {
		temp = sincos*(S12[i] + S21[i]);
		S1[i] = S11[i]*cos2 + S22[i]*sin2 + temp;
		S2[i] = S11[i]*sin2 + S22[i]*cos2 - temp;
		temp = sincos*(S22[i] - S11[i]);
		D1[i] =  S12[i]*cos2 - S21[i]*sin2 + temp;
		D2[i] = -S12[i]*sin2 + S21[i]*cos2 + temp;
	}

	/* join the unrotated part to the rotated part in the overlap */
	if (overlap == 1) {
		i = itmin;
		if (i >= 0 && i < nt) {
			S1[i] = 0.5*(S1[i] + S11[i]);
			D1[i] = 0.5*(D1[i] + S12[i]);
			D2[i] = 0.5*(D2[i] + S21[i]);
			S2[i] = 0.5*(S2[i] + S22[i]);
		}
	} else if (overlap > 1) {
		for (i=0; i<overlap; i++) {
			int k = itmin - taper + i;
			if (k < 0 || k >= nt) continue;
			c1 = cos((PI/2)*((float) i/(float) (overlap-1)));
			c2 = 1-c1;
			S1[k] = c1*S11[k] + c2*S1[k];
			D1[k] = c1*S12[k] + c2*D1[k];
			D2[k] = c1*S21[k] + c2*D2[k];
			S2[k] = c1*S22[k] + c2*S2[k];
		}
	}
	memcpy(S11, S1.data(), nt*sizeof(float));
	memcpy(S12, D1.data(), nt*sizeof(float));
	memcpy(S21, D2.data(), nt*sizeof(float));
	memcpy(S22, S2.data(), nt*sizeof(float));
}
