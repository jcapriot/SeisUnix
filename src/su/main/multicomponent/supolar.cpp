/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUPOLAR: $Revision: 1.9 $ ; $Date: 2015/08/07 22:21:04 $	*/

#include "su.h"
#include <vector>

/* window shape identifiers */
#define WBOXCAR 0
#define WBARTLETT 1
#define WHANNING 2
#define WWELSH 3
#include "segy.h"
#include "header.h"

/*********************** self documentation *****************************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                       ",
" SUPOLAR - POLarization analysis of three-component data               ",
"                                                                       ",
" supolar <stdin [optional parameters]                                  ",
"                                                                       ",
" Required parameters:                                                  ",
"    none                                                               ",
"                                                                       ",
" Optional parameters:                                                  ",
"    dt=(from header)  time sampling intervall in seconds               ",
"    wl=0.1            correlation window length in seconds             ",
"    win=boxcar        correlation window shape, choose \"boxcar\",     ",
"                      \"hanning\", \"bartlett\", or \"welsh\"          ",
"    file=polar        base of output file name(s)                      ",
"    rl=1              1 = rectilinearity evaluating 2 eigenvalues,     ",
"                      2, 3 = rectilinearity evaluating 3 eigenvalues   ",
"    rlq=1.0           contrast parameter for rectilinearity            ",
"    dir=1             1 = 3 components of direction of polarization    ",
"                          (the only three-component output file)       ",
"    tau=0             1 = global polarization parameter                ",
"    ellip=0           1 = principal, subprincipal, and transverse      ",
"                          ellipticities e21, e31, and e32              ",
"    pln=0             1 = planarity measure                            ",
"    f1=0              1 = flatness or oblateness coefficient           ",
"    l1=0              1 = linearity coefficient                        ",
"    amp=0             1 = amplitude parameters: instantaneous,         ",
"                          quadratic, and eigenresultant ir, qr, and er ",
"    theta=0           1, 2, 3 = incidence angle of principal axis      ",
"    phi=0             1, 2, 3 = horizontal azimuth of principal axis   ",
"    angle=rad         unit of angles theta and phi, choose \"rad\",    ",
"                      \"deg\", or \"gon\"                              ",
"    all=0             1, 2, 3 = set all output flags to that value     ",
"    verbose=0         1 = echo additional information                  ",
"                                                                       ",
"                                                                       ",
" Notes:                                                                ",
"    Three adjacent traces are considered as one three-component        ",
"    dataset.                                                           ",
"    Correct calculation of angles theta and phi requires the first of  ",
"    these traces to be the vertical component, followed by the two     ",
"    horizontal components (e.g. Z, N, E, or Z, inline, crossline).     ",
"    Significant signal energy on Z is necessary to resolve the 180 deg ",
"    ambiguity of phi (options phi=2,3 only).                           ",
"                                                                       ",
"    Each calculated polarization attribute is written into its own     ",
"    SU file. These files get the same base name (set with \"file=\")   ",
"    and the parameter flag as an extension (e.g. polar.rl).            ",
"                                                                       ",
"    In case of a tapered correlation window, the window length wl may  ",
"    have to be increased compared to the boxcar case, because of their ",
"    smaller effective widths (Bartlett, Hanning: 1/2, Welsh: 1/3).     ",
"                                                                       ",
" Range of values:                                                      ",
"    parameter     option  interval                                     ",
"    rl            1, 2    0.0 ... 1.0   (1.0: linear polarization)     ",
"    rl            3      -1.0 ... 1.0                                  ",
"    tau, l1       1       0.0 ... 1.0   (1.0: linear polarization)     ",
"    pln, f1       1       0.0 ... 1.0   (1.0: planar polarization)     ",
"    e21, e31, e32 1       0.0 ... 1.0   (0.0: linear polarization)     ",
"    theta         1      -pi/2... pi/2  rad                            ",
"    theta         2, 3    0.0 ... pi/2  rad                            ",
"    phi           1      -pi/2... pi/2  rad                            ",
"    phi           2      -pi  ... pi    rad   (see notes above)        ",
"    phi           3       0.0 ... 2 pi  rad   (see notes above)        ",
"                                                                       ",
"                                                                       ",
NULL};
#endif

/* 
 * Author: Nils Maercklin, 
 *         GeoForschungsZentrum (GFZ) Potsdam, Germany, 1998-2001.
 *         E-mail: nils@gfz-potsdam.de
 * 
 *
 * References:
 *    Jurkevics, A., 1988: Polarization analysis of three-component
 *       array data. Bulletin of the Seismological Society of America, 
 *       vol. 78, no. 5.
 *    Kanasewich, E. R., 1981: Time Sequence Analysis in Geophysics.
 *       The University of Alberta Press.
 *    Kanasewich, E. R., 1990: Seismic Noise Attenuation.
 *       Handbook of Geophysical Exploration, Pergamon Press, Oxford.
 *    Meyer, J. H. 1988: First Comparative Results of Integral and
 *       Instantaneous Polarization Attributes for Multicomponent Seismic
 *       Data. Institut Francais du Petrole.
 *    Press, W. H., Teukolsky, S. A., Vetterling, W. T., and Flannery, B. P.
 *       1996: Numerical Recipes in C - The Art of Scientific Computing.
 *       Cambridge University Press, Cambridge.
 *    Samson, J. C., 1973: Description of the Polarisation States of Vector
 *       Processes: Application to ULF Electromagnetic Fields.
 *       Geophysical Journal vol. 34, p. 403-419.
 *    Sheriff, R. E., 1991: Encyclopedic Dictionary of Exploration
 *       Geophysics. 3rd ed., Society of Exploration Geophysicists, Tulsa.
 *
 * Trace header fields accessed: ns, dt
 * Trace header fields modified: none
 */
/**************** end self doc *******************************************/

/* calculate time window weights for a smooth time window,        */
/* after Sheriff, 1991, p. 335; tapered windows are suggested by  */
/* Jurkevics, 1988                                                */

static void calc_window(float *w, int iwl, int iwin)
{
    int i,j;      /* loop indices within window */
    float m;      /* half of time window (iwl/2) */
    
    m = (float) iwl/2;
    
    switch (iwin) {
        case WBOXCAR :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 1.0;
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;
        case WBARTLETT :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 1.0 - (fabs((float) i) / m);
                
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;
        case WHANNING :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = 0.5 + 0.5*cos(PI*fabs((float) i) / m);
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;            
        case WWELSH :
            for (i=-iwl/2,j=0;i<iwl/2;i++,j++) {
                w[j] = fabs( fabs((float) i) / m - 1.0);
                w[j] *= w[j];            
            }
            if (j<iwl) w[iwl-1] = w[0];
            break;            
    }
}


/* covariance of two time sequences *data1 and *data2,   */
/* calculated in a window of length iwl,                 */
/* weighted with factors *w and starting at index istart */

static float covar(float *data1, float *data2, int istart, int iwl, float *w)
{
    int i, j;
    float cov=0.0;
    float mean1=0.0;
    float mean2=0.0;
    
    for (i=istart,j=0;i<(istart+iwl);i++,j++) {
        mean1 += data1[i] * w[j];
        mean2 += data2[i] * w[j];
    }
    mean1 = mean1 / (float) iwl;
    mean2 = mean2 / (float) iwl;

    for (i=istart,j=0;i<(istart+iwl);i++,j++) {
	/* modified by Mega G. Baitoff 6 May 2015 */
        /* cov += (data1[i]-mean1) * (data2[i]-mean2) * w[j]; */
        cov += (data1[i] * w[j] - mean1) * (data2[i] * w[j] - mean2);              
    }
    cov = cov / (float) iwl;
    return cov;
}   


/**********************************************************************/
/* Functions defining polarization attributes and angles              */
/**********************************************************************/

/* rectilinearity factor rl (3 different definitions) */

static float calc_rl(float *d, float rlq, int opt)
{   
    float rl;
    
    if (!(d[1]==0.0)) {
        switch (opt) {
            default: 
                /* case 1: */
                /* rl definition after Kanasewich, 1981 */
                rl = 1 - pow(fabs(d[2]/d[1]), rlq);
                break;
            case 2:
                /* rl definition after Jurkevics, 1988 */
                rl = 1 - pow( 0.5*(fabs(d[2]/d[1]) + fabs(d[3]/d[1]) ), rlq);
                break;
            case 3:
                /* rl definition after Meyer, 1988 */
                rl = 1 - pow(fabs(d[2]/d[1]) + fabs(d[3]/d[1]), rlq);
                break;
        }
        return rl;
    }
    else {
        return 0.0;
    }
}

/* vertical (incidence) angle theta (different definitions) */
/* assumes a right handed coordinate system starting        */
/* with the vertical component Z (Z=1, N,R=2, E,T=3)        */

static float calc_theta(float **v, int opt)
{
    float theta, horiz;
    
    switch (opt) {
        default:
            /* case 1: */
            /* definition after Kanasewich, 1981 */
            /* interval -pi/2 <= theta <= pi/2 */ 
            if (v[1][1]) {
                horiz = sqrt( v[2][1]*v[2][1] + v[3][1]*v[3][1] );
                theta = atan( horiz / v[1][1] );
            }
            else {
                theta = 0.0;
            }
            break;
        case 2:
        case 3:
            /* definition after Jurkevics, 1988 */
            /* interval 0.0 <= theta <= pi/2 */
            theta = acos( MIN(1.0, fabs(v[1][1])) );   /* (an eigenvector can be 1 + epsilon long: NaN for the program) */
            break;
    }
    return theta;            
}

/* horizontal azimuth phi (different definitions)    */
/* assumes a right handed coordinate system starting */
/* with the vertical component Z (Z=1, N,R=2, E,T=3) */

/* Note: The SIGN function is introduced to resolve  */
/* the 180 deg ambiguity by taking the positive      */
/* vertical component of v[1][1] (Jurkevics, 1988).  */
#define VSIGN ( (v[1][1]<0.0) ? -1.0 : 1.0 )

static float calc_phi(float **v, int opt)
{
    float phi;
    
    switch (opt) {
        default:
            /* case 1: */
            /* definition after Kanasewich, 1981 */
            /* interval -pi/2 <= phi <= pi/2 */ 
            if (v[2][1]) {
                phi = atan( v[3][1] / v[2][1] );
            }
            else {
                phi = (v[3][1]>0.0) ? 0.5*PI : -0.5*PI;
            }
            break;
        case 2:
        case 3:
            /* definitions after Jurkevics, 1988 */
            /* interval -pi <= phi <= pi */
            if (v[2][1]) {
                phi = atan2( v[3][1]*VSIGN, v[2][1]*VSIGN);
            }
            else {
                phi = (v[3][1]>0.0) ? 0.5*PI*VSIGN : -0.5*PI*VSIGN;
            }
            
            /* interval 0.0 <= phi <= 2*pi */
            if (phi<0.0 && opt==3) phi += 2.0*PI;
            break;
    }   
    return phi;
}
#undef VSIGN

/* global polarization parameter tau (Samson, 1973) */

static float calc_tau(float *d)
{
    float x1, x2, x3, x4, tau;
    
    if (d[1]) {   
        x1 = pow(( 1 - fabs(d[2]/d[1])), 2.);
        x2 = pow(( 1 - fabs(d[3]/d[1])), 2.);
        x3 = pow(( fabs(d[2]/d[1]) - fabs(d[3]/d[1]) ), 2.);
        x4 = pow(( 1 + fabs(d[2]/d[1]) + fabs(d[3]/d[1])), 2.);
        tau = sqrt( (x1+x2+x3) / (2*x4) );
        return tau;
    }
    else
        return 0.0;
}

/* ellipticities e_ik */

static float calc_ellip(float *d, int d1, int d2)
{
    float ellip;
    
    if (d[d2]) {
        ellip = sqrt( fabs(d[d1]/d[d2]) );
        return ellip;
    }
    else
        return 0.0;
}

/* planarity after Jurkevics, 1988 */

static float calc_plan(float *d)
{
    float pln;
    
    if (d[1]+d[2]) {
        pln = 1.0 - 2.0*d[3] / (d[1] + d[2]);
        return pln;
    }
    else
        return 0.0;
}

/* flatness coefficient f1 after Benhama et. al, 1988 */

static float calc_f1(float *d)
{
    float f1,x1,x2;
    
    x1 = 3.0 * calc_ellip(d,3,1);
    x2 = 1.0 + calc_ellip(d,2,1) + calc_ellip(d,3,1);
    f1 = 1.0 - x1 / x2; 
    return f1;
}

/* linearity coefficient l1 */

static float calc_l1(float *d)
{
    float l1,x1,x2;
    
    x1 = 3. * ( calc_ellip(d,2,1) + calc_ellip(d,3,1) );
    x2 = 2. * ( 1 + calc_ellip(d,2,1) + calc_ellip(d,3,1));
    l1 = 1. - x1 / x2;
    return l1;
}

/* direction cosines or directivity functions (3 components) */

static void calc_dir(float **data3c_dir, float **v, int it)
{
    int i;
    for (i=1;i<=3;i++) {
        data3c_dir[i][it]=v[i][1];
    }
}

/* amplitude parameters */

/* eigenresultant */

static float calc_er(float *d)
{
    float er;
    er = sqrt(fabs(d[1]));
    return er;
}

/* instantaneous and quadratic resultant (Meyer, 1988) */

static void ampparams(float **indata, float *data_ir, float *data_qr, int nt, int iwl)
{
    int i,it;
    float sqrsum;
    
    for (it=0;it<nt;it++) {
        /* instantaneous resultant */
        data_ir[it] = sqrt(indata[1][it]*indata[1][it]+ \
            indata[2][it]*indata[2][it]+indata[3][it]*indata[3][it]);
        
        /* quadratic resultant */
        if ( (it >= iwl/2) && (it < nt-iwl/2) ) {
            sqrsum=0.0;
            for (i=it-iwl/2;i<it-iwl/2+iwl;i++) {
                sqrsum += indata[1][i]*indata[1][i]+ \
                    indata[2][i]*indata[2][i]+indata[3][i]*indata[3][i];
            }
            data_qr[it] = sqrsum / (float) iwl;
        }
    }
}




/* The analysis of one three-component station of the program (the three traces are the vertical component, followed by
 * the two horizontal ones), for the correlation window of iwl samples and its shape iwin.
 *
 * data[3][nt]		the three components
 * rl, theta, phi, tau, ellip, pln, f1, l1, amp, dir	the flags of the attributes (0 = not): rl 1, 2 or 3 for the
 *			definitions of the rectilinearity (rlq is the contrast), theta and phi 1, 2 or 3
 * fangle		unit factor for theta and phi (1 radians, 180 / pi degrees, 200 / pi gon)
 * out[16][nt]		the attributes, zero where they are not computed (the caller zeros them, and the first and last iwl / 2
 *			samples of each are not):
 *			0 rl, 1 theta, 2 phi, 3 tau, 4 e21, 5 e31, 6 e32, 7 pln, 8 f1, 9 l1, 10 er, 11 ir, 12 qr,
 *			13, 14, 15 the components of the direction of polarization
 *
 * Differences from the program, which are its bugs: the horizontal azimuth phi of option 1 tests the component v[1][2]
 * (the vertical component of the second eigenvector) where it means v[2][1] (the north component of the first), and
 * the quadratic resultant sums over the window that goes from it - iwl / 2 for 1.5 iwl samples (reading past the end
 * of the trace for the last samples), where the window is iwl samples.
 */
void su_polar_station(int nt, int iwl, int iwin, float rlq, int rl, int theta, int phi, int tau, int ellip, int pln,
	int f1, int l1, int amp, int dir, float fangle, const float *data, float *out)
{
	int i, j, it;
	/* (the program's arrays are numbered from 1) */
	std::vector<float*> data3c(4), a(4), v(4), data3c_dir(4);
	std::vector<float> abuf(9), vbuf(9), d(4), w(iwl > 0 ? iwl : 1);
	std::vector<float> copy((size_t) 3*nt);

	for (i=1;i<=3;i++) {
		data3c[i] = &copy[(size_t) (i-1)*nt];
		memcpy(data3c[i], data+(size_t) (i-1)*nt, nt*sizeof(float));
		a[i] = &abuf[(i-1)*3] - 1;
		v[i] = &vbuf[(i-1)*3] - 1;
		data3c_dir[i] = out + (size_t) (12+i)*nt;
	}
	calc_window(w.data(), iwl, iwin);

	for (it=iwl/2;it<nt-iwl/2;it++) {
		/* covariance matrix */
		for (i=1;i<=3;i++) {
			for (j=i;j<=3;j++) {
				a[i][j]=a[j][i]=covar(data3c[i], data3c[j], it-iwl/2, iwl, w.data());
			}
		}
		/* compute eigenvalues and -vectors */
		eig_jacobi(a.data(),d.data(),v.data(),3);
		sort_eigenvalues(d.data(),v.data(),3);

		/* polarization parameters */
		if (rl) out[(size_t) 0*nt+it]=calc_rl(d.data(),rlq,rl);
		if (theta) out[(size_t) 1*nt+it]=calc_theta(v.data(), theta) * fangle;
		if (phi) out[(size_t) 2*nt+it]=calc_phi(v.data(), phi) * fangle;
		if (tau) out[(size_t) 3*nt+it]=calc_tau(d.data());
		if (ellip) {
			out[(size_t) 4*nt+it]=calc_ellip(d.data(),2,1);
			out[(size_t) 5*nt+it]=calc_ellip(d.data(),3,1);
			out[(size_t) 6*nt+it]=calc_ellip(d.data(),3,2);
		}
		if (pln) out[(size_t) 7*nt+it]=calc_plan(d.data());
		if (f1) out[(size_t) 8*nt+it]=calc_f1(d.data());
		if (l1) out[(size_t) 9*nt+it]=calc_l1(d.data());
		if (amp) out[(size_t) 10*nt+it]=calc_er(d.data());
		if (dir) calc_dir(data3c_dir.data(),v.data(),it);
	}

	/* compute amplitude parameters */
	if (amp) {
		std::vector<float*> rows(4);
		rows[1] = data3c[1]; rows[2] = data3c[2]; rows[3] = data3c[3];
		ampparams(rows.data(), out + (size_t) 11*nt, out + (size_t) 12*nt, nt, iwl);
	}
}
