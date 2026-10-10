/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                */

/* SUREMEL2DAN: $Revision: 1.4 $ ; $Date: 2015/08/11 22:52:55 $         */


/* 
  Elastic anisotropic 2D Fourier method modeling with REM time integration
*/

#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                        ",
" SUREMEL2DAN - Elastic anisotropic 2D Fourier method modeling with high ",
"               accuracy Rapid Expansion Method (REM) time integration   ",
"                                                                        ",
" suremel2dan [parameters]                                               ",
"                                                                        ",
" Required parameters:                                                   ",
"                                                                        ",
" nx=         number of grid points in horizontal direction              ",
" nz=         number of grid points in vertical direction                ",
" nt=         number of time samples                                     ",
" dx=         spatial increment in horizontal direction                  ",
" dz=         spatial increment in vertical direction                    ",
" dt=         time sample interval in seconds                            ",
" isx=        grid point # of horizontal source positions                ",
" isz=        grid point # of vertical source positions                  ",
" styp=       source types (pressure, shear, single forces)              ",
" samp=       amplitudes of sources                                      ",
" amode=      0: isotropic,  1: anisotropic                              ",
" vmax=       global maximum velocity (only if amode=1)                  ",
" vmin=       global minimum velocity (only if amode=1)                  ",
"                                                                        ",
" Optional parameters:                                                   ",
" fx=0.0      first horizontal coordinate                                ",
" fz=0.0      first vertical coordinate                                  ",
" irx=        horizontal grid point # of vertical receiver lines         ",
" irz=        vertical grid point # of horizontal receiver lines         ",
" rxtyp=      types of horizontal receiver lines                         ",
" rztyp=      types of vertical receivers lines                          ",
" sntyp=      types of snapshots                                         ",
"             0: P,  1: S,  2: UX,  3: UZ                                ",
" w=0.1       width of spatial source distribution (see notes)           ",
" sflag=2     source time function                                       ",
"             0: user supplied source function                           ",
"             1: impulse (spike at t=0)                                  ",
"             2: Ricker wavelet                                          ",
" fmax=       maximum frequency of Ricker (default) wavelet              ",
" amps=1.0    amplitudes of sources                                      ",
" prec=0      1: precompute Bessel coefficients b_k (see notes)          ",
"             2: use precomputed Bessel coefficients b_k                 ",
" vmaxu=      user-defined maximum velocity                              ",
" dtsnap=0.0  time interval in seconds of wave field snapshots           ",
" iabso=1     apply absorbing boundary conditions (0: none)              ",
" abso=0.1    damping parameter for absorbing boundaries                 ",
" nbwx=20     horizontal width of absorbing boundary                     ",
" nbwz=20     vertical width of absorbing boundary                       ",
" verbose=0   1: show parameters used                                    ",
"             2: print maximum amplitude at every expansion term         ",
"                                                                        ",
" c11file=c11       c11 filename                                         ",
" c13file=c13       c13 filename                                         ",
" c15file=c15       c15 filename                                         ",
" c33file=c33       c33 filename                                         ",
" c35file=c35       c35 filename                                         ",
" c55file=c55       c55 filename                                         ",
" vpfile=vp         P-velocity filename                                  ",
" vsfile=vs         S-velocity filename                                  ",
" densfile=dens     density filename                                     ",
"                                                                        ",
" sname=wavelet.su  user supplied source time function filename          ",
"                                                                        ",
" Basenames of seismogram and snapshot files:                            ",
" xsect=xsect_     x-direction section files basename                    ",
" zsect=zsect_     z-direction section files basename                    ",
" snap=snap_       snapshot files basename                               ",
"                                                                        ",
" jpfile=stderr        diagnostic output                                 ",
"                                                                        ",
" Notes:                                                                 ",
"  0. The combination of the Fourier method with REM time integration    ",
"     allows the computation of synthetic seismograms which are free     ",
"     of numerical grid dispersion. REM has no restriction on the        ",
"     time step size dt. The Fourier method requires at least two        ",
"     grid points per shortest wavelength.                               ",
"  1. nx and nz must be valid numbers for pfafft transform lengths.      ",
"     nx and nz must be odd numbers. For valid numbers see e.g.          ",
"     numbers in structure 'nctab' in source file                        ",
"     $CWPROOT/src/cwp/lib/pfafft.c.                                     ",
"  2. Velocities and densities are stored as plain C style files         ",
"     of floats where the fast dimension is along the z-direction.       ",
"  3. Units must be consistent, e.g. m, s and m/s.                       ",
"  4. A 20 grid points wide border at the sides and the bottom of        ",
"     the modeling grid is used for sponge boundary conditions           ",
"     (default: iabso=1).                                                ",
"     Source and receiver lines should be placed some (e.g. 10) grid     ",
"     points away from the absorbing boundaries in order to reduce       ",
"     reflections due to obliquely incident wavefronts.                  ",
"  5. Dominant frequency is about fmax/2 (sflag=2), absolute maximum     ",
"     is delayed by 3/fmax from beginning of wavelet.                    ",
"  6. If source is not single force (i.e. pressure or shear source)      ",
"     it should be not a spike in space; the parameter w determines      ",
"     at which distance (in grid points) from the source's center        ",
"     the Gaussian weight decays to 10 percent of its maximum.           ",
"     w=2 may be a reasonable choice; however, the waveform will be      ",
"     distorted.                                                         ",
"  7. Horizontal and vertical receiver line sections are written to      ",
"     separate files. Each file can hold more than one line.             ",
"  8. Parameter vmaxu may need to be chosen larger than the highest      ",
"     propagation velocity if the modeling run becomes unstable.         ",
"     This happens if the largest eigenvalue of the modeling             ",
"     operator L is larger than estimated from the largest velocity      ",
"     due to variations of the density.                                  ",
"  9. Bessel coefficients can be precomputed (prec=1) and stored on      ",
"     disk to save CPU time when several shots need to be run.           ",
"     In this case computation of Bessel coefficients can be skipped     ",
"     and read from disk file for reuse (prec=2).                        ",
"     For reuse of Bessel coefficients the user may need to define       ",
"     the overall maximum velocity (vmaxu).                              ",
" 10. If snapshots are not required, a spike source (sflag=1) may be     ",
"     applied and the resulting impulse response seismograms can be      ",
"     convolved later with a desired wavelet.                            ",
" 11. Output is written to SU style files.                               ", 
"     Basenames of seismogram and snapshot output files will be          ",
"     extended by the type of the data (p, s, ux, or uz).                ",
"     Additionally seismogram files will be consecutively numbered.      ",
"                                                                        ",
" Caveat:                                                                ",
"     Time sections and snapshots are kept entirely in memory during     ",
"     run time. Therefore, lots of time section and snapshots may        ",
"     eat up a large amount of memory.                                   ",
NULL};
#endif


/*
  Elastic anisotropic 2D Fourier method modeling with REM time integration

  Reference: 
  Kosloff, D., Fihlo A.Q., Tessmer, E. and Behle, A., 1989,
    Numerical solution of the acoustic and elastic wave equations by a
    new rapid expansion method, Geophysical Prospecting, 37, 383-394
  
 * Credits:
 *      University of Hamburg: Ekkehart Tessmer, July 2013
 */
/**************** end self doc ***********************************/

typedef struct {
	float *ax, *rkx1, *az, *rkz1;	/* the work arrays and wave numbers of difx and difz */
	int *ist;			/* the sponge boundary table of tbc2d */
} Remel2dCtx;


static void rk(float *ak, int n, float d, int ind);
static void bessel_jn(double x, int n, double *bes);
static void legcoeffs(int n, double *x, double *a);
static void null(int n, double *x0, double *p, double *pm, double eps, int *mxit);
static void legpol(int n, double x, double *p, double *pm, double *ps);

static void g2d(float w, int nx, int nz, float **g)
/*
  2-D Gaussian weights for spatial distribution of source 

  input:
    w:     width in grid points where Gaussian is one tenth of its maximum
    nx:    number of gridpoints in horizontal direction
    nz:    number of gridpoints in vertical direction

  output:
    g:     weights for spatial distribution of source
*/
{
  int i, i0, k, k0;
  float pi, sigma, fac, sum, r2, arg;

  pi = 4.*atan(1.);

  i0 = nx/2;
  k0 = nz/2;

  if (w < 0.1) w = 0.1; /* spatial spike */

  /* w is width where Gaussian is one tenth of its maximum */
  sigma = sqrt(-w*w/(2.*log(0.1)));

  fac = 1./(sigma*sqrt(2.*pi));

  sum = 0.;
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      r2 = (i-i0)*(i-i0)+(k-k0)*(k-k0);
      arg = r2/(2.*sigma*sigma);
      g[k][i] = 0.;
      if (arg < 80.) g[k][i] = fac*exp(-arg);
      sum += g[k][i];
    }
  }
  /* normalization */
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) g[k][i] /= sum;
  }
}
/**********************************************************************/
static void difx(Remel2dCtx *ctx, float **vin, float **vout, float **bulk, int nx, int nz, float dx, 
	  int iadd)
/* 
   compute 1st derivative along x-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dx:     grid spacing in x-direction
     iadd:   flag

   output:
     vout:   2D output array

     iadd=0: set result in vout
     iadd=1: add result to vout
     iadd=2: set result multiplied by bulk in vout
     iadd=3: add result multiplied by bulk to vout
*/
   
{
  int i, k;
  float  tmp;
  float *a, *rkx;
  int ieo, m1=-1, p1=1;
  int num;

  ieo = nz % 2;
  num = ((nz + 1) / 2) * 2;
  if (ctx->ax == NULL) ctx->ax = ealloc1float(2*nx);
  if (ctx->rkx1 == NULL) {
    ctx->rkx1 = ealloc1float(nx);
    rk(ctx->rkx1,nx,dx,1);
  }
  a = ctx->ax;
  rkx = ctx->rkx1;

  for (k=0; k<num; k+=2) {
	
    /* Load data from vin */
    for (i=0; i<nx; i++) {
      a[2*i] = vin[k][i];
      if (!ieo || !(k==num-2))
	a[2*i+1] = vin[k+1][i];
      else
	a[2*i+1] = 0.;
    }

    /* Perform derivative */
    pfacc(m1, nx, (complex*)a);
	
    for (i=0; i<nx; i++) {
      tmp      =  a[2*i]   * rkx[i];
      a[2*i]   = -a[2*i+1] * rkx[i];
      a[2*i+1] =  tmp;
    }
	
    pfacc(p1, nx, (complex*)a);

    if (iadd == 0) {

      /* Set data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i];
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] = a[2*i] * bulk[k][i];
	if (!ieo || !(k==num-2)) vout[k+1][i] = a[2*i+1] * bulk[k+1][i];
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (i=0; i<nx; i++) {
	vout[k][i] += a[2*i] * bulk[k][i];;
	if (!ieo || !(k==num-2)) vout[k+1][i] += a[2*i+1] * bulk[k+1][i];
      }
    }
  }
}
/**********************************************************************/
static void difz(Remel2dCtx *ctx, float **vin, float **vout, float **bulk, int nx, int nz, int nzpad, float dz,
	  int iadd, int iload)
/* 
   compute 1st derivative along z-direction 

   input:
     vin:    2D input array
     bulk:   2D array for multiplication with output array
     nx:     length of array in x-direction
     nz:     length of array in z-direction
     dz:     grid spacing in z-direction
     iadd:   flag

   output:
     vout:   2D output array

   iadd=0: set result in vout
   iadd=1: add result to vout
   iadd=2: set result multiplied by bulk in vout
   iadd=3: add result multiplied by bulk to vout

   iload=0: load nz rows and zero-pad, unload nzpad rows
   iload=1: load nzpad rows, unload nz rows
*/

{
  int i, k, nzunl;
  float  tmp;
  float *a, *rkz;
  int ieo, m1=-1, p1=1;
  int num;

  ieo = nx % 2;
  num = ((nx + 1) / 2) * 2;
  if (ctx->az == NULL) ctx->az = ealloc1float(2*nzpad);
  if (ctx->rkz1 == NULL) {
    ctx->rkz1 = ealloc1float(nzpad);
    rk(ctx->rkz1,nzpad,dz,1);
  }
  a = ctx->az;
  rkz = ctx->rkz1;

  for (i=0; i<num; i+=2) {

    if (iload == 0) {
      /* Load data from vin */
      for (k=0; k<nz; k++) {
	a[2*k] = vin[k][i];
	if (!ieo || !(i==num-2))
	  a[2*k+1] = vin[k][i+1];
	else
	  a[2*k+1] = 0.;
      }
      /* zero padding for free surface */
      for (k=nz; k<nzpad; k++) {
	a[2*k]   = 0.;
	a[2*k+1] = 0.;
      }      
      nzunl = nzpad;
    } else {
      for (k=0; k<nzpad; k++) {
	a[2*k] = vin[k][i];
	if (!ieo || !(i==num-2))
	  a[2*k+1] = vin[k][i+1];
	else
	  a[2*k+1] = 0.;
      }
      nzunl = nz;
    }

    /* Perform derivative */
    pfacc(m1, nzpad, (complex*)a);
	
    for (k=0; k<nzpad; k++) {
      tmp      =  a[2*k]   * rkz[k];
      a[2*k]   = -a[2*k+1] * rkz[k];
      a[2*k+1] =  tmp;
    }
	
    pfacc(p1, nzpad, (complex*)a);

    if (iadd == 0) {
	
      /* Set data into vout */
      for (k=0; k<nzunl; k++) {
	vout[k][i] = a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1];
      }
    } else if (iadd == 1) {

      /* Add data into vout */
      for (k=0; k<nzunl; k++) {
	vout[k][i] += a[2*k];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1];
      }
    } else if (iadd == 2) {

      /* Set data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] = a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[k][i+1];
      }
      if (nzunl == nzpad) {
	for (k=nz; k<nzpad; k++) {
	  vout[k][i] = a[2*k] * bulk[0][i];
	  if (!ieo || !(i==num-2)) vout[k][i+1] = a[2*k+1] * bulk[0][i+1];
	}
      }
    } else if (iadd == 3) {

      /* Add data multiplied by bulk into vout */
      for (k=0; k<nz; k++) {
	vout[k][i] += a[2*k] * bulk[k][i];
	if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[k][i+1];
      }
      if (nzunl == nzpad) {
	for (k=nz; k<nzpad; k++) {
	  vout[k][i] += a[2*k] * bulk[0][i];
	  if (!ieo || !(i==num-2)) vout[k][i+1] += a[2*k+1] * bulk[0][i+1];
	}
      }
    }
  }
}
 
/**********************************************************************/
static void rk(float *ak, int n, float d, int ind)
/*
  wave-numbers for Fourier transform

  input:
    n:      length of Fourier transform
    d:      grid spacing
    ind:    order of derivative

  output:
    ak:     array of wave numbers
*/
{
  int i,isign,n2;
  float pi, dn, c; 
  
  pi = 4.*atan(1.);
  isign = (int)(-pow((-1.),(float)ind));
  n2 = n/2;
  dn = 2.*pi/(n*d);
  
  for (i=0; i<n; i++) {
    c=dn*i;
    if (i > n2) c = -dn*(n-i);
    ak[i] = pow(c,(float)ind)*isign/n;
  }
}
/**********************************************************************/
static void bessel_jn(double x, int n, double *bes)
/*
  compute Bessel function values J0(x), J1(x), ...

  input:
    x:    argument
    n:    highest order of Bessel function

  output:
    bes:  function values of n orders of Bessel functions at argument x

  Notes:  1.  n must be larger than x
          2.  for accuracy of higher order Bessel function  choose n 
              larger (e.g. +50) than needed
*/
{
  int i, j;
  double sum;
  double eps=1.0e-130;

  if (x > (double)n)
    err("bessel_jn: n(%d) must be larger than x(%e)!\n",n,x);

  if (fabs(x) > eps) {

    bes[n-1] = 0.;
    bes[n-2] = 1.;
    
    sum = 0.;
    for (i=n-3; i>=0; i--) {
      bes[i] = (2*(i+1)/x)*bes[i+1] - bes[i+2];
      
      if (i%2 == 0) {
	sum = sum + 2. * bes[i];
	if (fabs(sum) > 1.e30) {
	  for (j=i; j<n; j++) bes[j] /= sum;
	  sum = 1.;
	}
      }
    }
    
    sum = bes[0];
    for (i=2; i<n; i+=2) sum += 2. * bes[i];
    for (i=0; i<n; i++) bes[i] /= sum;
  } else { /* if fabs(x) very small assume x=0. */
    bes[0] = 1.;
    for (i=1; i<n; i++) bes[i] = 0.;
  }
}

/*
  Routines for numerical convolution of Bessel functions with wavelet
  using Gauss-Legendre quadrature
*/

static void null(int n, double *x0, double *p, double *pm, double eps, 
	  void (*func)(int, double, double*, double*, double*));
static void legpol(int n, double x, double *p, double *pm, double *ps);
static void legcoeffs(int n, double *x, double *a);
static void bessel_jn(double x, int n, double *bes);

/*----------------------------------------------------------------------*/
static void legcoeffs(int n, double *x, double *a)
/*
  Abscissas and corresponding weights for Gauss-Legendre quadrature
  (abscissas are zeros of Legendre polynomial)

  only positive zeros need to be estimated since they come in +/- pairs;
  they are arranged in descending order;

  Input:
    n  polynomial degree

  Output:
    x  abscissas
    a  weights
*/
{
  int i, n2;
  double eps, con, p, pm, ps;

  eps = DBL_EPSILON;

  n2 = n/2;

  /* approximate location of zeros */
  con = 1.0 - 1.0/(8.0*n*n) + 1.0/(8.0*n*n*n);
  
  for (i=0; i<n2; i++)
    x[i] = con * cos(PI*(4.0*(i+1)-1.0)/(4.0*n+2.));
  
  /* accurate location of zeros using Newton-Raphson */
  for (i=0; i<n2; i++) {
    null(n, &x[i], &p, &pm, eps, legpol);
    a[i] = 2.0 * (1.0-x[i]*x[i])/(n*n*pm*pm);
  }
  
  if (n % 2) { /* if degree is odd */
    i = n2;
    x[i] = 0.;
    legpol(n, x[i], &p, &pm, &ps);
    a[i] = 2./(n*n*pm*pm);
  }

  /* generate negative zeros and corresponding weights */
  for (i=0; i<n2; i++) {
    x[n-1-i] = -x[i];
    a[n-1-i] =  a[i];
  }

  return;
}
/*----------------------------------------------------------------------*/
static void null(int n, double *x0, double *p, double *pm, double eps,
	  void (*func)(int, double, double*, double*, double*))
/*
  Newton-Raphson method

  Input:
    n     degree of function func
    x0    approximate extimate of zero
    eps   termination accuracy criterion epsilon
    func  function to be evaluated

  Output:
    p     value of polynomial of degree n
    pm    value of polynomial of degree n-1

*/
{
  int iter;
  double d, ps, x1;

  iter = 0;
  d = 2.*eps;
  while (d > eps) {
    iter++;
    if (iter > 10) {
      fprintf(stderr,"stop in null: too many iterations!\n");
      exit (1);
    }
    func(n, *x0, p, pm, &ps);
    x1 = *x0 - *p/ps;
    d = fabs(x1 - *x0);
    *x0 = x1;
  }
  return;
}
/*----------------------------------------------------------------------*/
static void legpol(int n, double x, double *p, double *pm, double *ps)
/*
  recursive computation of Legendre polynomials

  input:
    n   polynomial degree
    x   argument

  output:
    p   polynomial value (degree n)
    pm  polynomial value (degree n-1)
    ps  first derivative of p (degree n)
*/
{
  int i;
  double p0, p1, tmp;
  
  p0 = 1.0;
  p1 = x;
  for (i=1; i<n; i++) {
    tmp = ((2*i+1) * x * p1 - i * p0) / (i+1);
    p0 = p1;
    p1 = tmp;
  }

  *p = p1;
  *pm = p0;
  *ps = n * (x * p1 - p0) / (x * x - 1.);
  
  return;
}
/*----------------------------------------------------------------------*/
static void lintgr(float t, float r, float dt, int nt, int m, float *wave, double* sum)
/*
  integrate J_n(tR) h(t-tau) using Gauss-Legendre quadrature

  input:
    t      time
    r      parameter used in Bessel function argument
    dt     sample increment of wavelet stored in wave
    nt     number of wavelet samples 
    m      highest order of Bessel functions
    wave   array of wavelet sample values

  output:
    sum    quadrature result
*/
{
  int i, j, m2, mm, n;
  float *u, *xout, *yout;
  double arg, f, wmax;
  double *x, *a, *bes;

  n = t*r; /* number of quadrature points */
  n /= 2;
  if (n < 20) n = 20;

  m2 = m/2;
  mm = m + 50;  /* added indices for better accuracy of high J_n(x) */

  bes  = ealloc1double(mm+1);
  x    = ealloc1double(n);
  a    = ealloc1double(n);
  u    = ealloc1float(n);
  xout = ealloc1float(n);
  yout = ealloc1float(n);

  legcoeffs(n, x, a);

  for (j=0; j<m2; j++) sum[j] = 0.;

  for (i=0; i<n; i++) {
    /* map interval [-1,1] to [0,t] */
    u[i] = 0.5 * t * (x[i] + 1.);
    xout[i] = t - u[i];
  }

  wmax = 0.0;
  for (i=0; i<nt; i++)
    if (fabs(wave[i]) > wmax) wmax = fabs(wave[i]);
  wmax *= 1.0e-6;

  /* sinc interpolation of wavelet to quadrature points */
  ints8r(nt,dt,dt,wave,0.,0.,n,xout,yout);
  for (i=0; i<n; i++) {
    /* ignore very small contributions */
    if (fabs(yout[i]) < wmax) continue;

    f = yout[i] * a[i];

    arg = u[i]*r;
    bessel_jn(arg,mm,bes);

    for (j=0; j<m2; j++) sum[j] += bes[2*j+1] * f;
  }

  /* normalize by interval length */
  for (j=0; j<m2; j++) sum[j] *= 0.5 * t;
  
  free1float(yout);
  free1float(xout);
  free1double(bes);
  free1float(u);
  free1double(a);
  free1double(x);

  return;
} /* end function lintgr */
/**********************************************************************/
static float fwave(float t, float fmax)
/*
!
! --- Ricker wavelet ---
!
! t    : time
! fmax : maximum frequency of wavelet spectrum
!
! returns time sample at time t of a Ricker-wavelet with 
! maximum frequency fmax
!
*/
{
  float pi, pi2, agauss, tcut, s, res;

  pi = 4.*atan(1.);
  pi2 = 2.*pi;
  agauss = 0.5*fmax;
  tcut = 1.5/agauss;

  s = (t-tcut) * agauss;

  if (fabs(s) < 4.) {
    res = exp(-2.*s*s) * cos(pi2*s);
  } else {
    res = 0.;
  }

  return res;
} /* end function fwave */
/**********************************************************************/
static void damp(int nbw, float abso, float *bw)
{
/* 
    compute weights for aborbing boundary region

    input:
      nbw:    number of gridpoints of boundary width
      abso:   strength of absorption

    output:
      bw:     array of absorption weights

*/
  int i;
  float pi, delta;

  pi = 4. * atan(1.);
  delta = pi / nbw;
  
  for (i=0; i<nbw; i++) {
    bw[i] = 1.0 - abso * (1.0 + cos(i*delta)) * 0.5;
  }
 
  return;
}
/**********************************************************************/
static void tbc2d(Remel2dCtx *ctx, float *a, float *wbx, float *wbz, int nwbx, int nwbz,
           int nx, int nz)
/*
   2D absorbing boundary condition

   input:
   a:     array to be treated at boundaries
   wbx:   array of weights for absorbing boundaries in horizontal direction
   wbz:   array of weights for absorbing boundaries in vertical direction
   nbwx:  number of grid points of absorbing zone in x-direction
   nbwz:  number of grid points of absorbing zone in z-direction
   nx:    number of grid points of computational area in x-direction
   nz:    number of grid points of computational area in z-direction

   output:
   a:    array after absorbing boundaries treatment

   Note: The top absorbing region is placed below the bottom absorbing
         region. This is possible due to the cyclic nature of the
         discrete Fourier transform.

     (free) surface
   +----------------+
   |*              *|     * - absorbing region
   |*              *|
   |*              *|
   |*              *|
   |****************|
   |****************|
   +----------------+
         bottom
*/
{
  int *ist;
  int i, k;
  int ifl, ihv, istart, is, iz;
  int ind1, ind2, nwb, ng1, ng2, nz0;   
  float wb;
  
  if (ctx->ist == NULL) {
    ctx->ist = ealloc1int(nz);
    ist = ctx->ist;
    ng1 = nwbx;
    ng2 = nwbz;
    ihv = 0;
    if (ng1 > ng2) {
      ng1 = nwbz;
      ng2 = nwbx;
      ihv = 1;
    }
    ifl = -ng2;
    istart = nwbx;
    iz = 0;
    if (ihv == 1) ist[0] = nwbx;
    for (i=0; i<ng2; i++) {
      if (ifl >= 0) {
	ifl = ifl-ng2;
	iz++;
	if (ihv == 1) ist[iz] = nwbx-i;
	istart--;
      }
      ifl = ifl + ng1;
      if (ihv == 0) ist[i] = istart;
      /* printf("i= %d ist= %d\n",i,ist[i]); */
    }
  }
  ist = ctx->ist;
/*
   +------------------------+
   |   Treatment of sides   |
   +------------------------+
*/
  nz0 = nz - 2*nwbz;
  for (k=0; k<nz; k++) {
    nwb = nwbx;
    if (k >= nz0+nwbz) nwb = ist[nz-k-1];
    if (k >= nz0 && k < nz0+nwbz) nwb = ist[k-nz0];

    for (i=0; i<nwb; i++) {
      wb = wbx[i];
      ind1 = i + k*nx;
      a[ind1] *= wb;
      ind2 = (nx-i-1) + k*nx;
      a[ind2] *= wb;
    }
  }

/*
   +-------------------------+
   |   Treatment of bottom   |
   +-------------------------+
*/
  for (k=0; k<nwbz; k++) {
    is = ist[k];
    wb = wbz[nwbz-k-1];

    for (i=is; i<nx-is; i++) {
      ind1 = i + (k+nz0)*nx;
      ind2 = i + (nz-k-1)*nx;
      a[ind1] *= wb;
      a[ind2] *= wb;
    }
  }
}


/* Library version of SUREMEL2DAN
 *
 * The main program (the parameters, the files of the model and of the wavelet, the trace headers, the output files) is left to the
 * caller. What is here is the modelling of one shot, in a single call, with every array owned by the caller:
 *
 *	su_remel2dan:	the sections along horizontal and vertical receiver lines, and the snapshots
 *
 * The 2-D arrays are flat, [nz][nx] with x the fast axis (the program reads files with z the fast axis). The types of the
 * receivers and snapshots are 0 P, 1 S, 2 UX, 3 UZ. The sections are xsect[nsectx][nt][nx], zsect[nsectz][nt][nz] and the
 * snapshots snap[nsntyp][nsnap][nz][nx], nsnap = NINT(nt*dt/dtsnap). The model is the density and either (amode=0) the velocities vp
 * and vs, or (amode=1) the Voigt stiffnesses c11, c13, c15, c33, c35 and c55 with the global velocities vmax0 and vmin0.
 * The source types styp are 0 pressure, 1 shear, 2 single force in x, 3 single force in z.
 *
 * Differences from the program:
 *  - The Fourier differentiation and the sponge keep their tables in a structure of the call (the program had static variables),
 *    so the function can be called repeatedly and from several threads.
 *  - The Bessel coefficients are always computed (prec=0): the option of keeping them in the file b_k is not part of the library.
 *  - The pressure sections along a vertical receiver line (rztyp=0) were filled in where the shear ones (rztyp=1) were tested
 *    (so that they were never made, and the shear ones were overwritten by the pressure); they are made for rztyp=0.
 *  - The arrays that hold the source in the space domain start at 0 and are cleared for each source (the program did not set them,
 *    and the second and later sources included the boxes of the earlier ones again).
 *  - The odd-size check works (the program's was always false), and the dispersion check is only made if fmax is known.
 *  - The checks of the parameters return a code instead of stopping the program.
 */

#define XSECT(il,it,i) xsect[((size_t)(il)*nt+(it))*nx+(i)]
#define ZSECT(il,it,k) zsect[((size_t)(il)*nt+(it))*nz+(k)]
#define SNAP(il,it,k,i) snap[(((size_t)(il)*nsnap+(it))*nz+(k))*nx+(i)]

#define NXB 21
#define NZB 21

/* Returns 0, or: -1 invalid nx or nz (not a pfafft length, or even), -2 a source box beyond the model, -3 a parameter that is not
 * allowed, -4 vmax is zero, -5 fmax too high for the spatial sampling, -6 vmaxu is smaller than vmax, -7 a receiver line outside of
 * the model. */
int su_remel2dan(int amode, int nx, int nz, int nt, float dx, float dz, float dt,
	int nsourc, const int *isx, const int *isz, const int *styp, const float *samp, float w,
	int sflag, float fmax, int nwav_in, float dtwave_in, const float *wave_in, float vmaxu, float vmax0, float vmin0,
	float dtsnap, int nsntyp, const int *sntyp, int nsnap, int iabso, float abso, int nbwx, int nbwz,
	const float *dens, const float *vp, const float *vs,
	const float *c11in, const float *c13in, const float *c15in, const float *c33in, const float *c35in, const float *c55in,
	int nsectx, const int *irz, const int *rxtyp, int nsectz, const int *irx, const int *rztyp,
	float *xsect, float *zsect, float *snap)
{
	int i, k, l, il, it, indx, indz, m, m2, nzpad, nwav = 0, rc = 0, verbose = 0;
	float tmax, vmax, vmin, r, pi, dtwave = 0., theta, tmp1, tmp2, tmp3, t, amx, amz;
	float *bwx = NULL, *bwz = NULL, *wave = NULL;
	float **qux = NULL, **quz = NULL, **q2ux = NULL, **q2uz = NULL, **a1 = NULL, **a2 = NULL, **a3 = NULL;
	float **c11 = NULL, **c13 = NULL, **c15 = NULL, **c33 = NULL, **c35 = NULL, **c55 = NULL, **lambda = NULL, **twmy = NULL;
	float **rhoinv = NULL, **gbox = NULL;
	double **bessn = NULL, **bestr = NULL, *btemp = NULL;
	Remel2dCtx ctxs, *ctx = &ctxs;

	memset(ctx, 0, sizeof(Remel2dCtx));
	(void)amx; (void)amz; (void)verbose;
	if ((amode!=0 && amode!=1) || nt<1 || nx<1 || nz<1 || dx<=0.0 || dz<=0.0 || dt<=0.0 || dtsnap<0.0 || nsourc<1 || nbwx<1 || nbwz<1) return -3;
	if (nsnap!=((dtsnap>0.0) ? NINT(nt*dt/dtsnap) : 0)) return -3;
	if (nsnap>0 && nsntyp<1) return -3;
	if (sflag==2 && fmax<=0.0) return -3;
	for (l=0; l<nsourc; l++) if (styp[l]<0 || styp[l]>3) return -3;
	for (il=0; il<nsectx; il++) if (irz[il]<0 || irz[il]>=nz) return -7; else if (rxtyp[il]<0 || rxtyp[il]>3) return -3;
	for (il=0; il<nsectz; il++) if (irx[il]<0 || irx[il]>=nx) return -7; else if (rztyp[il]<0 || rztyp[il]>3) return -3;
	for (il=0; il<nsntyp; il++) if (sntyp[il]<0 || sntyp[il]>3) return -3;

	nzpad = nz;
	if (npfa(nx) != nx) return -1;
	if (!(nx % 2)) return -1;
	if (npfa(nz) != nz) return -1;
	if (!(nz % 2)) return -1;

	for (l=0; l<nsourc; l++) {
		if (isx[l] - (NXB+1)/2 < 0 || isx[l] + (NXB+1)/2 > nx-1) return -2;
		if (isz[l] - (NZB+1)/2 < 0 || isz[l] + (NZB+1)/2 > nz-1) return -2;
	}

	tmax = nt * dt;

	bwx    = ealloc1float(nbwx);
	bwz    = ealloc1float(nbwz);
	qux    = ealloc2float(nx,nz);
	quz    = ealloc2float(nx,nz);
	q2ux   = ealloc2float(nx,nz);
	q2uz   = ealloc2float(nx,nz);
	a1     = ealloc2float(nx,nz);
	a2     = ealloc2float(nx,nz);
	a3     = ealloc2float(nx,nz);
	rhoinv = ealloc2float(nx,nz);
	memset(a1[0], 0, sizeof(float)*(size_t)nx*nz);
	memset(a2[0], 0, sizeof(float)*(size_t)nx*nz);
	memset(a3[0], 0, sizeof(float)*(size_t)nx*nz);
	if (amode) {
		c11 =    ealloc2float(nx,nz);
		c13 =    ealloc2float(nx,nz);
		c15 =    ealloc2float(nx,nz);
		c33 =    ealloc2float(nx,nz);
		c35 =    ealloc2float(nx,nz);
		c55 =    ealloc2float(nx,nz);
	} else {
		lambda = ealloc2float(nx,nz);
		twmy   = ealloc2float(nx,nz);
	}
	gbox   = ealloc2float(NXB,NZB);

	if (sflag == 0) {
		if (nwav_in<1 || dtwave_in<=0.0) { rc = -3; goto done; }
		nwav = nwav_in;
		dtwave = dtwave_in;
		wave = ealloc1float(nwav);
		for (it=0; it<nwav; it++) wave[it] = wave_in[it];
	}

	if (sflag == 2) {
		dtwave = dt;
		nwav = 8./(dt*fmax);
		wave = ealloc1float(nwav>0 ? nwav : 1);
		for (it=0; it<nwav; it++) wave[it] = fwave((it+1)*dt,fmax);
	}

	/* --- the subsurface structure --- */
	vmin = 1.0e30;
	vmax = 0.;

	for (k=0; k<nz; k++)
		for (i=0; i<nx; i++) rhoinv[k][i] = 1./dens[(size_t)k*nx+i];

	if (amode) {
		for (k=0; k<nz; k++) {
			for (i=0; i<nx; i++) {
				size_t j = (size_t)k*nx+i;
				c11[k][i] = c11in[j];
				c13[k][i] = c13in[j];
				c15[k][i] = c15in[j];
				c33[k][i] = c33in[j];
				c35[k][i] = c35in[j];
				c55[k][i] = c55in[j];
			}
		}
		vmax = vmax0;
		vmin = vmin0;
	} else {
		for (k=0; k<nz; k++) {
			for (i=0; i<nx; i++) {
				size_t j = (size_t)k*nx+i;
				twmy[k][i] = 2.0*vs[j]*vs[j]/rhoinv[k][i];
				vmin = (vs[j] < vmin ? vs[j] : vmin);
				lambda[k][i] = vp[j]*vp[j]/rhoinv[k][i] - twmy[k][i];
				vmax = (vp[j] > vmax ? vp[j] : vmax);
			}
		}
	}

	if (vmax == 0.) { rc = -4; goto done; }
	/* (the program used fmax here even when it was not given, for a spike or a user wavelet: it is checked if it is) */
	if (fmax > 0.0 && (dx > dz ? dx : dz) > vmin/(2.*fmax)) { rc = -5; goto done; }

	/* change absorbing zone for no reflections at top */
	if (iabso) {
		if (amode) {
			for (i=0; i<nx; i++) {
				for (k=nz-nbwz; k<nz; k++) {
					rhoinv[k][i] = rhoinv[0][i];
					c11[k][i] = c11[0][i];
					c13[k][i] = c13[0][i];
					c15[k][i] = c15[0][i];
					c33[k][i] = c33[0][i];
					c35[k][i] = c35[0][i];
					c55[k][i] = c55[0][i];
				}
			}
		} else {
			for (i=0; i<nx; i++) {
				for (k=nz-nbwz; k<nz; k++) {
					rhoinv[k][i] = rhoinv[0][i];
					twmy[k][i] = twmy[0][i];
					lambda[k][i] = lambda[0][i];
				}
			}
		}
	}

	/* setup aborbing boundaries */
	if (iabso) {
		damp(nbwx,abso,bwx);
		damp(nbwz,abso,bwz);
	}

	/* prepare Bessel coefficients convolved with wavelet */
	if (vmaxu > 0.) {
		if (vmax > vmaxu) { rc = -6; goto done; }
		vmax = vmaxu;
	}

	/* largest EV of operator L */
	pi = 4.*atan(1.);
	r = 1.1 * pi * vmax * sqrt(1./(dx*dx) + 1./(dz*dz));

	m = tmax * r; /* highest index of terms in expansion */
	m2 = m/2; /* number of expansion terms */
	if (m2 < 2) { rc = -3; goto done; }

	bessn = ealloc2double(m2,nsnap>0 ? nsnap : 1);
	bestr = ealloc2double(m2,nt);
	btemp = ealloc1double(m+50);

	/* Bessel coefficients */
	for (it=0; it<nsnap; it++) {
		t = (it+1)*dtsnap;
		if (sflag == 1) {
			bessel_jn(t*r,m+50,btemp);
			for (l=0; l<m2; l++) bessn[it][l] = btemp[2*l+1];
		} else {
			lintgr(t,r,dtwave,nwav,m,wave,bessn[it]);
		}
	}

	for (it=0; it<nt; it++) {
		t = (it+1)*dt;
		if (sflag == 1) {
			bessel_jn(t*r,m+50,btemp);
			for (l=0; l<m2; l++) bestr[it][l] = btemp[2*l+1];
		} else {
			lintgr(t,r,dtwave,nwav,m,wave,bestr[it]);
		}
	}

	for (k=0; k<nz; k++) {
		for (i=0; i<nx; i++) qux[k][i] = quz[k][i] = 0.;
		for (i=0; i<nx; i++) rhoinv[k][i] *= 4./(r*r);
	}

	/* Gaussian source box */
	g2d(w, NXB, NZB, gbox);

  /* initialize first expansion term Q_1 */
  /* source distribution */ 
  for (l=0; l<nsourc; l++) { /* pressure sources */
    memset(a1[0], 0, sizeof(float)*(size_t)nx*nz);
    memset(a2[0], 0, sizeof(float)*(size_t)nx*nz);
    if (styp[l] == 0) {
      for (k=0; k<NZB; k++) {
	indz = isz[l]+k-NZB/2;
	for (i=0; i<NXB; i++) {
	  indx = isx[l]+i-NXB/2;
	  a1[indz][indx] = -0.5 * gbox[k][i] * samp[l];
	}
      }
      difx(ctx,a1,qux,NULL,nx,nz,dx,1);
      difz(ctx,a1,quz,NULL,nx,nz,nzpad,dz,1,0);
    } else if (styp[l] == 1) { /* shear sources */
      for (k=0; k<NZB; k++) {
	indz = isz[l]+k-NZB/2;
	for (i=0; i<NXB; i++) {
	  indx = isx[l]+i-NXB/2;
	  a1[indz][indx] = -0.5 * gbox[k][i] * samp[l];
	  a2[indz][indx] = -a1[indz][indx];
	}
      }
      difx(ctx,a2,quz,NULL,nx,nz,dx,1);
      difz(ctx,a1,qux,NULL,nx,nz,nzpad,dz,1,0);
    } else if (styp[l] == 2) { /* single forces (x-direction) */
      for (k=0; k<NZB; k++) {
	indz = isz[l]+k-NZB/2;
	for (i=0; i<NXB; i++) {
	  indx = isx[l]+i-NXB/2;
	  qux[indz][indx] += gbox[k][i] * samp[l];
	}
      }
    } else if (styp[l] == 3) { /* single forces (z-direction) */
      for (k=0; k<NZB; k++) {
	indz = isz[l]+k-NZB/2;
	for (i=0; i<NXB; i++) {
	  indx = isx[l]+i-NXB/2;
	  quz[indz][indx] += gbox[k][i] * samp[l];
	}
      }
    }
  }

  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      qux[k][i] *= 1./(r*dx*dz);
      quz[k][i] *= 1./(r*dx*dz);
    }
  }

  /* initialize second expansion term Q_3 */
  /* strains */
  /* E_xx */
  difx(ctx,qux,a1,NULL,nx,nz,dx,0);
  /* E_zz */
  difz(ctx,quz,a2,NULL,nx,nz,nzpad,dz,0,0);
  /* E_xz */
  difx(ctx,quz,a3,NULL,nx,nz,dx,0);
  
  /* shear snapshots (part I) */
  for (il=0; il<nsntyp; il++) {
    if (sntyp[il] == 1) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) = 2.*a3[k][i]*bessn[it][0];
	}
      }
    }
  }

  /* shear sections (part I) */
  for (il=0; il<nsectx; il++) {
    if (rxtyp[il] == 1) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) = 2.*a3[irz[il]][i]*bestr[it][0];
      }
    }
  }
  for (il=0; il<nsectz; il++) {
    if (rztyp[il] == 1) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) = 2.*a3[k][irx[il]]*bestr[it][0];
      }
    }
  }

  difz(ctx,qux,a3,NULL,nx,nz,nzpad,dz,1,0);

  /* shear snapshots (part II) */
  for (il=0; il<nsntyp; il++) {
    if (sntyp[il] == 1) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) -= a3[k][i]*bessn[it][0];
	}
      }
    }
  }
  
  /* shear sections (part II) */
  for (il=0; il<nsectx; il++) {
    if (rxtyp[il] == 1) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) -= a3[irz[il]][i]*bestr[it][0];
      }
    }
  }
  for (il=0; il<nsectz; il++) {
    if (rztyp[il] == 1) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) -= a3[k][irx[il]]*bestr[it][0];
      }
    }
  }

#if 0
  /* displacement divergence snapshots */
  for (il=0; il<nsntyp; il++) {
    if (sntyp[il] == 0) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) = 0.5*(a1[k][i]+a2[k][i])*bessn[it][0];
	}
      }
    }
  }
  
  /* displacement divergence sections */
  for (il=0; il<nsectx; il++) {
    if (rxtyp[il] == 0) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) = 0.5*(a1[irz[il]][i]+a2[irz[il]][i])*bestr[it][0];
      }
    }
  }
  for (il=0; il<nsectz; il++) {
    if (rztyp[il] == 1) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) = 0.5*(a1[k][irx[il]]+a2[k][irx[il]])*bestr[it][0];
      }
    }
  }
#endif

  /* --- stress-strain relation --- */
  if (amode) {
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	tmp1 = a1[k][i]*c11[k][i] + a2[k][i]*c13[k][i] + a3[k][i]*c15[k][i];
	tmp2 = a1[k][i]*c13[k][i] + a2[k][i]*c33[k][i] + a3[k][i]*c35[k][i];
	tmp3 = a1[k][i]*c15[k][i] + a2[k][i]*c35[k][i] + a3[k][i]*c55[k][i];
	a1[k][i] = tmp1;
	a2[k][i] = tmp2;
	a3[k][i] = tmp3;
      }
    }
  } else {
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	theta = a1[k][i] + a2[k][i];
	a1[k][i] = lambda[k][i]*theta + twmy[k][i]*a1[k][i];
	a2[k][i] = lambda[k][i]*theta + twmy[k][i]*a2[k][i];
	a3[k][i] = twmy[k][i]*a3[k][i]*0.5;
      }
    }
  }
  
  /* pressure snapshots */
  for (il=0; il<nsntyp; il++) {
    if (sntyp[il] == 0) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) = 0.5*(a1[k][i]+a2[k][i])*bessn[it][0];
	}
      }
    }
  }
  
  /* pressure sections */
  for (il=0; il<nsectx; il++) {
    if (rxtyp[il] == 0) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) = 0.5*(a1[irz[il]][i]+a2[irz[il]][i])*bestr[it][0];
      }
    }
  }
  for (il=0; il<nsectz; il++) {
    if (rztyp[il] == 0) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) = 0.5*(a1[k][irx[il]]+a2[k][irx[il]])*bestr[it][0];
      }
    }
  }
  
  /* --- accellerations --- */
  difx(ctx,a1,q2ux,rhoinv,nx,nz,dx,2);
  difz(ctx,a3,q2ux,rhoinv,nx,nz,nzpad,dz,3,0);
  
  difz(ctx,a2,q2uz,rhoinv,nx,nz,nzpad,dz,2,0);
  difx(ctx,a3,q2uz,rhoinv,nx,nz,dx,3);

  if (iabso) {
    tbc2d(ctx,qux[0],bwx,bwz,nbwx,nbwz,nx,nz);
    tbc2d(ctx,quz[0],bwx,bwz,nbwx,nbwz,nx,nz);
    tbc2d(ctx,q2ux[0],bwx,bwz,nbwx,nbwz,nx,nz);
    tbc2d(ctx,q2uz[0],bwx,bwz,nbwx,nbwz,nx,nz);
  }

  /* generation of second term */
  for (k=0; k<nz; k++) {
    for (i=0; i<nx; i++) {
      q2ux[k][i] += 3.*qux[k][i];
      q2uz[k][i] += 3.*quz[k][i];
    }
  }

  /* displacement snapshots */
  for (il=0; il<nsntyp; il++) {
    if (sntyp[il] == 2) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) = qux[k][i]*bessn[it][0];
	}
      }
    }
    if (sntyp[il] == 3) {
      for (it=0; it<nsnap; it++) {
	for (k=0; k<nz; k++) {
	  for (i=0; i<nx; i++)
	    SNAP(il,it,k,i) = quz[k][i]*bessn[it][0];
	}
      }
    }
  }
  
  /* displacement sections */
  for (il=0; il<nsectx; il++) {
    if (rxtyp[il] == 2) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) = qux[irz[il]][i]*bestr[it][0];
      }
    }
    if (rxtyp[il] == 3) {
      for (it=0; it<nt; it++) {
	for (i=0; i<nx; i++)
	  XSECT(il,it,i) = quz[irz[il]][i]*bestr[it][0];
      }
    }
  }

  for (il=0; il<nsectz; il++) {
    if (rztyp[il] == 2) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) = qux[k][irx[il]]*bestr[it][0];
      }
    }
    if (rztyp[il] == 3) {
      for (it=0; it<nt; it++) {
	for (k=0; k<nz; k++)
	  ZSECT(il,it,k) = quz[k][irx[il]]*bestr[it][0];
      }
    }
  }
  
  /* calculate all other odd expansion terms Q_2n+1 */
  for (l=1; l<m2; l++) {

    if (verbose == 2) {
      amx = 0.;
      amz = 0.;
      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  amx = (fabs(q2ux[k][i]) > amx ? fabs(q2ux[k][i]) : amx);
	  amz = (fabs(q2uz[k][i]) > amz ? fabs(q2uz[k][i]) : amz);
	}
      }
      fprintf(stderr,"l=%5d    amx=%e  amz=%e\n",l,amx,amz);
    }

    /* --- recursion relation --- */
    for (k=0; k<nz; k++) {
      for (i=0; i<nx; i++) {
	a1[k][i] = q2ux[k][i];
	a2[k][i] = q2uz[k][i];

	q2ux[k][i] = 2.*q2ux[k][i] - qux[k][i];
	q2uz[k][i] = 2.*q2uz[k][i] - quz[k][i];

	qux[k][i] = a1[k][i];
	quz[k][i] = a2[k][i];
      }
    }

    /* strains */
    /* E_xx */
    difx(ctx,qux,a1,NULL,nx,nz,dx,0);
    /* E_zz */
    difz(ctx,quz,a2,NULL,nx,nz,nzpad,dz,0,0);
    /* E_xz */
    difx(ctx,quz,a3,NULL,nx,nz,dx,0);
    
    /* shear snapshots (part I) */
    for (il=0; il<nsntyp; il++) {
      if (sntyp[il] == 1) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) += 2.*a3[k][i]*bessn[it][l];
	  }
	}
      }
    }
    
    /* shear sections (part I) */
    for (il=0; il<nsectx; il++) {
      if (rxtyp[il] == 1) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) += 2.*a3[irz[il]][i]*bestr[it][l];
	}
      }
    }
    for (il=0; il<nsectz; il++) {
      if (rztyp[il] == 1) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) += 2.*a3[k][irx[il]]*bestr[it][l];
	}
      }
    }
    
    difz(ctx,qux,a3,NULL,nx,nz,nzpad,dz,1,0);
    
    /* shear snapshots (part II) */
    for (il=0; il<nsntyp; il++) {
      if (sntyp[il] == 1) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) -= a3[k][i]*bessn[it][l];
	  }
	}
      }
    }
    
    /* shear sections (part II) */
    for (il=0; il<nsectx; il++) {
      if (rxtyp[il] == 1) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) -= a3[irz[il]][i]*bestr[it][l];
	}
      }
    }
    for (il=0; il<nsectz; il++) {
      if (rztyp[il] == 1) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) -= a3[k][irx[il]]*bestr[it][l];
	}
      }
    }

#if 0    
    /* displacment divergence snapshots */
    for (il=0; il<nsntyp; il++) {
      if (sntyp[il] == 0) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) += 0.5*(a1[k][i]+a2[k][i])*bessn[it][l];
	  }
	}
      }
    }
    
    /* displacment divergence sections */
    for (il=0; il<nsectx; il++) {
      if (rxtyp[il] == 0) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) += 0.5*(a1[irz[il]][i]+a2[irz[il]][i])*bestr[it][l];
	}
      }
    }
    for (il=0; il<nsectz; il++) {
      if (rztyp[il] == 1) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) += 0.5*(a1[k][irx[il]]+a2[k][irx[il]])*bestr[it][l];
	}
      }
    }
#endif

    /* --- stress-strain relation --- */
    if (amode) {
      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  tmp1 = a1[k][i]*c11[k][i] + a2[k][i]*c13[k][i] + a3[k][i]*c15[k][i];
	  tmp2 = a1[k][i]*c13[k][i] + a2[k][i]*c33[k][i] + a3[k][i]*c35[k][i];
	  tmp3 = a1[k][i]*c15[k][i] + a2[k][i]*c35[k][i] + a3[k][i]*c55[k][i];
	  a1[k][i] = tmp1;
	  a2[k][i] = tmp2;
	  a3[k][i] = tmp3;
	  }
      }
    } else {
      for (k=0; k<nz; k++) {
	for (i=0; i<nx; i++) {
	  theta = a1[k][i] + a2[k][i];
	  a1[k][i] = lambda[k][i]*theta + twmy[k][i]*a1[k][i];
	  a2[k][i] = lambda[k][i]*theta + twmy[k][i]*a2[k][i];
	  a3[k][i] = twmy[k][i]*a3[k][i]*0.5;
	}
      }
    }
    
    /* pressure snapshots */
    for (il=0; il<nsntyp; il++) {
      if (sntyp[il] == 0) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) += 0.5*(a1[k][i]+a2[k][i])*bessn[it][l];
	  }
	}
      }
    }
    
    /* pressure sections */
    for (il=0; il<nsectx; il++) {
      if (rxtyp[il] == 0) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) += 0.5*(a1[irz[il]][i]+a2[irz[il]][i])*bestr[it][l];
	}
      }
    }
    for (il=0; il<nsectz; il++) {
      if (rztyp[il] == 0) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) += 0.5*(a1[k][irx[il]]+a2[k][irx[il]])*bestr[it][l];
	}
      }
    }

    /* --- accellerations --- */
    difx(ctx,a1,q2ux,rhoinv,nx,nz,dx,3);
    difz(ctx,a3,q2ux,rhoinv,nx,nz,nzpad,dz,3,0);
    
    difz(ctx,a2,q2uz,rhoinv,nx,nz,nzpad,dz,3,0);
    difx(ctx,a3,q2uz,rhoinv,nx,nz,dx,3);

    if (iabso) {
      tbc2d(ctx,qux[0],bwx,bwz,nbwx,nbwz,nx,nz);
      tbc2d(ctx,quz[0],bwx,bwz,nbwx,nbwz,nx,nz);
      tbc2d(ctx,q2ux[0],bwx,bwz,nbwx,nbwz,nx,nz);
      tbc2d(ctx,q2uz[0],bwx,bwz,nbwx,nbwz,nx,nz);
    }

    /* displacement snapshots */
    for (il=0; il<nsntyp; il++) {
      if (sntyp[il] == 2) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) += qux[k][i]*bessn[it][l];
	  }
	}
      }
      if (sntyp[il] == 3) {
	for (it=0; it<nsnap; it++) {
	  for (k=0; k<nz; k++) {
	    for (i=0; i<nx; i++)
	      SNAP(il,it,k,i) += quz[k][i]*bessn[it][l];
	  }
	}
      }
    }
    
    /* displacement sections */
    for (il=0; il<nsectx; il++) {
      if (rxtyp[il] == 2) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) += qux[irz[il]][i]*bestr[it][l];
	}
      }
      if (rxtyp[il] == 3) {
	for (it=0; it<nt; it++) {
	  for (i=0; i<nx; i++)
	    XSECT(il,it,i) += quz[irz[il]][i]*bestr[it][l];
	}
      }
    }

    for (il=0; il<nsectz; il++) {
      if (rztyp[il] == 2) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) += qux[k][irx[il]]*bestr[it][l];
	}
      }
      if (rztyp[il] == 3) {
	for (it=0; it<nt; it++) {
	  for (k=0; k<nz; k++)
	    ZSECT(il,it,k) += quz[k][irx[il]]*bestr[it][l];
	}
      }
    }
  
  } /* end loop over expansion terms */

done:
	if (bwx) free1float(bwx);
	if (bwz) free1float(bwz);
	if (qux) free2float(qux);
	if (quz) free2float(quz);
	if (q2ux) free2float(q2ux);
	if (q2uz) free2float(q2uz);
	if (a1) free2float(a1);
	if (a2) free2float(a2);
	if (a3) free2float(a3);
	if (rhoinv) free2float(rhoinv);
	if (gbox) free2float(gbox);
	if (c11) { free2float(c11); free2float(c13); free2float(c15); free2float(c33); free2float(c35); free2float(c55); }
	if (lambda) { free2float(lambda); free2float(twmy); }
	if (wave) free1float(wave);
	if (bessn) free2double(bessn);
	if (bestr) free2double(bestr);
	if (btemp) free1double(btemp);
	if (ctx->ax) free1float(ctx->ax);
	if (ctx->rkx1) free1float(ctx->rkx1);
	if (ctx->az) free1float(ctx->az);
	if (ctx->rkz1) free1float(ctx->rkz1);
	if (ctx->ist) free1int(ctx->ist);
	return rc;
}
