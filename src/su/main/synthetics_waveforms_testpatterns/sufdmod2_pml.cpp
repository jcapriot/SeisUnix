/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUFDMOD2_PML: $Revision: 1.9 $ ; $Date: 2011/11/12 00:40:42 $	*/


#include "par.h"
#include "su.h"
#include "segy.h"

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
" 									",
" SUFDMOD2_PML - Finite-Difference MODeling (2nd order) for acoustic wave",
"    equation with PML absorbing boundary conditions.			",
" Caveat: experimental PML absorbing boundary condition version,	",
"may be buggy!								",
" 									",
" sufdmod2_pml <vfile >wfile nx= nz= tmax= xs= zs= [optional parameters]",
" 									",
" Required Parameters:							",
" <vfile		file containing velocity[nx][nz]		",
" >wfile		file containing waves[nx][nz] for time steps	",
" nx=			number of x samples (2nd dimension)		",
" nz=			number of z samples (1st dimension)		",
" xs=			x coordinates of source				",
" zs=			z coordinates of source				",
" tmax=			maximum time					",
" 									",
" Optional Parameters:							",
" nt=1+tmax/dt		number of time samples (dt determined for stability)",
" mt=1			number of time steps (dt) per output time step	",
" 									",
" dx=1.0		x sampling interval				",
" fx=0.0		first x sample					",
" dz=1.0		z sampling interval				",
" fz=0.0		first z sample					",
" 									",
" fmax = vmin/(10.0*h)	maximum frequency in source wavelet		",
" fpeak=0.5*fmax	peak frequency in ricker wavelet		",
" 									",
" dfile=		input file containing density[nx][nz]		",
" vsx=			x coordinate of vertical line of seismograms	",
" hsz=			z coordinate of horizontal line of seismograms	",
" vsfile=		output file for vertical line of seismograms[nz][nt]",
" hsfile=		output file for horizontal line of seismograms[nx][nt]",
" ssfile=		output file for source point seismograms[nt]	",
" verbose=0		=1 for diagnostic messages, =2 for more		",
" 									",
" abs=1,1,1,1		Absorbing boundary conditions on top,left,bottom,right",
" 			sides of the model. 				",
" 		=0,1,1,1 for free surface condition on the top		",
"                                                                       ",
" ...PML parameters....                                                 ",
" pml_max=1000.0        PML absorption parameter                        ",
" pml_thick=0           half-thickness of pml layer (0 = do not use PML)",
" 									",
" Notes:								",
" This program uses the traditional explicit second order differencing	",
" method. 								",
" 									",
" Two different absorbing boundary condition schemes are available. The ",
" first is a traditional absorbing boundary condition scheme created by ",
" Hale, 1990. The second is based on the perfectly matched layer (PML)	",
" method of Berenger, 1995.						",
" 									",
NULL};
#endif


/*
 * Authors:  CWP:Dave Hale
 *           CWP:modified for SU by John Stockwell, 1993.
 *           CWP:added frequency specification of wavelet: Craig Artley, 1993
 *           TAMU:added PML absorbing boundary condition: 
 *               Michael Holzrichter, 1998
 *           CWP/WesternGeco:corrected PML code to handle density variations:
 *               Greg Wimpey, 2006
 *
 * References: (Hale's absobing boundary conditions)
 * Clayton, R. W., and Engquist, B., 1977, Absorbing boundary conditions
 * for acoustic and elastic wave equations, Bull. Seism. Soc. Am., 6,
 *	1529-1540. 
 *
 * Clayton, R. W., and Engquist, B., 1980, Absorbing boundary conditions
 * for wave equation migration, Geophysics, 45, 895-904.
 *
 * Hale, D.,  1990, Adaptive absorbing boundaries for finite-difference
 * modeling of the wave equation migration, unpublished report from the
 * Center for Wave Phenomena, Colorado School of Mines.
 *
 * Richtmyer, R. D., and Morton, K. W., 1967, Difference methods for
 * initial-value problems, John Wiley & Sons, Inc, New York.
 *
 * Thomee, V., 1962, A stable difference scheme for the mixed boundary problem
 * for a hyperbolic, first-order system in two dimensions, J. Soc. Indust.
 * Appl. Math., 10, 229-245.
 *
 * Toldi, J. L., and Hale, D., 1982, Data-dependent absorbing side boundaries,
 * Stanford Exploration Project Report SEP-30, 111-121.
 *
 * References: (PML boundary conditions)
 * Jean-Pierre Berenger, ``A Perfectly Matched Layer for the Absorption of
 * Electromagnetic Waves,''  Journal of Computational Physics, vol. 114,
 * pp. 185-200.
 *
 * Hastings, Schneider, and Broschat, ``Application of the perfectly
 * matched layer (PML) absorbing boundary condition to elastic wave
 * propogation,''  Journal of the Accoustical Society of America,
 * November, 1996.
 *
 * Allen Taflove, ``Electromagnetic Modeling:  Finite Difference Time
 * Domain Methods'', Baltimore, Maryland: Johns Hopkins University Press,
 * 1995, chap. 7, pp. 181-195.
 *
 *
 * Trace header fields set: ns, delrt, tracl, tracr, offset, d1, d2,
 *                          sdepth, trid
 */
/**************** end self doc ********************************/


struct FdmodPml {
	float pml_max;
	int pml_thick;
	int pml_thickness;
	float **cax_b, **cax_r;
	float **cbx_b, **cbx_r;
	float **caz_b, **caz_r;
	float **cbz_b, **cbz_r;
	float **dax_b, **dax_r;
	float **dbx_b, **dbx_r;
	float **daz_b, **daz_r;
	float **dbz_b, **dbz_r;
	float **ux_b, **ux_r;
	float **uz_b, **uz_r;
	float **v_b, **v_r;
	float **w_b, **w_r;
	float dvv_0, dvv_1, dvv_2, dvv_3;
	float sigma, sigma_ex, sigma_ez, sigma_mx, sigma_mz;
};

static void pml_init (FdmodPml *P, int nx, int nz, float dx, float dz, float dt,
		      float **dvv, float **od);
static void pml_absorb (FdmodPml *P, int nx, float dx, int nz, float dz, float dt,
        float **dvv, float **od, float **pm, float **p, float **pp, int *abs);
static float **zalloc2 (size_t n1, size_t n2);

static void exsrc (int ns, float *xs, float *zs, float *vs, float (*xsd)[4], float (*zsd)[4],
	int nx, float dx, float fx,
	int nz, float dz, float fz,
	float dt, float t, float fmax, float **s)
/*****************************************************************************
exsrc - update source pressure function for an extended source
******************************************************************************
Input:
ns		number of x,z coordinates for extended source
xs		array[ns] of x coordinates of extended source
zs		array[ns] of z coordinates of extended source
nx		number of x samples
dx		x sampling interval
fx		first x sample
nz		number of z samples
dz		z sampling interval
fz		first z sample
dt		time step (ignored)
t		time at which to compute source function
fmax		maximum frequency

Output:
s		array[nx][nz] of source pressure at time t+dt
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 03/01/90
******************************************************************************/
{
	int ix,iz,ixv,izv;
	float sigma,tbias,ascale,tscale,ts,xn,zn,
		v,xv,zv,dxdv,dzdv,xvn,zvn,amp,dv,dist,distprev;

	/* zero source array */
	for (ix=0; ix<nx; ++ix)
		for (iz=0; iz<nz; ++iz)
			s[ix][iz] = 0.0 *dt ;
	
	/* compute time-dependent part of source function */
	sigma = 0.25/fmax;
	tbias = 3.0*sigma;
	ascale = -exp(0.5)/sigma;
	tscale = 0.5/(sigma*sigma);
	if (t>2.0*tbias) return;
	ts = ascale*(t-tbias)*exp(-tscale*(t-tbias)*(t-tbias));
	
	/* loop over extended source locations */
	for (v=vs[0],distprev=0.0,dv=1.0; dv!=0.0; distprev=dist,v+=dv) {
		
		/* determine x(v), z(v), dx/dv, and dz/dv along source */
		intcub(0,ns,vs,xsd,1,&v,&xv);
		intcub(0,ns,vs,zsd,1,&v,&zv);
		intcub(1,ns,vs,xsd,1,&v,&dxdv);
		intcub(1,ns,vs,zsd,1,&v,&dzdv);
		
		/* determine increment along extended source */
		if (dxdv==0.0)
			dv = dz/ABS(dzdv);
		else if (dzdv==0.0)
			dv = dx/ABS(dxdv);
		else
			dv = MIN(dz/ABS(dzdv),dx/ABS(dxdv));
		if (v+dv>vs[ns-1]) dv = vs[ns-1]-v;
		dist = dv*sqrt(dzdv*dzdv+dxdv*dxdv)/sqrt(dx*dx+dz*dz);
		
		/* determine source amplitude */
		amp = (dist+distprev)/2.0;
		
		/* let source contribute within limited distance */
		xvn = (xv-fx)/dx;
		zvn = (zv-fz)/dz;
		ixv = NINT(xvn); 
		izv = NINT(zvn);
		for (ix=MAX(0,ixv-3); ix<=MIN(nx-1,ixv+3); ++ix) {
			for (iz=MAX(0,izv-3); iz<=MIN(nz-1,izv+3); ++iz) {
				xn = ix-xvn;
				zn = iz-zvn;
				s[ix][iz] += ts*amp*exp(-xn*xn-zn*zn);
			}
		}
	}
}

/* prototype of subroutine used internally */

static void pml_absorb (FdmodPml *P, int nx, float dx, int nz, float dz, float dt,
        float **dvv, float **od, float **pm, float **p, float **pp,
        int *abs)
/**************************************************************************
  pml_absorb - uses the perfectly matched layer absorbing boundary condition.
***************************************************************************
Notes:
   The PML formulation is specialized to the acoustic case.


   Array        Size        Location of [0][0] with respect to p [0][0]
   -----   ---------------  -------------------------------------------
   P->ux_b    (nx+pml, pml+2)  (0, nz-1)
   P->uz_b         "               "
   P->dax_b        "               "
   P->dbx_b        "               "
   P->daz_b        "               "
   P->dbz_b        "               "
 
   P->ux_r    (pml+2, nz)      (nx-1, 0)
   P->uz_r         "               "
   P->dax_r        "               "
   P->dbx_r        "               "
   P->daz_r        "               "
   P->dbz_r        "               "

   P->v_b     (nx+pml, pml+3)  (0, nz-1.5)
   P->cax_b        "               "
   P->cbx_b        "               "

   P->w_b     (nx+pml, pml+2)  (-0.5, nz-1)
   P->caz_b        "               "
   P->cbz_b        "               "

   P->v_r     (pml+2, nz)      (nx-1, -0.5)
   P->cax_r        "               "
   P->cbx_r        "               "

   P->w_r     (pml+3, nz)      (nx-1.5, 0)
   P->caz_r        "               "
   P->cbz_r        "               "
***************************************************************************
References:
   Jean-Pierre Berenger, ``A Perfectly Matched Layer for the Absorption of
   Electromagnetic Waves,''  Journal of Computational Physics, vol. 114,
   pp. 185-200.

   Hastings, Schneider, and Broschat, ``Application of the perfectly
   matched layer (PML) absorbing boundary condition to elastic wave
   propogation,''  Journal of the Accoustical Society of America,
   November, 1996.

   Allen Taflove, ``Electromagnetic Modeling:  Finite Difference Time
   Domain Methods'', Baltimore, Maryland: Johns Hopkins University Press,
   1995, chap. 7, pp. 181-195.

   The PML ABC is implemented by extending the modeled region on
   the bottom and right sides and treating the modeled region as
   periodic.

   In the extended region, the differential equations of PML are
   modeled.  The extension is accomplished by using additional
   arrays which record the state in the extended regions.  The
   result is a nasty patchwork of arrays.  (It is possible to use
   the PML differential equations to model the absorbing and
   non-absorbing regions.  This greatly simplifies things at the
   expense of memory.)

   The size of the new arrays and the location of their (0,0)
   element in the coordinate space of the main p arrays are:
*************************************************************************
   Author:  TAMU: Michael Holzrichter, 1998
*************************************************************************/
{
        int ix, iz, jx, jz;

   /* Calculate v for bottom pad above and below main domain */

   for (ix=0, jz=P->pml_thickness+2; ix<nx; ++ix) {
      P->v_b[ix][ 0] = P->cax_b [ix][ 0] * P->v_b [ix][ 0] +
                      P->cbx_b [ix][ 0] * (P->ux_b [ix][   0] + P->uz_b [ix][   0]
                                       -((abs[2]!=0) ? p [ix][nz-2] : 0.0));

      P->v_b[ix][jz] = P->cax_b [ix][jz] * P->v_b [ix][jz] +
                      P->cbx_b [ix][jz] * (((abs[0]!=0) ? p [ix][   1] : 0.0)
                                       -P->ux_b [ix][jz-1] - P->uz_b [ix][jz-1]);
   }	

   /* Calculate v for bottom pad above and below right pad */

   for (ix=nx, jx=1, jz=P->pml_thickness+2; ix<nx+P->pml_thickness; ++ix, ++jx) {
      P->v_b  [ix][ 0] = P->cax_b [ix][ 0] * P->v_b [ix][ 0] +
                      P->cbx_b [ix][ 0] * (P->ux_b [ix][   0] + P->uz_b [ix][   0]
                                       -P->ux_r [jx][nz-2] - P->uz_r [jx][nz-2]);

      P->v_b  [ix][jz] = P->cax_b [ix][jz] * P->v_b [ix][jz] +
                      P->cbx_b [ix][jz] * (P->ux_r [jx][   1] + P->uz_r [jx][   1]
                                       -P->ux_b [ix][jz-1] - P->uz_b [ix][jz-1]);
   }


   /* Calculate v for main part of bottom pad */

   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=1; iz<P->pml_thickness+2; ++iz) {
         P->v_b [ix][iz] = P->cax_b [ix][iz] * P->v_b [ix][iz] +
                        P->cbx_b [ix][iz] * (P->ux_b [ix][iz  ] + P->uz_b [ix][iz  ]
                                         -P->ux_b [ix][iz-1] - P->uz_b [ix][iz-1]);
      }
   }


   /* Calculate w for left edge of bottom pad */

   for (iz=0, ix=nx+P->pml_thickness-1; iz<P->pml_thickness+2; ++iz) {
      P->w_b [ 0][iz] = P->caz_b [ 0][iz] * P->w_b [ 0][iz] +
                     P->cbz_b [ 0][iz] * (P->ux_b [ix][iz] + P->uz_b [ix][iz]
                                      -P->ux_b [ 0][iz] - P->uz_b [ 0][iz]);
   }


   /* Calculate w for main part of bottom pad */

   for (ix=1; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {
         P->w_b [ix][iz] = P->caz_b [ix][iz] * P->w_b [ix][iz] +
                        P->cbz_b [ix][iz] * (P->ux_b [ix-1][iz] + P->uz_b [ix-1][iz]
                                         -P->ux_b [ix  ][iz] - P->uz_b [ix  ][iz]);
      }
   }


   /* Calculate v along top and bottom edge of right pad */

   for (ix=0, jx=nx-1, jz=P->pml_thickness; ix<P->pml_thickness+2; ++ix, ++jx) {
      if (jx == nx+P->pml_thickness) jx = 0;

      P->v_r [ix][   0] = P->cax_r [ix][   0] * P->v_r [ix][   0] +
                       P->cbx_r [ix][   0] * (P->ux_r [ix][ 0] + P->uz_r [ix][ 0]
                                          -P->ux_b [jx][jz] - P->uz_b [jx][jz]);

      P->v_r [ix][nz-1] = P->cax_r [ix][nz-1] * P->v_r [ix][nz-1] +
                       P->cbx_r [ix][nz-1] * (P->ux_b [jx][   0] + P->uz_b [jx][   0]
                                          -P->ux_r [ix][nz-2] - P->uz_r [ix][nz-2]);
   }


   /* Calculate v in rest of right pad */

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=1; iz<nz-1; ++iz) {
         P->v_r [ix][iz] = P->cax_r [ix][iz] * P->v_r [ix][iz] +
                        P->cbx_r [ix][iz] * (P->ux_r [ix][iz  ] + P->uz_r [ix][iz  ]
                                         -P->ux_r [ix][iz-1] - P->uz_r [ix][iz-1]);
      }
   }


   /* Calculate w along left and right sides of right pad */

   for (iz=0, jx=P->pml_thickness+2; iz<nz; ++iz) {
      P->w_r [ 0][iz] = P->caz_r [ 0][iz] * P->w_r [ 0][iz] +
                     P->cbz_r [ 0][iz] * (((abs[3]!=0) ? p [nx-2][iz] : 0.0)
                                      -P->ux_r [   0][iz] - P->uz_r [   0][iz]);

      P->w_r [jx][iz] = P->caz_r [jx][iz] * P->w_r [jx][iz] +
                     P->cbz_r [jx][iz] * (P->ux_r [jx-1][iz] + P->uz_r [jx-1][iz]
                                      -((abs[1]!=0) ? p [1][iz] : 0.0));
   }


   /* Calculate w in main part of right pad */

   for (ix=1; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz; ++iz) {
         P->w_r [ix][iz] = P->caz_r [ix][iz] * P->w_r [ix][iz] +
                        P->cbz_r [ix][iz] * (P->ux_r [ix-1][iz] + P->uz_r [ix-1][iz]
                                         -P->ux_r [ix  ][iz] - P->uz_r [ix  ][iz]);
      }
   }


   /* Calculate ux and uz in bottom pad */

   for (ix=0; ix<nx+P->pml_thickness-1; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {
         P->ux_b [ix][iz] = P->dax_b [ix][iz] * P->ux_b [ix  ][iz] +
                         P->dbx_b [ix][iz] * (P->w_b [ix][iz  ] - P->w_b [ix+1][iz]);
      }
   }

   for (ix=nx+P->pml_thickness-1, iz=0; iz<P->pml_thickness+2; ++iz) {
      P->ux_b [ix][iz] = P->dax_b [ix][iz] * P->ux_b [ix  ][iz] + 
                      P->dbx_b [ix][iz] * (P->w_b [ix][iz  ] - P->w_b [   0][iz]);
   }


   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {
         P->uz_b [ix][iz] = P->daz_b [ix][iz] * P->uz_b [ix][iz  ] +
                         P->dbz_b [ix][iz] * (P->v_b [ix][iz+1] - P->v_b [ix][iz  ]);
      }
   }


   /* Calculate ux and uz in right pad */

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz; ++iz) {
         P->ux_r [ix][iz] = P->dax_r [ix][iz] * P->ux_r [ix  ][iz] +
                         P->dbx_r [ix][iz] * (P->w_r [ix  ][iz] - P->w_r [ix+1][iz]);
      }
   }

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz-1; ++iz) {
         P->uz_r [ix][iz] = P->daz_r [ix][iz] * P->uz_r [ix][iz  ] +
                         P->dbz_r [ix][iz] * (P->v_r [ix][iz+1] - P->v_r [ix][iz  ]);
      }
   }

   for (ix=0, iz=nz-1, jx=nx-1; ix<P->pml_thickness+2; ++ix, ++jx) {
      if (jx == nx+P->pml_thickness) jx = 0;

      P->uz_r [ix][ 0] = P->uz_b [jx][P->pml_thickness+1];
      P->uz_r [ix][iz] = P->uz_b [jx][              0];
   }


   /* Update top and bottom edge of main grid with new field values */

   for (ix=0, jz=P->pml_thickness+1; ix<nx; ++ix) {
      if (abs [0] != 0) pp [ix][0] = P->ux_b [ix][jz] + P->uz_b [ix][jz];
      if (abs [2] != 0) pp [ix][nz-1] = P->ux_b [ix][ 0] + P->uz_b [ix][ 0];
   }


   /* Update left and right edges of main grid with new field values */

   for (iz=1, jx=P->pml_thickness+1; iz<nz-1; ++iz) {
      if (abs [1] != 0) pp [0][iz] = P->ux_r [jx][iz] + P->uz_r [jx][iz];
      if (abs [3] != 0) pp [nx-1][iz] = P->ux_r [ 0][iz] + P->uz_r [ 0][iz];
   }
}


static void pml_init (FdmodPml *P, int nx, int nz, float dx, float dz, float dt,
		      float **dvv, float **od)
{
   int ix, iz;

   /* Allocate arrays for pad on right */

   P->cax_r = zalloc2 (nz, P->pml_thickness+2);
   P->cbx_r = zalloc2 (nz, P->pml_thickness+2);
   P->caz_r = zalloc2 (nz, P->pml_thickness+3);
   P->cbz_r = zalloc2 (nz, P->pml_thickness+3);
   P->dax_r = zalloc2 (nz, P->pml_thickness+2);
   P->dbx_r = zalloc2 (nz, P->pml_thickness+2);
   P->daz_r = zalloc2 (nz, P->pml_thickness+2);
   P->dbz_r = zalloc2 (nz, P->pml_thickness+2);

   P->ux_r  = zalloc2 (nz, P->pml_thickness+2);
   P->uz_r  = zalloc2 (nz, P->pml_thickness+2);
   P->v_r   = zalloc2 (nz, P->pml_thickness+2);
   P->w_r   = zalloc2 (nz, P->pml_thickness+3);


   /* Zero out arrays for pad on right */

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz; ++iz) {
         P->ux_r  [ix][iz] = P->uz_r  [ix][iz] = 0.0;
         P->v_r   [ix][iz] = P->w_r   [ix][iz] = 0.0;

         P->cax_r [ix][iz] = P->cbx_r [ix][iz] = 0.0;
         P->caz_r [ix][iz] = P->cbz_r [ix][iz] = 0.0;
         P->dax_r [ix][iz] = P->dbx_r [ix][iz] = 0.0;
         P->daz_r [ix][iz] = P->dbz_r [ix][iz] = 0.0;
      }
   }


   /* Zero out extra bit on right pad */

   for (ix=P->pml_thickness+2, iz=0; iz<nz; ++iz) {
      P->caz_r [ix][iz] = P->cbz_r [ix][iz] = 0.0;
      P->w_r   [ix][iz] = 0.0;
   }


   /* Allocate arrays for pad on bottom */

   P->cax_b = zalloc2 (P->pml_thickness+3, nx + P->pml_thickness);
   P->cbx_b = zalloc2 (P->pml_thickness+3, nx + P->pml_thickness);
   P->caz_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->cbz_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->dax_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->dbx_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->daz_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->dbz_b = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);

   P->ux_b  = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->uz_b  = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
   P->v_b   = zalloc2 (P->pml_thickness+3, nx + P->pml_thickness);
   P->w_b   = zalloc2 (P->pml_thickness+2, nx + P->pml_thickness);
	

   /* Zero out arrays for pad on bottom */

   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {
         P->ux_b  [ix][iz] = P->uz_b  [ix][iz] = 0.0;
         P->v_b   [ix][iz] = P->w_b   [ix][iz] = 0.0;

         P->cax_b [ix][iz] = P->cbx_b [ix][iz] = 0.0;
         P->caz_b [ix][iz] = P->cbz_b [ix][iz] = 0.0;
         P->dax_b [ix][iz] = P->dbx_b [ix][iz] = 0.0;
         P->daz_b [ix][iz] = P->dbz_b [ix][iz] = 0.0;
      }
   }


   /* Zero out extra bit on bottom pad */

   for (ix=0, iz=P->pml_thickness+2; ix<nx+P->pml_thickness; ++ix) {
      P->cax_b [ix][iz] = P->cbx_b [ix][iz] = 0.0;
      P->v_b   [ix][iz] = 0.0;
   }


   /* Initialize cax & cbx arrays */

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz; ++iz) {
         P->sigma_ez = 0.0;

         P->cax_r [ix][iz] = (2.0 - (P->sigma_ez * dt)) / (2.0 + (P->sigma_ez * dt));
         P->cbx_r [ix][iz] = (2.0 * dt / dz)         / (2.0 + (P->sigma_ez * dt));
      }
   }

   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+3; ++iz) {
         if ((iz == 0) || (iz == P->pml_thickness + 2)) {
            P->sigma_ez = 0.0;
         } else {
            P->sigma_ez = P->pml_max * 0.5 * (1.0 - cos (2*PI*(iz-0.5)/(P->pml_thickness+1)));
         }

         P->cax_b [ix][iz] = (2.0 - (P->sigma_ez * dt)) / (2.0 + (P->sigma_ez * dt));
         P->cbx_b [ix][iz] = (2.0 * dt / dz)         / (2.0 + (P->sigma_ez * dt));
      }
   }


   /* Initialize caz & cbz arrays */

   for (ix=0; ix<P->pml_thickness+3; ++ix) {
      for (iz=0; iz<nz; ++iz) {
         if ((ix == 0) || (ix == P->pml_thickness+2)) {
            P->sigma_ex = 0.0;
         } else {
            P->sigma_ex = P->pml_max * 0.5 * (1.0 - cos (2*PI*(ix-0.5)/(P->pml_thickness+1)));
         }

         P->caz_r [ix][iz] = (2.0 - (P->sigma_ex * dt)) / (2.0 + (P->sigma_ex * dt));
         P->cbz_r [ix][iz] = (2.0 * dt / dz)         / (2.0 + (P->sigma_ex * dt));
      }
   }


   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {
         if (ix == 0) {
            P->sigma_ex = P->pml_max * 0.5 * (1.0 - cos (2*PI*(0.5)/(P->pml_thickness+1)));
         } else if (ix < nx) {
            P->sigma_ex = 0.0;
         } else {
            P->sigma_ex = P->pml_max * 0.5 * (1.0 - cos (2*PI*(ix-nx+0.5)/(P->pml_thickness+1)));
         }

         P->caz_b [ix][iz] = (2.0 - (P->sigma_ex * dt)) / (2.0 + (P->sigma_ex * dt));
         P->cbz_b [ix][iz] = (2.0 * dt / dz)         / (2.0 + (P->sigma_ex * dt));
      }
   }


   /* Initialize right pad's dax dbx, daz, & dbz arrays */

   for (ix=0; ix<P->pml_thickness+2; ++ix) {
      for (iz=0; iz<nz; ++iz) {

         /* Determine P->sigma_mx and P->sigma_mz */

         if ((ix == 0) || (ix == P->pml_thickness+1)) {
            P->sigma_mx = 0.0;
         } else {
            P->sigma_mx = P->pml_max * 0.5 * (1.0 - cos (2*PI*(ix)/(P->pml_thickness+1)));
         }

         P->sigma_mz = 0.0;


         /* Determine velocity, interpolate */

	 if (od!=NULL) {
	    /*	    if (verbose) warn("initializing right pml with density");*/
	    if (ix == 0) {
	       P->dvv_0 = sqrt (dvv [nx-1][iz] * od [nx-1][iz]);
	    } else if (ix == P->pml_thickness+1) {
	       P->dvv_0 = sqrt (dvv [   0][iz] * od [   0][iz]);
	    } else {
	       P->dvv_0 = sqrt (dvv [nx-1][iz] * od [nx-1][iz]);
	       P->dvv_1 = sqrt (dvv [   0][iz] * od [   0][iz]);

	       P->dvv_0 = ((ix) * P->dvv_1 + (1+P->pml_thickness-ix)*P->dvv_0);
	       P->dvv_0 /= (P->pml_thickness+1);
	    }
	 }
	 else {
	    if (ix == 0) {
	       P->dvv_0 = sqrt (dvv [nx-1][iz]);
	    } else if (ix == P->pml_thickness+1) {
	       P->dvv_0 = sqrt (dvv [   0][iz]);
	    } else {
	       P->dvv_0 = sqrt (dvv [nx-1][iz]);
	       P->dvv_1 = sqrt (dvv [   0][iz]);

	       P->dvv_0 = ((ix) * P->dvv_1 + (1+P->pml_thickness-ix)*P->dvv_0);
	       P->dvv_0 /= (P->pml_thickness+1);
	    }
	 }

         P->dvv_0 = P->dvv_0 * P->dvv_0;


         P->dax_r [ix][iz] = (2.0 - (P->sigma_mx * dt)) / (2.0 + (P->sigma_mx * dt));
         P->dbx_r [ix][iz] = (2.0 * dt * P->dvv_0 / dx) / (2.0 + (P->sigma_mx * dt));

         P->daz_r [ix][iz] = (2.0 - (P->sigma_mz * dt)) / (2.0 + (P->sigma_mz * dt));
         P->dbz_r [ix][iz] = (2.0 * dt * P->dvv_0 / dz) / (2.0 + (P->sigma_mz * dt));
      }
   }


   /* Initialize bottom pad's dax, dbx, daz, & dbz arrays */

   for (ix=0; ix<nx+P->pml_thickness; ++ix) {
      for (iz=0; iz<P->pml_thickness+2; ++iz) {

         /* Determine P->sigma_mx and P->sigma_mz */

         if (ix < nx) {
            P->sigma_mx = 0.0;
         } else {
            P->sigma_mx = P->pml_max * 0.5 * (1.0 - cos (2*PI*(ix-nx+1)/(P->pml_thickness+1)));
         }

         if ((iz == 0) || (iz == P->pml_thickness+1)) {
            P->sigma_mz = 0.0;
         } else {
            P->sigma_mz = P->pml_max * 0.5 * (1.0 - cos (2*PI*(iz)/(P->pml_thickness+1)));
         }


         /* Determine velocity, interpolate */
	 if (od!=NULL) {
	    /*	    if (verbose) warn("initializing bottom pml with density");*/
	    if (ix < nx) {
	       if (iz == 0) {
		  P->dvv_0 = sqrt (dvv [ix][nz-1] * od [ix][nz-1]);
	       } else if (iz == P->pml_thickness+1) {
		  P->dvv_0 = sqrt (dvv [ix][   0] * od [ix][   0]);
	       } else {
		  P->dvv_0 = sqrt (dvv [ix][nz-1] * od [ix][nz-1]);
		  P->dvv_1 = sqrt (dvv [ix][   0] * od [ix][   0]);

		  P->dvv_0 = ((iz) * P->dvv_1 + (1+P->pml_thickness-iz)*P->dvv_0);
		  P->dvv_0 /= (P->pml_thickness+1);
	       }
	    } else {
	       if (iz == 0) {
		  P->dvv_0 = sqrt (dvv [nx-1][nz-1] * od [nx-1][nz-1]);
		  P->dvv_1 = sqrt (dvv [   0][nz-1] * od [   0][nz-1]);
	       } else if (iz == P->pml_thickness+1) {
		  P->dvv_0 = sqrt (dvv [nx-1][   0] * od [nx-1][   0]);
		  P->dvv_1 = sqrt (dvv [   0][   0] * od [   0][   0]);
	       } else {
		  P->dvv_2 = sqrt (dvv [nx-1][nz-1] * od [nx-1][nz-1]);
		  P->dvv_3 = sqrt (dvv [nx-1][   0] * od [nx-1][   0]);

		  P->dvv_0 = ((iz) * P->dvv_3 + (1+P->pml_thickness-iz)*P->dvv_2);
		  P->dvv_0 /= (P->pml_thickness+1);

		  P->dvv_2 = sqrt (dvv [0][nz-1] * od [0][nz-1]);
		  P->dvv_3 = sqrt (dvv [0][0] * od [0][0]);

		  P->dvv_1 = ((iz) * P->dvv_3 + (1+P->pml_thickness-iz)*P->dvv_2);
		  P->dvv_1 /= (P->pml_thickness+1);
	       }

	       P->dvv_0 = ((ix-nx+1) * P->dvv_1 + (nx+P->pml_thickness-ix)*P->dvv_0);
	       P->dvv_0 /= (P->pml_thickness+1);
	    }
	 }

	 else {
	    if (ix < nx) {
	       if (iz == 0) {
		  P->dvv_0 = sqrt (dvv [ix][nz-1]);
	       } else if (iz == P->pml_thickness+1) {
		  P->dvv_0 = sqrt (dvv [ix][   0]);
	       } else {
		  P->dvv_0 = sqrt (dvv [ix][nz-1]);
		  P->dvv_1 = sqrt (dvv [ix][   0]);

		  P->dvv_0 = ((iz) * P->dvv_1 + (1+P->pml_thickness-iz)*P->dvv_0);
		  P->dvv_0 /= (P->pml_thickness+1);
	       }
	    } else {
	       if (iz == 0) {
		  P->dvv_0 = sqrt (dvv [nx-1][nz-1]);
		  P->dvv_1 = sqrt (dvv [   0][nz-1]);
	       } else if (iz == P->pml_thickness+1) {
		  P->dvv_0 = sqrt (dvv [nx-1][   0]);
		  P->dvv_1 = sqrt (dvv [   0][   0]);
	       } else {
		  P->dvv_2 = sqrt (dvv [nx-1][nz-1]);
		  P->dvv_3 = sqrt (dvv [nx-1][   0]);

		  P->dvv_0 = ((iz) * P->dvv_3 + (1+P->pml_thickness-iz)*P->dvv_2);
		  P->dvv_0 /= (P->pml_thickness+1);

		  P->dvv_2 = sqrt (dvv [0][nz-1]);
		  P->dvv_3 = sqrt (dvv [0][0]);

		  P->dvv_1 = ((iz) * P->dvv_3 + (1+P->pml_thickness-iz)*P->dvv_2);
		  P->dvv_1 /= (P->pml_thickness+1);
	       }

	       P->dvv_0 = ((ix-nx+1) * P->dvv_1 + (nx+P->pml_thickness-ix)*P->dvv_0);
	       P->dvv_0 /= (P->pml_thickness+1);
	    }
	 }

         P->dvv_0 = P->dvv_0 * P->dvv_0;

         P->dax_b [ix][iz] = (2.0 - (P->sigma_mx * dt)) / (2.0 + (P->sigma_mx * dt));
         P->dbx_b [ix][iz] = (2.0 * dt * P->dvv_0 / dx) / (2.0 + (P->sigma_mx * dt));

         P->daz_b [ix][iz] = (2.0 - (P->sigma_mz * dt)) / (2.0 + (P->sigma_mz * dt));
         P->dbz_b [ix][iz] = (2.0 * dt * P->dvv_0 / dz) / (2.0 + (P->sigma_mz * dt));
      }
   }
}


/* Library version of SUFDMOD2_PML
 *
 * The main program is left to the caller, as for SUFDMOD2 (su_fdmod2_*), with the same point source (su_fdmod2_ptsrc, with a
 * strength of 1 and the ricker wavelet), and the same finite-difference star. What is here:
 *
 *	su_fdmod2pml_exsrc:	the source pressure of an extended source (a gaussian derivative wavelet; the splines are made by
 *				su_fdmod2_exsrc_setup)
 *	su_fdmod2_pml_init:	the coefficients and the state of the perfectly matched layer (pml_thick > 0), for a velocity model
 *	su_fdmod2_pml_tstep:	one time step, with the PML as the absorbing boundary condition
 *	su_fdmod2_pml_free:	free the layer
 *
 * The global variables of the program (the coefficients and the state of the layer) are in the structure that su_fdmod2_pml_init
 * makes, so that models do not share them. The arrays of the layer start at 0 (the program did not set all of them).
 * The 2-D arrays are flat, with z the fast axis: [nx][nz].
 */

typedef struct FdmodPml FdmodPml;

static float **zalloc2 (size_t n1, size_t n2)
{
	float **a = alloc2float(n1,n2);
	memset(a[0],0,sizeof(float)*n1*n2);
	return a;
}

static float **pml_rows(const float *base, int n, int m)
{
	int i;
	float **rows = (float**)ealloc1(n,sizeof(float*));
	for (i=0; i<n; ++i) rows[i] = (float*)base+(size_t)i*m;
	return rows;
}

void su_fdmod2pml_exsrc(int ns, const float *xs, const float *zs, const float *vs, const float *xsd, const float *zsd,
	int nx, float dx, float fx, int nz, float dz, float fz, float dt, float t, float fmax, float *s)
{
	float **rows = pml_rows(s,nx,nz);
	exsrc(ns,(float*)xs,(float*)zs,(float*)vs,(float(*)[4])xsd,(float(*)[4])zsd,nx,dx,fx,nz,dz,fz,dt,t,fmax,rows);
	free1(rows);
}

/* Returns the layer, or NULL for pml_thick<=0. od is NULL for a constant density of 1 */
FdmodPml *su_fdmod2_pml_init(int nx, int nz, float dx, float dz, float dt, int pml_thick, float pml_max,
	const float *dvv, const float *od)
{
	FdmodPml *P;
	float **rdvv, **rod = NULL;
	if (pml_thick<=0) return NULL;
	P = (FdmodPml*)calloc(1,sizeof(FdmodPml));
	P->pml_thick = pml_thick;
	P->pml_thickness = 2*pml_thick;
	P->pml_max = pml_max;
	rdvv = pml_rows(dvv,nx,nz);
	if (od!=NULL) rod = pml_rows(od,nx,nz);
	pml_init(P,nx,nz,dx,dz,dt,rdvv,rod);
	free1(rdvv);
	if (rod!=NULL) free1(rod);
	return P;
}

void su_fdmod2_pml_free(FdmodPml *P)
{
	if (P==NULL) return;
	free2float(P->cax_b); free2float(P->cax_r); free2float(P->cbx_b); free2float(P->cbx_r);
	free2float(P->caz_b); free2float(P->caz_r); free2float(P->cbz_b); free2float(P->cbz_r);
	free2float(P->dax_b); free2float(P->dax_r); free2float(P->dbx_b); free2float(P->dbx_r);
	free2float(P->daz_b); free2float(P->daz_r); free2float(P->dbz_b); free2float(P->dbz_r);
	free2float(P->ux_b); free2float(P->ux_r); free2float(P->uz_b); free2float(P->uz_r);
	free2float(P->v_b); free2float(P->v_r); free2float(P->w_b); free2float(P->w_r);
	free(P);
}

void su_fdmod2_pml_tstep(FdmodPml *P, int nx, float dx, int nz, float dz, float dt, const float *dvv, const float *od,
	const float *s, const float *pm, const float *p, float *pp, const int *abs)
{
	float **rdvv, **rod = NULL, **rpm, **rp, **rpp;
	su_fdmod2_star(nx,dx,nz,dz,dt,dvv,od,s,pm,p,pp);
	rdvv = pml_rows(dvv,nx,nz);
	if (od!=NULL) rod = pml_rows(od,nx,nz);
	rpm = pml_rows(pm,nx,nz); rp = pml_rows(p,nx,nz); rpp = pml_rows(pp,nx,nz);
	pml_absorb(P,nx,dx,nz,dz,dt,rdvv,rod,rpm,rp,rpp,(int*)abs);
	free1(rdvv); if (rod!=NULL) free1(rod); free1(rpm); free1(rp); free1(rpp);
}
