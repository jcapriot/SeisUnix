/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* su.h - include file for SU programs
 *
 * $Author: john $
 * $Source: /usr/local/cwp/src/su/include/RCS/su.h,v $
 * $Revision: 1.34 $ ; $Date: 2011/11/11 23:56:14 $
 */

#ifndef SU_H
#define SU_H

#include "par.h"

/* TYPEDEFS */
typedef union { /* storage for arbitrary type */
	char s[8];
	short h;
	unsigned short u;
	long l;
	unsigned long v;
	int i;
	unsigned int p;
	float f;
	double d;
	unsigned int U:16;
	unsigned int P:32;
} Value;


/* DEFINES */
#define CHECK_NT(label,nt) \
	if(nt > SU_NFLTS) err("%s=%d must not exceed %d",label,nt,SU_NFLTS)
#define NALLOC	(524288)
#define NFALLOC	(NALLOC/FSIZE)
#define NIALLOC	(NALLOC/ISIZE)
#define NDALLOC	(NALLOC/DSIZE)
#define LOWBYTE(w) ((w) & 0xFF)
#define HIGHBYTE(w) LOWBYTE((w) >>8)
#define LOWWORD(w) ((w) & 0xFFFF)
#define HIGHWORD(w) LOWWORD((w) >>16)
#define ISNEGCHAR(c) ((c) & 0x80)
#define SIGNEXTEND(c) (~0xFF | (int) (c))

/*	READ_OK  - read  permission for access(2)
 *	WRITE_OK - write permission for access(2)
 *	EXEC_OK  - exec  permission for access(2)
 *	FILE_OK  - file  existence  for access(2)
 *	Note: these are changed from the usual defines in file.h
 *	      because this include exists on some machines and
 *	      not others, often overlaps fcntl.h, etc.  Lint is
 *            happier with a fresh start.
 *	Note: Post-ANSI sometimes R_OK in unistd.h (this isn't
 *	      an ANSI file).
 */
#define		READ_OK		4
#define		WRITE_OK	2
#define		EXEC_OK		1
#define		FILE_OK		0

/* For plotting by keyword */
#define	 IS_DEPTH(str)	((  STREQ(str,"gelev")	|| \
			   STREQ(str,"selev")	|| \
			   STREQ(str,"sdepth")	|| \
			   STREQ(str,"gdel")	|| \
			   STREQ(str,"sdel")	|| \
			   STREQ(str,"swdep")	|| \
			   STREQ(str,"gwdep")	      )?cwp_true:cwp_false)

#define	 IS_COORD(str)	(( STREQ(str,"sx")   || \
			   STREQ(str,"sy")   || \
			   STREQ(str,"gx")   || \
			   STREQ(str,"gy")	   )?cwp_true:cwp_false)

/* FUNCTION PROTOTYPES */

/* valpkge */
int vtoi(cwp_String type, Value val);
long vtol(cwp_String type, Value val);
float vtof(cwp_String type, Value val);
double vtod(cwp_String type, Value val);
int valcmp(cwp_String type, Value val1, Value val2);
void printfval(cwp_String type, Value val);
void fprintfval(FILE *stream, cwp_String type, Value val);
void scanfval(cwp_String type, Value *valp);
void atoval(cwp_String type, cwp_String keyval, Value *valp);
void getparval(cwp_String name, cwp_String type, int n, Value *valp);
Value valtoabs(cwp_String type, Value val);

/* segy coordinate scalar utilities */
short elco_scalar(int ncoords, double c[]);
double from_segy_elco_multiplier(short segy_scalar);
double to_segy_elco_multiplier(short segy_scalar); /* reciprocal of from */


//BEGIN SU MAIN functions

// The programs of su/main are in the library as functions that do what the program does to a trace (or a gather)
// and that have no static state: the lookup tables and scratch arrays are arguments that the caller owns.

// filters
void su_bfhighpass(int zerophase, int npoles, float f3db, size_t nt, float *data_in, float *data_out);
void su_bflowpass(int zerophase, int npoles, float f3db, size_t nt, float *data_in, float *data_out);
void polygonalFilter(float *f, float *amps, int npoly, int nfft, float dt, float *filter, int *intfr);
void su_median_across(int n, int nwin, const float *rows, size_t stride, float *out, float *scratch);
void su_mix_across(int n, int nwin, const float *rows, size_t stride, const float *w, float *out);
void su_frac_filter(int nf, int nfft, float dt, float power, float sign, float phasefac, float scale, float *filt);
void su_phase_spectrum(int nf, float *ct, float a, float b, float c);
float su_kolmogoroff_log(int n, float pnoise, float *cx);
void su_kolmogoroff_fold(int n, float *cx);
void su_kolmogoroff_exp(int n, float rmax, float *cx);
void su_tvband_filter(const float *f, int nfft, int nfreq, float dt, float scale, float *filter);
void su_tvband_blend(int i0, int i1, const float *first, const float *second, float *out);

// transforms: suradon
float su_radon_gofx(int igopt, float offset, float intercept_off, float refdepth);
void su_radon_taupmute(int ip, int ipa, int ipb, int nt, int ltap, float *rt);
int su_radon_forward(int nx, int nt, const float *data, const float *g, float dt, int ntfft, int np, float pmin,
	float dp, const float *mutetime, int nk, float f1, float f2, float prewhite, float *ptau);
void su_radon_inverse(int np, int ntfft, int nx, const float *ptau, const float *g, float dt, float pmin, float dp,
	int ip1, float f1, float f2, float *out);
void su_radon_xinterp(int nt, int nx, int ninterp, const float *data, float freq1, float freq2, int lagc,
	int lent, int lenx, int xopt, float dt, float *out);

// datuming
int su_datumfd(int nt, float dt, int nx, int nsx, int nz, float dz, float dx, float x_0, int buff, int tap_len, int mx,
	int mz, const float *vel1, const float *vel2, const float *gx, const float *in, float *out);

int su_datumk2d_geometry(int nu0, const float *x0, const float *z0, int nu1, const float *x1, const float *z1,
	float dsmax, int nxi, float fxi, float dxi, int ns, float fs, float ds, float *szif, float *sz, float *nangl);
void su_datumk2d_tables(int ns, int nr, int nzt, float dr, float dzt, float fzt, float v0, float dvz, const float *sz,
	int nxt, float fxt, float dxt, float fs, float ds, float *tb, float *pb, float *sigb, float *angb, float *ttab);
int su_datumk2d_trace(const float *trace, int nt, float ft, float dt, float *dat, float aperx, int nx, float fx, float dx,
	int nz, float fz, float dz, int mtmax, float xm, float fmax, int nxi, float fxi, float dxi, float angmax,
	int ns, float fs, float ds, int nr, const float *tb, const float *pb, const float *angb, const float *sigb,
	const float *ttab, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt, int antiali, int sgn,
	const float *szif, const float *nangl);
int su_kdmdcr_trace(const float *trace, int nt, float ft, float dt, float sx, float gx, float *dat, float aperx,
	int nx, float fxin, float dxgo, int nzo, float fzo, float dzo, int nxi, float fxi, float dxi, int ns, float fs,
	float ds, const float *ttab, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt, const float *szif,
	float v0, float freq);
int su_kdmdcs_trace(const float *trace, int nt, float ft, float dt, float sx, float gx, float *dat, float aperx,
	int nx, float fxin, float dxso, int nzo, float fzo, float dzo, int nxi, float fxi, float dxi, int ns, float fs,
	float ds, const float *ttab, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt, const float *szif,
	float v0, int dc, int ngl, float freq);

// multicomponent
void su_polar_station(int nt, int iwl, int iwin, float rlq, int rl, int theta, int phi, int tau, int ellip, int pln,
	int f1, int l1, int amp, int dir, float fangle, const float *data, float *out);
void su_eipofi_station(int nt, int iwl, float pwr, int cubic, float *data, float *polar);
float su_alford_angle(const float *S11, const float *S12, const float *S21, const float *S22, float dt, int itmin,
	int itmax, float ang_inc, int taper, float *quality, float *inputquality, float *lag, float *xcorr);
void su_alford_rotate(float *S11, float *S12, float *S21, float *S22, float theta, int nt, int itmin, int taper);

// interp_extrap
int su_interp_traces(int nt, int nx, int ninterp, const float *data, float freq1, float freq2, int lagc, int lent,
	int lenx, int xopt, float dt, int iopt, int deriv, int linear, float *out);
void su_ocext_offset(float offset, float offextr, float fmax, float sdmo, int nx, float dx, int nt, float dt, float ft,
	const float *vrms_in, float *ptx);

// statics
void su_resstat(int nt, int ntraces, const float *data_in, const int *icmp, const int *irec, const int *ishot,
	int ncmp, int nrec, int nshot, int ntpick, int niter, int mode, float perc, int sub, float *tsstat,
	float *trstat);

// decon_shaping
int su_fxdecon(int nsp, int ntrc, const float *tdatain, float *tdataout, float dt, float taper, float fmin,
	float fmax, float twlen, int ntrw, int ntrf);

// amplitudes
void su_centsamp(float *rt, float *ct, float *mt, const float *time, int nt, float dt, int nvals_min);
float su_dipdivcor_scale(void);
void su_dipdivcor_table(int nt, int np, float dt, const float *tt, const float *vt, float *vs, float (*vind)[4],
	float *divcor, int trans, int norm);
void su_dipdivcor_filter(float k, float dpx, float dt, int np, int nw, int nt, const float *div,
	const complex *p, complex *q, complex *kq, complex *qq);
void su_gain_tpow_table(float *tpowfac, int nt, float tmin, float dt, float tpow, float tred);
void su_gain_epow_table(float *epowfac, int nt, float tmin, float dt, float epow, float etpow);
void su_gain(float *data, float tpow, float epow, float gpow,
	int agc, int gagc, int qbal, int pbal, int mbal, float scale, float bias,
	float trap, float clip, float qclip, int iwagc,
	int nt, int maxbal, float pclip, float nclip,
	const float *tpowfac, const float *epowfac,
	float *absdata, float *agcdata, float *d2, float *w, float *s);
void su_pgc_gain(int nt, const float *sum, int icount, int lw, float *g);

// operations (suop)
void su_op_saf(float *data, int nt, float *tmp);
void su_op_freq(float *data, int nt, float dt, float *tmp, float *tmp1);
void su_op_despike(float *data, int nt, int nw, float *tmp, float *tomed);

// stretching_moveout_resamp
void su_lintrp(const float *q, float *w, int *it, int lp, int lq);
void su_stretch(float *q, const float *p, const float *w, const int *it, int lq, int nw);
void su_ttoz_table(int nt, float dt, float ft, const float *v, int nz, float dz, float fz, float *tz, float *z);
void su_ztot_table(int nz, float dz, float fz, const float *v, int nt, float dt, float ft, float *zt, float *t);
void su_nmo_tables(int nt, float dt, float ft, float offset, const float *ovvt, float smute, int upward,
	int invert, int sscale, float *ttn, float *atn, float *tnt, float *at, int *itmute_out);
int su_nmo_a_tables(int nt, float dt, float ft, float offset, const float *ovvt, const float *oa1t, const float *oa2t,
	float smute, int upward, int invert, int sscale, float *ttn, float *atn, float *tnt, float *at, int *itmute_out);
void su_nmo(float *data, int nt, float dt, float ft, int itmute, int lmute, int sscale, int invert,
	const float *ttn, const float *atn, const float *tnt, const float *at, float *q);
void su_taupnmo_tables(int nt, float dt, float ft, float p, const float *vvt, float smute, float *ttn, float *atn,
	int *itmute_out);
void su_taupnmo(float *data, int nt, float dt, float ft, int itmute, int lmute, int sscale,
	const float *ttn, const float *atn, float *q);

// velocity analysis
int su_velan_accumulate(int nt, float dt, float ft, float offset, int nv, float dv, float fv, float anis1,
	float anis2, float smute, const float *data, float *num, float *den, float *nnz);
void su_velan_semblance(int nt, int ntout, int dtratio, int nsmooth, float pwr, const float *num, const float *den,
	const float *nnz, float *sem);
void su_relan_accumulate(int nz, float dz, float fz, float offset, int nr, float dr, float fr, float smute,
	const float *data, float *num, float *den, float *nnz);
void su_tivel_table(float a, float b, float e, float d, float maxangle, int nangle, int normalize, float *p, float *v,
	float *dv, float *d2v, float *vnmo);
float su_tivel_vnmoweak(float a, float b, float e, float d, float p);
void su_vel2df_accumulate(int nt, float dt, float ft, float offset, float tn_zero, int nv, float dv, float fv, int nvh,
	float dvh, float fvh, int nsmooth, const float *data, float *num, float *den, float *nnz);
void su_relanan_accumulate(int nz, float dz, float fz, float half_offset, float Z, int nr1, float dr1, float fr1,
	int nr2, float dr2, float fr2, int nsmooth, const float *data, float *num, float *den, float *nnz);
void su_velan_cc(int nt, float dt, float ft, int nx, const float *offsets, const float *data, int nv, float dv, float fv,
	float smute, int dtratio, int nsmooth, int normalize, int selective, float tau, float pwr, float *sem);

// synthetics
int su_goupillaudpo_tmax(int n, int l, int k);
int su_goupillaudpo(int n, const float *r, int l, int k, int tmax, int pV, float *x, float *out, int *odd);
int su_goupillaud_tmax(int n, int l, int k);
int su_goupillaud(int n, const float *rin, int l, int k, int tmax, int pV, float *out, int *odd);
int su_imp_nfft(int nt);
void su_imp2d_trace(int nt, float dt, int nfft, float c, float rs, float rg, const float *tout, float *rt, complex *ct,
	float *data);
void su_imp3d_trace(int nt, float dt, int nfft, float c, float rs, float rg, int dir, float rd, const float *tout, float *temp,
	float *rt, complex *ct, float *data);
int su_syncz_tables(int ninf, int ntr, float dx, const float *zint, const float *dip, const float *v, const float *rho,
	float *dipr, float *xl, float *xr, float *theta, float *m, float *b, float *d, float *k, float *w, float *meet,
	float *trcoefs, int *dorow, int *info);
void su_syncz_trace(int ninf, float x, int nt, float dt, float tdelay, const float *v, const float *m, const float *b,
	const float *d, const float *k, const float *w, const float *meet, const float *trcoefs, const int *dorow,
	const float *tout, float *data, float *trace);

// tapering
float su_taper_envelope(int tap_type, float f, float min, float max);
void su_taper_time(float t1, float t2, int tap_type, float dt, float *trace, int nt);
void su_ramp(float *data, int nt, int ntaper1, int ntaper2);

// windowing_sorting_muting
void su_mute_taper(int ntaper, float *taper);
void su_mute_above(float *data, int nt, float t, float tmin, float dt, int ntaper, const float *taper);
void su_mute_below(float *data, int nt, float t, float tmin, float dt, int ntaper, const float *taper);
void su_mute_line(float *data, int nt, float t, float tmin, float dt, int ntaper, const float *taper,
	float fval, float linvel, float tm0);
void su_mute_hyperbola(float *data, int nt, float t, float tmin, float dt, int ntaper, const float *taper,
	float fval, float linvel, float tm0);
void su_mute_polygon(float *data, int nt, float t, float tw, float dt, int ntaper, const float *taper);

// noise
// random number generators that keep their state in a su_rng (see suaddnoise.c)
typedef struct {
	int i, j;
	float c;
	float u[17];
} su_rng;
void su_rng_seed_uniform(su_rng *r, int seed);
float su_rng_uniform(su_rng *r);
void su_rng_seed_normal(su_rng *r, int seed);
float su_rng_normal(su_rng *r);

// attributes_parameter_estimation
void su_attr_differentiate(int n, float h, float *f);
void su_attr_unwrap_phase(int n, float w, float *phase);
void su_attr_envelope(int n, const float *re, const float *im, float *out);
void su_attr_phase(int n, const float *re, const float *im, float unwrap, float *out);
void su_attr_freq(int n, const float *re, const float *im, float dt, float unwrap, float *out);
void su_attr_normamp(int n, const float *re, const float *im, float *out);
void su_attr_fdenv(int n, const float *re, const float *im, float dt, float *out);
void su_attr_sdenv(int n, const float *re, const float *im, float dt, float *out);
void su_attr_bandwidth(int n, const float *re, const float *im, float dt, float *out, float *scratch);
void su_attr_q(int n, const float *re, const float *im, float dt, float unwrap, float *out, float *scratch);

// stacking
void su_stack_add(int n, int ncomp, const float *x, float *sum, int *nnz);
void su_stack_normalize(int n, int ncomp, float *sum, const int *nnz, float normpow);
void su_divstack_power(int nt, int ncomp, int ntwin, int peak, const float *x, float *intp);
void su_divstack_add(int nt, int ncomp, const float *x, const float *intp, float *sumdata, float *sumscale);
void su_divstack_finish(int nt, int ncomp, const float *sumdata, const float *sumscale, float *out);
void su_pws_accumulate(int nt, const float *data, const float *hdata, float *stdata, float *psct);
void su_pws_weights(int nt, const float *psct, int ntr, float pwr, int isl, float *psdata, float *scratch);
void su_pws_apply(int nt, const float *psdata, int ntr, float *stdata);
void su_stackup_add(int n, const float *x, double *sum, float *samplefold);
void su_stackup_finish(int n, const double *sum, const float *samplefold, float *out);

// transforms
void su_st_analytic(int len, float *h);
void su_st_row(int len, int n, const float *h, float *g);
void su_gabor_filter(float fcent, float dt, int nfft, float alpha, float band, float scale, float *filter);
float su_cwt_wavelet(int nwavelet, float xmin, float xcenter, float xmax, float sigma, float *waveletsum);
int su_cwt_filter(int nwavelet, const float *waveletsum, float scale, float dx, float width, float *filt);
void su_cwt_trace(int ns, int nconv, const float *convolution, float scale, float *rt);
void su_wfft_flatten(int nf, const float *ct, float w0, float w1, float w2, float *out);
void su_clogfft_spectrum(int nf, const float *ct, float *log_amp, float *phase);
void su_iclogfft_spectrum(int nf, const float *clog, int sym, float *ct);

// convolution_correlation
void su_acorfrac_spectrum(int nf, float *ct, float a, float b, int sym);

// synthetics_waveforms_testpatterns
void su_vibro_linear(float *data, int nt, float fs, float fe, float T, float dt, float phz);
void su_vibro_segments(float *data, int nt, const float *freq, const float *time, int isegm, float T, float dt, float phz);
void su_vibro_octave(float *data, int nt, float fs, float fe, float T, float dt, float swconst, float phz);
void su_vibro_hertz(float *data, int nt, float fs, float fe, float T, float dt, float swconst, float phz);
void su_vibro_tpower(float *data, int nt, float fs, float fe, float T, float dt, float swconst, float phz);

// synthetics
void su_addsinc_table(void);
void su_synlv(float *data,
	float xs, float zs, float xg, float zg,
	size_t nt, float dt, float ft,
	float v00, float dvdx, float dvdz,
	int ls, int er, int ob, Wavelet *w, int nr, Reflector *r, int lhd, int nhd, float *hd
);


void su_synlvcw(float *data,
	float xs, float zs, float xg, float zg,
	size_t nt, float dt, float ft,
	float v00, float dvdx, float dvdz, float gamma,
	int ls, int er, int ob, int sp,
	Wavelet *w, int nr, Reflector *r, int lhd, int nhd, float *hd
);


void su_synlvfti(float *trace, float v00, float dvdx, float dvdz,
	int ls, int er, int ob, Wavelet *w, int trans,
	int nitmax, float epst, int zeroff,
	int ntries, float epsx, float angxs,
	float xs, float zs, float xg, float zg,
	float a, float f, float l,
	int nr, Reflector *r, int nt, float dt, float ft,
	int lhd, int nhd, float *hd
);


int su_synvxz_offset(const float *vel, int nx, int nz, float dx, float dz, float fx,
	int nxb, int nxd, float xo, int nxm, float dxm, float fxm,
	int ls, Wavelet *w, int nr, Reflector *r, int nt, float dt, float ft,
	int lhd, int nhd, float *hd, float *out
);


int su_synvxzcs_shot(const float *vel, const float *pert, int nx, int nz, float dx, float dz, float fx,
	int nxb, int nxd, int nxc, int nzc, float xs, int nxg, float dxg, float gx0,
	int ls, Wavelet *w, int nr, Reflector *r, int nt, float dt, float ft,
	int lhd, int nhd, float *hd, float *out
);


void su_kdsyn2d_reference(int nr, int nzt, float dr, float dzt, float fzt, float dvz, float v0,
	float *tb, float *pb, float *sigb, float *cosb);
void su_kdsyn2d_residual(int ns, float fs, float ds, int nxt, float fxt, float dxt, int nzt, int nr, float dr,
	const float *tb, float *ttab);
void su_kdsyn2d_integrate(const float *mig, int nz, float dz, int nx, int m, float *migi);
int su_kdsyn2d_trace(float *trace, int nt, float ft, float dt, float xs, float xg,
	const float *mig, const float *migi, float aperx,
	int nx, float fx, float dx, int nz, float fz, float dz,
	int mzmax, int ls, float angmax, float v0, float fmax, Wavelet *w,
	const float *tb, const float *pb, const float *sigb, const float *cosb, int nr,
	const float *ttab, int ns, float fs, float ds, float *tsum, float *tt,
	int nxt, float fxt, float dxt, int nzt, float fzt, float dzt,
	int lhd, int nhd, float *hd);


int su_fdmod1_plan(const float *rv, int nz, float dz, float tmax, int nt_in, int styp, float freq,
	float *dt, int *nt, float *t0, int *ies);
int su_fdmod1_run(const float *rv, const float *rd, int nz, float dz, float dt, int nt, float t0, int ies,
	int isz, int irz, const int *abs, int styp, float freq, int td, int zd, int press,
	float *sismo, float *snaps);


void su_fdmod2_exsrc_setup(int ns, const float *xs, const float *zs, float *vs, float *xsd, float *zsd);
void su_fdmod2_ptsrc(float sstrength, float xs, float zs, int nx, float dx, float fx, int nz, float dz, float fz,
	float dt, float t, float fmax, float fpeak, int mono, float *s);
void su_fdmod2_exsrc(int ns, const float *xs, const float *zs, const float *vs, const float *xsd, const float *zsd,
	int nx, float dx, float fx, int nz, float dz, float fz,
	float dt, float t, float fpeak, int pwt, int mono, float *s);
void su_fdmod2_tstep(int nx, float dx, int nz, float dz, float dt, const float *dvv, const float *od, const float *s,
	const float *pm, const float *p, float *pp, const int *abs);


void su_fdmod2_star(int nx, float dx, int nz, float dz, float dt, const float *dvv, const float *od, const float *s,
	const float *pm, const float *p, float *pp);
typedef struct FdmodPml FdmodPml;
void su_fdmod2pml_exsrc(int ns, const float *xs, const float *zs, const float *vs, const float *xsd, const float *zsd,
	int nx, float dx, float fx, int nz, float dz, float fz, float dt, float t, float fmax, float *s);
FdmodPml *su_fdmod2_pml_init(int nx, int nz, float dx, float dz, float dt, int pml_thick, float pml_max,
	const float *dvv, const float *od);
void su_fdmod2_pml_free(FdmodPml *P);
void su_fdmod2_pml_tstep(FdmodPml *P, int nx, float dx, int nz, float dz, float dt, const float *dvv, const float *od,
	const float *s, const float *pm, const float *p, float *pp, const int *abs);


struct Ea2dfModel;
struct Ea2dfModel *su_ea2df_create(float dt, float ft, float lt, int nx, float dx, float fx, int nz, float dz, float fz,
	float sx, float sz, const char *stype, float sang, const char *wtype, float ts, float favg,
	int qsw, int aniso, int tsw, const int *bcin, float bc_a, float bc_r, float hsz, float vsx,
	const float *c11in, const float *c55in, const float *rhoin, const float *c13in, const float *c33in,
	const float *c15in, const float *c35in, const float *qin, int *err);
void su_ea2df_info(const struct Ea2dfModel *M, int *nxpadded, int *nzpadded, int *nt, int *wbc, int *done);
void su_ea2df_lines(const struct Ea2dfModel *M, const float **hu, const float **hw, const float **vu, const float **vw);
void su_ea2df_snapshot(const struct Ea2dfModel *M, int which, float *out);
int su_ea2df_step(struct Ea2dfModel *M);
void su_ea2df_free(struct Ea2dfModel *M);


int su_remac2d(int opflag, int nx, int nz, int nt, float dx, float dz, float dt,
	int nsourc, const int *isx, const int *isz, const float *amps_in, float w,
	int sflag, float fmax, int nwav_in, float dtwave_in, const float *wave_in, int fsflag, float vmaxu,
	float dtsnap, int iabso, float abso, int nbwx, int nbwz,
	const float *vel, const float *dens, int nsectx, const int *irz, int nsectz, const int *irx,
	float *sectx, float *sectz, int nsnap, float *snap);


int su_remel2dan(int amode, int nx, int nz, int nt, float dx, float dz, float dt,
	int nsourc, const int *isx, const int *isz, const int *styp, const float *samp, float w,
	int sflag, float fmax, int nwav_in, float dtwave_in, const float *wave_in, float vmaxu, float vmax0, float vmin0,
	float dtsnap, int nsntyp, const int *sntyp, int nsnap, int iabso, float abso, int nbwx, int nbwz,
	const float *dens, const float *vp, const float *vs,
	const float *c11in, const float *c13in, const float *c15in, const float *c33in, const float *c35in, const float *c55in,
	int nsectx, const int *irz, const int *rxtyp, int nsectz, const int *irx, const int *rztyp,
	float *xsect, float *zsect, float *snap);


struct FctModel;
struct FctModel *su_fctanismod_create(int nx, int nz, int nt, float dx, float dz, float dt, int sx, int sz,
	int receiverdepth, int vspnx, int impulse, int source, int isurf, int dofct, int fctxbeg, int fctzbeg, int fctxend, int fctzend,
	int forcex, int forcey, int forcez, int wavelet, int movebc, float fpeak, float eta0, float eta,
	float deta0dx, float deta0dz, float detadx, float detadz,
	const float *aa, const float *cc, const float *ff, const float *ll, const float *nn, const float *rho,
	const float *xzsource, int *err);
float su_fctanismod_vmax(const struct FctModel *M);
int su_fctanismod_done(const struct FctModel *M);
int su_fctanismod_step(struct FctModel *M);
void su_fctanismod_snapshot(const struct FctModel *M, int which, float *out);
void su_fctanismod_records(const struct FctModel *M, const float **refl, const float **vsp);
void su_fctanismod_free(struct FctModel *M);



// migration
int su_stolt(int nx, int nt, float dt, float dx, int ntmig, const float *tmig_in, const float *vmig_in, float vscale,
	float smig, float fmax, int lstaper, int lbtaper, const float *gather, float *out);
int su_migfd(int nx, int nt, int nz, float dz, float dt, float dx, int dip, const float *pin, const float *vin, float *cresult_out);
int su_migffd(int nx, int nt, int nz, float dz, float dt, float dx, const float *pin, const float *vin, float *cresult_out);
int su_migps(int nx, int nt, float dt, float dx, const float *ffil_in, int nxpad, int ltaper, int np, int ntflag,
	const float *vt, const float *data, float *out);
int su_migpspi(int nx, int nt, int nz, float dz, float dt, float dx, const float *pin, const float *vin, float *cresult_out);


void su_kdmig2d_reference(int nr, int nzt, float dxo, float dzt, float fzt, float dvz, float v0, float *tb, float *pb, float *cs0b, float *angb);
void su_kdmig2d_residual(int ns, float fs, float ds, int nxt, float fxt, float dxt, int nzt, int nr, float dxo, const float *tb, float *ttab);
int su_kdmig2d_trace(const float *data, int nt, float ft, float dt, float sx, float gx, float *mig, float aperx, int nxo, float fxo, float dxo, int nzo, float fzo, float dzo, int ls, int mtmax, float dxm, float fmax, float angmax, const float *tb, const float *pb, const float *cs0b, const float *angb, int nr, const float *ttab, int ns, float fs, float ds, float es, float offmax, float *tsum, float *tt, int nzt, float fzt, float dzt, int nxt, float fxt, float dxt, int npv, const float *tv, const float *cs, float *tvsum, float *cssum, float *mig1);


int su_ktmig2d(int ntr, int nt, float dt, float dx, float h, float angmax, int nfc, int fwidth, const float *data_in, const float *vel_in, float *mig_out);
void su_ints8_tables(void);


// well logs
float gassmann_rockT0(float h, float kmin, float mumin, float rmin, float phi, float a, float b, float c, float api, float g, float s, float so, float sg, float t, float p);
float gassmann_rockR0(float v1, float r1, float kmin, float mumin, float rmin, float phi, float a, float b, float c, float api, float g, float s, float so, float sg, float t, float p);
float gassmann_rockIp(float kmin, float mumin, float rmin, float phi, float a, float b, float c, float api, float g, float s, float so, float sg, float t, float p);
float gassmann_rockVp(float kmin, float mumin, float rmin, float phi, float a, float b, float c, float api, float g, float s, float so, float sg, float t, float p);
float gassmann_rockDensity(float rmin, float phi, float api, float g, float s, float so, float sg, float t, float p);
float gassmann_fluidVelocity(float api, float g, float s, float so, float sg, float t, float p);
float gassmann_fluidModulus(float api, float g, float s, float so, float sg, float t, float p);
float gassmann_fluidDensity(float api, float g, float s, float so, float sg, float t, float p);
float gassmann_brineVelocity(float s, float t, float p);
float gassmann_brineModulus(float s, float t, float p);
float gassmann_brineDensity(float s, float t, float p);
float gassmann_oilModulus(float api, float g, float t, float p);
float gassmann_oilVelocity(float api, float g, float t, float p);
float gassmann_oilDensity(float api, float g, float t, float p);
float gassmann_gasVelocity(float g, float t, float p);
float gassmann_gasModulus(float g, float t, float p);
float gassmann_gasDensity(float g, float t, float p);


// dip moveout
void su_dmofk_vrms(int ndmo, const float *tdmo_in, const float *vdmo_in, int nt, float dt, float ft, float *vrms);
int su_dmofk_nxfft(int nx, float dx, float offset);
void su_dmofk_offset(float offset, float fmax, float sdmo, int nx, float dx, int nt, float dt, float ft, const float *vrms_in, float *ptx);
int su_dmotx_tables(float dxcdp, float dt, float offmax, float tmute, float vrms, int nsmax, int lds, int ifds, float *ts, float *as, float *ds);
void su_dmotx_trace(int ns, const float *ts, const float *as, float offset, int icdp, int ncdp, float dxcdp, int nt, float dt, float ft, const float *p, float *q);
void su_dmotx_shape(int ncdp, int nt, int lds, int ifds, const float *ds_in, float *q);
int su_dmofkcw_table(int ntable, float gamma, float *zoh, float *boh);
void su_dmofkcw_stretch(float sdmo, float gamma, float *s1, float *s2);
void su_dmofkcw_offset(float offset, float fmax, int nx, float dx, int nt, float dt, float ft, const float *vrms_in, float *ptx, float gamma, const float *boh, const float *zoh, int ntable, float s1, float s2, float sign);
struct DmoVz;
struct DmoVz *su_dmovz_create(int nt, float dt, const float *vt_in, int nptab, int *err);
int su_dmovz_nxfft(int nx, float dx, float offset);
int su_dmovz_offset(const struct DmoVz *M, int nx, float dx, float offset, float fmax, float smute, float speed, float *q);
void su_dmovz_free(struct DmoVz *M);
struct DmoTivz;
struct DmoTivz *su_dmotivz_create(int nt, float dt, const float *vn_in, const float *v_in, const float *eta_in, int nptab, int *err);
int su_dmotivz_nxfft(int nx, float dx, float offset);
int su_dmotivz_offset(const struct DmoTivz *M, int nx, float dx, float offset, float fmax, float smute, float speed, float *q);
void su_dmotivz_free(struct DmoTivz *M);
int su_tihaledmo(int nx, int nt, float dt, float dx, float h, float v, int option, int ntpad, int nxpad, int np, float dp, const float *vnmo, float e, float d, const float *data, float *out);


#endif
