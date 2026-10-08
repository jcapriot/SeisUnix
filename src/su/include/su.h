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
#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

/* valpkge */
int vtoi(register cwp_String type, Value val);
long vtol(register cwp_String type, Value val);
float vtof(register cwp_String type, Value val);
double vtod(register cwp_String type, Value val);
int valcmp(register cwp_String type, Value val1, Value val2);
void printfval(register cwp_String type, Value val);
void fprintfval(FILE *stream, register cwp_String type, Value val);
void scanfval(register cwp_String type, Value *valp);
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

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif

#endif
