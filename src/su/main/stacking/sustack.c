/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* SUSTACK: $Revision: 1.27 $ ; $Date: 2011/11/16 23:14:54 $	*/

#include "su.h"
#include "segy.h"

#if 0 /* the program, which is not part of the library */
/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" SUSTACK - stack adjacent traces having the same key header word",
" 								",
"     sustack <stdin >stdout [Optional parameters]		",
" 							        ",
" Required parameters:						",
" 	none							",
" 							        ",
" Optional parameters: 						",
" 	key=cdp		header key word to stack on		",
" 	normpow=1.0	each sample is divided by the		",
"			normpow'th number of non-zero values	",
"			stacked (normpow=0 selects no division)	",
"	repeat=0	=1 repeats the stack trace nrepeat times",
"	nrepeat=10	repeats stack trace nrepeat times in	",
"	          	output file				",
" 	verbose=0	verbose = 1 echos information		",
" 							        ",
" Notes:							",
" ------							",
" The offset field is set to zero on the output traces, unless	",
" the user is stacking with key=offset. In that case, the value ",
" of the offset field is left unchanged. 		        ",
" 							        ",
" Sushw can be used afterwards if this is not acceptable.	",
"								",
" For VSP users:						",
" The stack trace appears ten times in the output file when	",
" setting repeat=1 and nrepeat=10. Corridor stacking can be	",
" achieved by properly muting the upgoing data with SUMUTE	",
" before stacking.						",
"								",
NULL};

/* Credits:
 *	SEP: Einar Kjartansson
 *	CWP: Jack K. Cohen, Dave Hale
 *	CENPET: Werner M. Heigl - added repeat trace functionality
 *
 * Note:
 *	The "valxxx" subroutines are in su/lib/valpkge.c.  In particular,
 *      "valcmp" shares the annoying attribute of "strcmp" that
 *		if (valcmp(type, val, valnew) {
 *			...
 *		}
 *	will be performed when val and valnew are different.
 *
 * Trace header fields accessed: ns
 * Trace header fields modified: nhs, tracl, offset
 */
/**************** end self doc ***********************************/


segy intrace, outtrace;

int
main(int argc, char **argv)
{
	cwp_String key;	/* header key word from segy.h		*/
	cwp_String type;/* ... its type				*/
	int indx;	/* ... its index			*/
	int nt;		/* number of data points on trace	*/
	int nsegy;	/* number of bytes in the segy		*/
	Value val;	/* value of key in current gather	*/
	Value valnew;	/* value of key in trace being treated	*/
	int fold;	/* number of traces stacked		*/
	int *nnz;	/* number of non-zero values stacked	*/
	float normpow;	/* divide by nnz[i]^normpow		*/
	int newtracl;	/* tracl for stacked traces		*/
	int repeat;	/* flag for stack trace repeating	*/
	int nrepeat;	/* no. of times stack trace is repeated	*/
	int verbose;	/* verbose flag				*/

	cwp_Bool is_offset=cwp_false;	/* stacking on offset?  */


	/* Initialize */
	initargs(argc, argv);
	requestdoc(1);


	/* Set parameters */
	if (!getparint   ("verbose", &verbose))	 verbose = 0;
	if (!getparfloat ("normpow", &normpow))	 normpow = 1.0;
	if (!getparstring("key", &key))		 key = "cdp";
	if (!getparint   ("repeat",&repeat))	 repeat = 0;
	if (!getparint   ("nrepeat",&nrepeat))	 nrepeat = 10;


        checkpars();

	if (STREQ(key,"offset"))		is_offset=cwp_true;

	type = hdtype(key);
	indx = getindex(key);

	/* Set up for first trace (must compare new key field each time) */
	nsegy = gettr(&intrace);
	nt = intrace.ns;
	memcpy( (void *) &outtrace, (const void *) &intrace, nsegy);
	nnz = ealloc1int(nt);
	{ register int i;
	  for (i = 0; i < nt; i++){
		if (intrace.data[i] != 0.0)  nnz[i]=1;
		else nnz[i] = 0;
	  }
        }
	if (intrace.nhs!=0)
		fold = intrace.nhs;
	else
		fold = 1;


	/* Loop over traces */
	newtracl = 1;
	gethval(&intrace, indx, &val);
	while (nsegy) {		     /* While previous trace non-empty */
		nsegy = gettr(&intrace);
		gethval(&intrace, indx, &valnew);
		if (valcmp(type, val, valnew) || !nsegy) {	
			/* Either val and valnew differ, indicating a  */
			/* new gather or nsegy is zero, indicating the */
		        /* end of the traces.                          */
			if (verbose) {
				fprintf(stderr, "val=");
				fprintfval(stderr, type, val);
				fprintf(stderr, "\tfold=%d\n", fold);
			}

			/* Add header info and output stack */
			outtrace.nhs = fold;
			outtrace.tracl = newtracl++;
			if(!is_offset) outtrace.offset = 0;
			if (normpow && fold != 1) {
			        register int i;
				for (i = 0; i < nt; ++i) {
				    float nnzi = nnz[i];
				    if (nnzi)
					outtrace.data[i] /= pow(nnzi, normpow);
				}
			}
			if (repeat) {
				register int i;
				for (i=0;i<nrepeat;i++)
					puttr(&outtrace);
			} else	puttr(&outtrace);

			/* Set up for next gather */
			memcpy( (void *) &outtrace,
					(const void *) &intrace, nsegy);
			{ register int i;
	  		  for (i = 0; i < nt; i++){
				if (intrace.data[i] != 0.0)  nnz[i]=1;
				else nnz[i] = 0;
			  }
			}

			if (intrace.nhs!=0)
				fold = intrace.nhs;
			else
				fold = 1;

			val = valnew;

		} else {	/* still in same gather */
			register int i;
			for (i = 0; i < nt; ++i) {
				float datum = intrace.data[i];
				if (!(datum == 0.0))  ++nnz[i];
				outtrace.data[i] += datum;
			}
			++fold;
		}
	}


	return(CWP_Exit());
}
#endif

/* Library version of SUSTACK
 *
 * The main program is not part of the library. What it does with the traces of a gather is two functions. A trace has n samples of
 * ncomp floats each (1 for a real trace, 2 for a complex one, the real and imaginary part next to each other).
 *
 *  su_stack_add()        adds a trace to the sums sum[n * ncomp], and counts the samples that were not 0 in nnz[n]
 *  su_stack_normalize()  divides each sample by the number of values that were not 0 to the power normpow
 *                        (normpow 0 is no division, a sample that was always 0 is left)
 */
void su_stack_add(int n, int ncomp, const float *x, float *sum, int *nnz)
{
	int i, c, nonzero;

	for (i = 0; i < n; ++i) {
		nonzero = 0;
		for (c = 0; c < ncomp; ++c) {
			float datum = x[i*ncomp+c];
			if (datum != 0.0) nonzero = 1;
			sum[i*ncomp+c] += datum;
		}
		if (nonzero) ++nnz[i];
	}
}

void su_stack_normalize(int n, int ncomp, float *sum, const int *nnz, float normpow)
{
	int i, c;

	for (i = 0; i < n; ++i) {
		if (nnz[i]) {
			float divisor = pow(nnz[i], normpow);
			for (c = 0; c < ncomp; ++c) sum[i*ncomp+c] /= divisor;
		}
	}
}
