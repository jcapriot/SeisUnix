/* Copyright (c) Colorado School of Mines, 2022.*/
/* All rights reserved.                       */

/* SUSEMBASE: $Revision: 1.1 $ ; $Date: 2023/03/02 23:05:59 $      */

#include <stdio.h>
#include <string.h>

#include "su.h"
#include <vector>
#include "segy.h"

segy tr;
segy trd;

/*********************** self documentation **********************/
#if 0 /* the sdoc of the program, which is not part of the library */
char *sdoc[] = {
"                                                                            ",
" SUSEMBASE - Semblance Based Structure-Oriented Interpolate and Smooth.     ",
"                                                                            ",
" This program performs semblance scanning to derive (apparent) dips in the  ",
" direction of cdp numbers (or other key values). It then uses those dips to ",
" align trace sample values while performing interpolation and/or smoothing. ",
"                                                                            ",
" susembase    <in.su >out.su  (other parameters)                            ",
"                                                                            ",
" Parameter overview:                                                        ",
"                                                                            ",
" lockey=cdp   Trace key that indicates trace location. Input traces must be ",
"              in increasing-OR-EQUAL order of these values. This key value  ",
"              is reset on output. See outputinc and outputalign parameters. ",
"        Note: The input locations must be in increasing-OR-EQUAL order.     ",
"              So, traces can have duplicate lockeys provided the numrolling ",
"              parameter is set large enough. Whether these duplicates have  ",
"              serious effects on the output traces varies with situation.   ",
"                                                                            ",
" groupkey=    Trace key that indicates the end of a set of lockey values.   ",  
"              If specified, traces must be ordered by lockey and groupkey.  ",  
"              If a 3D stacked dataset is input, this parameter should often ",  
"              be either igi or igc, depending on which direction is being   ",  
"              processed in this run (see program subincvs).                 ",  
"        Note: The keys for susort are listed in reverse order with respect  ",  
"              to lockey and groupkey here.                                  ",  
"                                                                            ",
" outputinc=1.0    Increment for output locations. Example: if a cdp stack   ",
"                  is input and cdp numbers increment by 1 then setting this ",
"                  parameter to 0.5 will output about twice as many traces.  ",
"            Note: Since lockey is (almost always) integer, it cannot hold   ",
"                  decimal fractions, so if this outputinc parameter is      ",
"                  between 0 and 1, the lockey output values are inversed:   ",
"                   (input_values - outputalign) / outputinc + outputshift   ",
"                  If outputinc is greater than 1 then the lockey output     ",
"                  values are:    (input_values - outputalign) + outputshift ",
"                                                                            ",
" outputalign=0.0  Alignment for output locations. Example: if a cdp stack   ",
"                  is input and cdp numbers increment by 1 and outputinc=1.0 ",
"                  then setting this parameter to 0.25 will create output    ",
"                  traces at input cdp numbers plus 0.25                     ",
"                                                                            ",
" outputshift=0.0  Add this value to output location numbers. If lockey is   ",
"                  cdp and first input cdp number is small, traces may be    ",
"                  extrapolated to negative cdp numbers. This may be an      ",
"                  issue for some situations, so you can add this constant.  ",
"                                                                            ",
" locmaxscan=3 Maximum Relative lockey to include in scanning for dips.      ",
"              Must be greater than or equal to 0. This program scans traces ",
"              from -locmaxscan to locmaxscan (inclusive) at each location.  ",
"           =0 Do not scan for dips. You must specify parameter locmaxsum    ",
"              and you cannot specify scanstep, scanmax, scanmin, dipext.    ",
"                                                                            ",
" scanstep=0.5   Time Step to Scan (milliseconds per 1 lockey unit).         ",
" scanmax =12.0  Maximum Time to Scan (milliseconds per 1 lockey unit).      ",
" scanmin =-scanmax  Minimum Time to Scan (milliseconds per 1 lockey unit).  ",
"                                                                            ",
" dipext=12. Dip Extender (milliseconds).                                    ",
"            At each location the semblance scanning produces an estimate    ",
"            of the dip at each sample time. This parameter performs a       ",
"            weighted-average of those estimated dips symmetrically from     ",
"            -dipext to +dipext. The weights used are the sum of the         ",
"            squared amplitudes of seismic values along the estimated dip    ",
"            in the locmaxscan range of input traces. In cruder terms: the   ",
"            power along the dip value chozen by the semblance scanning.     ",
"            Effectively, this causes the dip values from large horizons     ",
"            to be trusted more and extended up and down nearby times.       ",
"                                                                            ",
" dipuse=1  Estimate and use dips at input locations.                        ",
"       =0  Estimate and use dips at output locations.                       ",
"                                                                            ",
" locmaxsum=locmaxscan  Maximum Relative lockey to include in sumopt options.",
"                       Must be greater than 0. This program uses the input  ",
"                       traces from -locmaxsum to locmaxsum (inclusive) at   ",
"                       each output location.                                ",
"                                                                            ",
" sumopt=stack All options here only include input traces within -locmaxsum  ",
"              to +locmaxsum range (inclusive) around each output location.  ",
"              This default just adds them all together and divides by total.",
"              Note that correct dips from the semblance scanning will cause ",
"              these traces to be locally flat during stacking and therefore ",
"              the output will not have as much smearing as traces stacked   ",
"              without dip compensation.                                     ",
"       =sum   Just like option stack but do not divide by total.            ",
"              (This option is intended to help with parameter setting).     ",
"       =one   Find the input trace with lockey value nearest the output     ",
"              location, apply the dips, and output it. Effectively, this    ",
"              is interpolation using just the nearest trace.                ",
"       =two   Find the two input traces with lockey value nearest output,   ",
"              apply the dips, sum both traces and divide by 2, and output.  ",
"       =twolinear Find the two input traces with lockey values nearest      ",
"                  output, apply the dips, stack them with linear weight     ",
"                  based on their locations, and output. (Often one of the   ",
"                  two traces will get all the weight, of course).           ",
"       =dadstack  Detect-Amplitude-Differences. You should expect to have   ",
"                  to do considerable testing of the other dad parameters    ",
"                  for this option. Precise setting of the dad parameters    ",
"                  may result in the output traces being sharper at rapid    ",
"                  amplitude changes (such as locations where horizons start ",
"                  or terminate. These dad options search for rapid changes  ",
"                  in amplitude ACROSS LOCATIONS at each sample time and     ",
"                  only stacks the sample values from -locmaxsum TO the      ",
"                  detected location -or- FROM the detected location to the  ",
"                  +locmaxsum location. Which side is stacked depends on     ",
"                  WHERE detected location is relative to output location.   ",
"             ***  If a rapid change in amplitude is not detected at this    ",
"                  sample time then ALL samples are stacked.                 ",
"       =dadsum  Just like option dadstack but do not divide by total.       ",
"                (This option is intended to help with parameter setting).   ",
"                                                                            ",
" dadscale=1.0 Detect-Amplitude-Differences scaler.                          ",
"             This parameter can only be set if sumopt=dadstack or dadsum.   ",
"             When a rapid change in amplitude is detected this parameter    ",
"             causes the range of locations stacked at this sample to        ",
"             increase or decrease.                                          ",
"       Note: Usual range for this value is 0.9 to 1.1                       ",
"                                                                            ",
" dadnoise=0.0  Detect-Amplitude-Differences Noise Absolute Amplitude.       ",
"              This parameter can only be set if sumopt=dadstack or dadsum.  ",
"              If the average absolute amplitude at the sample time across   ",
"              the range of input traces is less than this value, this       ",
"              sample time is not considered to have detected a rapid change ",
"              in amplitude (so all locations are stacked for output sample).",
"                                                                            ",
" dadext=12.   Detect-Amplitude-Differences extender (milliseconds).         ",
"              This parameter can only be set if sumopt=dadstack or dadsum.  ",
"              At each output location the amplitude difference detector     ",
"              estimates if-and-where a rapid difference is located in the   ",
"              locmaxsum range of input traces. This parameter performs a    ",
"              weighted-average of those estimated rapid difference locations",
"              symmetrically in time from -dadext to +dadext. The weights    ",
"              used are the sum of absolute amplitudes of the seismic values ",
"              at the same samples in the locmaxsum range of input traces.   ",
"              Effectively, this means the detected-amplitude-difference     ",
"              locations at large horizons are trusted more, and thus they   ",
"              are extended up and down in time near those larger horizons.  ",
"                                                                            ",
" outdip=   Do not output this file.                                         ",
"       =(file name) Output the dip values chosen by the semblance scanning  ",
"           into trace samples in this file. The trace header values are     ",
"           copied from the seismic traces. If dipuse=1 then they are the    ",
"           input headers, if dipuse=0 then they are the output headers.     ",
"     Note: These dips are after dipext processing. But even with dipext=0   ",
"           they will not have scanstep granularity because the scanned      ",
"           maximum semblance dip undergoes a 3-point interpolation to get   ",
"           a better approximation of the dip at the true peak in semblance. ",
"     Note: It is technically more correct to call the dips made and used    ",
"           by this program APPARANT dips. The dip values are put into the   ",
"           trace samples and are in milliseconds per 1 location unit.       ",
"           That is, if the lockey=cdp then the dips are ms per cdp. There   ",
"           is no dip azimuth output. The dips are in the direction from one ",
"           location to next (if lockey=cdp this is one cdp toward next cdp).",
"                                                                            ",
" outpow=   Do not output this file.                                         ",
"       =(file name) Output the average power along the dip direction chosen ",
"           by the semblance scanning. (Same header values as for outdip).   ",
"                                                                            ",
" outsem=   Do not output this file.                                         ",
"       =(file name) Output the semblance value computed and chosen by the   ",
"           semblance scanning. (Same header values as for outdip).          ",
"                                                                            ",
" indip=    Name of trace file which contains (apparant) dip values.         ",
"           Dip scanning is disabled if you use this parameter. Thus you     ",
"           cannot specify parameters scanstep, scanmax, scanmin, dipext,    ",
"           and locmaxscan - and you must specify parameter locmaxsum.       ",
"           This file is usually initially made by specifying outdip above.  ",
"           If dipuse=1 then this file must have the same order as the       ",
"           order as the seismic traces you are inputting in this run.       ",
"           If dipuse=0 then this file must have the same INPUT order as     ",
"           the OUTPUT order of this run.                                    ",
"     Note: The intention here is to allow you to output dips from one run   ",
"           of this program and use them in a second run of this program.    ",
"           The first run can input smoothed, denoised seismic to get dips,  ",
"           and the second run can input the un-smoothed, noisy seismic and  ",
"           use the dips from the first run. As well, you can also actually  ",
"           smooth and denoise the DIP VALUES in this file before input here.",
"                                                                            ",
"                                                                            ",
" bypass=0  Treat all output locations the same.                             ",
"       =1  If an input trace exists exactly at the output location, COPY    ",
"           the input trace without any change except to reset the lockey    ",
"           to the output location value. Process all other output locations ",
"           as indicated by the other options you have specified.            ",
"       =-1 Only process and output traces where an input trace is exactly   ",
"           on an output location. Typically you use this option when you    ",
"           just want to smooth the input, without creating any more output  ",
"           locations (so set outputalign=0 and usually set outputinc=1).    ",
"  CAUTION: Options 1 and -1 only retain the first trace with duplicate      ",
"           input location value. For instance if lockey=offset and          ",
"           groupkey=cdp then only the first trace with the same offset      ",
"           within cdps will be copied or processed to output.               ",
"                                                                            ",
" numrolling=2*locmaxscan+2*locmaxsum+5. This is maximum number of traces    ",
"              to store internally. Traces are read-in as needed to keep the ",
"              output locations centred in this storage (a rolling buffer).  ",
"          *** For certain setups, this default is too large. For instance,  ",
"              when lockey=offset and groupkey=cdp this value can be set to  ",
"              the maximum cdp fold (but must always be set greater than 4). ",
"          *** For certain setups, this default is too small. For instance,  ",
"              when duplicate input lockey values exist.                     ",
"        Note: This program expects input locations to have approximately    ",
"              the same interval between them and to be distributed nicely.  ",
"              The input locations must be in increasing-OR-EQUAL order.     ",
"              So, traces can have duplicate lockey values provided this     ",
"              parameter is set large enough. Whether these duplicates have  ",
"              serious effects on the output traces varies with situation.   ",
"                                                                            ",
" Note: This program enforces non-hierarchical parameter error checking.     ",
"       This means if you disable a main option, all parameters related      ",
"       to it must be removed from the command line. That is, for novice     ",
"       users this program is sufficiently complex that it refuses to guess  ",
"       which parameter was specified in error.                              ",
"                                                                            ",
" ------------------------------------------------------------------------   ", 
"                                                                            ",
NULL};
#endif

/* Credits:                                                                  */
/* Author: Andre Latour, Feb 2023                                            */ 
/*                                                                           */
/* Trace key usually involved: cdp   (user options: maybe igi,igc others)    */
/*                                                                           */
/**************** end self doc ***********************************************/

/* The scanning and the smoothing of the dips of susembase (and susemgroup, which has them as well), as library functions
 * of arrays. The programs keep their traces in a rolling buffer and their locations in the headers: these are the
 * functions of the dips (scanit and wsmooth of the programs), the rest is done with the traces in memory.
 *
 * scanit had its working buffers passed by the caller, for the numrolling traces; they are made here.
 */

/* Semblance Scanner                                                          */
/*                                                                            */
/* Semblance herein is computed as:                                           */
/*                                                                            */
/*  (Square of (Sum of amplitudes in locminscan to locmaxscan range)) /       */
/*             (Sum of squared amplitudes in that range))                     */
/*             divided by the number in that range.                           */
/*                                                                            */
/* Scans from smin to smin*sstep*(jsmax-1) are done for each sample time and  */
/* the maximum semblance value computed as above is found. That maximum is    */
/* then 3-point interpolated using the power at that maximum and the power    */
/* at its 2 neighbours (i.e. the values of their denominators above).         */
/* Note that this 3 point interpolation improves the result                   */
/* and also has the benefit of reducing the sstep granularity.                */
/*   Note this semblance is computed over just 1 time sample, but see dipext. */
/*                                                                            */
/*                                                                            */
/*                                                                            */
/* data       = The input traces, [numrolling][lendip].                        */
/* aloc       = Location value for each input trace.                          */
/* numrolling = Number of traces in atrace and locations in aloc.             */
/* scanloc    = The location value where the scan takes place.                */
/* locminscan = The minimum relative location to include in scan (inclusive). */
/* locmaxscan = The maximum relative location to include in scan (inclusive). */
/* smin       = Minimum relative sample to scan.                              */
/* sstep      = Fractional sample interval to scan.                           */
/* jsmax      = Number of intervals to scan (from smin, increment by sstep).  */
/*                                                                            */
/* adip       = Output apparant dip (in samples per 1 location difference).   */
/* apow       = Output power along the chosen dip.                            */
/* asem       = Output semblance of the chosen dip.                           */
/*                                                                            */
/* lendip     = Length of traces and adip,apow,asem,awrk,awrq.                */
/* numout     = Number of traces in scanloc+locminscan to scanloc+locminscan. */
/*              0 returns all adip=0 and asem=0 and apow=0                    */
/*              1 returns all adip=0 and asem=1 and apow=trace amp squared    */


void su_sem_scan(const float *data, const int *aloc, int numrolling, double scanloc,
            double locminscan, double locmaxscan, float smin, float sstep, int jsmax,
            float *adip, float *apow, float *asem, int lendip, int *numout)
{

  int i=0;
  int j=0;
  int m=0;
  int n=0;
  float alin=0.;
  float aval=0.;
  float sval=0.;
  int numinrange=0;
  std::vector<float> awrk(jsmax), awrq(jsmax);
  std::vector<int> arng(numrolling), ifori((size_t) numrolling*jsmax);
  std::vector<float> afori((size_t) numrolling*jsmax);

  float thesem=0.;
  int thej=0; 
  float sl=0.;
  float sr=0.;
  float sp=0.;

  for (n=0; n<lendip; n++) {
    adip[n] = 0.;
    apow[n] = 0.;
    asem[n] = 0.;
  }

  for (i=0; i<numrolling; i++) {
    if(aloc[i] >= scanloc+locminscan && aloc[i] <= scanloc+locmaxscan) {
      arng[numinrange] = i;
      numinrange++;
    }
  }

  *numout = 0;
  if(numinrange<1) return; 
  *numout = 1;

/* If only 1 trace in range, then dip is presumed 0 (as initialized above)    */
/* And set power and semblance values correctly.                              */

  if(numinrange<2) {
    for (n=0; n<lendip; n++) {
      apow[n] = data[(size_t) arng[0]*lendip + (n)] * data[(size_t) arng[0]*lendip + (n)];
      asem[n] = 1.;
    }
    return; 
  }

/* Precompute the whole and fractional samples for linear interpolation so    */
/* it does not get repeated for every trace sample.                           */
/* (It is possible the compiler is smart enough to do this, but ????)         */
/* Note to future programmers: you could also use ints8r (or similar)         */
/*      to get more accurate seismic values than linear interpolation.        */
/*      But that would require passing a buffer to hold the shifted traces.   */

  for(j=0; j<jsmax; j++) {
    sval = smin + j*sstep;
    for (i=0; i<numinrange; i++) {
      m = ((int) (sval*(aloc[arng[i]]-scanloc)));
      if(sval*(aloc[arng[i]]-scanloc) > 0.) m++;
      alin = m - (sval*(aloc[arng[i]]-scanloc));
      ifori[(size_t) i*jsmax + j] = m;
      afori[(size_t) i*jsmax + j] = alin;
    }
  }

  for (n=0; n<lendip; n++) {                 
    for(j=0; j<jsmax; j++) {
      awrk[j] = 0.;                                                                
      awrq[j] = 1.e-30;  /* prevents divide by 0 when all samples = 0.0       */
      for (i=0; i<numinrange; i++) {
        m = ifori[(size_t) i*jsmax + j];
        if(n+m > 0 && n+m < lendip) { /* totals near top,bottom != numinrange */
          alin = afori[(size_t) i*jsmax + j];        
          aval = alin * data[(size_t) arng[i]*lendip + (n+m-1)] + (1.-alin) * data[(size_t) arng[i]*lendip + (n+m)];
          awrk[j] += aval;
          awrq[j] += aval*aval;
        }
      } 
    } /* end of  for(j=0; j<jsmax; j++) { */ 

/* Keep biggest semblance value and which scan number this is (thej).         */

    thesem = 0.;
    thej = jsmax/2; /* in case all zeroes (like a mute zone) */
    for(j=0; j<jsmax; j++) {
      if(awrk[j]*awrk[j]/awrq[j] > thesem) {
        thesem  = awrk[j]*awrk[j]/awrq[j];
        thej    = j;
      }
    } /* end of  for(j=0; j<jsmax; j++) { */

 /* At this point in the code, the dip values have sstep granularity.         */
 /* So perform 3 point interpolation around the maximum semblance value.      */
 /* This is done because switching from one sstep value to another might      */
 /* cause some effects if a reflection has slowly changing actual dip.        */
 /* Also, it allows users easier testing of what sstep is good enough (if the */
 /* 3 point interpolated dip of a big sstep gives as good results as the      */
 /* 3 point interpolated dip of a small sstep, use big and save CPU time).    */
           
    if(thesem < 1.e-30) {   /* if all zeroes (mute zones).                    */
      adip[n] = 0.0;  
      asem[n] = 0.0;  
      apow[n] = 0.0;
    }
    else if(thej>0 && thej<jsmax-1) {
      sl = awrk[thej-1]*awrk[thej-1]/awrq[thej-1];
      sr = awrk[thej+1]*awrk[thej+1]/awrq[thej+1];
      sp = 0.5 * (sl-sr) / (sl - 2.0*thesem + sr + 1.e-30);
      adip[n] = smin + thej*sstep + sp*sstep;     /* dip at interpolated peak */
      asem[n] = (thesem - 0.25 * (sl-sr) * sp)/numinrange; /* semblnc of peak */
      apow[n] = awrq[thej]/numinrange; /* power along dip of max semblnc pick */
    }
    else {
      adip[n] = smin + thej*sstep;
      asem[n] = thesem/numinrange;
      apow[n] = awrq[thej]/numinrange;
    }

  } /* end of  for (n=0; n<lendip; n++) { */ 

  return;
}


/* Weighted smoothing                                                         */
/*                                                                            */
/* val = input values                                                         */
/* wgt = input weights                                                        */
/* len = length of val,wgt,sval.                                              */
/* minsd = lower range to include in sum for i-th value (usually = -maxsd)    */
/* maxsd = higher range to include in sum for i-th value (usually positive)   */
/*                                                                            */
/* sval = output smoothed values                                              */
/*                                                                            */
/* Notes:                                                                     */
/* 1. Weights must be >= 0 or you could get divide by 0. (not checked herein) */
/* 2. Could do this by add and subtract to running totals as i increments,    */
/*    which would be faster for big minsd,maxsd. But do not expect that here. */
/*    (And the compiler-optimizer will unroll the inner loop).                */

void su_sem_wsmooth (const float *val, const float *wgt, int len, int minsd, int maxsd, float *sval) {

  int i=0;
  int k=0;
  float twgt=0.;
  float tval=0.;

/* A trace shorter than the smoother has no interior to smooth (the program reads outside of its arrays). */
  if(len-maxsd <= 0-minsd) {
    for(i=0; i<len; i++) sval[i] = val[i];
    return;
  }

  for(i=0-minsd; i<len-maxsd; i++) { 
    twgt = 1.e-30;
    tval = 0.;
    for(k=minsd; k<=maxsd; k++) {
      twgt += wgt[i+k];
      tval += wgt[i+k] * val[i+k];
    }
    sval[i] = tval / twgt;
  }

  for(i=0; i<0-minsd; i++) sval[i] = sval[0-minsd]; 
  for(i=len-maxsd; i<len; i++) sval[i] = sval[len-maxsd-1];

  return;
}
