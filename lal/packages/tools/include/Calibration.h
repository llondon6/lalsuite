/**** <lalVerbatim file="CalibrationHV">
 * Author: P. R. Brady, J. D. E. Creighton
 * $Id$
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * \section{Header \texttt{Calibration.h}}
 *
 * \subsection*{Synopsis}
 * \begin{verbatim}
 * #include <lal/Calibration.h>
 * \end{verbatim}
 *
 **** </lalLaTeX> */

#ifndef _CALIBRATION_H
#define _CALIBRATION_H

#include <lal/LALDatatypes.h>
#include <lal/BandPassTimeSeries.h>

#ifdef  __cplusplus
extern "C" {
#pragma } /** to match the previous brace **/
#endif

/**** <lalLaTeX>
 *
 * \subsection*{Error conditions}
 *
 **** </lalLaTeX> */
/**** <lalErrTable> */
#define CALIBRATIONH_ENULL 001
#define CALIBRATIONH_ESIZE 002
#define CALIBRATIONH_ESZMM 004
#define CALIBRATIONH_EZERO 010
#define CALIBRATIONH_ETIME 020
#define CALIBRATIONH_EUNIT 040

#define CALIBRATIONH_MSGENULL "Null pointer"
#define CALIBRATIONH_MSGESIZE "Invalid size"
#define CALIBRATIONH_MSGESZMM "Size mismatch"
#define CALIBRATIONH_MSGEZERO "Zero factor"
#define CALIBRATIONH_MSGETIME "Time out of range"
#define CALIBRATIONH_MSGEUNIT "Incompatible units"
/**** </lalErrTable> */

/**** <lalLaTeX>
 *
 * \subsection*{Structures}
 * \idx[Type]{CalibrationType}
 * \idx[Type]{CalibrationRecord}
 * \idx[Type]{CalibrationFunctions}
 * \idx[Type]{CalibrationUpdateParams}
 *
 * \subsubsection*{Type \texttt{CalibrationType}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef enum
{
  CalibrationAmplitude = 001,
  CalibrationOffset    = 002,
  CalibrationDelay     = 004,
  CalibrationTransfer  = 010,
  CalibrationZPG       = 020
}
CalibrationType;
/**** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 * \subsubsection*{Type \texttt{CalFactors}}
 *
 **** </lalLaTeX> */  
/**** <lalVerbatim> */
typedef struct 
tagCalFactors 
{
  COMPLEX16 alpha;    
  COMPLEX16 alphabeta;
  COMPLEX16 beta;
  COMPLEX16 exc;    
  COMPLEX16 asq;    
  COMPLEX16 darm;   
} 
CalFactors;
/**** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 * \subsubsection*{Type \texttt{UpdateFactorsParams}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagUpdateFactorsParams
{
   REAL8 lineFrequency;
   COMPLEX16 openloop;
   COMPLEX16 digital;
   REAL4TimeSeries *darmCtrl;
   REAL4TimeSeries *asQ;
   REAL4TimeSeries *exc;
}
UpdateFactorsParams;
/**** </lalVerbatim> */

/**** <lalLaTeX>
 *
 * Document \verb+CalibrationType+
 *
 * \subsubsection*{Type \texttt{CalibrationRecord}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationRecord
{
  CHAR                     name[LALNameLength];
  LIGOTimeGPS              epoch;
  REAL8                    duration;
  CHAR                     reference[LALNameLength];
  LALUnit                  units;
  UINT4                    type;
  REAL8                    conversion;
  REAL8                    offset;
  REAL8                    delay;
  COMPLEX8FrequencySeries *transfer;
  REAL8Vector             *zeros;
  REAL8Vector             *poles;
  REAL8                    gain;
}
CalibrationRecord;
/**** </lalVerbatim> */


/**** <lalLaTeX>
 *
 * Document \verb+CalibrationRecord+
 *
 * \subsubsection*{Type \texttt{CalibrationFunctions}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationFunctions
{
  COMPLEX8FrequencySeries *responseFunction;
  COMPLEX8FrequencySeries *sensingFunction;
}
CalibrationFunctions;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * The type \texttt{CalibrationFunctions} contains two calibration functions,
 * the sensing function $C(f)$ and the response function $R(f)$.  While the
 * response function is the function that is most often wanted, the sensing
 * function is needed in updating calibration from one epoch to another.
 *
 * \subsubsection*{Type \texttt{CalibrationUpdateParams}}
 *
 **** </lalLaTeX> */
/**** <lalVerbatim> */
typedef struct
tagCalibrationUpdateParams
{
  LIGOTimeGPS epoch;
  LIGOTimeGPS duration;
  CHAR *ifo;
  COMPLEX8 alpha;
  COMPLEX8 alphabeta;
  COMPLEX8TimeSeries *openLoopFactor;
  COMPLEX8TimeSeries *sensingFactor;
}
CalibrationUpdateParams;
/**** </lalVerbatim> */
/**** <lalLaTeX>
 * The type \texttt{CalibrationUpdateParams} contains two time series
 * representing an overall gain factor for the open-loop gain function $H(f)$
 * and the sensing function $C(f)$.  These transfer functions are known to
 * change (by an overall factor) with time, and these two factors can be
 * tracked using the injected calibration lines.  The factors are stored
 * in this structure as (very-slowly varying) time series, and are to be
 * used in updating the calibration functions described previously.
 * (The response function can be computed from the open-loop gain and the
 * sensing function.  It is simply $R(f)=[1+H(f)]/C(f)$.)  In addition, this
 * structure contains the present epoch and the duration of the data to be 
 * calibrated to identify the particular set of
 * factors (from those recorded in the time series) to use.  
 *
 * \vfill{\footnotesize\input{CalibrationHV}}
 * \newpage\input{ComputeTransferC}
 *
 **** </lalLaTeX> */

typedef
struct StrainOutTag {
  REAL8TimeSeries h;         /* timeseries containing h(t) */
  COMPLEX16TimeSeries alpha; /* alpha time series */
  COMPLEX16TimeSeries beta;  /* beta time series */
} StrainOut;
                                                                                                                               
typedef
struct StrainInTag {
  REAL4TimeSeries AS_Q ;   /* timeseries containing ASQ */
  REAL4TimeSeries DARM ;   /* timeseries containing DARM_CTRL */
  REAL4TimeSeries EXC ;    /* timeseries containing the excitation */
  COMPLEX16 Do;            /* digital filter at cal line frequency */
  COMPLEX16 Go;            /* OLG at cal line frequency */
  REAL8 f;                 /* calibration line frequency */
  REAL4 To;                /* factors integration time */
  REAL8IIRFilter *Cinv;    /* Filters for inverse of sensing function */
  INT4 CinvUSF;            /* Upsampling factor for sensing function */
  INT4 CinvDelay;          /* Overall inverse sensing function delay */
  REAL8IIRFilter *AA;      /* Filters for analog actuation function */
  INT4 AADelay;            /* Overall analog actuation function delay */
  REAL8IIRFilter *AX;      /* Digital filters for x arm actuation function */
  REAL8IIRFilter *AY;      /* Digital filters for y arm actuation function */
  REAL8IIRFilter *G;       /* Digital servo filters */
  INT4 NCinv;              /* Numbers of filters of each type */
  INT4 NAA;
  INT4 NAX;
  INT4 NAY;
  INT4 NG;
} StrainIn;

typedef
struct MyIIRFilter {
  INT4 yOrder;
  INT4 xOrder;
  REAL8 a[20];
  REAL8 b[20];
  REAL8 yhist[20];
  REAL8 xhist[20];
} MyIIRFilter;



void LALComputeTransfer( LALStatus *status, CalibrationRecord *calrec );


void
LALUpdateCalibration(
    LALStatus               *status,
    CalibrationFunctions    *output,
    CalibrationFunctions    *input,
    CalibrationUpdateParams *params
    );

void
LALResponseConvert(
    LALStatus               *status,
    COMPLEX8FrequencySeries *output,
    COMPLEX8FrequencySeries *input
    );

void LALComputeCalibrationFactors(
    LALStatus              *status,
    CalFactors             *output,    
    UpdateFactorsParams    *input
    );

void LALComputeStrain(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    );

void XLALGetFactors(
    LALStatus              *status,
    StrainOut              *output,    
    StrainIn               *input
    );

int XLALhROverAlpha(REAL8TimeSeries *hR, StrainOut *output);
int XLALhCTimesBeta(REAL8TimeSeries *hC, StrainOut *output);
int XLALUpsamplehR(REAL8TimeSeries *uphR, REAL8TimeSeries *hR, int up_factor);

void XLALReadFilters(
     LALStatus *status, 
     REAL8IIRFilter **Cinv, 
     REAL8IIRFilter **G,
     REAL8IIRFilter **AA,
     REAL8IIRFilter **AX,
     REAL8IIRFilter **AY,
     char *fname,
     int GPS,
     int *USF,
     int *NC, 
     int *NG, 
     int *NAA, 
     int *NAX, 
     int *NAY, 
     int *CD, 
     int *AAD);

void XLALWriteFilters(
     LALStatus *status, 
     REAL8IIRFilter *Cinv, 
     REAL8IIRFilter *G,
     REAL8IIRFilter *AA,
     REAL8IIRFilter *AX,
     REAL8IIRFilter *AY,
     char *fname,
     int GPS,
     int USF,
     int NC, int NG, int NAA, int NAX, int NAY, int CD, int AAD);


#ifdef  __cplusplus
#pragma { /** to match the next brace **/
}
#endif

#endif /* _CALIBRATION_H */
