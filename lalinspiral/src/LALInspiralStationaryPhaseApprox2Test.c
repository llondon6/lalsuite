/*
*  Copyright (C) 2007 Jolien Creighton, B.S. Sathyaprakash, Thomas Cokelaer
*
*  This program is free software; you can redistribute it and/or modify
*  it under the terms of the GNU General Public License as published by
*  the Free Software Foundation; either version 2 of the License, or
*  (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with with program; see the file COPYING. If not, write to the
*  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
*  MA  02111-1307  USA
*/

/**** <lalVerbatim file="LALInspiralStationaryPhaseApprox2CV">
 * Author: B.S. Sathyaprakash
 **** </lalVerbatim> */

/**** <lalLaTeX>
 *
 *
 * \subsection{Module \texttt{LALInspiralStationaryPhaseApprox2.c}}
 * %% A one-line description of the function(s) defined in this module.
 * This module computes the usual stationary phase approximation to the
 * Fourier transform of a chirp waveform
 * given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 *
 * \subsubsection*{Prototypes}
 * \input{LALInspiralStationaryPhaseApprox2CP}
 * \idx{LALInspiralStationaryPhaseApprox2()}
 * \begin{itemize}
 * \item {\tt signalvec:} Output containing the inspiral waveform.
 * \item {\tt params:} Input containing binary chirp parameters.
 * \end{itemize}
 *
 * \subsubsection*{Description}
 *
 * %% A description of the data analysis task performed by this function;
 * %% this is the main place to document the module.
 * Computes the Fourier transform of the chirp signal in the stationary
 * phase approximation and returns the real and imagninary parts of the
 * Fourier domain signal in the convention of fftw. For a signal vector
 * of length {\tt n=signalvec->length} ({\tt n} even):
 * \begin{itemize}
 * \item {\tt signalvec->data[0]} is the {\it real} 0th frequency component of the Fourier transform.
 * \item {\tt signalvec->data[n/2]} is the {\it real} Nyquist frequency component of the Fourier transform.
 * \item {\tt signalvec->data[k]} and {\tt signalvec->data[n-k],} for {\tt k=1,\ldots, n/2-1,} are
 * the real and imaginary parts of the Fourier transform at a frequency $k\Delta f=k/T,$ $T$ being
 * the duration of the signal and $\Delta f=1/T$ is the frequency resolution.
 * \end{itemize}
 *
 * \subsubsection*{Algorithm}
 *
 * %% A description of the method used to perform the calculation.
 *
 * The standard SPA is given by Eq.~(\ref{eq:InspiralFourierPhase:f2}).
 * We define a variable function pointer {\tt LALInspiralTaylorF2Phasing} and point
 * it to one of the {\texttt static} functions defined within this function
 * that explicitly calculates the Fourier phase at the PN order chosen by the user.
 * The reference points are chosen so that on inverse Fourier transforming
 * the time-domain waveform will
 * \begin{itemize}
 * \item be padded with zeroes in the first {\tt params->nStartPad} bins,
 * \item begin with a phase shift of {\tt params->nStartPhase} radians,
 * \item have an amplitude of ${\tt n} v^2.$
 * \end{itemize}
 * \subsubsection*{Uses}
 * \begin{verbatim}
   LALInspiralSetup
   LALInspiralChooseModel
   LALInspiralTaylorF2Phasing[0234567]PN
 * \end{verbatim}
 *
 * %% List of any external functions called by this function.
 * \begin{verbatim}
 * None
 * \end{verbatim}
 * \subsubsection*{Notes}
 *
 * %% Any relevant notes.
 *
 * If it is required to compare the output of this module with a time domain
 * signal one should use an inverse Fourier transform routine that packs data
 * in the same way as fftw. Moreover, one should divide the resulting inverse
 * Fourier transform by a factor ${\tt n}/2$ to be consistent with the
 * amplitude used in time-domain signal models.
 *
 * \vfill{\footnotesize\input{LALInspiralStationaryPhaseApprox2CV}}
 *
 **** </lalLaTeX> */

#include "LALInspiral.h"

#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

// NB: PROTOTYPE RESIDES IN LALINSPIRAL.H

NRCSID (LALINSPIRALSTATIONARYPHASEAPPROX2C, "$Id$");

/*  <lalVerbatim file="LALInspiralStationaryPhaseApprox2CP"> */
void
LALInspiralStationaryPhaseApprox2Test (
   LALStatus        *status,
   REAL4Vector      *signalvec,
   InspiralTemplate *params
   )
{ /* </lalVerbatim>  */
   REAL8 Oneby3, UNUSED h1, UNUSED h2, pimmc, f, v, df, shft, phi, amp0, amp, psif, psi;
   INT4 n, nby2, i, f0, fn;
   expnCoeffs ak;
   expnFunc func;
   //void (*LALInspiralTaylorF2PhasingTest)(InspiralTemplate *, REAL8 , REAL8 *) = NULL;

   INITSTATUS (status, "LALInspiralStationaryPhaseApprox2", LALINSPIRALSTATIONARYPHASEAPPROX2C);
   ATTATCHSTATUSPTR(status);
   ASSERT (signalvec,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->data,  status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (params, status, LALINSPIRALH_ENULL, LALINSPIRALH_MSGENULL);
   ASSERT (signalvec->length>2,  status, LALINSPIRALH_ECHOICE, LALINSPIRALH_MSGECHOICE);

   n = signalvec->length;
   nby2 = n/2;
   memset( &ak, 0, sizeof( ak ) );
   LALInspiralSetup(status->statusPtr, &ak, params);
   CHECKSTATUSPTR(status);
   LALInspiralChooseModel(status->statusPtr, &func, &ak, params);
   CHECKSTATUSPTR(status);

   Oneby3 = 1.L/3.L;
   df = params->tSampling/signalvec->length;
   pimmc = LAL_PI * params->totalMass * LAL_MTSUN_SI;
   f0 = params->fLower;
   fn = (params->fCutoff < ak.fn) ? params->fCutoff : ak.fn;
   v = pow(pimmc*f0, Oneby3);

   /* If we want to pad with zeroes in the beginning then the instant of
    * coalescence will be the chirp time + the duration for which padding
    * is needed. Thus, in the equation below nStartPad occurs with a +ve sign.
    * This code doesn't support non-zero start-time. i.e. params->startTime
    * should be necessarily zero.
    */
   shft = 2.L*LAL_PI * (ak.tn + params->nStartPad/params->tSampling + params->startTime);
   phi =  params->startPhase + LAL_PI/4.L;
   amp0 = params->signalAmplitude * ak.totalmass * pow(LAL_PI/12.L, 0.5L) * df;
/*
   Compute the standard stationary phase approximation.
*/
   h1 = signalvec->data[0] = 0.L;
   h2 = signalvec->data[nby2] = 0.L;
   for (i=1; i<nby2; i++) {
      f = i * df;
      if (f < f0 || f > fn)
      {
	      /*
	       * All frequency components below f0 and above fn are set to zero
	       */
	      signalvec->data[i] = 0.;
	      signalvec->data[n-i] = 0.;
      }
      else
      {
	      v = pow(pimmc * f, Oneby3);
	      LALInspiralTaylorF2PhasingTest(params, f, &psif);
	      psi = shft * f + phi + psif;
	      /*
		 dEnergybyFlux computes 1/(dv/dt) while what we need is 1/(dF/dt):
		 dF/dt=(dF/dv)(dv/dt)=-dEnergybyFlux/(dF/dv)=-dEnergybyFlux 3v^2/(pi m^2)
		 Note that our energy is defined as E=-eta v^2/2 and NOT -eta m v^2/2.
		 This is the reason why there is an extra m in the last equation above
		 amp = amp0 * pow(-dEnergybyFlux(v)/v^2, 0.5) * v^2;
		     = amp0 * pow(-dEnergybyFlux(v), 0.5) * v;
	      */
	      amp = amp0 * pow(-func.dEnergy(v,&ak)/func.flux(v,&ak),0.5L) * v;
	      signalvec->data[i] = (REAL4) (amp * cos(psi));
	      signalvec->data[n-i] = (REAL4) (-amp * sin(psi));

      }
      /*
	 printf ("%e %e \n", v, psif);
	 printf ("%e %e %e %e %e\n", f, pow(h1,2.)+pow(h2,2.), h2, psi, psif);
	 printf ("&\n");
       */

   }
   params->fFinal = fn;
   DETATCHSTATUSPTR(status);
   RETURN(status);
}

void LALInspiralTaylorF2PhasingTest(
                          InspiralTemplate *params,
                          REAL8 f,
                          REAL8 *psif)
{
    
    // TOTAL PHASE
    REAL8 phaseParams[10] = {0.0};
    
    // FILL PHASE PARAMETERS
    TaylorF2fillPhaseParams(params, phaseParams, 3, 0.0);
    
    switch (params->order)
    {
        case LAL_PNORDER_NEWTONIAN:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
            );
        case LAL_PNORDER_HALF:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
            );
            break;
        case LAL_PNORDER_ONE:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
            );
            break;
        case LAL_PNORDER_ONE_POINT_FIVE:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
             + phaseParams[3]*pow(f, -2.0/3.0)
            );
            break;
        case LAL_PNORDER_TWO:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
             + phaseParams[3]*pow(f, -2.0/3.0)
             + phaseParams[4]*pow(f, -1.0/3.0)
            );
            break;
        case LAL_PNORDER_TWO_POINT_FIVE:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
             + phaseParams[3]*pow(f, -2.0/3.0)
             + phaseParams[4]*pow(f, -1.0/3.0)
            + phaseParams[5]
            + phaseParams[6]*log(f)
            );
            break;
        case LAL_PNORDER_THREE:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
             + phaseParams[3]*pow(f, -2.0/3.0)
             + phaseParams[4]*pow(f, -1.0/3.0)
            + phaseParams[5]
            + phaseParams[6]*log(f)
            + phaseParams[7]*pow(f, 1.0/3.0)
            + phaseParams[8]*pow(f, 1.0/3.0)*log(f)
            );
            break;
        case LAL_PNORDER_THREE_POINT_FIVE:
            *psif =  
            (  phaseParams[0]*pow(f, -5.0/3.0)
             + phaseParams[1]*pow(f, -4.0/3.0)
             + phaseParams[2]*pow(f, -3.0/3.0)
             + phaseParams[3]*pow(f, -2.0/3.0)
             + phaseParams[4]*pow(f, -1.0/3.0)
            + phaseParams[5]
            + phaseParams[6]*log(f)
            + phaseParams[7]*pow(f, 1.0/3.0)
            + phaseParams[8]*pow(f, 1.0/3.0)*log(f)
            + phaseParams[9]*pow(f, 2.0/3.0)
            );
            break;
        default:
            printf("INVALID PN ORDER!");
            exit(-1);
    }
    
    return;
}

void TaylorF2fillPhaseParams(
                                         InspiralTemplate *params,
                                         REAL8 *phaseParams,
                                         INT4    testParam,
                                         REAL8 testParamValue)
{
    
    // SYSTEM DEPENDENT PARAMETER - DUMMIES, NEED TO GET FROM PARAMETER STRUCTURES
    REAL8 mtot = params->totalMass;
    REAL8 eta = params->eta;
    
    REAL8 twopimtot = 2*LAL_TWOPI*mtot*LAL_MTSUN_SI;
    REAL8 comprefac = 3.0/(256.0*eta);
    
    // POPULATE INDIVIDUAL PHASE PARAMETERS
    // SEE arXiv:gr-qc/0411146
    /*
    phaseParams[0] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-5.0/3.0)* 1.0; //phi0
    phaseParams[1] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-4.0/3.0)* 0.0; //phi1
    phaseParams[2] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-3.0/3.0)* 20.0/9.0*(743.0/336.0 + 11.0/4.0*pow(eta,2.0)); //phi2
    phaseParams[3] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-2.0/3.0)* -16.0*LAL_PI; //phi3
    phaseParams[4] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-1.0/3.0)* 10.0*(3058673.0/1016064.0 + 5429.0/1008.0*eta + 617.0/144.0*pow(eta,2.0)); // phi4
    phaseParams[5] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-0.0/3.0)* LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5
    phaseParams[6] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,-0.0/3.0)* 3.0*LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5l
    phaseParams[7] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,1.0/3.0)* ((11583231236531.0/4694215680.0 - 640.0/3.0*pow(LAL_PI, 2.0) - 6848.0/21.0*LAL_GAMMA) + eta*(-15335597827.0/3048192.0 + 2255.0/12.0*pow(LAL_PI, 2.0) - 1760.0/3.0*LAL_THETA + 12320.0/9.0*LAL_LAMBDA) + 76055.0/1728.0*pow(eta, 2.0) - 127825.0/1296.0*pow(eta, 3.0) + -6848.0/21.0*log(4*pow(LAL_PI*mtot, 1.0/3.0))); //phi6 + some factors originally coming from phi6l
    phaseParams[8] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,1.0/3.0)* -6848.0/21.0; //phi6l
    phaseParams[9] = 3.0/(128.0*eta)*pow(LAL_PI*LAL_MTSUN_SI*mtot,2.0/3.0)* LAL_PI*(77096675.0/254016.0 + 378515.0/1512.0*eta - 74045.0/756.0*pow(eta, 2.0)); //phi7
     */
    
    
    // SEE arXiv:1005.0304
    phaseParams[0] = comprefac*pow(twopimtot, -5.0/3.0); //phi0
    phaseParams[1] = comprefac*pow(twopimtot,-4.0/3.0)* 0.0; //phi1
    phaseParams[2] = comprefac*pow(twopimtot,-3.0/3.0)* (3715.0/756.0 + 55.0/9.0*eta); //phi2
    phaseParams[3] = comprefac*pow(twopimtot,-2.0/3.0)* -16.0*LAL_PI; //phi3
    phaseParams[4] = comprefac*pow(twopimtot,-1.0/3.0)* (15293365.0/508032.0 + 27145.0/504.0*eta + 3085.0/72.0*pow(eta, 2.0)); // phi4
    phaseParams[5] = comprefac*pow(twopimtot,-0.0/3.0)* LAL_PI*(38645.0/756.0 - 65.0/9.0*eta)*(1.0+log(twopimtot*pow(6.0, 1.5))); //phi5
    phaseParams[6] = comprefac*pow(twopimtot,-0.0/3.0)* LAL_PI*(38645.0/756.0 - 65.0/9.0*eta); //phi5l
    phaseParams[7] = comprefac*pow(twopimtot,1.0/3.0)* ((11583231236531.0/4694215680.0 - 640.0/3.0*pow(LAL_PI, 2.0) - 6848.0/21.0*LAL_GAMMA) + eta*(-15737765635.0/3048192.0 + 2255.0/12.0*pow(LAL_PI, 2.0)) + 76055.0/1728.0*pow(eta, 2.0) - 127825.0/1296.0*pow(eta, 3.0) + -6848.0/63.0*log(64.0*twopimtot)); //phi6
    phaseParams[8] = comprefac*pow(twopimtot,1.0/3.0)* -6848.0/63.0; //phi6l
    phaseParams[9] = comprefac*pow(twopimtot,2.0/3.0)* LAL_PI*(77096675.0/254016.0 + 378515.0/1512.0*eta - 74045.0/756.0*pow(eta, 2.0)); //phi7
    
    if      (testParam==0){phaseParams[0] += testParamValue;}
    else if (testParam==1){phaseParams[1] += testParamValue;}
    else if (testParam==2){phaseParams[2] += testParamValue;}
    else if (testParam==3){phaseParams[3] += testParamValue;}
    else if (testParam==4){phaseParams[4] += testParamValue;}
    else if (testParam==5){phaseParams[5] += testParamValue;}
    else if (testParam==6){phaseParams[6] += testParamValue;}
    else if (testParam==7){phaseParams[7] += testParamValue;}
    else if (testParam==8){phaseParams[8] += testParamValue;}
    else if (testParam==9){phaseParams[9] += testParamValue;}
    
    return;
}
