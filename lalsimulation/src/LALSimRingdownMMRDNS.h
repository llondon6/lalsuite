#ifndef _LALSIM_RINGDOWN_MMRDNS_H
#define _LALSIM_RINGDOWN_MMRDNS_H


/* ************************************************************  */
/*
 * Copyright (C) 2016 Lionel London
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

 #if defined(__cplusplus)
 extern "C" {
 #elif 0
 } /* so that editors will match preceding brace */
 #endif

/* Include the desired Libs */
#include <stdbool.h>
#include <math.h>
#include <complex.h>
/* LAL specific libs  */

#include <lal/LALDatatypes.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include <lal/LALSimIMR.h>
#include <lal/LALConfig.h>
#include <lal/LALConstants.h>
#include <lal/Date.h>
#include <lal/FrequencySeries.h>
#include <lal/Units.h>
#include <lal/LALSimInspiral.h>

/* ---------------------------------------- */
/* General model methods and parameters     */
/* ---------------------------------------- */

/* Mode included within the model */
static const int MMRDNS_MODES[9][3] =  { {2,2,0},
                                         {2,2,1},
                                         {3,3,0},
                                         {3,3,1},
                                         {4,4,0},
                                         {5,5,0},
                                         {2,1,0},
                                         {3,2,0},
                                         {4,3,0}
                                       };

/* --> Domain mapping for final spin, Kappa */
static double KAPPA( double jf, int l, int m );
/* --> Fits for complex QNM frequencies */
static double CW07102016( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n );      /* Overtone Number */
/* --> Fits for complex QNM separation constants; needed to calculate spheroidal harmonics */
static double SC07102016( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n );      /* Overtone Number */
/* --> Fits for spheroidal harmonic normalization constants; needed to calculate spheroidal harmonics */
static double CC09102016( double kappa,  /* Domain mapping for remnant BH's spin */ (Dimensionless)
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n );      /* Overtone Number */
/* --> Final Mass fit arXiv:1404.3197 */
static double MF07102016( double eta );
/* --> Final Spin fit arXiv:1404.3197 */
static double JF07102016( double eta );

/* ------------------------------------------------ # */
/* Angular parameter functions
/* ------------------------------------------------ # */
static double K1( int m, int s );
static double K2( int m, int s );
static complex double ALPHA( int m, int s, int p );
static complex double BETA( int m, int s, int p, complex double aw, complex double A_lm, complex double w_lm );
static complex double GAMMA( int m, int s, int p, complex double aw );

/*
* Spheroical Harmonic Functions (Leaver's Formulation circa 1986/85)
*/
static double XLALSpinWeightedSpheroidalHarmonic( double jf,           /* Spin of remnant */
                   int l, int m, int n, /* QNM indeces */
                   double theta,        /* polar angle */
                   double phi,          /* azimuthal angle */
                   bool norm           /* boolean toggle for normalization  */
                 );

/* Gnerate Time domain ringdown waveform  */
int XLALSimRingdownMMRDNSTD(
        REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        REAL8 phiRef,                   /**< orbital phase at reference pt. */
        REAL8 inclination,              /**< inclination angle */
        REAL8 deltaT,                   /**< sampling interval (s) */
        REAL8 m1,                       /**< mass of companion 1 (kg) */
        REAL8 m2,                       /**< mass of companion 2 (kg) */
        REAL8 r,                        /**< distance of source (m) */
        const LALSimInspiralTestGRParam *extraParams /**< linked list containing the extra testing GR parameters */
        );


/* ************************************************************  */
#endif	/* of #ifndef _LALSIM_RINGDOWN_MMRDNS_H */
