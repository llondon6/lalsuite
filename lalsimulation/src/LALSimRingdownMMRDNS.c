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


/* .................... */
/* HEADER SECTION       */
/* .................... */
#ifdef __GNUC__
#define UNUSED __attribute__ ((unused))
#else
#define UNUSED
#endif

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_min.h>
#include <gsl/gsl_spline.h>
#include <lal/Units.h>
#include <lal/SeqFactories.h>
#include <lal/Date.h>
#include <lal/StringInput.h>
#include <lal/Sequence.h>
#include <lal/FileIO.h>

#include <lal/XLALError.h>

#include <lal/LALConstants.h>
#include <lal/LALStdio.h>
#include <lal/LALSimSphHarmSeries.h>
#include <lal/LALStdlib.h>
#include <lal/LALSimInspiral.h>
#include "check_series_macros.h"

#include <lal/TimeSeries.h>
#include <lal/FrequencySeries.h>

#include "LALSimRingdownMMRDNS.h"
#include "LALSimRingdownCW.h"
#include <lal/LALSimInspiralTestGRParams.h>

/*
* -------------------------------------------------------------------------------- *
* Low level models: QNM Frequencies, Separation Constants and Spheroidal Harmonics
* -------------------------------------------------------------------------------- *
*/

/*
* QNM Separation Constants: Note that name encodes date of writing
*/
static double complex SC07102016( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless) */
                          int l,        /* Polar eigenvalue */
                          int input_m,  /* Azimuthal eigenvalue */
                          int n ) {     /* Overtone Number */

  /* Predefine powers to increase efficiency */
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;
  double kappa7 = kappa6 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken */
  int m = abs(input_m);

  /**/
  complex double j = _Complex_I;

  /* Initialize the answer */
  double complex ans;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit. */
    ans = 0.55262405 + 6.54272463*cexp(0.24443847*j)*kappa + 5.94664565*cexp(3.88409012*j)*kappa2 + 5.39298183*cexp(1.01651284*j)*kappa3 + 3.58701474*cexp(4.53395559*j)*kappa4 + 1.36858235*cexp(1.57079633*j)*kappa5 + 0.18520700*cexp(4.71238898*j)*kappa6 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.55229247 + 7.94074969*cexp(0.64081239*j)*kappa + 12.55567057*cexp(4.41980669*j)*kappa2 + 13.68518711*cexp(1.48039237*j)*kappa3 + 10.43884041*cexp(4.72599435*j)*kappa4 + 4.20731453*cexp(1.57079633*j)*kappa5 + 0.76232588*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.  */
      ans = 8.18542769*cexp(6.27603422*j) + 1.55192720*cexp(1.79088081*j)*kappa + 8.94654695*cexp(5.18681710*j)*kappa2 + 28.66050158*cexp(1.63658858*j)*kappa3 + 60.77789497*cexp(4.72114050*j)*kappa4 + 72.13239907*cexp(1.57079633*j)*kappa5 + 45.38115278*cexp(4.71238898*j)*kappa6 + 11.84706755*cexp(1.57079633*j)*kappa7;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 13.05294185 + 9.23462388*cexp(0.14179514*j)*kappa + 7.09045393*cexp(3.69184561*j)*(kappa2) + 6.46711175*cexp(0.89254551*j)*(kappa3) + 4.96905278*cexp(4.43853588*j)*(kappa4) + 2.62299932*cexp(1.57079633*j)*(kappa5) + 0.58168681*cexp(4.71238898*j)*(kappa6);

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit. */
      ans = 3.10089518*cexp(6.25822093*j) + 2.69208437*cexp(1.95853947*j)*kappa + 16.58575360*cexp(4.98423605*j)*kappa2 + 57.84090876*cexp(1.63720921*j)*kappa3 + 118.21761290*cexp(4.72674943*j)*kappa4 + 135.93985738*cexp(1.57079633*j)*kappa5 + 82.81742189*cexp(4.71238898*j)*kappa6 + 20.85173245*cexp(1.57079633*j)*kappa7;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 5.70465254 + 7.94433155*cexp(0.18039136*j)*kappa + 6.55099749*cexp(3.77926384*j)*kappa2 + 6.31422768*cexp(0.93863733*j)*kappa3 + 4.81214531*cexp(4.46906976*j)*kappa4 + 2.38927043*cexp(1.57079633*j)*kappa5 + 0.48077965*cexp(4.71238898*j)*kappa6;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 5.70318420 + 8.94926548*cexp(0.49834140*j)*kappa + 12.70528736*cexp(4.31772419*j)*kappa2 + 15.63533560*cexp(1.39390017*j)*kappa3 + 14.19057659*cexp(4.66913674*j)*kappa4 + 7.33238119*cexp(1.57079633*j)*kappa5 + 1.53701758*cexp(4.71238898*j)*kappa6;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 15.28866348 + 0.75297352*cexp(0.22048290*j)*kappa + 3.64936150*cexp(0.61644055*j)*kappa2 + 8.02530641*cexp(4.82756576*j)*kappa3 + 12.47205664*cexp(1.67334685*j)*kappa4 + 10.30282199*cexp(4.71238898*j)*kappa5 + 3.52885679*cexp(1.57079633*j)*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 22.52292196 + 10.44137664*cexp(0.11607502*j)*kappa + 7.79707643*cexp(3.61247422*j)*kappa2 + 6.59989026*cexp(0.83792606*j)*kappa3 + 4.90367451*cexp(4.40545635*j)*kappa4 + 2.59913853*cexp(1.57079633*j)*kappa5 + 0.58985077*cexp(4.71238898*j)*kappa6;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  /* If m<0, then take the conjugate */
  if ( input_m < 0 ) {
    /**/
    ans = conj( ans );
  }

  return ans;

} /* END of SC07102016 */


/*
* Spheroidal Harmonic Normalization Constants: Note that name encodes date of writing
*/
UNUSED static double CC09102016( double kappa,  /* Domain mapping for remnant BH's spin (Dimensionless)*/
                          int l,        /* Polar eigenvalue*/
                          int input_m,  /* Azimuthal eigenvalue*/
                          int n ) {     /* Overtone Number*/

  /* Predefine powers to increase efficiency*/
  double kappa2 = kappa  * kappa;
  double kappa3 = kappa2 * kappa;
  double kappa4 = kappa3 * kappa;
  double kappa5 = kappa4 * kappa;
  double kappa6 = kappa5 * kappa;

  /* NOTE that |m| will be used to determine the fit to use, and if input_m < 0, then a conjugate will be taken*/
  int m = abs(input_m);

  /* Initialize the answer*/
  double complex ans;

  /* Use If-Else ladder to determine which mode function to evaluate */
  if ( 2==l && 2==m && 0==n  ){

    /* Fit for (l,m,n) == (2,2,0). This is a zero-damped mode in the extremal Kerr limit.*/
    ans = 7.86366171 - 3.61447483*kappa + 3.48996689*kappa2 - 2.29347705*kappa3 + 0.74425069*kappa4 ;

    } else if ( 2==l && 2==m && 1==n ) {

      /* Fit for (l,m,n) == (2,2,1). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 7.86298703 - 3.59872285*kappa + 2.88459437*kappa2 - 0.92740734*kappa3 - 0.04445478*kappa4;

    } else if ( 3==l && 2==m && 0==n ) {

      /* Fit for (l,m,n) == (3,2,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.74845717 - 0.08157463*kappa + 1.03748092*kappa2 - 3.27926931*kappa3 + 7.24584503*kappa4 - 7.41316799*kappa5 + 3.06056035*kappa6;

    } else if ( 4==l && 4==m && 0==n ) {

      /* Fit for (l,m,n) == (4,4,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 1.75389888 + 1.00111258*kappa + 1.55498487*kappa2 - 1.22344804*kappa3 + 1.64621074*kappa4;

    } else if ( 2==l && 1==m && 0==n ) {

      /* Fit for (l,m,n) == (2,1,0). This is NOT a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.04393302 - 0.06877527*kappa + 0.87671129*kappa2 - 3.92206769*kappa3 + 8.59631959*kappa4 - 8.52199526*kappa5 + 3.31150324*kappa6;

    } else if ( 3==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (3,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51631915 + 0.16499714*kappa + 1.30114387*kappa2 - 0.83622153*kappa3 + 0.82020713*kappa4;

    } else if ( 3==l && 3==m && 1==n ) {

      /* Fit for (l,m,n) == (3,3,1). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 3.51530809 + 0.19285707*kappa + 0.96814190*kappa2 - 0.00547882*kappa3 + 0.24982172*kappa4;

    } else if ( 4==l && 3==m && 0==n ) {

      /* Fit for (l,m,n) == (4,3,0). This is a zero-damped mode in the extremal Kerr limit.*/
      ans = 0.39542385 - 0.09918352*kappa + 1.52850262*kappa2 - 5.09932727*kappa3 + 10.95647104*kappa4 - 10.99914124*kappa5 + 4.52212985*kappa6;

    } else if ( 5==l && 5==m && 0==n ) {

      /* Fit for (l,m,n) == (5,5,0). This is a zero-damped mode in the extremal Kerr limit. */
      ans = 0.91349889 + 0.89568178*kappa + 2.54404526*kappa2 - 2.82437113*kappa3 + 3.28143852*kappa4;

    } else {

      /**/
      ans  = 0.0;

  } /* END of IF-ELSE Train for QNM cases */

  return ans;

} /* END of CC09102016 */


/*
* Final Spin (Dimensionless)
*/
UNUSED static double JFNS07102016( double eta ) {
  /* Implement Final Spin Fit from arXiv:1404.3197 */
  return eta*( 3.4339 - eta*( 3.7988 + eta*(5.7733 - 6.3780*eta) ) );
}

/*
* Final Mass (Dimensionless: Relative to M=1 initially)
*/
UNUSED static double MFNS07102016( double eta ) {
  /* Implement Final Mass Fit from arXiv:1404.3197 */
  return 1.0 - eta*(0.046297 - eta*(0.71006 + eta*(1.5028 - eta*(4.0124 - eta*0.28448))));
}


/* ------------------------------------------------
          Angular parameter functions
 ------------------------------------------------ */
UNUSED static double K1( int m, int s ){
  return 0.5*abs(m-s);
}
UNUSED static double K2( int m, int s ){
  return 0.5*abs(m+s);
}
UNUSED static complex double ALPHA( int m, int s, int p ){
  /**/
  double k1 = 0.5*abs(m-s);
  return -2.0*(p+1.0)*(p+2.0*k1+1.0);
}
UNUSED static complex double BETA( int m, int s, int p, complex double aw, complex double A_lm ){
  /**/
  double k1 = K1(m,s);
  double k2 = K2(m,s);
  return  p*(p-1.0)+2.0*p*(k1+k2+1.0-2.0*aw) - ( 2.0*aw*(2.0*k1+s+1.0)-(k1+k2) * (k1+k2+1.0) ) - ( aw*aw + s*(s+1.0) + A_lm);
}
UNUSED static complex double GAMMA( int m, int s, int p, complex double aw ){
  /**/
  double k1 = K1(m,s);
  double k2 = K2(m,s);
  return 2.0*aw*(p+k1+k2+s);
}

/*
* Spheroical Harmonic Functions (Leaver's Formulation circa 1986/85)
*/
UNUSED static complex double XLALSpinWeightedSpheroidalHarmonic( double jf,           /* Spin of remnant */
                   int l, int m, int n, /* QNM indeces */
                   double theta,        /* polar angle */
                   double phi          /* azimuthal angle */
                 ) {

  /* Set spin weight */
  const int s = -2;

  /* Use tabulated cw and sc values from the core package*/
  double complex cw, sc, aw;
  double kappa = KAPPA(jf,l,m);

  cw = CW07102016( kappa, l, m, n );
  sc = SC07102016( kappa, l, m, n );

  /* Define dimensionless deformation parameter */
  aw = jf*cw;

  /* ------------------------------------------------
      Calculate the angular eighenfunction
   ------------------------------------------------ */

  /* Variable map for theta */
  double u = cos(theta);

  /* the non-sum part of eq 18 */
  complex double X = 1.0;
  X = X * cexp(aw*u) * pow(1.0+u, K1(m,s) );
  X = X * pow(1.0-u,K2(m,s));

  /* NOTE that here we apply the normalization constant */
  X = X / CC09102016(kappa,l,m,n);

  /* initial series values */
  double a0 = 1.0;

  double a1 = -BETA( m, s, 0, aw, sc )/ALPHA( m, s, 0 );

  /* the sum part */
  complex double Y = a0;
  complex double S;
  complex double dY;
  double xx;
  complex double a2;
  int done = 0;
  int k, j;
  Y = Y + a1*(1.0+u);
  k = 1;
  int kmax = 2e3;
  double et = 1e-8;
  while ( ! done ) {
    k += 1;
    j = k-1;
    a2 = -1.0*( BETA( m, s, j, aw, sc )*a1 + GAMMA(m,s,j,aw)*a0 ) / ALPHA(m,s,j);
    dY = pow(a2*(1.0+u),k);
    Y += dY;
    xx = cabs( dY );

    done = (k>=l) && ( (xx<et && k>30) || k>kmax );
    done = done || xx<et;
    a0 = a1;
    a1 = a2;
  }

  /* together now */
  S = X*Y*cexp( _Complex_I * m * phi );

  /* Use the same sign convention as spherical harmonics
  e.g http://en.wikipedia.org/wiki/Spin-weighted_spherical_harmonics#Calculating */
  double C = 1.0;
  C = C * pow(-1,fmax(-m,-s)) * pow(-1,l);
  S = C * S;

  /**/
  return S;

} /* END of Spheroical Harmonic function */

/**/
int XLALSimRingdownMMRDNSTD(
        UNUSED REAL8TimeSeries **hplus,        /**< OUTPUT h_+ vector */
        UNUSED REAL8TimeSeries **hcross,       /**< OUTPUT h_x vector */
        UNUSED REAL8 phiRef,                   /**< orbital phase at reference pt. */
        UNUSED REAL8 inclination,              /**< inclination angle */
        UNUSED REAL8 deltaT,                   /**< sampling interval (s) */
        UNUSED REAL8 m1,                       /**< mass of companion 1 (kg) */
        UNUSED REAL8 m2,                       /**< mass of companion 2 (kg) */
        UNUSED REAL8 r,                        /**< distance of source (m) */
        UNUSED const LALSimInspiralTestGRParam *nonGRparams /**< linked list containing the extra testing GR parameters */
        ){

        /* Declarations */
        double T0 = 0;
        double kappa = KAPPA( 0.68, 2, 2 );
        double w22 = creal( CW07102016(kappa,2,2,0) );
        double tau22 = cimag( CW07102016(kappa,2,2,0) );
        int waveform_length = 1000;


        /* Create the plus and cross polarization vectors */
        *hplus = XLALCreateREAL8TimeSeries("hplus", T0, 0.0, deltaT, &lalStrainUnit, waveform_length);
        *hcross = XLALCreateREAL8TimeSeries("hcross", T0, 0.0, deltaT, &lalStrainUnit, waveform_length);


  return 0;
}


/* XLALSimRingdownSingleModeTD: Time domain waveformgenerator for single QNM without angular dependence (i.e. this function generates multipole moments only ). In  */
int XLALSimRingdownGenerateSingleModeTD(
  UNUSED COMPLEX16TimeSeries **hk,  /**< OUTPUT vector for QNM time series */
  UNUSED REAL8 T0,                            /**< Start time  */
  UNUSED REAL8 deltaT,                        /**< sampling interval (s) */
  UNUSED REAL8 Nsamples,                      /**< Number of time domain samples */
  UNUSED complex double Ak,                   /**< COMPLEX QNM Strain Amplitude */
  UNUSED complex double CWk                   /**< COMPLEX QNM Frequency */
){

  /*  */
  *hk = XLALCreateCOMPLEX16TimeSeries( "", t0, 0.0, deltaT, &lalStrainUnit, length);
  memset( (*hk)->data->data, 0, length*sizeof(COMPLEX16) ) ;

  /**/
  COMPLEX16 cwdt = 0;

  /**/
  cwdt = 0; /* CWk*deltaT; */

  /**/
  for ( k=0, k < Nsamples, ++k ){

    /**/
    COMPLEX16 h;
    h = Ak * cexp(-I * cwdt * k);
    (*hk)->data->data[k] = h;

  }

  return 0;

}
