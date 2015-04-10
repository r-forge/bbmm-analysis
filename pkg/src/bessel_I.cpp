/**
 * The code in this file is based on code from the GNU Scientific library,
 * subject to the following license:
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "bessel.h"

double bessel_I0_scaled(double x) {
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=exp(-ax)*(1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
         +y*(0.2659732+y*(0.360768e-1+y*0.45813e-2))))));
   } else {
      y=3.75/ax;
      ans=(0.39894228+y*(0.1328592e-1
         +y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2
         +y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1
         +y*0.392377e-2)))))))) / sqrt(ax);
   }
   return ans;
}

double bessel_I1_scaled(double x) {
   double ax,ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y=x/3.75,y=y*y;
      ans=exp(-ax)*(ax*(0.5+y*(0.87890594+y*(0.51498869+y*(0.15084934
         +y*(0.2658733e-1+y*(0.301532e-2+y*0.32411e-3)))))));
   } else {
      y=3.75/ax;
      ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
         -y*0.420059e-2));
      ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
         +y*(0.163801e-2+y*(-0.1031555e-1+y*ans))));
      ans /= sqrt(ax);
   }
   return x < 0.0 ? -ans : ans;
}

double bessel_In_scaled(int n, const double x) {
  const double ax = fabs(x);

  n = abs(n);  /* I(-n, z) = I(n, z) */

  /* CHECK_POINTER(result) */

  if(n == 0) {
    return bessel_I0_scaled(x);
  } else if(n == 1) {
    return bessel_I1_scaled(x);
  } else if(x == 0.0) {
    return 0.0;
  } else if(x*x < 10.0*(n+1.0)/M_E) {
    double ex   = exp(-ax);
    double t = ex * bessel_IJ_taylor(n, ax, 1, 50, DBL_EPSILON);
    if(x < 0.0 && n%2) t = -t;
    return t;
  } else { // if(n < 150 && ax < 1e7) {
    double rat = bessel_I_CF1_ser(n, ax);
    double Ikp1 = rat * sqrt(DBL_MIN);
    double Ik   = sqrt(DBL_MIN);
    double Ikm1;
    int k;
    for(k=n; k >= 1; k--) {
      Ikm1 = Ikp1 + 2.0*k/ax * Ik;
      Ikp1 = Ik;
      Ik   = Ikm1;
    }
    double t = bessel_I0_scaled(ax) * (sqrt(DBL_MIN) / Ik);
    if(x < 0.0 && n%2) t = -t;
    return t;
//  } else if( GSL_MIN( 0.29/(n*n), 0.5/(n*n + x*x) ) < 0.5*rootn(DBL_EPSILON, 3)) {
//    int stat_as = gsl_sf_bessel_Inu_scaled_asymp_unif_e((double)n, ax, result);
//    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
//    return stat_as;
//  } else {
//    const int nhi = 2 + (int) (1.2 / rootn(DBL_EPSILON, 6));
//    gsl_sf_result r_Ikp1;
//    gsl_sf_result r_Ik;
//    int stat_a1 = gsl_sf_bessel_Inu_scaled_asymp_unif_e(nhi+1.0,     ax, &r_Ikp1);
//    int stat_a2 = gsl_sf_bessel_Inu_scaled_asymp_unif_e((double)nhi, ax, &r_Ik);
//    double Ikp1 = r_Ikp1.val;
//    double Ik   = r_Ik.val;
//    double Ikm1;
//    int k;
//    for(k=nhi; k > n; k--) {
//      Ikm1 = Ikp1 + 2.0*k/ax * Ik;
//      Ikp1 = Ik;
//      Ik   = Ikm1;
//    }
//    result->val = Ik;
//    result->err = Ik * (r_Ikp1.err/r_Ikp1.val + r_Ik.err/r_Ik.val);
//    if(x < 0.0 && GSL_IS_ODD(n)) result->val = -result->val;
//    return GSL_ERROR_SELECT_2(stat_a1, stat_a2);
  }
}
