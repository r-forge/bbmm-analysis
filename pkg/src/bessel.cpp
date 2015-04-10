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

#include <math.h>
#include <float.h>
#include "bessel.h"

double bessel_IJ_taylor(const int n, const double x,
                        const int sign,
                        const int kmax,
                        const double threshold) {
  if(x == 0.0) {
  	return (n == 0.0) * 1.0; // conditional return value without branching
  } else {
    int k;
  	/* The original implementation defines prefactor := (x/2)^n / Gamma(n+1).
  	 * For integer n, this equals (x/2)^n / n!
  	 */
    double prefactor = pow(0.5*x, n);
    for (k=1; k<=n; k++) { prefactor /= k; } // divide by n!.
    // Future improvement: limit the number of divisions to lose less precision
    
    double sum = 1.0;

    /* Evaluate the sum.
     * [Abramowitz+Stegun, 9.1.10]
     * [Abramowitz+Stegun, 9.6.7]
     */
    {
      const double y = sign * 0.25 * x*x;
      double term = 1.0;

      for(k=1; k<=kmax; k++) {
        term *= y/((n+k)*k);
        sum += term;
        if(fabs(term/sum) < threshold) break;
      }
    }

	return prefactor * sum;
  }
}

/* Evaluate the continued fraction CF1 for I_{nu+1}/I_nu
 * using Gautschi (Euler) equivalent series.
 */
double bessel_I_CF1_ser(const int n, const double x) {
  const int maxk = 20000;
  double tk   = 1.0;
  double sum  = 1.0;
  double rhok = 0.0;
  int k;

  for(k=1; k<maxk; k++) {
    double ak = 0.25 * (x/(n+k)) * x/(n+k+1);
    rhok = -ak*(1.0 + rhok)/(1.0 + ak*(1.0 + rhok));
    tk  *= rhok;
    sum += tk;
    if(fabs(tk/sum) < DBL_EPSILON) break;
  }

  return x/(2*(n+1)) * sum;
}


