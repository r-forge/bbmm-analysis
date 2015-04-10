/**
 * Various utility function for internal use in the package
 */

#ifndef __MA_UTIL_H__
#define __MA_UTIL_H__

#include <R.h>
#include <Rinternals.h>
#include <math.h>
#include "vec2d.h"

// d2 is the squared distance between mu and p in the next function
inline double pdf_norm2d(double d2, double var)
{
	// straightforward implementation of the expression for the pdf of N(mu, sigma^2)
    return exp(d2 / (-2.0 * var)) / (2 * M_PI * var);
}

inline double pdf_norm2d(Vec2D<double> p, Vec2D<double> mu, double var)
{
	// Call the squared-distance based function for the PDF
    return pdf_norm2d((p-mu).norm2(), var);
}

inline double cdf_norm2d(Vec2D<double> p, double hCellSize, Vec2D<double> mu, double var)
{
	double s2V = sqrt(2 * var);
	// P(X inside cell of size cellSize centered at p) =
	// 	P(p.x - cs/2 <= X.x <= p.x + cs/2) P(p.y - cs/2 <= X.y <= p.y + cs/2)
	return (erf((p.x + hCellSize - mu.x)/s2V) - erf((p.x - hCellSize - mu.x)/s2V))
			* (erf((p.y + hCellSize - mu.y)/s2V) - erf((p.y - hCellSize - mu.y)/s2V));
}

double marcumQ(const double alpha, const double beta);
inline double cdf_rice(double x, double nu, double sigma) {
	return 1.0 - marcumQ(nu/sigma, x/sigma);
}

/**
 * Functions relating to the integration of diffusion coefficients
 */

/**
 * construct an object suitable for passing to integrate() given a diffusion function.
 */
SEXP initIntegrate(SEXP diff);

/**
 * Do the actual integration over the interval [l,h],
 * given the result of a call to initIntegrate()
 */
double integrate(SEXP call, double l, double h);

#endif // __MA_UTIL_H__

