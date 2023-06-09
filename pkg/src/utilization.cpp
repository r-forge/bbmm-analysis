#include "vec2d.h"
#include <math.h>
#include <algorithm>

#include "encounter.h"
#include "util.h"
#include <Rinternals.h>
#include <Rdefines.h>

#define Dx(i) (data[4*(i)    ])
#define Dy(i) (data[4*(i) + 1])
#define Dv(i) (data[4*(i) + 2])
#define Dt(i) (data[4*(i) + 3])

#define Rx(i) (res[4*(i)    ])
#define Ry(i) (res[4*(i) + 1])
#define Rv(i) (res[4*(i) + 2])
#define Rw(i) (res[4*(i) + 3])

extern "C" {

SEXP UDTimesteps(SEXP bData, SEXP diff, SEXP timeSteps) {
	SEXP result, dims,
			dInt, dIntL, dIntH; // call to integrate, pointers to lower and uppper bounds
	diff = VECTOR_ELT(diff, 0);

	PROTECT(dInt = initIntegrate(diff));
	
	// Extract the number of relocations in the input data
	PROTECT(dims = GET_DIM(bData));
	int nLoc = INTEGER(dims)[1];
	UNPROTECT(1);    //dims

	// Create the result matrix
	PROTECT(result = allocMatrix(REALSXP, 4, length(timeSteps)));

	// Extract the data from certain objects
	double *ts = REAL(timeSteps);
	double *data = REAL(bData);
	double *res = REAL(result);
	
	off_t i = 1;
	double diff0T = integrate(dInt, Dt(i-1), Dt(i));
	for (off_t k = 0; k < length(timeSteps); k++) {
		double t = ts[k];
		
		/**
		 * The diffusion coefficient is allowed to vary over time. The computations
		 * are derived in Sijben, Stef. "Computational Movement Analysis Using 
		 * Brownian Bridges." (2013). Master's thesis, TU Eindhoven.
		 */
		while (t > Dt(i) && i < (nLoc-1)) {
			i++;
			// Set the limits of integration, then compute the integral
			diff0T = integrate(dInt, Dt(i-1), Dt(i));
		}
		double alpha;
		if (diff0T > 0) {
			alpha = integrate(dInt, Dt(i-1), t) / diff0T;
		} else {
			// diffusion coefficient is zero everywhere, fall back to "old" method
			alpha = (t - Dt(i-1)) / (Dt(i) - Dt(i-1));
		}
		
		if (t < Dt(0)) {
			Rx(k) = Dx(0);
			Ry(k) = Dy(0);
			Rv(k) = Dv(0);
			Rw(k) = 0.5 * std::max(0.0, ts[k+1] - Dt(0)); // ignore the part before Dt(0), as we have no data there
		} else if (t > Dt(nLoc - 1)) {
			// We know here that i == length(timeSteps) - 1
			Rx(k) = Dx(i);
			Ry(k) = Dy(i);
			Rv(k) = Dv(i);
			Rw(k) = 0.5 * std::max(0.0, Dt(i) - ts[k-1]); // ignore the part after Dt(i), as we have no data there	
		} else {
			Rx(k) = (1.0-alpha) * Dx(i-1) + alpha * Dx(i);
			Ry(k) = (1.0-alpha) * Dy(i-1) + alpha * Dy(i);
			Rv(k) = alpha * (1.0-alpha) * diff0T
					+ (1.0-alpha)*(1.0-alpha) * Dv(i-1)
					+ alpha*alpha * Dv(i);
		
			// Set the weights correctly	
			if (k == 0) {
				Rw(k) = 0.5 * (ts[1] - t);
			} else if (k == length(timeSteps) - 1) {
				Rw(k) = 0.5 * (t - ts[k-1]);
			} else {
		 		Rw(k) = 0.5 * (ts[k+1] - ts[k-1]);
			}
		}

	}
	
	UNPROTECT(2); // dInt, result
	return result;
} // UDTimesteps


/**
 * Computes the spatial distribution of a trajectory
 */
void utilizationDistribution(double *result, int *resultSize,
		int *nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		double *tsX, double *tsY, double *tsVar, double *tsW,
		double *xc, double *yc, int *nrow, int *ncol) {

	double halfCellSize = 0.5 * (xc[1] - xc[0]);
	for (off_t x = 0; x < *ncol; x++) {
		for (off_t y = 0; y < *nrow; y++) {
			Vec2D<double> pos(xc[x], yc[y]);
			double *value = result + (x * *nrow + y);
			*value = 0.0;
			for (off_t k = 0; k < *nsteps; k++) {
				*value += tsW[k] * cdf_norm2d(pos, halfCellSize,
						Vec2D<double>(tsX[k], tsY[k]), tsVar[k]);
			}
			
			// TODO: is there a simple way to normalize the result?
		}
	}
}

} // extern "C"

