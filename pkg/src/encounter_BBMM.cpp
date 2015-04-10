#include "vec2d.h"
#include <math.h>
#include <algorithm>

#include "encounter.h"
#include "util.h"

extern "C" {

/**
 * Compute the expected duration of encounters between two bursts using the BBMM.
 */
void encounterBBMM(double *result, double *threshold, 
		int *nsteps, BBMM_timestep<double> *data) {
	size_t n = (size_t)*nsteps;
	
	// find the intersection of the time intervals described by the trajectories
	// iterate over *nAlpha time intervals
	for (size_t i = 0; i < n; i++) {
		*result += data[i].weight
				* cdf_rice(*threshold, data[i].mu.norm(), sqrt(data[i].var));
	}
}

/**
 * Computes the spatial distribution of expected encounter durations
 */
//void encounterUD(double *result, double *threshold, 
//		int *nloc1, BBMM_measurement<double> *data1,
//		int *nloc2, BBMM_measurement<double> *data2,
//		int *nAlpha,
//		double *xc, double *yc, int *nrow, int *ncol) {
		
void encounterUD(double *result, int *resultSize,
		double *threshold, int *nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		double *tsX,  double *tsY,  double *tsVar, double *tsW,
		double *ts2X, double *ts2Y, double *ts2SD, // and the distribution of the other trajectory
		double *xc, double *yc, int *nrow, int *ncol) {
	
	for (off_t x = 0; x < *ncol; x++) {
		//Rprintf("\n%d ", x);
		for (off_t y = 0; y < *nrow; y++) {
			//Rprintf("(%d,%d)\n", x, y);
			Vec2D<double> pos(xc[x], yc[y]);
			double *value = result + (x * *nrow + y);
			*value = 0.0;
			for (off_t k = 0; k < *nsteps; k++) {
				*value += tsW[k] * pdf_norm2d(pos, Vec2D<double>(tsX[k], tsY[k]), tsVar[k])
						* cdf_rice(*threshold, (Vec2D<double>(ts2X[k], ts2Y[k])-pos).norm(), ts2SD[k]);
				
				//Rprintf("(%d, %d) %d/%d %10f %10f %10f %10f %10f %10f\n", x, y, k, *nsteps, *value, tsX[k], tsY[k], tsVar[k], tsW[k]);
			}
			
			// TODO: is there a simple way to normalize the result?
		}
		//Rprintf("\n\n");
	}
}


void riceTest(double *result, int *resultSize,
		double *nu, double *sd, double *x) {
	for (off_t i = 0; i < *resultSize; i++) {
		result[i] = cdf_rice(x[i], *nu, *sd);
	}
}
} // extern "C"

