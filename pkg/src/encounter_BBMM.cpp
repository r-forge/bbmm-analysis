#include "vec2d.h"
#include <math.h>
#include <algorithm>

#include "encounter.h"

extern "C" {

/**
 * Compute the expected duration of encounters between two bursts using the BBMM.
 */
void encounterBBMM(double *result, double *threshold, 
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2,
		double *timestepSize) {
	size_t n1 = (size_t)*nloc1;
	size_t n2 = (size_t)*nloc2;
	
	// find the intersection of the time intervals described by the trajectories
	double startTime = std::max(data1[0].t,    data2[0].t);
	double endTime   = std::min(data1[n1-1].t, data2[n2-1].t);

	size_t i = 1;
	size_t j = 1;
	
	*result = 0.0;

	// iterate over *nAlpha time intervals
	for (double t = startTime; t <= endTime; t += *timestepSize) {
		double alpha1 = (t-data1[i-1].t) / (data1[i].t - data1[i-1].t);
		while (alpha1 > 1.0 && i < n1) {
			i++;
			alpha1 = (t-data1[i-1].t) / (data1[i].t - data1[i-1].t);
		}
		
		double alpha2 = (t-data2[j-1].t) / (data2[j].t - data2[j-1].t);
		while (alpha2 > 1.0 && j < n2) {
			j++;
			alpha2 = (t-data2[j-1].t) / (data2[j].t - data2[j-1].t);
		}
		
		Vec2D<double> meanDistance = 
			((1.0-alpha1) * data1[i-1].mu + alpha1 * data1[i].mu)
			-
			((1.0-alpha2) * data2[j-1].mu + alpha2 * data2[j].mu);
			
		double varDistance =
			(data1[i].t - data1[i-1].t) * alpha1 * (1.0-alpha1) * data1[i].vDiff // TODO: check whether V11 at i or i-1
				 + (1.0-alpha1)*(1.0-alpha1) * data1[i-1].vMeas + alpha1*alpha1 * data1[i].vMeas
			+
			(data2[j].t - data2[j-1].t) * alpha2 * (1.0-alpha2) * data2[j].vDiff // TODO: check whether V11 at i or i-1
				 + (1.0-alpha2)*(1.0-alpha2) * data2[j-1].vMeas + alpha2*alpha2 * data2[j].vMeas;
		
		*result += cdf_rice(*threshold, meanDistance.norm(), sqrt(varDistance));
	}
	
	// up to now, we assumed each interval had unit length, scale the result
	*result *= *timestepSize;
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

