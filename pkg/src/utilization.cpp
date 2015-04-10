#include "vec2d.h"
#include <math.h>
#include <algorithm>

#include "encounter.h"

extern "C" {

void UDTimesteps(BBMM_timestep<double> *result,
		int *nloc, BBMM_measurement<double> *data,
		int *nsteps, double *timeLimits) {
	int n = *nloc;
	
	double startTime = timeLimits[0];
	double endTime   = timeLimits[1];

	off_t i = 1;
	
	double timeStep = (endTime-startTime) / (*nsteps-1);
	for (off_t k = 0; k < *nsteps; k++) {
		double t = startTime + k * timeStep;
		
		double alpha = (t-data[i-1].t) / (data[i].t - data[i-1].t);
		while (alpha > 1.0 && i < (*nloc-1)) {
			i++;
			alpha = (t-data[i-1].t) / (data[i].t - data[i-1].t);
		}
		
		result[k].mu = ((1.0-alpha) * data[i-1].mu + alpha * data[i].mu);
		result[k].var = (data[i].t-data[i-1].t) * alpha * (1.0-alpha) * data[i-1].vDiff
				+ (1.0-alpha)*(1.0-alpha) * data[i-1].vMeas + alpha*alpha * data[i].vMeas;
		if (k == 0) {
			if (t < data[0].t) {
				result[0].mu  = data[0].mu;
				result[0].var = data[0].vMeas;
			}
			result[0].weight = timeStep / 2.0 + t - data[0].t;
		} else if (k == *nsteps - 1) {
			if (t > data[n-1].t) {
				result[*nsteps - 1].mu  = data[n-1].mu;
				result[*nsteps - 1].var = data[n-1].vMeas;
			}
			result[*nsteps - 1].weight = timeStep / 2.0 + data[n-1].t - t;
		} else {
	 		result[k].weight = timeStep;
		}
	}
}

/**
 * Computes the spatial distribution of a trajectory
 */
void utilizationDistribution(double *result, int *resultSize,
		int *nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		double *tsX, double *tsY, double *tsVar, double *tsW,
		double *xc, double *yc, int *nrow, int *ncol) {

	for (off_t x = 0; x < *ncol; x++) {
		for (off_t y = 0; y < *nrow; y++) {
			//Rprintf("(%d,%d)\n", x, y);
			Vec2D<double> pos(xc[x], yc[y]);
			double *value = result + (x * *nrow + y);
			*value = 0.0;
			for (off_t k = 0; k < *nsteps; k++) {
				*value += tsW[k] * pdf_norm2d(pos, Vec2D<double>(tsX[k], tsY[k]), tsVar[k]);
			}
			
			// TODO: is there a simple way to normalize the result?
		}
	}
}

} // extern "C"

