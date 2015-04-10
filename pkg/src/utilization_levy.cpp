//#include "vec2d.h"
//#include <math.h>
//#include <limits>
//#include <algorithm>

//#include "encounter.h"

//extern "C" {

//void UDLevySteps(LevyMM_timestep<double> *result,
//		int *nloc, BBMM_measurement<double> *data,
//		int *nsteps, double *timeLimits) {
//	int n = *nloc;
//	
//	double startTime = timeLimits[0];
//	double endTime   = timeLimits[1];

//	off_t i = 1;
//	
//	double timeStep = (endTime-startTime) / (*nsteps-1);
//	for (off_t k = 0; k < *nsteps; k++) {
//		double t = startTime + k * timeStep;
//		
//		double alpha = (t-data[i-1].t) / (data[i].t - data[i-1].t);
//		while (alpha > 1.0 && i < (*nloc-1)) {
//			i++;
//			alpha = (t-data[i-1].t) / (data[i].t - data[i-1].t);
//		}
//		
//		double T = data[i].t-data[i-1].t;
//		result[k].scaleFactor = pow(2 * PI, 2/3) / (
//		  ((data[i].mu-data[i-1].mu).norm2() + pow(T*data[i-1].vDiff, 2))
//		* T * alpha * (1-alpha));
//		
//		// Store the coefficients of the polynomial describing the distribution.
//		// If distr. is degenerate, store the params of measured loc. distr.
//		if (alpha <= 0) {
//			result[k].c3 = data[i-1].mu;
//			result[k].c2 = data[i-1].vMeas;
//			result[k].scaleFactor = std::numeric_limits<double>::infinity();
//		} else if (alpha >= 1.0) {
//			result[k].c3 = data[i-1].mu;
//			result[k].c2 = data[i-1].vMeas;
//			result[k].scaleFactor = std::numeric_limits<double>::infinity();
//		} else {
//			double a2t = data[i-1].mu.norm2() + pow((t-data[i-1].t)*data[i-1].vDiff, 2);
//			double b2t = data[i].mu.norm2()   + pow((data[i].t-t)  *data[i-1].vDiff, 2);
//		
//			result[k].c3 = -2.0 * (data[i-1].mu + data[i].mu);
//			result[k].c2 = (a2t + b2t + 4 * data[i-1].mu.dot(data[i].mu));
//			result[k].c1 = -2.0 * (b2t * data[i-1].mu + a2t * data[i].mu);
//			result[k].c0 = a2t * b2t;
//		}
//		result[k].weight = timeStep;
//	}
//}

///**
// * Computes the spatial distribution of a trajectory
// */
//void UDLevy(double *result, int *resultSize,
//		int *nsteps,
//		LevyMM_timestep<double> *ts,
//		double *xc, double *yc, int *nrow, int *ncol) {

//	for (off_t x = 0; x < *ncol; x++) {
//		for (off_t y = 0; y < *nrow; y++) {
//			//Rprintf("(%d,%d)\n", x, y);
//			Vec2D<double> pos(xc[x], yc[y]);
//			double p2 = pos.norm2();
//			
//			double *value = result + (x * *nrow + y);
//			*value = 0.0;
//			for (off_t k = 0; k < *nsteps; k++) {
//				double v = ts[k].scaleFactor * (
//					  p2 * p2
//					+ p2 * ts[k].c3.dot(pos)
//					+ p2 * ts[k].c2
//					+ ts[k].c1.dot(pos)
//					+ ts[k].c0
//				);
//				v = pow(v, -1.5);

//				*value += ts[k].weight * v;
//			}
//			//Rprintf("(%f, %f) %E %E\n", pos.x, pos.y, p2, *value);
//		}
//	}
//}

//} // extern "C"

