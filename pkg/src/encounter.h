#ifndef __ENCOUNTER_H__
#define __ENCOUNTER_H__

#include <R.h>
#include "vec2d.h"

template<class T>
struct BBMM_measurement {
	Vec2D<T> mu; // the mean location
	T vDiff; // the diffusion coefficient for the bridge preceding this measurement // TODO: check this
	T vMeas; // the variance of the measured location
	T t; // the time of the measurement
};

template<class T>
struct BBMM_timestep {
	Vec2D<T> mu; // the mean location
	T var; // The variance at this time step
	T weight; // the weight of this time step in the total UD
};

extern "C" {

void UDTimesteps(BBMM_timestep<double> *result,
		int *nloc, BBMM_measurement<double> *data,
		int *nsteps, double *timeLimits);

/**
 * Compute a utilization distribution
 */
void getUD(double *result, int *resultSize,
		int *nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		double *tsX, double *tsY, double *tsVar, double *tsW,
		double *xc, double *yc, int *nrow, int *ncol);

void encounterLinear(double *result, double *threshold,
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2,
		double *timestepSize);
void encounterLinearIntervals(double *result, double *threshold,
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2,
		double *timestepSize, double *encounterIntervals);

/**
 * Compute the expected duration of encounters between two bursts using the BBMM.
 */
void encounterBBMM(double *result, double *threshold, 
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2,
		double *timestepSize);

/**
 * Computes the spatial distribution of expected encounter durations.
 *
 * Computes the expected duration of encounters occuring between two bursts
 * of location measurements as a distribution over a spatial grid defined by
 * the list of coordinates xc and yc.
 * In particular, for each x from xc, y from yc this function evaluates the
 * expected duration that the first trajectory is at (x,y) while simultaneously
 * the second trajectory is at most a distance *threshold away.
 */
void encounterUD(double *result, int *resultSize,
		double *threshold, int *nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		double *tsX,  double *tsY,  double *tsVar, double *tsW,
		double *ts2X, double *ts2Y, double *ts2SD, // and the distribution of the other trajectory
		double *xc, double *yc, int *nrow, int *ncol);

void diffusion(double *xys, double *Tr, 
	 int *nloc, double *Lr, double *sigma, int *nsig);

}; // extern "C"

inline double pdf_norm2d(Vec2D<double> p, Vec2D<double> mu, double var)
{
	// straightforward implementation of the expression for the pdf of N(mu, sigma^2)
    return exp((p-mu).norm2() / (-2.0 * var)) / (2 * PI * var);
}

double marcumQ(const double alpha, const double beta);
inline double cdf_rice(double x, double nu, double sigma) {
	return 1.0 - marcumQ(nu/sigma, x/sigma);
}
  
#endif // __ENCOUNTER_H__

