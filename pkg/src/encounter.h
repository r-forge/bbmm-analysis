#ifndef __ENCOUNTER_H__
#define __ENCOUNTER_H__

#include <R.h>
#include <Rinternals.h>
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

template<class T>
struct LevyMM_timestep {
	Vec2D<T> c3; // coefficients of a polynomial of degree 4,
	T        c2; // some of which are vectors.
	Vec2D<T> c1; // The coefficient for x^4 is always 1.
	T        c0;
	T        scaleFactor; // The whole polynomial is scaled by this factor.
	T weight; // the weight of this time step in the total UD
};

extern "C" {

/**
 * Compute the mean location and variance at every timestep for a call to getUD().
 */
SEXP UDTimesteps(SEXP bData, SEXP diff, SEXP timeSteps);

/**
 * Compute a utilization distribution
 */
void utilizationDistribution(double *result, int *resultSize,
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
  
#endif // __ENCOUNTER_H__

