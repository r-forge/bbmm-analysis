#include "vec2d.h"
#include <math.h>
#include <algorithm>
#include <limits>

#include "encounter.h"

// returns the number of (distinct) roots found
template<class T>
size_t quadraticRealRoots(T a, T b, T c, T &x1, T &x2) {
	double delta=0.0;
	const T eps = std::numeric_limits<T>::epsilon();
	// equation is: a*x^2 + b*x + c = 0; a cannot be zero.
	if ( fabs(a) <= eps ) { return false; }

	delta=b*b-4.0*a*c;
	if ( delta<0.0 ) { return 0; }
	if ( delta>=0.0 && delta <= eps )
	{ x1=(-b/(2.0*a)); x2=x1; return 1; }

	x1 = (-b - sqrt(delta)) / (2.0*a);
	x2 = (-b + sqrt(delta)) / (2.0*a);
	return 2;
}

// Returns the duration of encounters of two linearly moving points for t in [0,1]
template<class T>
T encounterFractionSegment(const T threshold,
		const Vec2D<T> &start1, const Vec2D<T> &end1,
		const Vec2D<T> &start2, const Vec2D<T> &end2,
		T *encStart = NULL, T *encEnd = NULL) {

	Vec2D<T> d0 = start2-start1; // the initial distance between the objects
	Vec2D<T> dr = (end2-end1) - d0; // the relative displacement of the objects

	// solve a quadratic equation of the form:
	// |dr|^2 t^2 + 2(d0.dr) t + |d0|^2 = threshold^2
	// note that this is a convex parabola and thus the encounter must occur between the roots, if any
	T root1, root2;
	size_t roots = quadraticRealRoots<T>(
			dr.norm2(),
			2 * d0.dot(dr),
			d0.norm2() - threshold*threshold,
			root1, root2
	);

	// clip the roots to the allowed range
	root1 = std::min(1.0, std::max(0.0, root1));
	root2 = std::min(1.0, std::max(0.0, root2));

	if (roots <= 1 || root2 <= root1) {
		// no encounter occurs
		if (encStart != NULL) { *encStart = std::numeric_limits<double>::quiet_NaN(); }
		if (encEnd   != NULL) { *encEnd   = std::numeric_limits<double>::quiet_NaN(); }
		
		return 0.0;
	} else {
		// pass the extent of the encounter intervals back if requested
		if (encStart != NULL) { *encStart = root1; }
		if (encEnd   != NULL) { *encEnd   = root2; }
	
		// encounter occurs in [root1,root2]
		return root2 - root1;
	}
}

extern "C" {

void encounterLinear(double *result, double *threshold,
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2) {
	encounterLinearIntervals(result, threshold, nloc1, data1, nloc2, data2, NULL);
}

void encounterLinearIntervals(double *result, double *threshold, 
		int *nloc1, BBMM_measurement<double> *data1,
		int *nloc2, BBMM_measurement<double> *data2,
		double *encounterIntervals) {
	size_t n1 = (size_t)*nloc1;
	size_t n2 = (size_t)*nloc2;

	*result = 0.0;
	
	size_t i = 1; // these point to the end times for the interval currently being processed
	size_t j = 1;

	// first make sure both intervals overlap
	while (i < n1 && data1[i].t < data2[0].t) {
		i++;
	}
	while (j < n2 && data2[j].t < data1[0].t) {
		j++;
	}

	// calculate the starting positions of the first interval
	double startTime = std::max(data1[0].t, data2[0].t);
	
	double alpha1 = (startTime - data1[i-1].t) / (data1[i].t - data1[i-1].t);
	Vec2D<double> start1 = (1-alpha1) * data1[i-1].mu + alpha1 * data1[i].mu;
	
	double alpha2 = (startTime - data2[j-1].t) / (data2[j].t - data2[j-1].t);
	Vec2D<double> start2 = (1-alpha2) * data2[j-1].mu + alpha2 * data2[j].mu;

	double currEncIntv[2];
	if (encounterIntervals != NULL) {
		for (size_t k = 0; k < 2*(n1 + n2); k++) {
			encounterIntervals[k] =  std::numeric_limits<double>::quiet_NaN();
		}
	}
	
	
	for (size_t k = 0; i < n1 && j < n2; k++) {
		double endTime = std::min(data1[i].t, data2[j].t);
		
		alpha1 = (endTime - data1[i-1].t) / (data1[i].t - data1[i-1].t);
		Vec2D<double> end1 = (1-alpha1) * data1[i-1].mu + alpha1 * data1[i].mu;
		
		alpha2 = (endTime - data2[j-1].t) / (data2[j].t - data2[j-1].t);
		Vec2D<double> end2 = (1-alpha2) * data2[j-1].mu + alpha2 * data2[j].mu;
		
		*result += (endTime - startTime) *
				encounterFractionSegment<double>(*threshold, 
					start1, end1, start2, end2, currEncIntv, currEncIntv+1);
		if (encounterIntervals != NULL) {
			encounterIntervals[2*k]     = startTime + currEncIntv[0] * (endTime-startTime);
			encounterIntervals[2*k + 1] = startTime + currEncIntv[1] * (endTime-startTime);
		}
				
		if (data1[i].t < data2[j].t) {
			i++;
		} else {
			j++;
		}
		
		startTime = endTime;
		start1 = end1;
		start2 = end2;
	}
}

} // extern "C"
