#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <limits>
#include <R.h>

#include "encounter.h"
#include "util.h"

/**********************************************************************
 * Maximisation of the likelihood for the Brownian bridge
 *********************************************************************/

template<class T>
struct BBMM_dcBridge {
	T d2; // squared distance between mean and observation
	T c1, c0; // coefficients for the linear function mapping diffusion coefficient to variance
};

extern "C" {
   
void diffusion_static(double *result,
		BBMM_dcBridge<double> *bridges, int *nbridges,
	 	double *sigma, int *nsig)
{
    size_t n = *nbridges;
    size_t ns = *nsig;

    double llMax=-std::numeric_limits<double>::infinity();
    
    // iterate over all possible values of sigma and compute the likelihood
    for (size_t r = 0; r < ns; r++) {
		double ll = 0.0;
		
		// calculate the probability for every other location measurement,
		// given a bridge between its neighbors
		for (size_t i=0; i < n; i++) {
			ll += log(pdf_norm2d(bridges[i].d2, bridges[i].c1 * sigma[r] + bridges[i].c0));
		}
		
		if (ll > llMax) {
			llMax = ll;
			*result = sigma[r];
		}
    }
}

} // extern "C"
