#include <math.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>
#include <R.h>

#include "encounter.h"

#define D_X(i)  (xys[3*(i)  ])
#define D_Y(i)  (xys[3*(i)+1])
#define D_S2(i) (xys[3*(i)+2])

/* *********************************************************************
   Maximisation of the likelihood for the Brownian bridge
   
   *********************************************************************/

extern "C" {

double norm2d(double x1, double y1, double moyx, double moyy,
	    double var)
{
	// straightforward implementation of the expression for the pdf of N(mu, sigma^2)
    double cste;
    cste = (1 / (2.0 * 3.141592653589793238 * var));
    return cste * exp( (-1.0 / (2.0 * var)) * (((x1 - moyx) * (x1 - moyx))+((y1 - moyy) * (y1 - moyy))));
}
   
void diffusion(double *xys, double *Tr, 
	 int *nloc, double *Lr, double *sigma, int *nsig)
{
    int i, nlo, ns, r;
    nlo = *nloc;
    ns = *nsig;
    
    double ai,muX, muY,sigmai;
    
    // iterate over all possible values of sigma and compute the likelihood
    for (r = 0; r < ns; r++) {
		Lr[r] = 0;
		
		// calculate the probability for every other location measurement,
		// given a bridge between its neighbors
		for (i=1; i+1 < nlo; i+=2) {
			ai = (Tr[i] - Tr[i-1])/(Tr[i+1] - Tr[i-1]);
	
			muX = (1.0-ai) * D_X(i-1) + ai * D_X(i+1);
			muY = (1.0-ai) * D_Y(i-1) + ai * D_Y(i+1);
	
			sigmai = ((Tr[i+1]-Tr[i-1]) * ai * (1-ai) * sigma[r]) + ((1.0-ai)*(1.0-ai) * D_S2(i-1)) + (ai*ai * D_S2(i+1));
	
			Lr[r] += log(norm2d(D_X(i), D_Y(i), muX, muY, sigmai));
		}
    }
}

} // extern "C"
