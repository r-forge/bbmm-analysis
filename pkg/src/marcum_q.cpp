#include "encounter.h"
#include "bessel.h"

/**
 * Compute Q_1 (alpha, beta)
 *
 * Uses the Bessel series approach:
 *
 *     Q_M (a, b) = e^(-(a*a + b*b)/2) + sum_{k=1-M}^{inf} (a/b)^k I_k(ab)
 *
 * When beta>=alpha, use a modified form of the series:
 *
 * 1 - Q_M (a, b) = e^(-(a*a + b*b)/2) + sum_{k=M}^{inf}   (b/a)^k I_k(ab)
 *
 * This performs best when alpha < beta, according to:
 *     Computation of Rice and Noncentral Chi-Squared Probabilities
 *     Robert T. Short
 *
 * To obtain I_k(alpha*beta) at each term of the sum, use the identity
 *   I_{n+1} (x) = - (2n/x) I_n (x) + I_{n-1} (x)
 * and evaluate the sum backwards for numerical stability.
 */
double marcumQ(const double alpha, const double beta) {
	// TODO: the summation seems to fail for some very large beta

	if (beta == 0) return 1.0;

	const double ab = alpha * beta;
	const double two_over_ab = 2.0/ab;

	// Decide what form of the sum we are going to evaluate
	int startIter = (alpha >= beta); // start at 0 if alpha < beta, 1 otherwise
	// the scaling factor for I_k in each iteration, either alpha/beta or its reciprocal
	double scaleFactor = pow(alpha/beta, startIter*-2 + 1);


	/* starting values */
	int n = 60; // TODO: determine the best possible value of n based on the function arguments
	double Inp1 = bessel_In_scaled(n+1, ab);
	double In   = bessel_In_scaled(n, ab);
	double Inm1;

	double result = 0.0;
	for(; n>=startIter; n--) {
		result *= scaleFactor;
		result += In;
		Inm1 = Inp1 + n * two_over_ab * In;
		Inp1 = In;
		In   = Inm1;
	}

	// If startIter != 0, the sum is still off by a factor (alpha/beta)^-startIter.
	result *= pow(scaleFactor, startIter);
	// Also apply the factor outside the sum. Since we used an exponentially scaled
	// Bessel routine returning e^(-ab)I_k(ab), we have a different scaling constant here.
	result *= exp(-0.5 * (alpha-beta) * (alpha-beta));

	// depending on the approach we used, we now have either Q(a,b) or 1-Q(a,b)
	return (1-startIter) * result  +  startIter * (1-result);
}

