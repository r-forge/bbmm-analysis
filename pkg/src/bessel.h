
double bessel_I0_scaled(double x);
double bessel_I1_scaled(double x);
double bessel_In_scaled(int n, double x);

double bessel_IJ_taylor(const int n, const double x,
                        const int sign,
                        const int kmax,
                        const double threshold);

/* Evaluate the continued fraction CF1 for I_{nu+1}/I_nu
 * using Gautschi (Euler) equivalent series.
 */
double bessel_I_CF1_ser(const int n, const double x);
