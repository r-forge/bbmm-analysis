/**
 * Evaluate the PDF of a bivariate circular normal distribution
 */
double4 pdf_norm2d(const double4 dist2, const double4 var) {
	return exp((double4)(-0.5) * dist2 / var)
				/ (var * (double4)(2 * M_PI));
}

double4 cdf_rice(const double4 x, const double4 nu, const double4 sd) {
	return (double4)(1.0 - marcumQ(nu.s0/sd.s0, x.s0/sd.s0),
	                 1.0 - marcumQ(nu.s1/sd.s1, x.s1/sd.s1),
	                 1.0 - marcumQ(nu.s2/sd.s2, x.s2/sd.s2),
	                 1.0 - marcumQ(nu.s3/sd.s3, x.s3/sd.s3));
}
