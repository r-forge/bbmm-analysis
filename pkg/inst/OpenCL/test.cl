
//__kernel void besselITest(__global double *result, const int resultSize, const int n, __global double *x) {
//	int idx = get_global_id(0);

//	if (idx < resultSize) {
//		result[idx] = bessel_In_scaled(n, x[idx]);
//	}
//}

//__kernel void riceTest(__global double *result, const int resultSize,
//		const double nu, const double sd, __global double *x) {
//	int idx = get_global_id(0);

//	if (idx < resultSize) {
//		result[idx] = cdf_rice(x[idx], nu, sd);
//	}
//}
