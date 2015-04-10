/**
 * Calculates the utilization distribution given a specific trajectory.
 * tsX and tsY contain the mean location at every time step. tsVar stores the
 * variance at these times and tsW allows to assign a weight to each time step.
 *
 * The UD is calculated on a grid for which the lines are given by xc and yc.
 *
 * Vectorized implementation. Should operate faster especially on CPUs,
 * where it takes advantage of SSE/AVX instructions
 *
 * On GPUs there is not as much advantage, but it effectively unrolls the for
 * loop a few times, which may benefit the performance too.
 */
__kernel void utilizationDistribution(__global double *result, const int resultSize,
		const int nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		__global double4 *tsX,  __global double4 *tsY,  __global double4 *tsVar, __global double4 *tsW,
		__global double *xc, __global double *yc, const int nxc, const int nyc) {
	int idx = get_global_id(0);
		
	while (idx < resultSize) {
		const double4 posX = (double4)(xc[idx / nxc]);
		const double4 posY = (double4)(yc[idx % nxc]);
	
		double density = 0.0;
		for (int i = 0; i < (nsteps+3) / 4; i++) {
			const double4 dx = tsX[i] - posX;
			const double4 dy = tsY[i] - posY;
			
			// taking the dot product with the weights computes a weighted sum of the density
			density += dot(pdf_norm2d((dx*dx + dy*dy), tsVar[i]), tsW[i]);
		}
		result[idx] = density;
		idx += get_global_size(0);
	}
}

__kernel void encounterUD (__global double *result, const int resultSize,
		const double threshold, const int nsteps,
		// mean X coord, mean Y coord, variance and weight for each time step
		__global double4 *tsX, __global double4 *tsY, __global double4 *tsVar, __global double4 *tsW,
		__global double4 *ts2X, __global double4 *ts2Y, __global double4 *ts2SD,
		__global double *xc, __global double *yc, const int nrow, const int ncol) {	
	int idx = get_global_id(0);
		
	while (idx < resultSize) {
		const double4 posX = (double4)(xc[idx / nrow]);
		const double4 posY = (double4)(yc[idx % nrow]);
	
		double density = 0.0;
		for (int i = 0; i < (nsteps+3) / 4; i++) {
			const double4 dx = tsX[i] - posX;
			const double4 dy = tsY[i] - posY;
			
			const double4 meanDist = hypot(ts2X[i] - posX, ts2Y[i] - posY);
			
			density += dot(pdf_norm2d((dx*dx + dy*dy), tsVar[i]) 
					* cdf_rice((double4)threshold, meanDist, ts2SD[i]), tsW[i]);
		}
		result[idx] = density;
		idx += get_global_size(0);
	}
}

