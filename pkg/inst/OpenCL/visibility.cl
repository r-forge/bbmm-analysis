
__kernel void DEMVisibility(__global double *result, const int resultSize,
		__global double *DEM, __global int2 *dPos, const int nPos) {
	int2 pos = (int2)(get_global_id(0), get_global_id(1));
	int2 size = (int2)(get_global_size(0), get_global_size(1));
	
	double hPos = DEM[pos.x * size.y + pos.y] + 0.5;
	
	for (int i = 0; i < nPos; i++) {
		int2 dp = dPos[i];
		int2 target = pos + dp;
		
		// be careful about target being outside the DEM; return NaN if it is
		double val = NAN;
		if (all(target >= 0 && target < size)) {
			double dh = DEM[target.x * size.y + target.y] + 0.5 - hPos;
			val = 1;
			
			
			// Walk along the x axis and investigate crossings with vertical edges
			if (abs(dp.x) >= 2) {
				for (int x = pos.x + sign(dp.x); x != target.x && val; x += sign(dp.x)) {
					double hLoS = hPos + (x-pos.x)/dp.x * dh;
					
					double y = pos.y + (x-pos.x)/dp.x * dp.y;
					int yi = (int)y;
					double hDEM = DEM[x * size.y + yi]
							+ (y-yi) * (DEM[x * size.y + yi + 1] - DEM[x * size.y + yi]);
					
					if (hLoS < hDEM) {
						val = 0;
					}
				}
			}
			
			// Walk along the y axis and investigate crossings with horizontal edges
			if (abs(dp.y) >= 2) {
				for (int y = pos.y + sign(dp.y); y != target.y && val; y += sign(dp.y)) {
					double hLoS = hPos + (y-pos.y)/dp.y * dh;
					
					double x = pos.x + (y-pos.y)/dp.y * dp.x;
					int xi = (int)x;
					double hDEM = DEM[xi * size.y + y] 
							+ (x-xi) * (DEM[(xi+1) * size.y + y] - DEM[xi * size.y + y]);
					
					if (hLoS < hDEM) {
						val = 0;
					}
				}
			}
		}
		
		result[i + nPos*(pos.x + size.x*pos.y)] = val;
	}
}

