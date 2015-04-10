//#include <R.h>
//#include <Rdefines.h>
//#include <cmath>
//#include "encounter.h"

//#define Z(x,y,d) (DEM[(d)+ DEMsize[0]*((x)+DEMsize[1]*(y))])
//#define sgn(x) ((x) < 0 ? -1 : (x)>0)
//#define cdf_norm(x, var) erf((x)/sqrt(2*(var))) // careful: not really the normal CDF

//inline double cellProb (double dx, double dy, double var, double halfCellSize) {
//	return 0.25
//			* (cdf_norm(dx+halfCellSize, var) - cdf_norm(dx-halfCellSize, var))
//			* (cdf_norm(dy+halfCellSize, var) - cdf_norm(dy-halfCellSize, var));
//}

//extern "C" {

//SEXP encounterProbability(SEXP visibility, SEXP p1, SEXP p2, SEXP s_xc, SEXP s_yc) {
//	int *DEMsize = INTEGER(GET_DIM(visibility));
//	double *xc = REAL(s_xc);
//	double *yc = REAL(s_yc);
//	double halfCellSize = 0.5 * (xc[1]-xc[0]);

//	SEXP result;
//	PROTECT(result = allocVector(REALSXP, 1));
//	
//	int i,j,t;
//	double prob = 0.0;
//	for (i=0; i < DEMsize[0]; i++) {
//		double dx = xc[i] - REAL(p1)[0];
//		for (j=0; j < DEMsize[1]; j++) {
//			double dy = yc[j] - REAL(p1)[1];
//			
//			//Rprintf("(%f, %f)\n", xc[i], yc[j]);
//			
//			SEXP s_visData = VECTOR_ELT(visibility, i + DEMsize[0]*j);
//			double *visData = REAL(s_visData);
//		
//			double p = 0.0;
//			for (t=0; t<LENGTH(s_visData); t+=3) {
//				//Rprintf("\t(%f, %f): %f\n", visData[t], visData[t+1], visData[t+2]);
//				double dxt = visData[t]   - REAL(p2)[0];
//				double dyt = visData[t+1] - REAL(p2)[1];
//				p += visData[t+2] * cellProb(dxt, dyt, REAL(p2)[2], halfCellSize);
//			}
//			prob += p * cellProb(dx, dy, REAL(p1)[2], halfCellSize);
//		}
//	}

//	REAL(result)[0] = prob; // integrating twice over an area element with area cellSize^2
//	UNPROTECT(1);
//	return result;
//}

//SEXP visibility(SEXP s_DEM)
//{
//	int *DEMsize = INTEGER(GET_DIM(s_DEM));
//	double *DEM = REAL(s_DEM);
//	SEXP result, s_visible;
//	
//	PROTECT(result = allocVector(VECSXP, DEMsize[1]*DEMsize[2]));
//	
//	int i,j,d; // indexing row, column and DEM instance
//	int ti, tj; // coordinates of viewing target
//	int *visible = (int*)malloc(sizeof(int) * DEMsize[1] * DEMsize[2]);
//	for (i = 0; i < DEMsize[1]; i++) {
//		for (j = 0; j < DEMsize[2]; j++) {
//			for (int t = 0; t < DEMsize[1]*DEMsize[2]; t++) { visible[t] = 0; }
//			int nVisible = 0;
//			
//			for (d = 0; d < DEMsize[0]; d++) {
//				double hPos = Z(i,j,d);
//			
//				for (ti = 0; ti < DEMsize[1]; ti++) {
//					int di = ti-i;
//					for (tj = 0; tj < DEMsize[2]; tj++) {
//						int dj = tj-j;
//						double dh = Z(ti,tj,d) - hPos;
//						
//						bool vis = true;
//						// Walk along the x axis and investigate crossings with vertical edges
//						for (int x = i + sgn(di); x != ti && vis; x += sgn(di)) {
//							double hLoS = hPos + (x-i)/di * dh;
//				
//							double y = j + (x-i)/(double)di * dj;
//							int yi = (int)y;
//							double hDEM = Z(x,yi,d)
//									+ (y-yi) * (Z(x,yi+1,d) - Z(x,yi,d));
//				
//							if (hLoS < hDEM) { vis = false; }
//						}
//						
//						// Walk along the y axis and investigate crossings with horizontal edges
//						for (int y = j + sgn(dj); y != tj && vis; y += sgn(dj)) {
//							double hLoS = hPos + (y-j)/dj * dh;
//				
//							double x = i + (y-j)/(double)dj * di;
//							int xi = (int)x;
//							double hDEM = Z(xi,y,d)
//									+ (x-xi) * (Z(xi+1,y,d) - Z(xi,y,d));
//				
//							if (hLoS < hDEM) { vis = false; }
//						}
//						
//						if (vis && visible[tj + DEMsize[2]*ti] == 0) { nVisible++; }
//						visible[tj + DEMsize[2]*ti] += vis;
//					}
//				}
//			}
//			
//			// collect all points with nonzero probability of being visible
//			PROTECT(s_visible = allocVector(INTSXP, nVisible * 3));
//			nVisible = 0;
//			for (ti = 0; ti < DEMsize[1]; ti++) {
//				for (tj = 0; tj < DEMsize[2]; tj++) {
//					if (visible[tj + DEMsize[2]*ti]) {
//						INTEGER(s_visible)[nVisible++] = ti;
//						INTEGER(s_visible)[nVisible++] = tj;
//						INTEGER(s_visible)[nVisible++] = visible[tj + DEMsize[2]*ti];
//					}
//				}
//			}
//			SET_VECTOR_ELT(result, i + DEMsize[1]*j, s_visible);
//		}
//	}
//	free(visible);
//	
//	UNPROTECT(1 + DEMsize[1]*DEMsize[2]);
//	return result;
//}

//} // extern "C"

