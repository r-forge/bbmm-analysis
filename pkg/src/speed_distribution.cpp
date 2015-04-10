#include "vec2d.h"
#include <math.h>
#include <algorithm>
#include <R.h>

#include "encounter.h"
#include "bessel.h"

#define BB_SCALE 18.42068

extern "C" {

void speedDistribution(double *wSum, double *weight,
		BBMM_measurement<double> *data, int *nloc,
		double *xc, double *yc, int *nrow, int *ncol,
		double *timeStep, double *dts, double *dtf) {

	ptrdiff_t iloc = -1; // the index of the measurement directly before t
	BBMM_measurement<double> di  = data[0];
	BBMM_measurement<double> di1 = data[0];
	
	ptrdiff_t xmin = 0;
	ptrdiff_t ymin = 0;
	for (double t = data[0].t - *dts; t <= data[*nloc-1].t - *dtf; t += *timeStep) {
		while (iloc < *nloc-1 && data[iloc+1].t < t) {
			iloc++;
			di  = di1;
			di1 = data[iloc+1];
		}
		
		double alpha = (t-di.t) / (di1.t-di.t);
		Vec2D<double> pMu = (1-alpha) * di.mu + alpha * di1.mu;
		double pVar = (1-alpha)*(1-alpha)*di.vMeas + alpha*alpha*di1.vMeas
				+ (di1.t-di.t)*alpha*(1-alpha)*di.vDiff;
		
		// find the indices in the coordinate lists where we start computing
		double dMax = sqrt(BB_SCALE*pVar);
		while (xmin < *nrow-1 && xc[xmin]    < pMu.x - dMax) { xmin++; }
		while (xmin > 0       && xc[xmin-1] >= pMu.x - dMax) { xmin--; }
		while (ymin < *ncol-1 && yc[ymin]    < pMu.y - dMax) { ymin++; }
		while (ymin > 0       && yc[ymin-1] >= pMu.y - dMax) { ymin--; }
		for (ptrdiff_t x = xmin; x < *nrow && xc[x] <= pMu.x + dMax; x++) {
			for (ptrdiff_t y = ymin; y < *ncol && yc[y] <= pMu.y + dMax; y++) {
				Vec2D<double> pos(xc[x], yc[y]);
				double d2 = (pos-pMu).norm2();
				
				if (d2 <= BB_SCALE*pVar) {
					double w = pdf_norm2d(d2, pVar);
					
					double ts = t + *dts;
					double tf = t + *dtf;
					
					Vec2D<double> sMu, fMu;
					double sVar, fVar;
					if (ts == t) {
						sMu  = pos;
						sVar = 0.0;
					} else if (ts < di.t) {
						// the start of the time interval is in the previous bridge
						BBMM_measurement<double> dim1  = data[iloc-1];
						double sAlpha = (di1.t - t) / (di1.t-di.t);
						double sBeta  = (di.t - ts) / (di.t-dim1.t);

						if (sAlpha <= 0.0 || sBeta >= 1.0) {
							BBMM_measurement<double> dim2  = data[iloc-2];
							BBMM_measurement<double> dim1  = data[iloc-1];
						
							sAlpha = (ts - dim2.t) / (dim1.t - dim2.t);
						
							sMu  = (1-sAlpha) * dim2.mu + sAlpha * dim1.mu;
							sVar = (1-sAlpha) * (1-sAlpha) * dim2.vMeas
									+ sAlpha * sAlpha * dim1.vMeas
									+ (dim1.t - dim2.t) * (1-sAlpha) * sAlpha * dim2.vDiff;
						} else {
							double v = ((1-sAlpha)/sAlpha)*((1-sAlpha)/sAlpha) * di1.vMeas
									+ (di1.t-di.t) * (1-sAlpha) / sAlpha * di.vDiff;
						
							sMu = sBeta * dim1.mu + (1-sBeta) * (
									  v * di.mu
									+ di.vMeas * (pos - (1-sAlpha)*di1.mu) / sAlpha
							) / (v + di.vMeas);
							sVar = sBeta*sBeta * dim1.vMeas
									+ (di.t-dim1.t) * sBeta * (1-sBeta) * dim1.vDiff
									+ (1-sBeta)*(1-sBeta) * v * di.vMeas / (v + di.vMeas);
						}
					} else {
						// The start of the time interval is in the same bridge as t
						double sAlpha = (ts - di.t) / (t - di.t);
						
						sMu = (1-sAlpha) * di.mu + sAlpha * pos;
						sVar = (t-di.t) * sAlpha * (1-sAlpha) * di.vDiff
								+ (1-sAlpha)*(1-sAlpha) * di.vMeas;
					}
					
					if (tf == t) {
						fMu  = pos;
						fVar = 0.0;
					} else if ((iloc < *nloc-2 && tf > data[iloc+2].t)) {
						
					} else if (tf > di1.t) {
						// the end of the time interval is in the next bridge
						BBMM_measurement<double> di2  = data[iloc+2];
						double fAlpha = (t - di.t) / (di1.t-di.t);
						double fBeta  = (tf - di1.t) / (di2.t-di1.t);
						
						if (fAlpha <= 0.0 || fBeta >= 1.0) {
							BBMM_measurement<double> di2  = data[iloc+2];
							BBMM_measurement<double> di3  = data[iloc+3];
						
							fAlpha = (tf - di2.t) / (di3.t - di2.t);
						
							fMu  = (1-fAlpha) * di2.mu + fAlpha * di3.mu;
							fVar = (1-fAlpha) * (1-fAlpha) * di2.vMeas
									+ fAlpha * fAlpha * di3.vMeas
									+ (di3.t - di2.t) * (1-fAlpha) * fAlpha * di2.vDiff;
						} else {
							double v = ((1-fAlpha)/fAlpha)*((1-fAlpha)/fAlpha) * di.vMeas
									+ (di1.t-di.t) * (1-fAlpha) / fAlpha * di.vDiff;
						
							fMu = fBeta * di2.mu + (1-fBeta) * (
									  v * di1.mu
									+ di1.vMeas * (pos - (1-fAlpha)*di.mu) / fAlpha
							) / (v + di1.vMeas);
							fVar = fBeta*fBeta * di2.vMeas
									+ (di2.t-di1.t) * fBeta * (1-fBeta) * di1.vDiff
									+ (1-fBeta)*(1-fBeta) * v * di1.vMeas / (v + di1.vMeas);
						}
					} else {
						// The end of the time interval is in the same bridge as t
						double fAlpha = *dtf / (di1.t - t);
						
						fMu = (1-fAlpha) * pos + fAlpha * di1.mu;
						fVar = (di1.t-t) * fAlpha * (1-fAlpha) * di.vDiff
								+ (fAlpha)*(fAlpha) * di1.vMeas;
					}
					
					double nu = (fMu-sMu).norm() / (*dtf-*dts);
					double sigma = sqrt(sVar + fVar) / (*dtf-*dts);
					
					double lx = -nu*nu / (2*sigma*sigma);
					double avgSpeed = sigma * sqrt(M_PI / 2)
							* ((1-lx) * bessel_I0_scaled(-lx/2) - lx * bessel_I1_scaled(-lx/2));
							
					wSum[x + *nrow*y]   += w * avgSpeed;
					weight[x + *nrow*y] += w;
				}
			}
		}
	}
}
	
} // extern "C"

