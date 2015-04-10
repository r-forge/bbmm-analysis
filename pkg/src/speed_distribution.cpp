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
					// determine the position distribution at the start of the time interval
					if (ts == t) {
						sMu  = pos;
						sVar = 0.0;
					} else if (iloc >= 2 && ts < data[iloc-1].t) {
							// start of interval at least two bridges away;
							// location independent of fixed position at t
							BBMM_measurement<double> dS = data[iloc-2];
							BBMM_measurement<double> dF = data[iloc-1];
							ptrdiff_t i = iloc - 2;
							while (i > 0 && ts < dS.t) {
								dF = dS;
								dS = data[--i];
							}
						
							double sAlpha = (ts - dS.t) / (dF.t - dS.t);
						
							sMu  = (1-sAlpha) * dS.mu + sAlpha * dF.mu;
							sVar = (1-sAlpha) * (1-sAlpha) * dS.vMeas
									+ sAlpha * sAlpha * dF.vMeas
									+ (dF.t - dS.t) * (1-sAlpha) * sAlpha * dS.vDiff;
					} else {
						double sAlpha = (di1.t - t) / (di1.t-di.t);						
						// See section 3.1.1 of MSc thesis for derivation of these quantities
						Vec2D<double> muI = di.mu;
						double varI = di.vMeas;
						if (sAlpha > 0) {
							Vec2D<double> muQ = (sAlpha-1)/sAlpha * di1.mu + pos/sAlpha;
							double varQ = ((1-sAlpha)/sAlpha)*((1-sAlpha)/sAlpha) * di1.vMeas
										+ (di1.t-di.t) * (1-sAlpha) / sAlpha * di.vDiff;
						
							muI = (muQ * di.vMeas + di.mu * varQ)
									/ (varQ + di.vMeas);
							varI = varQ * di.vMeas / (varQ + di.vMeas);
						}
						
						if (iloc >= 1 && ts < di.t) {
							BBMM_measurement<double> dim1  = data[iloc-1];
							double sBeta  = (di.t - ts) / (di.t-dim1.t);
						
							// the start of the time interval is in the previous bridge
							sMu = (1-sBeta) * muI + sBeta * dim1.mu;
							sVar = (1-sBeta)*(1-sBeta) * varI + sBeta*sBeta * dim1.vMeas
									+ (di.t-dim1.t)*sBeta*(1-sBeta) * dim1.vDiff;
						} else {
							// The start of the time interval is in the same bridge as t
							
							// *dts is a negative number; negate the sign of the denominator
							double sBeta = *dts / (di.t - t);
						
							sMu = (1-sBeta) * pos + sBeta * muI;
							sVar = sBeta*sBeta * varI
									+ (t-di.t) * sBeta * (1-sBeta) * di.vDiff;
						}
					}
					
					// determine the position distribution at the end of the time interval
					if (tf == t) {
						fMu  = pos;
						fVar = 0.0;
					} else if (iloc < *nloc-3 && tf > data[iloc+2].t) {
							// end of interval at least two bridges away;
							// location independent of fixed position at t
							BBMM_measurement<double> dS = data[iloc+2];
							BBMM_measurement<double> dF = data[iloc+3];
							ptrdiff_t i = iloc + 3;
							while (i < *nloc-1 && tf > dF.t) {
								dS = dF;
								dF = data[++i];
							}
						
							double fAlpha = (tf - dS.t) / (dF.t - dS.t);
						
							fMu  = (1-fAlpha) * dS.mu + fAlpha * dF.mu;
							fVar = (1-fAlpha) * (1-fAlpha) * dS.vMeas
									+ fAlpha * fAlpha * dF.vMeas
									+ (dF.t - dS.t) * (1-fAlpha) * fAlpha * dS.vDiff;
					} else {
						double fAlpha = (t - di.t) / (di1.t-di.t);
						
						// See section 3.1.1 of MSc thesis for derivation of these quantities
						Vec2D<double> muI1 = di1.mu;
						double varI1 = di1.vMeas;
						if (fAlpha > 0) {
							Vec2D<double> muQ = (fAlpha-1)/fAlpha * di.mu + pos/fAlpha;
							double varQ = ((1-fAlpha)/fAlpha)*((1-fAlpha)/fAlpha) * di.vMeas
										+ (di1.t-di.t) * (1-fAlpha) / fAlpha * di.vDiff;
						
							muI1 = (muQ * di1.vMeas + di1.mu * varQ)
									/ (varQ + di1.vMeas);
							varI1 = varQ * di1.vMeas / (varQ + di1.vMeas);
						}
						
						if (iloc < *nloc-2 && tf > di1.t) {
							BBMM_measurement<double> di2  = data[iloc+2];
							double fBeta  = (tf - di1.t) / (di2.t-di1.t);
						
							// the end of the time interval is in the next bridge
							fMu = (1-fBeta) * muI1 + fBeta * di2.mu;
							fVar = (1-fBeta)*(1-fBeta) * varI1 + fBeta*fBeta * di2.vMeas
									+ (di2.t-di1.t)*fBeta*(1-fBeta) * di1.vDiff;
						} else {
							// The end of the time interval is in the same bridge as t
							double fBeta = *dtf / (di1.t - t);
						
							fMu = (1-fBeta) * pos + fBeta * muI1;
							fVar = (fBeta)*(fBeta) * varI1
									+ (di1.t-t) * fBeta * (1-fBeta) * di.vDiff;
						}
					}
					
					// get parameters for Rician speed distribution
					double nu = (fMu-sMu).norm() / (*dtf-*dts);
					double sigma = sqrt(sVar + fVar) / (*dtf-*dts);
					
					// compute mean of Rice(nu, sigma)
					double lx = -nu*nu / (2*sigma*sigma);
					double avgSpeed = sigma * sqrt(M_PI / 2)
							* ((1-lx) * bessel_I0_scaled(-lx/2) - lx * bessel_I1_scaled(-lx/2));
					
					// update weighted average of speed at pos
					wSum[x + *nrow*y]   += w * avgSpeed;
					weight[x + *nrow*y] += w;
				}
			}
		}
	}
}
	
} // extern "C"

