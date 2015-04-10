#include "vec2d.h"
#include <math.h>
#include <limits>
#include <algorithm>
#include <R.h>

#include "encounter.h"
#include "util.h"
#include "bessel.h"
#include <Rinternals.h>
#include <Rdefines.h>

#define BB_SCALE 18.42068 // TODO: explain

extern "C" {

SEXP speedDistribution(SEXP burst, SEXP rasterCoords,
		SEXP timeStep, SEXP dt) {
	SEXP result, rSum, rWeight, xc, yc, dInt;
	
	double tStep = REAL(timeStep)[0], dts = REAL(dt)[0], dtf = REAL(dt)[1];
	
	xc = VECTOR_ELT(rasterCoords, 0);
	yc = VECTOR_ELT(rasterCoords, 1);
	double *xcd = REAL(xc), *ycd = REAL(yc);
	
	// allocate output data structure: a list of two matrices
	PROTECT(result = allocVector(VECSXP, 2));
	PROTECT(rSum = allocMatrix(REALSXP, length(xc), length(yc)));
	SET_VECTOR_ELT(result, 0, rSum);
	PROTECT(rWeight = allocMatrix(REALSXP, length(xc), length(yc)));
	SET_VECTOR_ELT(result, 1, rWeight);
	
	double *rs = REAL(rSum);
	double *rw = REAL(rWeight);
	UNPROTECT(2); // rWeight, rSum
	
	// Zero the result matrices
	for (ptrdiff_t i = 0; i < length(xc) * length(yc); i++) {
			rs[i] = 0.0;
			rw[i] = 0.0;
	}
	
	// get some data in an easier to handle form
	Vec2D<double> *dMu = (Vec2D<double>*)REAL(VECTOR_ELT(burst, 0));
	double *ts  = REAL(VECTOR_ELT(burst, 1));
	double *var = REAL(VECTOR_ELT(burst, 2));
	ptrdiff_t nloc = length(VECTOR_ELT(burst, 1)); // the number of relocations
	
	// construct a call to integrate the diffusion coefficient
	// burst@diff is the 4th element of burst, it is itself a one-element list.
	PROTECT(dInt = initIntegrate(VECTOR_ELT(VECTOR_ELT(burst, 3), 0)));
	
	ptrdiff_t iloc = -1; // the index of the measurement directly before t
	
	// precompute integral of diffusion coefficient over all measurement intervals
	// intDiff[i] holds the integral of diff coeff over [ts[i], ts[i+1]]
	double intDiff[nloc-1];
	for (off_t i = 0; i < nloc-1; i++) {
		intDiff[i] = integrate(dInt, ts[i], ts[i+1]);
	}
	
	ptrdiff_t xmin = 0;
	ptrdiff_t ymin = 0;
	for (double t = ts[0] - dts; t <= ts[nloc-1] - dtf; t += tStep) {
		while (iloc < nloc-1 && ts[iloc+1] < t) {
			iloc++;
		}
		
		double alpha = integrate(dInt, ts[iloc], t) / intDiff[iloc];
		Vec2D<double> pMu = (1-alpha) * dMu[iloc] + alpha * dMu[iloc+1];
		double pVar = (1-alpha)*(1-alpha)*var[iloc] + alpha*alpha*var[iloc+1]
				+ alpha*(1-alpha)*intDiff[iloc];
				
		// TODO: pull all calls to integrate out of the innner loop and put them somewhere around here.
		double sAlpha = std::numeric_limits<double>::quiet_NaN(),
			sBeta = sAlpha, sBetaD = sAlpha,
			fAlpha = sAlpha, fBeta = sAlpha, fBetaD = sAlpha;
		
		// find the indices in the coordinate lists where we start computing
		double dMax = sqrt(BB_SCALE*pVar);
		while (xmin < length(xc)-1 && xcd[xmin]    < pMu.x - dMax) { xmin++; }
		while (xmin > 0            && xcd[xmin-1] >= pMu.x - dMax) { xmin--; }
		while (ymin < length(yc)-1 && ycd[ymin]    < pMu.y - dMax) { ymin++; }
		while (ymin > 0            && ycd[ymin-1] >= pMu.y - dMax) { ymin--; }
		
		for (ptrdiff_t x = xmin; x < length(xc) && xcd[x] <= pMu.x + dMax; x++) {
			for (ptrdiff_t y = ymin; y < length(yc) && ycd[y] <= pMu.y + dMax; y++) {
				Vec2D<double> pos(xcd[x], ycd[y]);
				double d2 = (pos-pMu).norm2();
				
				if (d2 > BB_SCALE*pVar) {
					// We are outside the radius around the mean that we consider
					continue;
				}
				
				double w = pdf_norm2d(d2, pVar);
				
				double tStart  = t + dts;
				double tFinish = t + dtf;
				
				Vec2D<double> sMu, fMu;
				double sVar, fVar;
				// determine the position distribution at the start of the time interval
				if (dts == 0.0) {
					sMu  = pos;
					sVar = 0.0;
				} else if (iloc >= 2 && tStart < ts[iloc-1]) {
					// start of interval at least two bridges away;
					// location independent of fixed position at t
					ptrdiff_t i = iloc - 2;
					while (i > 0 && tStart < ts[i]) {
						i--;
					}
					
					if (sAlpha != sAlpha) { // sAlpha is NaN
						sAlpha = integrate(dInt, ts[i], tStart) / intDiff[i];
					}
				
					sMu  = (1-sAlpha) * dMu[i] + sAlpha * dMu[i+1];
					sVar = (1-sAlpha) * (1-sAlpha) * var[i]
							+ sAlpha * sAlpha * var[i+1]
							+ (1-sAlpha) * sAlpha * intDiff[i];
				} else {
					if (sAlpha != sAlpha) { // sAlpha is NaN
						sAlpha = integrate(dInt, t, ts[iloc+1]) / intDiff[iloc];
					}
					// See section 3.1.1 of MSc thesis for derivation of these quantities
					Vec2D<double> muI = dMu[iloc];
					double varI = var[iloc];
					if (sAlpha > 0) {
						Vec2D<double> muQ = (sAlpha-1)/sAlpha * dMu[iloc+1] + pos/sAlpha;
						double varQ = ((1-sAlpha)/sAlpha)*((1-sAlpha)/sAlpha) * var[iloc+1]
									+ (1-sAlpha) / sAlpha * intDiff[iloc];
					
						muI = (muQ * var[iloc] + dMu[iloc] * varQ)
								/ (varQ + var[iloc]);
						varI = varQ * var[iloc] / (varQ + var[iloc]);
					}
					
					if (iloc >= 1 && tStart < ts[iloc]) {
						if (sBeta != sBeta) { // sBeta is NaN
							sBeta  = integrate(dInt, tStart, ts[iloc])
									/ intDiff[iloc-1];
						}
					
						// the start of the time interval is in the previous bridge
						sMu = (1-sBeta) * muI + sBeta * dMu[iloc-1];
						sVar = (1-sBeta)*(1-sBeta) * varI + sBeta*sBeta * var[iloc-1]
								+ sBeta*(1-sBeta) * intDiff[iloc-1];
					} else {
						// The start of the time interval is in the same bridge as t
						
						// *dts is a negative number; negate the sign of the denominator
						if (sBeta != sBeta) { // sBeta is NaN
							sBetaD = integrate(dInt, ts[iloc], t);
							sBeta = integrate(dInt, tStart, t) / sBetaD;
						}
					
						sMu = (1-sBeta) * pos + sBeta * muI;
						sVar = sBeta*sBeta * varI
								+ sBeta * (1-sBeta) * sBetaD;
					}
				}
				
				// determine the position distribution at the end of the time interval
				if (tFinish == t) {
					fMu  = pos;
					fVar = 0.0;
				} else if (iloc < nloc-3 && tFinish > ts[iloc+2]) {
						// end of interval at least two bridges away;
						// location independent of fixed position at t
						ptrdiff_t i = iloc + 2;
						while (i < nloc-1 && tFinish > ts[i+1]) {
							i++;
						}
						
						if (fAlpha != fAlpha) { // fAlpha is NaN
							fAlpha = integrate(dInt, ts[i], tFinish) / intDiff[i];
						}
					
						fMu  = (1-fAlpha) * dMu[i] + fAlpha * dMu[i+1];
						fVar = (1-fAlpha) * (1-fAlpha) * var[i]
								+ fAlpha * fAlpha * var[i+1]
								+ (1-fAlpha) * fAlpha * intDiff[i];
				} else {
					if (fAlpha != fAlpha) { // fAlpha is NaN
						fAlpha = integrate(dInt, ts[iloc], t) / intDiff[iloc];
					}
					
					// See section 3.1.1 of MSc thesis for derivation of these quantities
					Vec2D<double> muI = dMu[iloc+1];
					double varI = var[iloc+1];
					if (fAlpha > 0) {
						Vec2D<double> muQ = (fAlpha-1)/fAlpha * dMu[iloc] + pos/fAlpha;
						double varQ = ((1-fAlpha)/fAlpha)*((1-fAlpha)/fAlpha) * var[iloc]
									+ (1-fAlpha) / fAlpha * intDiff[iloc];
					
						muI = (muQ * var[iloc+1] + dMu[iloc+1] * varQ)
								/ (varQ + var[iloc+1]);
						varI = varQ * var[iloc+1] / (varQ + var[iloc+1]);
					}
					
					if (iloc < nloc-2 && tFinish > ts[iloc+1]) {
						if (fBeta != fBeta) { // fBeta is NaN
							fBeta  = integrate(dInt, ts[iloc+1], tFinish) / intDiff[iloc+1];
						}
					
						// the end of the time interval is in the next bridge
						fMu = (1-fBeta) * muI + fBeta * dMu[iloc+2];
						fVar = (1-fBeta)*(1-fBeta) * varI + fBeta*fBeta * var[iloc+2]
								+ fBeta*(1-fBeta) * intDiff[iloc+1];
					} else {
						// The end of the time interval is in the same bridge as t
						if (fBeta != fBeta) { // fBeta is NaN
							fBetaD = integrate(dInt, t, ts[iloc+1]);
							fBeta  = integrate(dInt, t, tFinish) / fBetaD;
						}
					
						fMu = (1-fBeta) * pos + fBeta * muI;
						fVar = (fBeta)*(fBeta) * varI
								+ fBeta * (1-fBeta) * fBetaD;
					}
				}
				
				// get parameters for Rician speed distribution
				double nu = (fMu-sMu).norm() / (dtf-dts);
				double sigma = sqrt(sVar + fVar) / (dtf-dts);
				
				// compute mean of Rice(nu, sigma)
				double lx = -nu*nu / (2*sigma*sigma);
				double avgSpeed = sigma * sqrt(M_PI / 2)
						* ((1-lx) * bessel_I0_scaled(-lx/2) - lx * bessel_I1_scaled(-lx/2));
				
				// update weighted average of speed at pos
				rs[x + y * length(xc)] += w * avgSpeed;
				rw[x + y * length(xc)] += w;
			}
		}
	} // for t
	
	
	UNPROTECT(2); // dInt, result
	return result;
} // speedDistribution
	
} // extern "C"

