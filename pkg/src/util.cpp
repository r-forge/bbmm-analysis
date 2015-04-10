#include <R.h>
#include <Rinternals.h>
#include "util.h"

SEXP initIntegrate(SEXP diff) {
	if (!isFunction(diff)) {
		error("The diffusion coefficient must be a function");
	}

	// Construct a function call to the R function integrate(f, l, h)
	// and set or allocate arguments
	SEXP dInt;
	PROTECT(dInt = allocList(4));
	SET_TYPEOF(dInt, LANGSXP);
	SETCAR   (dInt, install("integrate"));
	SETCADR  (dInt, diff);
	SETCADDR (dInt, allocVector(REALSXP, 1));
	SETCADDDR(dInt, allocVector(REALSXP, 1));
	
	UNPROTECT(1);
	return dInt;
}

double integrate(SEXP call, double l, double h) {
	// Set the limits of integration in the call
	REAL(CADDR (call))[0] = l;
	REAL(CADDDR(call))[0] = h;
	
	// Call the R function
	SEXP iVal = eval(call, R_GlobalEnv);
	
	// TODO: check whether the message says "OK"
	
	// Extract the value of the result
	return REAL(VECTOR_ELT(iVal, 0))[0];
}
