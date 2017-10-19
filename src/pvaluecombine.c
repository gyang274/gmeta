#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include <R.h> 
#include <Rmath.h>
#include <Rdefines.h>
#include <Rinternals.h>

//using namespace std;

double nFactorial(double n);
double nChoosek(double n, double k);
double pConvolveUniform( double q, double n );

SEXP pvaluecombine( SEXP RpVec, SEXP Rmethod ) {
	int k = length(RpVec);
	const char * method = CHAR(STRING_ELT(Rmethod, 0));
	
	SEXP Rcmbdpvalue = PROTECT(allocVector(REALSXP, 1));
	memset(REAL(Rcmbdpvalue), 0.0, sizeof(double));
	
	double * cmbdpvalue = REAL(Rcmbdpvalue);
	if (!strcmp(method, "fisher")) {
		for (int i=0; i<k; i++) {
			*cmbdpvalue += log(REAL(RpVec)[i]);
		}
		*cmbdpvalue = 1 - pchisq(-2 * *cmbdpvalue, 2*k, 1, 0);
	} else if (!strcmp(method, "normal") || !strcmp(method, "stouffer")) {
		for (int i=0; i<k; i++) {
			*cmbdpvalue += qnorm(REAL(RpVec)[i], 0.0, 1.0, 1, 0);
		}
		*cmbdpvalue = *cmbdpvalue / sqrt(k);
		*cmbdpvalue = pnorm(*cmbdpvalue, 0.0, 1.0, 1, 0);
	} else if (!strcmp(method, "min") || !strcmp(method, "tippett")) {
		*cmbdpvalue = REAL(RpVec)[0];
		for (int i=1; i<k; i++) {
			*cmbdpvalue = fmin2(*cmbdpvalue, REAL(RpVec)[i]);
		}
		*cmbdpvalue = 1 - pow(1-*cmbdpvalue, k);
	} else if (!strcmp(method, "max")) {
		*cmbdpvalue = REAL(RpVec)[0];
		for (int i=1; i<k; i++) {
			*cmbdpvalue = fmax2(*cmbdpvalue, REAL(RpVec)[i]);
		}
		*cmbdpvalue = pow(*cmbdpvalue, k);
	} else if (!strcmp(method, "sum")) {
		for (int i=0; i<k; i++) {
			*cmbdpvalue += REAL(RpVec)[i];
		}
		if (k <= 30) {
			*cmbdpvalue = pConvolveUniform(*cmbdpvalue, (double)k);
		} else {
			*cmbdpvalue = pnorm(*cmbdpvalue, (double)k/2.0, sqrt((double)k/12.0), 1, 0);
		}
	} else {
		*cmbdpvalue = 3.1415926;
	}
	// return
	UNPROTECT(1);
	return(Rcmbdpvalue);
}

double pConvolveUniform( double q, double n ) {
	double p = 0.0;
	for(int i = 0; i < n; i++) {
		p += (double)((-1) * ((i % 2)*2-1) * pow((0<(q-(double)i))?(q-(double)i):0, n) / (nFactorial(n-i) * nFactorial(i)));
	}
	p = p>1.0?1.0:p;
	return(p);
}

double nChoosek( double n, double k ) {
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    
    int i; double rslt = n;
    for( i=2; i <= k; ++i ) {
        rslt *= (n-i+1);
        rslt /= i;
    }
    return rslt;
}

double nFactorial(double n) {
    if( n == 0 ) return 1;
    if( n  > 0 ) return n * nFactorial(n-1) ;
	return 0;
}
