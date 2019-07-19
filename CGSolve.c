#include "LAPACK.h"


void CGSolve(int* N, double** A, double* b, double* threshold, int* iterMax, int* iteration, double* x, double* alpha, double* beta, double** r){
	*iteration=100;

	char TRANS='T';
	double COEFF=1;
	int INC=1;

	dgemv_(&TRANS, N, N, &COEFF, &A[0][0], N, &b[0], &INC, &COEFF, &x[0], &INC);
}

	
