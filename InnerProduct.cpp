#include "LAPACK.h"
#include "matrix_util.h"
#include "InnerProduct.h"

double TwoNorm(int N, double* p){
	//calculate p \cdot p
	double norm=0;
	int i;
	for(i=0;i<N;i++){
		norm+=p[i]*p[i];
	}
	return norm;
}

double ANorm(int N, double* p, double** A){
	//calculate p \cdot A*p

	char TRANS='T';
	double COEFF=1;
	int INC=1;

	int i;
	double* Ap=alloc_dvector(N);
	for(i=0;i<N;i++){
		Ap[i]=0;
	}
	dgemv_(&TRANS, &N, &N, &COEFF, &A[0][0], &N, &p[0], &INC, &COEFF, &Ap[0], &INC);
	return InnerProduct(N, p, Ap);
}

double InnerProduct(int N, double* p, double* q){
	//calculate p \cdot q
	double prod=0;
	int i;
	for(i=0;i<N;i++){
		prod+=p[i]*q[i];
	}
	return prod;
}
