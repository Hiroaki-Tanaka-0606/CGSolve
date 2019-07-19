#include "LAPACK.h"
#include "InnerProduct.h"
#include "CGSolve.h"
#include "matrix_util.h"

#include <stdio.h>

void CGSolve(int* N, double** A, double* b, double* threshold, int* iterMax, int* iteration, double* x, double* alpha, double* beta, double** r, bool* verbose){
	//solve linear equation A*x=b by CG method

	double* p=alloc_dvector(*N);
	//constant for dgemv_
	char TRANS='T';
	double COEFF=1;
	int INC=1;

	int i,j;
	//initial value:
	//x[0]=zero vector
	//r[0]=b
	//p[0]=r[0]
	for(j=0;j<*N;j++){
		r[0][j]=b[j];
		p[j]=b[j];
		x[j]=0;
	}
	double rNorm=TwoNorm(*N,r[0]);
	for(*iteration=0; *iteration<*iterMax-1 && rNorm>*threshold;*iteration=i+1){
		i=*iteration;
		//alpha[i]=(r[i],r[i])/(p[i],Ap[i])
		alpha[i]=rNorm/ANorm(*N,p,A);
		//x[i+1]=x[i]+alpha[i]p[i]
		for(j=0;j<*N;j++){
			x[j]=x[j]+alpha[i]*p[j];
		}
		//r[i+1]=r[i]-alpha[i]Ap[i]
		for(j=0;j<*N;j++){
			r[i+1][j]=r[i][j];
		}
		double minusAlpha=-alpha[i];
		dgemv_(&TRANS, N, N, &minusAlpha, &A[0][0], N, &p[0], &INC, &COEFF, &r[i+1][0],&INC);
		//beta[i]=(r[i+1],r[i+1])/(r[i],r[i])
		double nextrNorm=TwoNorm(*N,r[i+1]);
		beta[i]=nextrNorm/rNorm;
		//p[i+1]=r[i+1]+beta[i]p[i];
		for(j=0;j<*N;j++){
			p[j]=r[i+1][j]+beta[i]*p[j];
		}
		rNorm=nextrNorm;
		if(*verbose){
			printf("i=%4d, rNorm=%10.5e\n",i,rNorm);
		}
	}
}

void ShiftedCGSolve(int* N, double* alpha, double* beta, double** r, int* iteration, double* s, double* x){
	int i,j;
	double* pi=alloc_dvector(*iteration+1);
	double* p=alloc_dvector(*N);

	//initial value
	//pi[0]=1
	pi[0]=1;
	//p[0]=r[0]
	//x[0]=zero vector
	for(j=0;j<*N;j++){
		p[j]=r[0][j];
		x[j]=0;
	}
	//pi[1]=(1+alpha[0]s)pi[0]
	pi[1]=(1+alpha[0]*(*s))*pi[0];
	//x[1]=x[0]+pi[0]/pi[1]alpha[0]p[0]
	for(j=0;j<*N;j++){
		x[j]=x[j]+pi[0]/pi[1]*alpha[0]*p[j];
	}
	//p[1]=r[1]/pi[1]+beta[0]*(1/pi[1])^2*p[0]
	for(j=0;j<*N;j++){
		p[j]=r[1][j]/pi[1]+beta[0]/(pi[1]*pi[1])*p[j];
	}
	for(i=1;i<*iteration;i++){
		//pi[i+1]=(1+alpha[i]s)pi[i]-(alpha[i]beta[i-1]/alpha[i-1])*(pi[i-1]-pi[i])
		pi[i+1]=(1+alpha[i]*(*s))*pi[i]-(alpha[i]*beta[i-1]/alpha[i-1])*(pi[i-1]-pi[i]);
		//x[i+1]=x[i]+pi[i]/pi[i+1]alpha[i]p[i]
		for(j=0;j<*N;j++){
			x[j]=x[j]+pi[i]/pi[i+1]*alpha[i]*p[j];
		}
		if(i!=*iteration-1){
			//p[i+1]=r[i+1]/pi[i+1]+beta[i]*(pi[i]/pi[i+1])^2*p[i]
			for(j=0;j<*N;j++){
				p[j]=r[i+1][j]/pi[i+1]+beta[i]*(pi[i]/pi[i+1])*(pi[i]/pi[i+1])*p[j];
			}
		}
	}
}
