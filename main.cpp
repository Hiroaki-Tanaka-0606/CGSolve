#include <stdio.h>
#include "matrix_util.h"
#include "CGSolve.h"
#include <time.h>

extern "C"{
	void dgesv_(int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
}

const char* format="%12.5f ";

int main(){
	int N=160;
	double threshold=1e-10;
	int iterMax=1000;
	double** A=alloc_dmatrix(N,N);
	double* b=alloc_dvector(N);
	double* x=alloc_dvector(N);
	double* alpha=alloc_dvector(iterMax);
	double* beta=alloc_dvector(iterMax);
	double** r=alloc_dmatrix(iterMax,N);
	int iteration=0;

	
	int i,j;
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			if(i==j){
				A[i][j]=2;
			}else if(i==j-1 || i==j+1){
				A[i][j]=1;
			}else{
				A[i][j]=0;
			}
		}
		if(i==0){
			b[i]=1;
		}else{
			b[i]=0;
			}
		//b[i]=i+1;
		x[i]=0;
	}

	//print matrix A
	printf("A: %d * %d matrix\n",N,N);
	for(i=0;i<N;i++){
		for(j=0;j<N;j++){
			printf(format,A[i][j]);
		}
		printf("\n");
	}

	//print vector b
	printf("b: %d * 1 vector\n",N);
	printf("transposed: ");
	for(i=0;i<N;i++){
		printf(format,b[i]);
	}
	printf("\n\n");

	//solve Ax=b by CG method
	printf("Solve A*x=b by CG method\n");
	bool verbose=true;
	CGSolve(&N,A,b,&threshold,&iterMax,&iteration,x,alpha,beta,r,&verbose);

	printf("Calculation end with %d iterations\n",iteration);
	if(iteration==iterMax){
		printf("Solution is not converged\n");
	}
	printf("Solution: x(transposed)=");
	for(i=0;i<N;i++){
		printf(format,x[i]);
	}
	printf("\n\n");

	//solve Ax=b by LAPACK 
	printf("Solve A*x=b by LAPACK dgesv\n");
	int NRHS=1;
	int INFO;
	int* IPIV=alloc_ivector(N);
	dgesv_(&N, &NRHS, &A[0][0], &N, &IPIV[0], &b[0], &N, &INFO);
	if(INFO==0){
		printf("Solution: x(transposed)=");
		for(i=0;i<N;i++){
			printf(format,b[i]);
		}
		printf("\n");
	}else{
		printf("Not Successful\n");
	}



	//solve (A+sI)x=b by shifted CG method
	printf("Solve (A+sI)x=b by shifted CG method\n");
	double s=0;
	double ds=0.0001;
	int count=100000;
	int t;
	printf("s=%.5f+%.5ft, t=1 to %d\n",s,ds,count);
	clock_t start=clock();
	for(t=0;t<count;t++){
		s+=ds;
		ShiftedCGSolve(&N,alpha,beta,r,&iteration,&s,x);
		/*
		printf("s=%.3f, Solution: x(transposed)=",s);
		for(i=0;i<N;i++){
			printf(format,x[i]);
		}
		printf("\n");*/
	}
	clock_t end=clock();
	double time=static_cast<double>(end-start)/CLOCKS_PER_SEC;
	printf("Calculation end\n");
	printf("Time: ");
	printf(format,time);
	printf(" sec.\n");

	//solve (A+sI)x=b by LAPACk
	printf("Solve (A+sI)x=b by LAPACK dgesv\n");
	s=0;
	printf("s=%.5f+%.5ft, t=1 to %d\n",s,ds,count);
  start=clock();
	for(t=0;t<count;t++){
		s+=ds;
		
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				if(i==j){
					A[i][j]=2+s;
				}else if(i==j-1 || i==j+1){
					A[i][j]=1;
				}else{
					A[i][j]=0;
				}
			}
			if(i==0){
				b[i]=1;
			}else{
				b[i]=0;
			}
			//b[i]=i+1;
		}
		
		dgesv_(&N, &NRHS, &A[0][0], &N, &IPIV[0], &b[0], &N, &INFO);
		/*
		if(INFO==0){
			printf("s=%.3f, Solution: x(transposed)=",s);
			for(i=0;i<N;i++){
				printf(format,b[i]);
			}
			printf("\n");
		}else{
			printf("Not Successful\n");
			}*/
	}
	
  end=clock();
	time=static_cast<double>(end-start)/CLOCKS_PER_SEC;
	printf("Calculation end\n");
	printf("Time: ");
	printf(format,time);
	printf(" sec.\n");
}
