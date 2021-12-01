// calcPSF photoemission structure factor calculation

#include <iostream>
#include <complex>
#include <chrono>
#include <omp.h>

using namespace std;
// test program of BLAS and OMP parallelization
extern "C"{
	double ddot_(
							 int* N,
							 double* CX,
							 int* INCX,
							 double* CY,
							 int* INCY);
	// in case of complex<double>, the result is the first argument
	void zdotu_(
												 complex<double>* result,
												 int* N,
												 complex<double>* CX,
												 int* INCX,
												 complex<double>* CY,
												 int* INCY);
}

complex<double> zdot(int* N, complex<double>* X, complex<double>* Y){
	int INC=1;
  complex<double> result;
	zdotu_(&result, N, &X[0], &INC, &Y[0], &INC);
	return result;
}

complex<double> zdot2(int* N, complex<double>* X, complex<double>* Y){
	complex<double> ret=complex<double>(0, 0);
	int i;
	for(i=0; i<*N; i++){
		ret+=X[i]*Y[i];
	}
	return ret;
}

double ddot(int* N, double* X, double* Y){
	int INC=1;
	return ddot_(N, &X[0], &INC, &Y[0], &INC);
}

void PSF_test(){
	int N=128;
	complex<double> matrix[N][N];
	complex<double> vec[N];
	int i,j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			matrix[i][j]=complex<double>(i+2*j, 0);
		}
		vec[i]=complex<double>(i*i, 0);
	}
	/*
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%4.1f ", matrix[i][j].real());
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++){
		printf("%4.1f ", vec[i].real());
	}
	printf("\n");*/

	
	complex<double> dot_product;
	int tMax=10000;
	int t;
	chrono::system_clock::time_point start=chrono::system_clock::now();
	for(t=0; t<tMax; t++){
		for(i=0; i<N; i++){
			dot_product=zdot(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product.real());
		}
	}
	chrono::system_clock::time_point end=chrono::system_clock::now();
	double duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Time1: %.3f [ms]\n", duration);

	complex<double> dot_product2;
	start=chrono::system_clock::now();
#pragma omp parallel shared(N, matrix, vec) private(dot_product2)
#pragma omp for
	for(i=0; i<N; i++){
		for(t=0; t<tMax; t++){
			dot_product2=zdot(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product2.real());
		}
	}
	end=chrono::system_clock::now();
	duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Time2: %.3f [ms]\n", duration);

	return;

	/*
	start=clock();
	for(t=0; t<tMax; t++){
		for(i=0; i<N; i++){
			dot_product=zdot2(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product.real());
		}
	}
  end=clock();
	duration=(double)(end-start)/CLOCKS_PER_SEC*1000;
	printf("Time3: %.3f [ms]\n", duration);
	*/
	/*
	double dmatrix[N][N];
	double dvec[N];
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			dmatrix[i][j]=i+2*j;
		}
		dvec[i]=i*i;
	}
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%4.1f ", dmatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++){
		printf("%4.1f ", dvec[i]);
	}
	printf("\n");
	
	double ddot_product;
	for(i=0; i<N; i++){
		ddot_product=ddot(&N, dmatrix[i], dvec);
		printf("M[%d] * V = %4.1f\n", i, ddot_product);
		}*/
	
}
