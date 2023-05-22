#include <complex>
using namespace std;

complex<double>** alloc_zmatrix(int n){
	complex<double>* buffer=new complex<double>[n*n];
	complex<double>** mat=new complex<double>*[n];
	for(int i=0; i<n; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

complex<double>** alloc_zmatrix(int m, int n){
	complex<double>* buffer=new complex<double>[m*n];
	complex<double>** mat=new complex<double>*[m];
	for(int i=0; i<m; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

complex<double>*** alloc_zpmatrix(int n){
	complex<double>** buffer=new complex<double>*[n*n];
	complex<double>*** mat=new complex<double>**[n];
	for(int i=0; i<n; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

void delete_zpmatrix(complex<double>*** mat){
	delete[] mat[0];
	delete[] mat;
}

void delete_zmatrix(complex<double>** mat){
	delete[] mat[0];
	delete[] mat;
}


double** alloc_dmatrix(int n){
	double* buffer=new double[n*n];
	double** ret=new double*[n];
	for(int i=0; i<n; i++){
		ret[i]=&buffer[i*n];
	}
	return ret;
}

double** alloc_dmatrix(int m, int n){
	double* buffer=new double[m*n];
	double** ret=new double*[m];
	for(int i=0; i<m; i++){
		ret[i]=&buffer[i*n];
	}
	return ret;
}

void delete_dmatrix(double** mat){
	delete[] mat[0];
	delete[] mat;
}
