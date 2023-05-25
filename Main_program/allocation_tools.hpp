#include <complex>
using namespace std;

complex<double>** alloc_zmatrix(int n);
complex<double>** alloc_zmatrix(int m, int n);
void delete_zmatrix(complex<double>** mat);

complex<double>*** alloc_zpmatrix(int n);
void delete_zpmatrix(complex<double>*** mat);

double** alloc_dmatrix(int n);
double** alloc_dmatrix(int m, int n);
void delete_dmatrix(double** mat);

int** alloc_imatrix(int m, int n);

double*** alloc_dcube(int l, int m, int n);
