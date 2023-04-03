#include <complex>

complex<double>**** convert_LCAO(int sk, int sb, int sm, int sd, double* LCAO_raw);

void outer_product(double* a, double* b, double* c);

double inner_product(double* a, double* b);

void operator_coefficient(char* polariation, double theta, double phi, complex<double>* Y_coeff);

void spherical_harmonics(double* r, complex<double>* Ylm);

double Gaunt(int lp, int mp, int l, int m);

double sp_bessel(int l, double x);

double interpolate_potential(double* r, int* count, double* cube, double* rec_cell);

double Fourier_expansion_z(double* V_buffer, int* V_count, double* g, double* atom_cell_buffer, complex<double>* Vg);
