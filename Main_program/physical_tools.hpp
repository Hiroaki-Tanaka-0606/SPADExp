#include <complex>

complex<double>**** convert_LCAO(int sk, int sb, int sm, int sd, double* LCAO_raw);

void outer_product(double* a, double* b, double* c);

double inner_product(double* a, double* b);

void operator_coefficient(char* polariation, double theta, double phi, complex<double>* Y_coeff);

void spherical_harmonics(double* r, complex<double>* Ylm);

double Gaunt(int lp, int mp, int l, int m);

double sp_bessel(int l, double x);

double interpolate_potential(double* r, int* count, double* cube, double* rec_cell);

void Fourier_expansion_z(double* V_buffer, int* V_count, double* g, double* atom_cell_buffer, complex<double>* Vg, double* Vg_abs);

void solve_final_state(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int V00_index, complex<double>** Vgg_buffer, double* g_vec_buffer, complex<double>* FP_loc_buffer);
