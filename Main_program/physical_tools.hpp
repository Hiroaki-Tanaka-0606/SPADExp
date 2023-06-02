#include <complex>

complex<double>**** convert_LCAO(int sk, int sb, int sm, int sd, double* LCAO_raw);

void outer_product(double* a, double* b, double* c);

double inner_product(double* a, double* b);

void operator_coefficient(char* polarization, double theta, double phi, complex<double>* Y_coeff);

void electric_field_vector(char* polarization, double theta, double phi, complex<double>* e_vec);

void spherical_harmonics(double* r, complex<double>* Ylm);

double Gaunt(int lp, int mp, int l, int m);

double sp_bessel(int l, double x);

double interpolate_potential(double* r, int* count, double* cube, double* atom_cell_buffer);

void Fourier_expansion_z(double* V_buffer, int* V_count, double* g, double* atom_cell_buffer, complex<double>* Vg, double* Vg_abs);

void solve_final_state_Numerov(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int z_start, int V00_index, complex<double>** Vgg_buffer, double* g_vec_buffer, complex<double>* FP_loc_buffer);

void solve_final_state_Matrix(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int z_start, int V00_index, complex<double>** Vgg_buffer, double* g_vec_buffer, complex<double>* FP_loc_buffer, complex<double>** left_matrix, complex<double>** right_matrix, int calc_mode, complex<double>* loc_edge);

void solve_final_state_from_bulk(double Ekin, double* k_para, double kz, int g_count, int z_count, int bulk_count, double dz, int z_start, int V00_index, complex<double>** Vgg_buffer, double* g_vec_buffer, complex<double>* bulk_z_buffer, complex<double>* FP_loc_buffer, complex<double>** left_matrix, complex<double>** right_matrix, complex<double>* bulk_coefs);

void solve_final_state_real_space(double Ekin, double* k_para, double kz, double* atom_cell_buffer, int* VKS_count, double* VKS_buffer, complex<double>* FP_loc_buffer);

complex<double> interpolate_fgz(double z, complex<double>* fgz, double dz, int z_count);

double spherical_harmonic_theta(int l, int m, double theta);

void solve_nonlocal_wfn(double Ekin, int l, int r_count, double* r_arr, double* V_loc, int VPS_l_length, int* VPS_l, double* VPS_E_buffer, double** VPS_nonloc_buffer, double* psi);

void calc_bulk_potential(complex<double>* Vgg, double dz_slab, int slab_count, double z_offset, double bulk_height, double dz_bulk, int z_count_bulk, int bulk_count, complex<double>* Vgg_average, double* Vgg_stdev);

complex<double> Fourier_expansion_1D(complex<double>* Vgg_bulk, double g, double dz, int count);

int solve_final_states_bulk(double Ekin, double* k_para, double gz, int g_count, double* g_vec_buffer, complex<double>* Vgg_buffer, int kz_count, double* dispersion_kz, double* dispersion_buffer, complex<double>*** final_states_pointer, double** kz_pointer);

double* interpolate_wfn(int wfn_length, double* wfn, double* r, int wfn_length_reduced, double dr);

void calc_bulk_dispersion(double* k_para, int kz_count, double* kz, int g_count, double* g_vec_buffer, complex<double>* Vgg_buffer, double* eigen_buffer);
