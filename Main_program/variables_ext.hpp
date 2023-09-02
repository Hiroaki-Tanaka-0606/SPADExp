// calcPSF variables declaration
// description are in variables.hpp
#include <cstdio>

// constants
extern int buffer_length;
extern int Calculation_length;
extern int Log_file_length;
extern int Log_length;
extern char* Log_buffer;
extern int Solution_length;
extern int Potential_length;
extern int Potential_file_length;
extern double Data_read_error;
extern int HDF5_file_length;
extern int PA_state_length;
extern int Ph_label_length;
extern double au_ang;
extern double Eh;

// variables
/// blocks
extern bool Ct_block_appeared;
extern bool TF_block_appeared;
extern bool Rg_block_appeared;
extern bool At_block_appeared;
extern bool SC_block_appeared;
extern bool Oc_block_appeared;
extern bool PA_block_appeared;
extern bool Ex_block_appeared;
extern bool Or_block_appeared;
extern bool Ph_block_appeared;
/// Ct block
extern char* Calculation;
extern bool Calculation_set;
extern char* Log_file;
extern bool Log_file_set;
extern bool Console_log;
extern bool Console_log_set;
extern char* Output_file;
extern bool Output_file_set;
/// Rg block
extern int Radial_grid_count;
extern double* Radial_grid_intervals;
extern int* Radial_grid_points;
extern double* x_coordinates;
extern int x_count;
/// TF block
extern bool TF_test;
extern bool TF_test_set;
extern double Initial_diff_offset;
extern bool Initial_diff_offset_set;
extern double Initial_diff_delta;
extern bool Initial_diff_delta_set;
extern int Initial_diff_size;
extern bool Initial_diff_size_set;
extern double Initial_diff_min;
extern bool Initial_diff_min_set;
extern double Initial_diff_max;
extern bool Initial_diff_max_set;
extern double TF_threshold;
extern bool TF_threshold_set;
extern char* TF_solution;
extern double TF_solution_set;
extern double* TF_phi;
extern double* TF_phi_diff;
/// At block
extern int n_min;
extern bool n_min_set;
extern int n_max;
extern bool n_max_set;
extern int n_single;
extern int n_single_set;
extern int l_min;
extern bool l_min_set;
extern int l_max;
extern bool l_max_set;
extern int l_single;
extern int l_single_set;
extern int Z_min;
extern bool Z_min_set;
extern int Z_max;
extern bool Z_max_set;
extern int Z_single;
extern bool Z_single_set;
extern char* At_potential;
extern bool At_potential_set;
extern char* At_potential_file;
extern bool At_potential_file_set;
extern char* At_solution;
extern bool At_solution_set;
extern double Bisection_step;
extern double Bisection_step_set;
extern double At_radius_factor;
extern bool At_radius_factor_set;
extern double At_E_threshold;
extern bool At_E_threshold_set;
extern double At_initial_diff;
extern double At_initial_diff2;
extern double* At_p_x;
extern double* At_p_diff_x;
extern double* At_v_x;
extern double At_bisection_threshold;
extern int At_min_iteration;
extern int At_max_iteration;
/// SC block
extern double SC_mix_weight;
extern bool SC_mix_weight_set;
extern double SC_criterion_a;
extern bool SC_criterion_a_set;
extern double SC_criterion_b;
extern bool SC_criterion_b_set;
extern double* SC_sigma_x;
extern int SC_orbital_count;
extern int* SC_n_list;
extern int* SC_l_list;
extern double* SC_eigen_list;
extern double** SC_p_x;
/// Oc block
extern int Occupation_count;
extern int** At_occupation;
/// Ex block
extern int Ex_energy_count;
extern double* Ex_energies;
/// Or block
extern int Ph_orbital_count;
extern char** Ph_orbital_labels;
extern int* Ph_l_list;
extern double* Ph_binding_energies;
extern int Ph_skip_points;
extern bool Ph_skip_points_set;
extern int Ph_calc_points;
extern bool Ph_calc_points_set;
/// PA block
extern char* PA_input_file;
extern bool PA_input_file_set;
extern double PA_E_min;
extern bool PA_E_min_set;
extern double PA_E_max;
extern bool PA_E_max_set;
extern double PA_dE;
extern bool PA_dE_set;
extern double PA_E_pixel;
extern bool PA_E_pixel_set;
extern char* PA_initial_state;
extern bool PA_initial_state_set;
extern char* PA_final_state;
extern bool PA_final_state_set;
extern double PA_final_state_step;
extern bool PA_final_state_step_set;
extern char* PA_polarization;
extern bool PA_polarization_set;
extern double PA_theta;
extern double PA_phi;
extern bool PA_theta_set;
extern bool PA_phi_set;
extern int PA_ext_up;
extern int PA_ext_ri;
extern int PA_ext_dn;
extern int PA_ext_le;
extern bool PA_ext_set;
extern char* PA_AO_file;
extern bool PA_AO_file_set;
extern bool PA_weighting;
extern bool PA_weighting_set;
extern double* PA_weighting_axis;
extern bool PA_weighting_axis_set;
extern char* PA_weighting_shape;
extern bool PA_weighting_shape_set;
extern double PA_weighting_origin;
extern bool PA_weighting_origin_set;
extern double PA_weighting_width; 
extern bool PA_weighting_width_set;
extern bool PA_use_angstrom;
extern bool PA_use_angstrom_set;
extern char* PA_output_data;
extern bool PA_output_data_set;
extern double PA_sigma_max;
extern double PA_decay_max;
extern bool PA_reflection;
extern bool PA_reflection_set;
extern double PA_reflection_coef;
extern double PA_reflection_coef_set;
extern bool PA_include_neg_depth;
extern bool PA_include_neg_depth_set;
extern double PA_excitation_energy;
extern bool PA_excitation_energy_set;
extern double PA_FPFS_energy_step;
extern bool PA_FPFS_energy_step_set;
extern int PA_FPFS_range;
extern double PA_FPFS_kRange;
extern bool PA_FPFS_kRange_set;
extern bool PA_FPFS;
extern char* PA_FPFS_file;
extern bool PA_FPFS_file_set;
extern char* PA_VPS_file;
extern bool PA_VPS_file_set;
extern int PA_Lebedev_order_ave;
extern int PA_Lebedev_order_int;
extern int PA_lp_max;
extern double PA_zero_threshold;
extern int PA_theta_points;
extern bool PA_ignore_core;
extern bool PA_ignore_core_set;
extern bool PA_add_nonorth_term;
extern bool PA_add_nonorth_term_set;
extern bool PA_calc_all_nonloc;
extern bool PA_calc_all_nonloc_set;
//extern bool PA_calc_all_loc;
//extern bool PA_calc_all_loc_set;
extern bool PA_orth_correction;
extern bool PA_orth_correction_set;
extern bool PA_ignore_nonlocal;
extern bool PA_ignore_nonlocal_set;
extern bool PA_FPFS_bulk_include_nonlocal;
extern bool PA_FPFS_bulk_include_nonlocal_set;
extern bool PA_FPFS_Numerov;
extern bool PA_FPFS_Numerov_set;
extern double PA_FPFS_bulk_volume;
//extern double PA_FPFS_margin;
//extern bool PA_FPFS_edge_smoothing;
//extern bool PA_FPFS_edge_smoothing_set;
//extern int PA_FPFS_smoothing_E;
//extern bool PA_FPFS_smoothing_E_set;
//extern int PA_FPFS_smoothing_k;
//extern bool PA_FPFS_smoothing_k_set;
extern bool PA_interpolate_wfn;
extern bool PA_interpolate_wfn_set;
extern double PA_interpolate_wfn_coef;
extern bool PA_interpolate_wfn_coef_set;
extern double PA_FPFS_bulk_min_ang;
extern double PA_FPFS_bulk_max_ang;
extern int PA_FPFS_bulk_count;
extern bool PA_FPFS_bulk_set;
extern int PA_FPFS_bulk_kz_steps;
extern int PA_FPFS_bulk_kappaz_steps_left;
extern int PA_FPFS_bulk_kappaz_steps_right;
extern int PA_FPFS_bulk_buffer_size;
extern double PA_FPFS_negligible_gap_size;
extern bool PA_FPFS_negligible_gap_size_set;
//extern bool PA_FPFS_perturbation;
//extern bool PA_FPFS_perturbation_set;
//extern bool PA_calc_complex_dispersion;
//extern bool PA_calc_complex_dispersion_set;
extern double PA_FPFS_real_eigenvalue_criterion;
extern bool PA_FPFS_real_eigenvalue_criterion_set;
extern double PA_FPFS_dispersion_group_criterion_left;
extern double PA_FPFS_dispersion_group_criterion_right;
extern bool PA_FPFS_dispersion_group_criterion_left_set;
extern bool PA_FPFS_dispersion_group_criterion_right_set;
extern double PA_FPFS_gap_upper_limit;
extern bool PA_FPFS_gap_upper_limit_set;
extern int PA_FPFS_kz_margin_index_size;
extern bool PA_FPFS_kz_margin_index_size_set;
extern double PA_FPFS_cspace_size;
extern bool PA_FPFS_cspace_size_set;
extern double PA_FPFS_kz_exclude_criterion;
extern bool PA_FPFS_kz_exclude_criterion_set;
extern double PA_FPFS_nonloc_offset;
extern bool PA_FPFS_nonloc_offset_set;

/// files
extern FILE* Log_file_obj;
extern FILE* Output_file_obj;

