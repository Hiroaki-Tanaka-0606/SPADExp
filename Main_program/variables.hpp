// calcPSF variables
#include <cstdio>

// constants
int buffer_length;        // size of the buffer used to load one line from the input file
int Calculation_length;   // size of the char* Calculation
int Log_file_length;      // size of the char* Log_file
int Log_length;           // max length of the one log message
char* Log_buffer;         // buffer for the log output
int Solution_length;      // size of the char* TF_solution, At_solution
int Potential_length;     // size of the char* At_potential
int Potential_file_length;// size of the char* At_potential_file
double Data_read_error;   // tolerable error in file reading (such as x_coordinates)
int HDF5_file_length;     // size of the char* PA_input_file, PA_AO_file
int PA_state_length;      // size of the char* PA_initial_state, PA_final_state
int Ph_label_length;      // size of the char* Ph_orbital_labels[i]
double au_ang;            // 1 bohr (Ang)
double Eh;                // 1 Hartree (eV)

// variables
/// blocks
bool Ct_block_appeared;        // true: the block "Control" appeared, false: not
bool TF_block_appeared;        // similar to Ct_block_appeared (description is ignored hereafter)
bool Rg_block_appeared;
bool At_block_appeared;
bool SC_block_appeared;
bool Oc_block_appeared;
bool PA_block_appeared;
bool Ex_block_appeared;
bool Or_block_appeared;
bool Ph_block_appeared;
/// Ct block
char* Calculation;             // name of the calculation
bool Calculation_set;          // true: variable Calculation has a value in the input file
char* Log_file;                // name of the output log file
bool Log_file_set;             // similar to Calculation_set (description is ignored hereafter)
bool Console_log;              // true(default): output the calculation log on the console
bool Console_log_set;
char* Output_file;             // name of the output data file
bool Output_file_set;
/// Rg block
int Radial_grid_count;         // number of different radial grid configurations
double* Radial_grid_intervals; // interval of each radial grid
int* Radial_grid_points;       // number of points in each radial grid
double* x_coordinates;         // x coordinates of the radial grid
int x_count;                   // number of points in x_coordinates
/// TF block
bool TF_test;         // true: calculate TF potentials for some initial values
bool TF_test_set;
double Initial_diff_offset;    // test-mode: initial value for the first calculation
bool Initial_diff_offset_set;
double Initial_diff_delta;     // test-mode: delta for the initial value
bool Initial_diff_delta_set;
int Initial_diff_size;         // test-mode: how many points 
bool Initial_diff_size_set;
double Initial_diff_min;       // actual-mode: minimum initial value
bool Initial_diff_min_set;
double Initial_diff_max;       // actual-mode: maximum initial value
bool Initial_diff_max_set;
double TF_threshold;           // actual-mode: threshold for the convergence
bool TF_threshold_set;
char* TF_solution;             // solution: RK1 (or Euler) or RK4 (default)
double TF_solution_set;
double* TF_phi;                // Thomas-Fermi potential
double* TF_phi_diff;           // Differential of TF_phi;
/// At block
/// Either of (min and max) and (single) is allowed
int n_min;                     // minimum of principal quantum number (n)
bool n_min_set;
int n_max;                     // maximum of principal quantum number (n)
bool n_max_set;
int n_single;                  // the principal quantum number (n)
bool n_single_set;
int l_min;                     // azimutal quantum number (l)
bool l_min_set;
int l_max;
bool l_max_set;
int l_single;
bool l_single_set;
int Z_min;                     // atomic number (Z)
bool Z_min_set;
int Z_max;
bool Z_max_set;
int Z_single;
bool Z_single_set;
char* At_potential;            // Potential: H-like (-Z/r), Thomas-Fermi (-Z/r phi(x)), file
bool At_potential_set;
char* At_potential_file;       // Potential_file: necessary if Potential is Thomas-Fermi or file
bool At_potential_file_set;
char* At_solution;             // solution: RK1 (Euler), RK4, Numerov
bool At_solution_set;
double Bisection_step;         // initial step of the bisection method
bool Bisection_step_set;
double At_radius_factor;       // matching_radius*radius_factor determines the last point to be calculated
bool At_radius_factor_set;
double At_E_threshold;
bool At_E_threshold_set;
double At_initial_diff;        // initial differntial of wfn p'(0) (fixed in initialize())
double At_initial_diff2;       // initial differential of wfn p'(inf)
double* At_p_x;                // scaled radial wave function
double* At_p_diff_x;           // differential of At_p_x;
double* At_v_x;                // scaled potential energy
double At_bisection_threshold;
int At_min_iteration;
int At_max_iteration;
/// SC block
double SC_mix_weight;          // mixing weight for the next iteration
bool SC_mix_weight_set;
double SC_criterion_a;         // SCF criterion alpha
bool SC_criterion_a_set;
double SC_criterion_b;         // SCF criterion beta
bool SC_criterion_b_set;
double* SC_sigma_x;            // number density, integrated in the angles
int SC_orbital_count;          // number of occupied orbitals
int* SC_n_list;                // list of n and l
int* SC_l_list;
double* SC_eigen_list;         // list of eigenvalues
double** SC_p_x;               // list of P(x)
/// Oc block
int Occupation_count;
int** At_occupation;           // list of occupation numbers
/// Ex block
int Ex_energy_count;           // number of excitation energies
double* Ex_energies;           // excitation energies
/// Or block
int Ph_orbital_count;          // number of orbital to calculate phase shift
char** Ph_orbital_labels;      // orbital label (like '1s')
int* Ph_l_list;                // list of l
double* Ph_binding_energies;   // binding energies of atomic states
/// Ph block
int Ph_skip_points;            // number of skipped points in the phase analysis
bool Ph_skip_points_set;
int Ph_calc_points;            // number of used point in the phase analysis
bool Ph_calc_points_set;
/// PA block
char* PA_input_file;           // input file (.hdf5), output of the postproc.o
bool PA_input_file_set;
double PA_E_min;               // energy to start the calculation of the dispersion
bool PA_E_min_set;
double PA_E_max;               // energy to end the calculation of the dispersion
bool PA_E_max_set;
double PA_dE;                  // width of Gauss broadening of the dispersion
bool PA_dE_set;
double PA_E_pixel;             // delta E of the calculation of the dispersion
bool PA_E_pixel_set;
char* PA_initial_state;        // initial state ("PAO" or "AO")
bool PA_initial_state_set;
char* PA_final_state;          // final state ("PW", "Calc", "FP_PAO", or "FP_AO")
bool PA_final_state_set;
double PA_final_state_step;    // calculation step of final states (Bohr^-1)
bool PA_final_state_step_set;
char* PA_polarization;         // polarization ("Linear", "RCircular", or "LCircular")
bool PA_polarization_set;
double PA_theta;               // polarization angle
double PA_phi;
bool PA_theta_set;
bool PA_phi_set;
int PA_ext_up;                 // extend length to the up, right, down, left
int PA_ext_ri;
int PA_ext_dn;
int PA_ext_le;
bool PA_ext_set;
char* PA_AO_file;              // PAO and AO file
bool PA_AO_file_set;
bool PA_weighting;             // weighting for the surface states calculation
bool PA_weighting_set;
double* PA_weighting_axis;     // axis perpendicular to the surface
bool PA_weighting_axis_set;
char* PA_weighting_shape;      // weighting shape ("Rect" rectangular, "Exp" exponential, "Tri" triangular, "Sqrt" sqrt of triangle)
bool PA_weighting_shape_set;
double PA_weighting_origin;    // the origin is on weighting axis, with distance of this value from (0, 0, 0), negative value means the origin is on the opposite direction from the weighting axis
bool PA_weighting_origin_set;
double PA_weighting_width;     // width of the rectangle or the decrease coefficient lambda, negative value means the inverted weighting function
bool PA_weighting_width_set;
bool PA_use_angstrom;          // use angstrom as the unit of Weighting_axis, Weighting_origin, and Weighting_width
bool PA_use_angstrom_set;
char* PA_output_data;          // what to output: Band, PAD (default), PAD_real, PAD_imag
bool PA_output_data_set;
double PA_sigma_max;           // max sigma of the envelope Gauss function
double PA_decay_max;           // max sigma of the exponential weighting function
bool PA_reflection;            // reflection correction
bool PA_reflection_set;
double PA_reflection_coef;     // coefficient for reflection correction
bool PA_reflection_coef_set;
bool PA_include_neg_depth;     // whether the atoms with negative depth (above surface) are included in the intensity calculations
bool PA_include_neg_depth_set; 
double PA_excitation_energy;   // Excitation energy for first-principles final states (FPFSs) in unit of eV
bool PA_excitation_energy_set;
double PA_FPFS_energy_step;    // energy step to calculate FPFSs
bool PA_FPFS_energy_step_set;
int PA_FPFS_range;             // index range used to determine the in-plane search area
double PA_FPFS_kRange;         // determine the k space included in the FPFS calculation with !Numerov
bool PA_FPFS_kRange_set;
bool PA_FPFS;                  // whether final_state is (FP_PAO or FP_AO)
char* PA_FPFS_file;            // export destinaion
bool PA_FPFS_file_set;
char* PA_VPS_file;             // VPS database file for FPFS
bool PA_VPS_file_set;
int PA_Lebedev_order_ave;  // Lebedev order for averaging
int PA_Lebedev_order_int;  // Lebedev order for integration
int PA_lp_max; // max value of l^\prime in FPFS calculations
double PA_zero_threshold; // threshold to judge the value is approximately equal to zero
int PA_theta_points;      // number of points in the theta integration
bool PA_ignore_core;      // ignore the core (<rc) in the radial integration
bool PA_ignore_core_set;
bool PA_add_nonorth_term; // include the \tau_i \cdot e < psi_F | psi_I> term
bool PA_add_nonorth_term_set;
bool PA_calc_all_nonloc;     // calculate the nonlocal wave functions for all atoms which can have zero weights
bool PA_calc_all_nonloc_set;
//bool PA_calc_all_loc;        // calculate the local wave functions for the whole range
//bool PA_calc_all_loc_set;
bool PA_orth_correction;  // apply the orthgonality correction
bool PA_orth_correction_set;
bool PA_ignore_nonlocal;  // ignore the nonlocal part of the pseudopotential
bool PA_ignore_nonlocal_set;
bool PA_FPFS_bulk_include_nonlocal; // bulk band calculations include the nonlocal part
bool PA_FPFS_bulk_include_nonlocal_set; 
bool PA_FPFS_Numerov;  // perform the nonlocal wfn calculation using the Numerov method
bool PA_FPFS_Numerov_set;
double PA_FPFS_bulk_volume;
//double PA_FPFS_margin;   // margin to determine the calculation range of local wfn if PA_calc_all_wfn=false
//bool PA_FPFS_edge_smoothing; // perform the box smoothing of the edge values
//bool PA_FPFS_edge_smoothing_set;
//int PA_FPFS_smoothing_E; // smoothing width (half-side)
//bool PA_FPFS_smoothing_E_set;
//int PA_FPFS_smoothing_k; // smoothing width (half-side)
//bool PA_FPFS_smoothing_k_set;
bool PA_interpolate_wfn; // reduce the number of points in radial wave functions
bool PA_interpolate_wfn_set;
double PA_interpolate_wfn_coef; // dz*coef becomes delta r
bool PA_interpolate_wfn_coef_set;
double PA_FPFS_bulk_min_ang; // bulk range (Ang)
double PA_FPFS_bulk_max_ang; // 
int PA_FPFS_bulk_count;      // number of bulk layers in the bulk range
bool PA_FPFS_bulk_set;
int PA_FPFS_bulk_kz_steps; // steps to calculate bulk band dispersion along kz
int PA_FPFS_bulk_kappaz_steps_left; // 2pi/c/steps is the delta kappaz 
int PA_FPFS_bulk_kappaz_steps_right;
double PA_FPFS_nonloc_offset; // radial offset
bool PA_FPFS_nonloc_offset_set;

double PA_FPFS_negligible_gap_size;
bool PA_FPFS_negligible_gap_size_set;
int PA_FPFS_bulk_buffer_size;
//bool PA_FPFS_perturbation; // Apply perturbation method in surface calculations
//bool PA_FPFS_perturbation_set;
//bool PA_calc_complex_dispersion;
//bool PA_calc_complex_dispersion_set;
double PA_FPFS_real_eigenvalue_criterion;
bool PA_FPFS_real_eigenvalue_criterion_set;
double PA_FPFS_dispersion_group_criterion_left;
bool PA_FPFS_dispersion_group_criterion_left_set;
double PA_FPFS_dispersion_group_criterion_right;
bool PA_FPFS_dispersion_group_criterion_right_set;
double PA_FPFS_gap_upper_limit;
bool PA_FPFS_gap_upper_limit_set;
int PA_FPFS_kz_margin_index_size;
bool PA_FPFS_kz_margin_index_size_set;
double PA_FPFS_cspace_size;
bool PA_FPFS_cspace_size_set;
double PA_FPFS_kz_exclude_criterion;
bool PA_FPFS_kz_exclude_criterion_set;

/// files
FILE* Log_file_obj;            // log file object
FILE* Output_file_obj;         // data file object
