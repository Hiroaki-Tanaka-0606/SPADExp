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
int HDF5_file_length;     // size of the char* PS_input_file, PS_AO_file
int PS_state_length;      // size of the char* PS_initial_state, PS_final_state
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
bool PS_block_appeared;
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
/// PS block
char* PS_input_file;           // input file (.hdf5), output of the postproc.o
bool PS_input_file_set;
double PS_E_min;               // energy to start the calculation of the dispersion
bool PS_E_min_set;
double PS_E_max;               // energy to end the calculation of the dispersion
bool PS_E_max_set;
double PS_dE;                  // width of Gauss broadening of the dispersion
bool PS_dE_set;
double PS_E_pixel;             // delta E of the calculation of the dispersion
bool PS_E_pixel_set;
char* PS_initial_state;        // initial state ("PAO" or "AO")
bool PS_initial_state_set;
char* PS_final_state;          // final state ("PW" or "Calc")
bool PS_final_state_set;
char* PS_polarization;         // polarization ("Linear", "RCircular", or "LCircular")
bool PS_polarization_set;
double PS_theta;               // polarization angle
double PS_phi;
bool PS_theta_set;
bool PS_phi_set;
int PS_ext_up;                 // extend length to the up, right, down, left
int PS_ext_ri;
int PS_ext_dn;
int PS_ext_le;
bool PS_ext_set;
char* PS_AO_file;              // PAO and AO file
bool PS_AO_file_set;
double PS_sigma_max;           // max sigma of the envelope Gauss function


/// files
FILE* Log_file_obj;            // log file object
FILE* Output_file_obj;         // data file object
