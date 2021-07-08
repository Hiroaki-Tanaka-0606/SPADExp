// calcPSF variables
#include <cstdio>

// constants
int buffer_length;        // size of the buffer used to load one line from the input file
int Calculation_length;   // size of the char* Calculation
int Log_file_length;      // size of the char* Log_file
int Log_length;           // max length of the one log message
char* Log_buffer;         // buffer for the log output

// variables
/// blocks
bool Ct_block_appeared;        // true: the block "Control" appeared, false: not
bool TF_block_appeared;        // similar to Ct_block_appeared (description is ignored hereafter)
bool Rg_block_appeared;
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
double* TF_phi;                // Thomas-Fermi potential
double* TF_phi_diff;           // Differential of TF_phi;

/// files
FILE* Log_file_obj;            // log file object
FILE* Output_file_obj;         // data file object
