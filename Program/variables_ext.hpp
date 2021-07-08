// calcPSF variables declaration
// description are in variables.hpp
#include <cstdio>

// constants
extern int buffer_length;
extern int Calculation_length;
extern int Log_file_length;
extern int Log_length;
extern char* Log_buffer;

// variables
/// blocks
extern bool Ct_block_appeared;
extern bool TF_block_appeared;
extern bool Rg_block_appeared;
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
extern double* TF_phi;
extern double* TF_phi_diff;

/// files
extern FILE* Log_file_obj;
extern FILE* Output_file_obj;

