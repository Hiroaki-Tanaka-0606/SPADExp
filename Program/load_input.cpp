// calcPSF load input

#include <string>
#include <iostream>
#include <cstring>

#include "variables_ext.hpp"
#include "parse_inputs.hpp"

using namespace std;

void output_error(int line_number, char* description){
	printf("Error in line %d: %s\n", line_number, description);
}

int load_input(){
	// buffer variables
	string input_line_s;
	char* input_line_c=new char[buffer_length+1];
	int actual_line_length;
	char* keyword_buffer=new char[buffer_length+1];
	int sscanf_status;
	int parse_status;
	char* value_buffer=new char[buffer_length+1];

	// status variables
	int status=1;            // 1: OK, 0: error
	int line_number=0;       // number of the line being loaded
	bool in_the_block=false; // true when the line is in the block, false when not
	string* block_name;      // name of the block which the line is being in
	int Radial_grid_index=0; // current index of the radial grid configuration

	while(getline(cin, input_line_s)){
		// cout << input_line_s << endl;

		//increase line number
		line_number++;
		
		// read one line from cin (stdin), copy to input_line_c
		actual_line_length=input_line_s.length();
		if(actual_line_length>buffer_length){
		  printf("Warning: line %d is longer than %d characters, the last part is ignored\n", line_number, buffer_length);
			actual_line_length=buffer_length;
		}
		input_line_s.copy(input_line_c, actual_line_length);
		input_line_c[actual_line_length]='\0';
		// cout << input_line_c << endl;

		// extract the first word
		sscanf_status=sscanf(input_line_c, "%s", keyword_buffer);
		/// sscanf_status<1 means nothing is in the line
		if(sscanf_status<1){
			printf("Line %d is a blank, ignoring\n", line_number);
			continue;
		}
		
		// distinguish the type of the line
		/// comment: 
		if(keyword_buffer[0]=='!' || keyword_buffer[0]== '#'){
			printf("Line %d is a comment, ignoring\n", line_number);
			continue;
		}
		/// start of the block
		if(keyword_buffer[0]=='&'){
			// check whether the line is already in a block or not
			if(in_the_block){
				output_error(line_number, (char*)"already in a block");
				status=0;
				goto FINALIZATION;
			}
			in_the_block=true;
			
			// extract the block name
			int keyword_buffer_length=strlen(keyword_buffer);
			block_name=new string(&keyword_buffer[1], &keyword_buffer[keyword_buffer_length]);

			// check whether the block name is valid and the block appears for the first time
			if(*block_name==string("Control")){
				// Ct: Control
				if(Ct_block_appeared){
					output_error(line_number, (char*)"block 'Control' already appeared"); status=0; goto FINALIZATION;
				}
				Ct_block_appeared=true;
			}else if(*block_name==string("Thomas-Fermi")){
				// TF: Thomas-Fermi
				if(TF_block_appeared){
					output_error(line_number, (char*)"block 'Thomas-Fermi' already appeared"); status=0; goto FINALIZATION;
				}
				TF_block_appeared=true;
			}else if(*block_name==string("Atomic-wfn")){
				// At: Atomic-wfn (wavefunction)
				if(At_block_appeared){
					output_error(line_number, (char*)"block 'Atomic-wfn' already appeared"); status=0; goto FINALIZATION;
				}
				At_block_appeared=true;
			}else if(*block_name==string("Radial-grid")){
				// Rg: Radial grid
				if(Rg_block_appeared){
					output_error(line_number, (char*)"block 'Radial-grid' already appeared"); status=0; goto FINALIZATION;
				}
				Rg_block_appeared=true;
				/// read the number of rows
				sscanf_status=sscanf(input_line_c, "%*s %d", &Radial_grid_count);
				if(sscanf_status<1 || Radial_grid_count<1){
					output_error(line_number, (char*)"invalid number of radial grids"); status=0; goto FINALIZATION;
				}
				Radial_grid_intervals=new double[Radial_grid_count];
				Radial_grid_points=new int[Radial_grid_count];
			}else{
				// none of the above
				output_error(line_number, (char*)"invalid block name"); status=0; goto FINALIZATION;
			}
			printf("Block %s started at line %d\n", block_name->c_str(), line_number);
			continue;
		}
		/// end of the block
		if(keyword_buffer[0]=='/' && strlen(keyword_buffer)==1){
			printf("Block %s ended at line %d\n", block_name->c_str(), line_number);
			in_the_block=false;
			if(*block_name==string("Radial-grid") && Radial_grid_index!=Radial_grid_count){
				output_error(line_number, (char*)"inequal number of radial grids"); status=0; goto FINALIZATION;
			}
			continue;
		}
		/// something invalid, out of the block
		if(!in_the_block){
			output_error(line_number, (char*)"something invalid appears out of the block");
			status=0;
			goto FINALIZATION;
		}
		
		/// special block
		if(*block_name==string("Radial-grid")){
			// Rg: (interval) (points)
			if(Radial_grid_index>=Radial_grid_count){
				output_error(line_number, (char*)"too much configurations in the Radial-grid block"); status=0; goto FINALIZATION;
			}
			sscanf_status=sscanf(input_line_c,"%lf %d", &Radial_grid_intervals[Radial_grid_index], &Radial_grid_points[Radial_grid_index]);
			if(sscanf_status<2 || Radial_grid_intervals[Radial_grid_index]<0 || Radial_grid_points[Radial_grid_index]<=0){
				output_error(line_number, (char*)"invalid value(s) of a radial grid"); status=0; goto FINALIZATION;
			}
			Radial_grid_index++;
			continue;
		}

		/// usual block, "(keyword) (value)"
		/* template
			/// Variable: type
			if(strcmp(keyword_buffer, "Variable")==0){
				if(Variable_set){
					output_error(line_number, (char*)"keyword Variable already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_type(input_line_c, (char*)"Variable", (arguments,) value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Variable"); status=0; goto FINALIZATION;
				}
			  Variable_set=true; continue;
			}*/
			
		if(*block_name==string("Control")){
			// Ct block
			/// Calculation: char[]
			if(strcmp(keyword_buffer, "Calculation")==0){
				if(Calculation_set){
					output_error(line_number, (char*)"keyword Calculation already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Calculation", Calculation, Calculation_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calculation"); status=0; goto FINALIZATION;
				}
				Calculation_set=true;
				continue;
			}
			/// Log_file: char[]
			if(strcmp(keyword_buffer, "Log_file")==0){
				if(Log_file_set){
					output_error(line_number, (char*)"keyword Log_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Log_file", Log_file, Log_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Log_file"); status=0; goto FINALIZATION;
				}
				Log_file_set=true;
				continue;
			}
			/// Console_log: bool
			if(strcmp(keyword_buffer, "Console_log")==0){
				if(Console_log_set){
					output_error(line_number, (char*)"keyword Console_log already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &Console_log, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Console_log"); status=0; goto FINALIZATION;
				}
				Console_log_set=true;
				continue;
			}
			/// Output_file: char[]
			if(strcmp(keyword_buffer, "Output_file")==0){
				if(Output_file_set){
					output_error(line_number, (char*)"keyword Output_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Output_file", Output_file, Log_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Output_file"); status=0; goto FINALIZATION;
				}
				Output_file_set=true;
				continue;
			}
		}else if(*block_name==string("Thomas-Fermi")){
			// TF block
			/// Calculation_test (TF_test): bool
			if(strcmp(keyword_buffer, "Calculation_test")==0){
				if(TF_test_set){
					output_error(line_number, (char*)"keyword Calculation_test already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &TF_test, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calculation_test"); status=0; goto FINALIZATION;
				}
			  TF_test_set=true; continue;
			}
			/// Initial_diff_offset: double
			if(strcmp(keyword_buffer, "Initial_diff_offset")==0){
				if(Initial_diff_offset_set){
					output_error(line_number, (char*)"keyword Initial_diff_offset already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &Initial_diff_offset);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_diff_offset"); status=0; goto FINALIZATION;
				}
			  Initial_diff_offset_set=true; continue;
			}
			/// Initial_diff_delta: double
			if(strcmp(keyword_buffer, "Initial_diff_delta")==0){
				if(Initial_diff_delta_set){
					output_error(line_number, (char*)"keyword Initial_diff_delta already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &Initial_diff_delta);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_diff_delta"); status=0; goto FINALIZATION;
				}
			  Initial_diff_delta_set=true; continue;
			}
			/// Initial_diff_size: int
			if(strcmp(keyword_buffer, "Initial_diff_size")==0){
				if(Initial_diff_size_set){
					output_error(line_number, (char*)"keyword Initial_diff_size already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Initial_diff_size);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_diff_size"); status=0; goto FINALIZATION;
				}
			  Initial_diff_size_set=true; continue;
			}
			/// Initial_diff_min: double
			if(strcmp(keyword_buffer, "Initial_diff_min")==0){
				if(Initial_diff_min_set){
					output_error(line_number, (char*)"keyword Initial_diff_min already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &Initial_diff_min);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_diff_min"); status=0; goto FINALIZATION;
				}
			  Initial_diff_min_set=true; continue;
			}
			/// Initial_diff_max: double
			if(strcmp(keyword_buffer, "Initial_diff_max")==0){
				if(Initial_diff_max_set){
					output_error(line_number, (char*)"keyword Initial_diff_max already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &Initial_diff_max);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_diff_max"); status=0; goto FINALIZATION;
				}
			  Initial_diff_max_set=true; continue;
			}
			/// Threshold(TF_threshold): double
			if(strcmp(keyword_buffer, "Threshold")==0){
				if(TF_threshold_set){
					output_error(line_number, (char*)"keyword Threshold already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &TF_threshold);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Threshold"); status=0; goto FINALIZATION;
				}
			  TF_threshold_set=true; continue;
			}
			/// Solution: char[]
			if(strcmp(keyword_buffer, "Solution")==0){
				if(TF_solution_set){
					output_error(line_number, (char*)"keyword Solution already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Solution", TF_solution, Solution_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Solution"); status=0; goto FINALIZATION;
				}
				TF_solution_set=true;
				continue;
			}
		}else if(*block_name==string("Atomic-wfn")){
			// At block
			/// n_min: int
			if(strcmp(keyword_buffer, "n_min")==0){
				if(n_min_set){
					output_error(line_number, (char*)"keyword n_min already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &n_min);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of n_min"); status=0; goto FINALIZATION;
				}
			  n_min_set=true; continue;
			}
			/// n_max: int
			if(strcmp(keyword_buffer, "n_max")==0){
				if(n_max_set){
					output_error(line_number, (char*)"keyword n_max already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &n_max);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of n_max"); status=0; goto FINALIZATION;
				}
			  n_max_set=true; continue;
			}
			/// n (n_single): int
			if(strcmp(keyword_buffer, "n")==0){
				if(n_single_set){
					output_error(line_number, (char*)"keyword n already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &n_single);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of n"); status=0; goto FINALIZATION;
				}
			  n_single_set=true; continue;
			}
			/// l_min: int
			if(strcmp(keyword_buffer, "l_min")==0){
				if(l_min_set){
					output_error(line_number, (char*)"keyword l_min already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &l_min);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of l_min"); status=0; goto FINALIZATION;
				}
			  l_min_set=true; continue;
			}
			/// l_max: int
			if(strcmp(keyword_buffer, "l_max")==0){
				if(l_max_set){
					output_error(line_number, (char*)"keyword l_max already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &l_max);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of l_max"); status=0; goto FINALIZATION;
				}
			  l_max_set=true; continue;
			}
			/// l (l_single): int
			if(strcmp(keyword_buffer, "l")==0){
				if(l_single_set){
					output_error(line_number, (char*)"keyword l already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &l_single);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of l"); status=0; goto FINALIZATION;
				}
			  l_single_set=true; continue;
			}
			/// Z_min: int
			if(strcmp(keyword_buffer, "Z_min")==0){
				if(Z_min_set){
					output_error(line_number, (char*)"keyword Z_min already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Z_min);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Z_min"); status=0; goto FINALIZATION;
				}
			  Z_min_set=true; continue;
			}
			/// Z_max: int
			if(strcmp(keyword_buffer, "Z_max")==0){
				if(Z_max_set){
					output_error(line_number, (char*)"keyword Z_max already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Z_max);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Z_max"); status=0; goto FINALIZATION;
				}
			  Z_max_set=true; continue;
			}
			/// Z (Z_single): int
			if(strcmp(keyword_buffer, "Z")==0){
				if(Z_single_set){
					output_error(line_number, (char*)"keyword Z already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Z_single);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Z"); status=0; goto FINALIZATION;
				}
			  Z_single_set=true; continue;
			}
			/// Potential (At_potential): char[]
			if(strcmp(keyword_buffer, "Potential")==0){
				if(At_potential_set){
					output_error(line_number, (char*)"keyword Potential already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Potential", At_potential, Potential_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Potential"); status=0; goto FINALIZATION;
				}
				At_potential_set=true;
				continue;
			}
			/// Potential_file (At_potential_file): char[]
			if(strcmp(keyword_buffer, "Potential_file")==0){
				if(At_potential_file_set){
					output_error(line_number, (char*)"keyword Potential_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Potential_file", At_potential_file, Potential_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Potential_file"); status=0; goto FINALIZATION;
				}
				At_potential_file_set=true;
				continue;
			}
			
			/// Solution (At_solution): char[]
			if(strcmp(keyword_buffer, "Solution")==0){
				if(At_solution_set){
					output_error(line_number, (char*)"keyword Solution already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Solution", At_solution, Solution_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Solution"); status=0; goto FINALIZATION;
				}
				At_solution_set=true;
				continue;
			}
			/// Bisection_step: double
			if(strcmp(keyword_buffer, "Bisection_step")==0){
				if(Bisection_step_set){
					output_error(line_number, (char*)"keyword Bisection_step already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &Bisection_step);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Bisection_step"); status=0; goto FINALIZATION;
				}
			  Bisection_step_set=true; continue;
			}
			/// E_threshold (At_E_threshold): double
			if(strcmp(keyword_buffer, "E_threshold")==0){
				if(At_E_threshold_set){
					output_error(line_number, (char*)"keyword E_threshold already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &At_E_threshold);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of E_threshold"); status=0; goto FINALIZATION;
				}
			  At_E_threshold_set=true; continue;
			}
			/// Radius_factor (At_radius_factor): double
			if(strcmp(keyword_buffer, "Radius_factor")==0){
				if(At_radius_factor_set){
					output_error(line_number, (char*)"keyword Radius_factor already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &At_radius_factor);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Radius_factor"); status=0; goto FINALIZATION;
				}
			  At_radius_factor_set=true; continue;
			}
		}else{
			// none of the above, which should not happen
			output_error(line_number, (char*)"unexpected block name error");
			status=0;
			goto FINALIZATION;
		}
		output_error(line_number, (char*)"invalid keyword");
		status=0;
		goto FINALIZATION;
	}

 FINALIZATION:
	delete input_line_c;
	delete keyword_buffer;
	return status;
}
