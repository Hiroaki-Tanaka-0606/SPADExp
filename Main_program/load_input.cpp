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
	int Oc_orbital_index=0;  // current index of the occupation configuration
	int Ex_energy_index=0;   // current index of the excitation energy configuration
	int Ph_orbital_index=0;  // current index of the orbital configuration

	int i;

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
			}else if(*block_name==string("SCF-atom")){
				// SC: SCF-atom
				if(SC_block_appeared){
					output_error(line_number, (char*)"block 'SCF-atom' already appeared"); status=0; goto FINALIZATION;
				}
				SC_block_appeared=true;
			}else if(*block_name==string("Occupation")){
				// Oc: Occupation of an atom
				if(Oc_block_appeared){
					output_error(line_number, (char*)"block 'Occupation' already appeared"); status=0; goto FINALIZATION;
				}
				Oc_block_appeared=true;
				/// read the number of rows
				sscanf_status=sscanf(input_line_c, "%*s %d", &Occupation_count);
				if(sscanf_status<1 || Occupation_count<1){
					output_error(line_number, (char*)"invalid number of Occupation_count"); status=0; goto FINALIZATION;
				}
				int* At_occupation_alloc=new int[Occupation_count*Occupation_count];
				for(int i=0; i<Occupation_count*Occupation_count; i++){
					At_occupation_alloc[i]=0;
				}
			  At_occupation=new int*[Occupation_count];
				for(int i=0; i<Occupation_count; i++){
					At_occupation[i]=&At_occupation_alloc[i*Occupation_count];
				}
			}else if(*block_name==string("Excitation-energy")){
				// Ex: Excitation energies
				if(Ex_block_appeared){
					output_error(line_number, (char*)"block 'Excitation-energy' already appeared"); status=0; goto FINALIZATION;
				}
				Ex_block_appeared=true;
				/// read the number of rows
				sscanf_status=sscanf(input_line_c, "%*s %d", &Ex_energy_count);
				if(sscanf_status<1 || Ex_energy_count<1){
					output_error(line_number, (char*)"invalid number of Ex_energy_count"); status=0; goto FINALIZATION;
				}
				Ex_energies=new double[Ex_energy_count];
			}else if(*block_name==string("Orbital")){
				// Or: orbitals to calculate phase shift
				if(Or_block_appeared){
					output_error(line_number, (char*)"block 'Orbital' already appeared"); status=0; goto FINALIZATION;
				}
				Or_block_appeared=true;
				/// read the number of rows
				sscanf_status=sscanf(input_line_c, "%*s %d", &Ph_orbital_count);
				if(sscanf_status<1 || Ph_orbital_count<1){
					output_error(line_number, (char*)"invalid number of Ph_orbital_count"); status=0; goto FINALIZATION;
				}
				Ph_orbital_labels=new char*[Ph_orbital_count];
				for(i=0; i<Ph_orbital_count; i++){
					Ph_orbital_labels[i]=new char[Ph_label_length+1];
				}
				Ph_l_list=new int[Ph_orbital_count];
				Ph_binding_energies=new double[Ph_orbital_count];
			}else if(*block_name==string("Phase-shift")){
				// Ph: Phase-shift
				if(Ph_block_appeared){
					output_error(line_number, (char*)"block 'Phase-shift' already appeared"); status=0; goto FINALIZATION;
				}
				Ph_block_appeared=true;
			}else if(*block_name==string("PAD")){
				// PA: Photoemission angular distribution (PAD)
				if(PA_block_appeared){
					output_error(line_number, (char*)"block 'PAD' already appeared"); status=0; goto FINALIZATION;
				}
				PA_block_appeared=true;
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
			if(*block_name==string("Occupation") && Oc_orbital_index!=Occupation_count){
				output_error(line_number, (char*)"inequal number of orbital indices"); status=0; goto FINALIZATION;
			}
			if(*block_name==string("Excitation-energy") && Ex_energy_index!=Ex_energy_count){
				output_error(line_number, (char*)"inequal number of excitation energies"); status=0; goto FINALIZATION;
			}
			if(*block_name==string("Orbital") && Ph_orbital_index!=Ph_orbital_count){
				output_error(line_number, (char*)"inequal number of orbitals"); status=0; goto FINALIZATION;
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
		if(*block_name==string("Occupation")){
			// Oc: (s orbital) (p orbital) (d orbtital) ...
			if(Oc_orbital_index>=Occupation_count){
				output_error(line_number, (char*)"too much configurations in the Occupation block"); status=0; goto FINALIZATION;
			}
			int number_of_targets=Oc_orbital_index+1;
			char* sscanf_template=new char[buffer_length+1];
			int i, j;
			for(i=0; i<number_of_targets; i++){
				sscanf_template[0]='\0';
				for(j=0; j<=i-1; j++){
					sprintf(sscanf_template, "%s %%*d", sscanf_template);
				}
				sprintf(sscanf_template, "%s %%d", sscanf_template);
				sscanf_status=sscanf(input_line_c, sscanf_template, &At_occupation[Oc_orbital_index][i]);
			
				if(sscanf_status<1 || At_occupation[Oc_orbital_index][i]<0){
					output_error(line_number, (char*)"invalid value(s) of a occupation number"); status=0; goto FINALIZATION;
				}
			}
			Oc_orbital_index++;
			delete[] sscanf_template;
			continue;
		}
		if(*block_name==string("Excitation-energy")){
			// Ex: (energy (eV))
			// Ex_energies[i] is in Eh
			if(Ex_energy_index>=Ex_energy_count){
				output_error(line_number, (char*)"too much configurations in the Excitation-energy block"); status=0; goto FINALIZATION;
			}
			double energy_eV;
			sscanf_status=sscanf(input_line_c,"%lf", &energy_eV);
			if(sscanf_status<1){
				output_error(line_number, (char*)"invalid value of an excitation energy"); status=0; goto FINALIZATION;
			}
			Ex_energies[Ex_energy_index]=energy_eV/Eh;
			Ex_energy_index++;
			continue;
		}
		if(*block_name==string("Orbital")){
			// Or: (label) (binding energy (eV)
			// Ph_binding_energies[i] is in Eh
			if(Ph_orbital_index>=Ph_orbital_count){
				output_error(line_number, (char*)"too much configurations in the Orbital block"); status=0; goto FINALIZATION;
			}
			double energy_eV;
			sscanf_status=sscanf(input_line_c,"%s %lf", Ph_orbital_labels[Ph_orbital_index], &energy_eV);
			if(sscanf_status<2){
				output_error(line_number, (char*)"invalid value(s) of an orbital"); status=0; goto FINALIZATION;
			}
			Ph_binding_energies[Ph_orbital_index]=energy_eV/Eh;
			if(strlen(Ph_orbital_labels[Ph_orbital_index])>Ph_label_length){
				output_error(line_number, (char*)"Error: too long orbital label"); status=0; goto FINALIZATION;
			}
			if(Ph_orbital_labels[Ph_orbital_index][1]=='s'){
				Ph_l_list[Ph_orbital_index]=0;
			}else if(Ph_orbital_labels[Ph_orbital_index][1]=='p'){
				Ph_l_list[Ph_orbital_index]=1;
			}else if(Ph_orbital_labels[Ph_orbital_index][1]=='d'){
				Ph_l_list[Ph_orbital_index]=2;
			}else if(Ph_orbital_labels[Ph_orbital_index][1]=='f'){
				Ph_l_list[Ph_orbital_index]=3;
			}else{
				output_error(line_number, (char*)"Error: invalid orbital"); status=0; goto FINALIZATION;
			}
			Ph_orbital_index++;
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
			
		}else if(*block_name==string("SCF-atom")){
			// SC block
			/// Mix_weight (SC_mix_weight): double
			if(strcmp(keyword_buffer, "Mix_weight")==0){
				if(SC_mix_weight_set){
					output_error(line_number, (char*)"keyword Mix_weight already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &SC_mix_weight);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Mix_weight"); status=0; goto FINALIZATION;
				}
			  SC_mix_weight_set=true; continue;
			}
			/// Criterion_a (SC_criterion_a): double
			if(strcmp(keyword_buffer, "Criterion_a")==0){
				if(SC_criterion_a_set){
					output_error(line_number, (char*)"keyword Criterion_a already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &SC_criterion_a);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of SC_criterion_a"); status=0; goto FINALIZATION;
				}
			  SC_criterion_a_set=true; continue;
			}
			/// Criterion_b (SC_criterion_b): double
			if(strcmp(keyword_buffer, "Criterion_b")==0){
				if(SC_criterion_b_set){
					output_error(line_number, (char*)"keyword Criterion_b already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &SC_criterion_b);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of SC_criterion_b"); status=0; goto FINALIZATION;
				}
			  SC_criterion_b_set=true; continue;
			}
		}else if(*block_name==string("PAD")){
			// PA block
			/// Input_file (PA_input_file): char[]
			if(strcmp(keyword_buffer, "Input_file")==0){
				if(PA_input_file_set){
					output_error(line_number, (char*)"keyword Input_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Input_file", PA_input_file, HDF5_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Input_file"); status=0; goto FINALIZATION;
				}
				PA_input_file_set=true; continue;
			}
			/// E_min (PA_E_min): double
			if(strcmp(keyword_buffer, "E_min")==0){
				if(PA_E_min_set){
					output_error(line_number, (char*)"keyword E_min already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_E_min);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of E_min"); status=0; goto FINALIZATION;
				}
			  PA_E_min_set=true; continue;
			}
			/// E_max (PA_E_max): double
			if(strcmp(keyword_buffer, "E_max")==0){
				if(PA_E_max_set){
					output_error(line_number, (char*)"keyword E_max already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_E_max);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of E_max"); status=0; goto FINALIZATION;
				}
			  PA_E_max_set=true; continue;
			}
			/// dE (PA_dE): double
			if(strcmp(keyword_buffer, "dE")==0){
				if(PA_dE_set){
					output_error(line_number, (char*)"keyword dE already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_dE);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of dE"); status=0; goto FINALIZATION;
				}
			  PA_dE_set=true; continue;
			}
			/// E_pixel (PA_E_pixel): double
			if(strcmp(keyword_buffer, "E_pixel")==0){
				if(PA_E_pixel_set){
					output_error(line_number, (char*)"keyword E_pixel already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_E_pixel);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of E_pixel"); status=0; goto FINALIZATION;
				}
			  PA_E_pixel_set=true; continue;
			}
			/// Initial_state (PA_initial_state): char[]
			if(strcmp(keyword_buffer, "Initial_state")==0){
				if(PA_initial_state_set){
					output_error(line_number, (char*)"keyword Initial_state already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Initial_state", PA_initial_state, PA_state_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Initial_state"); status=0; goto FINALIZATION;
				}
				PA_initial_state_set=true; continue;
			}
			/// Final_state (PA_final_state): char[]
			if(strcmp(keyword_buffer, "Final_state")==0){
				if(PA_final_state_set){
					output_error(line_number, (char*)"keyword Final_state already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Final_state", PA_final_state, PA_state_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Final_state"); status=0; goto FINALIZATION;
				}
				PA_final_state_set=true; continue;
			}
			/// Final_state_step (PA_final_state_step): double
			if(strcmp(keyword_buffer, "Final_state_step")==0){
				if(PA_final_state_step_set){
					output_error(line_number, (char*)"keyword Final_state_step already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_final_state_step);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Final_state_step"); status=0; goto FINALIZATION;
				}
			  PA_final_state_step_set=true; continue;
			}			
			/// Polarization (PA_polarization): char[]
			if(strcmp(keyword_buffer, "Polarization")==0){
				if(PA_polarization_set){
					output_error(line_number, (char*)"keyword Polarization already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Polarization", PA_polarization, PA_state_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Polarization"); status=0; goto FINALIZATION;
				}
				PA_polarization_set=true; continue;
			}
			/// Theta (PA_theta): double
			if(strcmp(keyword_buffer, "Theta")==0){
				if(PA_theta_set){
					output_error(line_number, (char*)"keyword Theta already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_theta);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Theta"); status=0; goto FINALIZATION;
				}
			  PA_theta_set=true; continue;
			}
			/// Phi (PA_phi): double
			if(strcmp(keyword_buffer, "Phi")==0){
				if(PA_phi_set){
					output_error(line_number, (char*)"keyword Phi already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_phi);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Phi"); status=0; goto FINALIZATION;
				}
			  PA_phi_set=true; continue;
			}
			/// Atomic_orbitals_file (PA_AO_file): char[]
			if(strcmp(keyword_buffer, "Atomic_orbitals_file")==0){
				if(PA_AO_file_set){
					output_error(line_number, (char*)"keyword Atomic_orbitals_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Atomic_orbitals_file", PA_AO_file, HDF5_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Atomic_orbitals_file"); status=0; goto FINALIZATION;
				}
				PA_AO_file_set=true; continue;
			}
			/// Extend (PA_ext_(up|ri|dn|le)): int[4]
			/// if dimension==1, only right and left are used
			if(strcmp(keyword_buffer, "Extend")==0){
				if(PA_ext_set){
					output_error(line_number, (char*)"keyword Extend already appeared"); status=0; goto FINALIZATION;
				}
				sscanf_status=sscanf(input_line_c, "%*s %d %d %d %d", &PA_ext_up, &PA_ext_ri, &PA_ext_dn, &PA_ext_le);
				if(sscanf_status!=4){
					output_error(line_number, (char*)"invalid value of Extend"); status=0; goto FINALIZATION;
				}
				PA_ext_set=true; continue;
			}
			/// Weighting (PA_weighting): bool
			if(strcmp(keyword_buffer, "Weighting")==0){
				if(PA_weighting_set){
					output_error(line_number, (char*)"keyword Weighting already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_weighting, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Weighting"); status=0; goto FINALIZATION;
				}
			  PA_weighting_set=true; continue;
			}
			/// Weighting_axis (PA_weighting_axis): double[3]
			if(strcmp(keyword_buffer, "Weighting_axis")==0){
				if(PA_weighting_axis_set){
					output_error(line_number, (char*)"keyword Weighting_axis already appeared"); status=0; goto FINALIZATION;
				}
				sscanf_status=sscanf(input_line_c, "%*s %lf %lf %lf", &PA_weighting_axis[0], &PA_weighting_axis[1], &PA_weighting_axis[2]);
				if(sscanf_status!=3){
					output_error(line_number, (char*)"invalid value of Weighting_axis"); status=0; goto FINALIZATION;
				}
				PA_weighting_axis_set=true; continue;
			}
			/// Weighting_shape (PA_weighting_shape): char
			if(strcmp(keyword_buffer, "Weighting_shape")==0){
				if(PA_weighting_shape_set){
					output_error(line_number, (char*)"keyword Weighting_shape already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Weighting_shape", PA_weighting_shape, PA_state_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of PA_weighting_shape"); status=0; goto FINALIZATION;
				}
				PA_weighting_shape_set=true; continue;
			}
			/// Weighting_origin (PA_weighting_origin): double
			if(strcmp(keyword_buffer, "Weighting_origin")==0){
				if(PA_weighting_origin_set){
					output_error(line_number, (char*)"keyword Weighting_origin already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_weighting_origin);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Weighting_origin"); status=0; goto FINALIZATION;
				}
			  PA_weighting_origin_set=true; continue;
			}
			/// Weighting_width (PA_weighting_width): double
			if(strcmp(keyword_buffer, "Weighting_width")==0){
				if(PA_weighting_width_set){
					output_error(line_number, (char*)"keyword Weighting_width already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_weighting_width);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Weighting_width"); status=0; goto FINALIZATION;
				}
			  PA_weighting_width_set=true; continue;
			}
			/// Use_angstrom: bool
			if(strcmp(keyword_buffer, "Use_angstrom")==0){
				if(PA_use_angstrom_set){
					output_error(line_number, (char*)"keyword Use_angstrom already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_use_angstrom, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Use_angstrom"); status=0; goto FINALIZATION;
				}
			  PA_use_angstrom_set=true; continue;
			}
			/// Output_data: char
			if(strcmp(keyword_buffer, "Output_data")==0){
				if(PA_output_data_set){
					output_error(line_number, (char*)"keyword Output_data already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"Output_data", PA_output_data, PA_state_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Output_data"); status=0; goto FINALIZATION;
				}
				PA_output_data_set=true; continue;
			}
			/// Reflection: bool
			if(strcmp(keyword_buffer, "Reflection")==0){
				if(PA_reflection_set){
					output_error(line_number, (char*)"keyword Reflection already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_reflection, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Reflection"); status=0; goto FINALIZATION;
				}
			  PA_reflection_set=true; continue;
			}
			/// Reflection_coef: double
			if(strcmp(keyword_buffer, "Reflection_coef")==0){
				if(PA_reflection_coef_set){
					output_error(line_number, (char*)"keyword Reflection_coef already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_reflection_coef);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Reflection_coef"); status=0; goto FINALIZATION;
				}
			  PA_reflection_coef_set=true; continue;
			}
			/// Include_neg_depth: bool
			if(strcmp(keyword_buffer, "Include_neg_depth")==0){
				if(PA_include_neg_depth_set){
					output_error(line_number, (char*)"keyword Include_neg_depth already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_include_neg_depth, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Include_neg_depth"); status=0; goto FINALIZATION;
				}
				PA_include_neg_depth_set=true; continue;
			}
			/// Excitation_energy: double
			if(strcmp(keyword_buffer, "Excitation_energy")==0){
				if(PA_excitation_energy_set){
					output_error(line_number, (char*)"keyword Excitation_energy already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_excitation_energy);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Excitation_energy"); status=0; goto FINALIZATION;
				}
			  PA_excitation_energy_set=true; continue;
			}
			/// FPFS_energy_step: double
			if(strcmp(keyword_buffer, "FPFS_energy_step")==0){
				if(PA_FPFS_energy_step_set){
					output_error(line_number, (char*)"keyword FPFS_energy_step already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_FPFS_energy_step);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_energy_step"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_energy_step_set=true; continue;
			}
			/// FPFS_kRange: double
			if(strcmp(keyword_buffer, "FPFS_kRange")==0){
				if(PA_FPFS_kRange_set){
					output_error(line_number, (char*)"keyword FPFS_kRange already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_FPFS_kRange);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_kRange"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_kRange_set=true; continue;
			}
			/// FPFS_file: char
			if(strcmp(keyword_buffer, "FPFS_file")==0){
				if(PA_FPFS_file_set){
					output_error(line_number, (char*)"keyword FPFS_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"FPFS_file", PA_FPFS_file, HDF5_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_file"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_file_set=true; continue;
			}
			/// VPS_file: char
			if(strcmp(keyword_buffer, "VPS_file")==0){
				if(PA_VPS_file_set){
					output_error(line_number, (char*)"keyword VPS_file already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_char(input_line_c, (char*)"VPS_file", PA_VPS_file, HDF5_file_length, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of VPS_file"); status=0; goto FINALIZATION;
				}
			  PA_VPS_file_set=true; continue;
			}
			/// Ignore_core: bool
			if(strcmp(keyword_buffer, "Ignore_core")==0){
				if(PA_ignore_core_set){
					output_error(line_number, (char*)"keyword Ignore_core already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_ignore_core, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Ignore_core"); status=0; goto FINALIZATION;
				}
			  PA_ignore_core_set=true; continue;
			}
			/// Add_nonorth_term: bool
			if(strcmp(keyword_buffer, "Add_nonorth_term")==0){
				if(PA_add_nonorth_term_set){
					output_error(line_number, (char*)"keyword Add_nonorth_term already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_add_nonorth_term, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Add_nonorth_term"); status=0; goto FINALIZATION;
				}
			  PA_add_nonorth_term_set=true; continue;
			}
			/// Calc_all_nonloc: bool
			if(strcmp(keyword_buffer, "Calc_all_nonloc")==0){
				if(PA_calc_all_nonloc_set){
					output_error(line_number, (char*)"keyword Calc_all_nonloc already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_calc_all_nonloc, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calc_all_nonloc"); status=0; goto FINALIZATION;
				}
			  PA_calc_all_nonloc_set=true; continue;
			}
			/// Calc_all_loc: bool
			if(strcmp(keyword_buffer, "Calc_all_loc")==0){
				if(PA_calc_all_loc_set){
					output_error(line_number, (char*)"keyword Calc_all_loc already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_calc_all_loc, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calc_all_loc"); status=0; goto FINALIZATION;
				}
			  PA_calc_all_loc_set=true; continue;
			}
			/// Orth_correction: bool
			if(strcmp(keyword_buffer, "Orth_correction")==0){
				if(PA_orth_correction_set){
					output_error(line_number, (char*)"keyword Orth_correction already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_orth_correction, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Orth_correction"); status=0; goto FINALIZATION;
				}
			  PA_orth_correction_set=true; continue;
			}
			/// Ignore_nonlocal: bool
			if(strcmp(keyword_buffer, "Ignore_nonlocal")==0){
				if(PA_ignore_nonlocal_set){
					output_error(line_number, (char*)"keyword Ignore_nonlocal already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_ignore_nonlocal, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Ignore_nonlocal"); status=0; goto FINALIZATION;
				}
			  PA_ignore_nonlocal_set=true; continue;
			}
			/// FPFS_Numerov: bool
			if(strcmp(keyword_buffer, "FPFS_Numerov")==0){
				if(PA_FPFS_Numerov_set){
					output_error(line_number, (char*)"keyword FPFS_Numerov already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_FPFS_Numerov, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_Numerov"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_Numerov_set=true; continue;
			}
			/// FPFS_edge_smoothing: bool
			if(strcmp(keyword_buffer, "FPFS_edge_smoothing")==0){
				if(PA_FPFS_edge_smoothing_set){
					output_error(line_number, (char*)"keyword FPFS_edge_smoothing already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_FPFS_edge_smoothing, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_edge_smoothing"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_edge_smoothing_set=true; continue;
			}
			/// FPFS_smoothing_E: int
			if(strcmp(keyword_buffer, "FPFS_smoothing_E")==0){
				if(PA_FPFS_smoothing_E_set){
					output_error(line_number, (char*)"keyword FPFS_smoothing_E already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &PA_FPFS_smoothing_E);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_smoothing_E"); status=0; goto FINALIZATION;
				}
				PA_FPFS_smoothing_E_set=true; continue;
			}
			/// FPFS_smoothing_k: int
			if(strcmp(keyword_buffer, "FPFS_smoothing_k")==0){
				if(PA_FPFS_smoothing_k_set){
					output_error(line_number, (char*)"keyword FPFS_smoothing_k already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &PA_FPFS_smoothing_k);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_smoothing_k"); status=0; goto FINALIZATION;
				}
				PA_FPFS_smoothing_k_set=true; continue;
			}
			/// FPFS_bulk: double(min, ang) double(max, ang) int(count)
			if(strcmp(keyword_buffer, "FPFS_bulk")==0){
				if(PA_FPFS_bulk_set){
					output_error(line_number, (char*)"keyword FPFS_bulk already appeared"); status=0; goto FINALIZATION;
				}
				sscanf_status=sscanf(input_line_c, "%*s %lf %lf %d", &PA_FPFS_bulk_min_ang, &PA_FPFS_bulk_max_ang, &PA_FPFS_bulk_count);
				if(sscanf_status!=3){
					output_error(line_number, (char*)"invalid value of FPFS_bulk"); status=0; goto FINALIZATION;
				}
				PA_FPFS_bulk_set=true; continue;
			}
			/// FPFS_perturbation: bool
			if(strcmp(keyword_buffer, "FPFS_perturbation")==0){
				if(PA_FPFS_perturbation_set){
					output_error(line_number, (char*)"keyword FPFS_perturbation already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_FPFS_perturbation, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of FPFS_perturbation"); status=0; goto FINALIZATION;
				}
			  PA_FPFS_perturbation_set=true; continue;
			}
			/// Calc_complex_dispersion: bool
			if(strcmp(keyword_buffer, "Calc_complex_dispersion")==0){
				if(PA_calc_complex_dispersion_set){
					output_error(line_number, (char*)"keyword Calc_complex_dispersion already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_calc_complex_dispersion, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calc_complex_dispersion"); status=0; goto FINALIZATION;
				}
			  PA_calc_complex_dispersion_set=true; continue;
			}
			/// interpolate_wfn: bool
			if(strcmp(keyword_buffer, "Interpolate_wfn")==0){
				if(PA_interpolate_wfn_set){
					output_error(line_number, (char*)"keyword Interpolate_wfn already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_bool(input_line_c, &PA_interpolate_wfn, value_buffer);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Interpolate_wfn"); status=0; goto FINALIZATION;
				}
			  PA_interpolate_wfn_set=true; continue;
			}
			/// Interpolate_wfn_coef: double
			if(strcmp(keyword_buffer, "Interpolate_wfn_coef")==0){
				if(PA_interpolate_wfn_coef_set){
					output_error(line_number, (char*)"keyword Interpolate_wfn_coef already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_double(input_line_c, &PA_interpolate_wfn_coef);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Interpolate_wfn_coef"); status=0; goto FINALIZATION;
				}
				PA_interpolate_wfn_coef_set=true; continue;
			}
		}else if(*block_name==string("Phase-shift")){
			// Ph block
			/// Skip_points (Ph_skip_points): int
			if(strcmp(keyword_buffer, "Skip_points")==0){
				if(Ph_skip_points_set){
					output_error(line_number, (char*)"keyword Skip_points already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Ph_skip_points);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Skip_points"); status=0; goto FINALIZATION;
				}
				Ph_skip_points_set=true; continue;
			}
			/// Calc_points (Ph_calc_points): int
			if(strcmp(keyword_buffer, "Calc_points")==0){
				if(Ph_calc_points_set){
					output_error(line_number, (char*)"keyword Calc_points already appeared"); status=0; goto FINALIZATION;
				}
				parse_status=parse_int(input_line_c, &Ph_calc_points);
				if(parse_status==0){
					output_error(line_number, (char*)"invalid value of Calc_points"); status=0; goto FINALIZATION;
				}
				Ph_calc_points_set=true; continue;
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
	delete[] input_line_c;
	delete[] keyword_buffer;
	return status;
}
