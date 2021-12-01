// calcPSF setup (initialization and finalization)

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "variables_ext.hpp"
#include "log.hpp"

using namespace std;

void initialize(){
	// constants
	buffer_length=1024;
	Calculation_length=64;
	Log_file_length=1024;
	Log_length=256;
	Log_buffer=new char[Log_length+1];
	Solution_length=16;
	Potential_length=32;
	Potential_file_length=1024;
	Data_read_error=0.01;
	HDF5_file_length=1024;
	PS_state_length=64;
	
	// variables
	/// blocks
	Ct_block_appeared=false;
	TF_block_appeared=false;
	At_block_appeared=false;
	/// Ct block
	Calculation=new char[Calculation_length+1];
	Calculation_set=false;
	Log_file=new char[Log_file_length+1];
	Log_file_set=false;
	Console_log=true;
	Console_log_set=false;
	Output_file=new char[Log_file_length+1];
	Output_file_set=false;
	/// TF block
  TF_test=false;
	TF_test_set=false;
	Initial_diff_offset_set=false;
	Initial_diff_delta_set=false;
	Initial_diff_size_set=false;
	Initial_diff_min=-1.49; // -1.52 for RK1
	Initial_diff_max=-1.48; // -1.51 for RK1
	Initial_diff_min_set=false;
	Initial_diff_max_set=false;
	TF_threshold=1e-5;
	TF_threshold_set=false;
	TF_solution=new char[Solution_length+1];
	strcpy(TF_solution, "RK4");
	TF_solution_set=false;
	/// At block
	n_min_set=false;
	n_max_set=false;
	n_single_set=false;
	l_min_set=false;
	l_max_set=false;
	l_single_set=false;
	Z_min_set=false;
	Z_max_set=false;
	Z_single_set=false;
	At_potential_set=false;
	At_potential=new char[Potential_length+1];
	At_potential_file_set=false;
	At_potential_file=new char[Potential_file_length+1];
	At_solution_set=false;
	At_solution=new char[Solution_length+1];
	strcpy(At_solution, "Numerov");
	At_radius_factor_set=false;
	At_radius_factor=8;
	At_E_threshold_set=false;
	At_E_threshold=1e-5;
	Bisection_step_set=false;
	Bisection_step=1e-3;
	At_initial_diff=1;
	At_initial_diff2=-1e-10;
	At_bisection_threshold=1.01;
	At_min_iteration=100;
	At_max_iteration=1000;
	/// PS block
	PS_input_file_set=false;
	PS_input_file=new char[HDF5_file_length+1];
	PS_E_min_set=false;
	PS_E_max_set=false;
	PS_dE_set=false;
	PS_E_pixel_set=false;
	PS_initial_state_set=false;
	PS_initial_state=new char[PS_state_length+1];
	PS_final_state_set=false;
	PS_final_state=new char[PS_state_length+1];
	PS_polarization_set=false;
	PS_polarization=new char[PS_state_length+1];
	PS_theta_set=false;
	PS_phi_set=false;
	PS_AO_file_set=false;
	PS_AO_file=new char[HDF5_file_length+1];		
}

void finalize(){
	delete Calculation;
	delete Log_file;
	delete Output_file;
	delete TF_solution;
	delete At_solution;
	if(Log_file_set && Log_file_obj!=NULL){
		fclose(Log_file_obj);
	}
	if(Output_file_set && Output_file_obj!=NULL){
		fclose(Output_file_obj);
	}
	if(Radial_grid_count>0){
		delete Radial_grid_intervals;
		delete Radial_grid_points;
	}
	delete PS_input_file;
}

void initialize_radial_grid(){
	Radial_grid_count=11;
	Radial_grid_intervals=new double[11];
	Radial_grid_intervals[0]=0.0025;
	Radial_grid_intervals[1]=0.005;
	Radial_grid_intervals[2]=0.01;
	Radial_grid_intervals[3]=0.02;
	Radial_grid_intervals[4]=0.04;
	Radial_grid_intervals[5]=0.08;
	Radial_grid_intervals[6]=0.16;
	Radial_grid_intervals[7]=0.32;
	Radial_grid_intervals[8]=0.64;
	Radial_grid_intervals[9]=1.28;
	Radial_grid_intervals[10]=2.56;
	Radial_grid_points=new int[11];
	Radial_grid_points[0]=40;
	Radial_grid_points[1]=40;
	Radial_grid_points[2]=40;
	Radial_grid_points[3]=40;
	Radial_grid_points[4]=40;
	Radial_grid_points[5]=40;
	Radial_grid_points[6]=40;
	Radial_grid_points[7]=40;
	Radial_grid_points[8]=40;
	Radial_grid_points[9]=40;
	Radial_grid_points[10]=40;
}

void generate_radial_grid(){
	int i, j, index;
	int points_count=1;
	for(i=0; i<Radial_grid_count; i++){
		points_count+=Radial_grid_points[i];
	}
	x_coordinates=new double[points_count];
	x_coordinates[0]=0;
	index=1;
	double x_current=0;
	for(i=0; i<Radial_grid_count; i++){
		for(j=0; j<Radial_grid_points[i]; j++){
			x_coordinates[index]=x_coordinates[index-1]+Radial_grid_intervals[i];
			index++;
		}
	}
	x_count=points_count;
}

void setup_radial_grid(){
	char* sprintf_buffer=new char[Log_length+1];
	write_log((char*)"----Radial-grid block----");
	// Radial grid
	if(Rg_block_appeared==false){
		write_log((char*)"Default radial grid is used");
		initialize_radial_grid();
	}
	generate_radial_grid();
	int i;
	int last_index=0;
	for(i=0; i<Radial_grid_count; i++){
		sprintf(sprintf_buffer, "    Block %3d: x = %8.3f to %8.3f, %3d points, interval %8f", (i+1), x_coordinates[last_index], x_coordinates[last_index+Radial_grid_points[i]], Radial_grid_points[i], Radial_grid_intervals[i]);
		write_log(sprintf_buffer);
		last_index+=Radial_grid_points[i];
	}
	delete sprintf_buffer;
}

int setup_potential(int Z, double mu){
	int status=1;
	int i;
	FILE* Potential_file_obj;
	string input_line_s;
	char* input_line_c=new char[buffer_length+1];
	int actual_line_length;
	int sscanf_status;
	int current_index=0;
	double current_x;
	double x_error;
	
	if(strcmp(At_potential, "H-like")==0){
		// hydrogren-like potential: v(x)=-Z/(mu x)
		for(i=0; i<x_count; i++){
			if(i==0){
				At_v_x[i]=0;
			}else{
				At_v_x[i]=-Z/(mu*x_coordinates[i]);
			}
		}
	}else if(strcmp(At_potential, "Thomas-Fermi")==0 || strcmp(At_potential, "file")==0){
		// Thomas-Fermi potential: v(x)=-Z/(mu x) * phi(x)
		// phi(x) is from the file
		// or
		// file: v(x) is from the file
		ifstream Potential_file_obj(At_potential_file);
		if(!Potential_file_obj.is_open()){
			write_log((char*)"Error: could not open the potential file");
			status=0; goto FINALIZATION;
		}
		
		while(getline(Potential_file_obj, input_line_s)){
			// read one line from cin (stdin), copy to input_line_c
			actual_line_length=input_line_s.length();
			if(actual_line_length>buffer_length){
				actual_line_length=buffer_length;
			}
			input_line_s.copy(input_line_c, actual_line_length);
			input_line_c[actual_line_length]='\0';
			// parse
			sscanf_status=sscanf(input_line_c, "%lf %lf", &current_x, &At_v_x[current_index]);
			if(sscanf_status==2){
				// parse succeeded
				if(current_index!=0){
					// error check
					x_error=abs(1.0-current_x/x_coordinates[current_index]);
					if(x_error>Data_read_error){
						// error is too large
						write_log((char*)"Error: x coordinate mismatch");
						status=0; goto FINALIZATION;
					}
					if(strcmp(At_potential, "Thomas-Fermi")==0){
						// only for Thomas-Fermi potential
						if(At_v_x[current_index]>1.0/Z){
							At_v_x[current_index]*=-Z/(mu*x_coordinates[current_index]);
						}else{
							At_v_x[current_index]=-1.0/(mu*x_coordinates[current_index]);
						}
					}
				}else{
					At_v_x[current_index]=0;
				}
				current_index++;
				if(current_index==x_count){
					break;
				}
			}
		}
	}
	goto FINALIZATION;
	
 FINALIZATION:
	delete input_line_c;
	return status;
}
