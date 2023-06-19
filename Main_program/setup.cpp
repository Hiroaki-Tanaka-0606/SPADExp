// calcPSF setup (initialization and finalization)

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <H5Cpp.h>
#include "variables_ext.hpp"
#include "log.hpp"
#include "HDF5_tools.hpp"

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
	Data_read_error=0.0001;
	HDF5_file_length=1024;
	PA_state_length=64;
	au_ang=0.529177; // (Ang)
	Eh=27.2114; // (eV)
	Ph_label_length=8;
	
	// variables
	/// blocks
	Ct_block_appeared=false;
	TF_block_appeared=false;
	At_block_appeared=false;
	Rg_block_appeared=false;
	SC_block_appeared=false;
	Oc_block_appeared=false;
	PA_block_appeared=false;
	Ex_block_appeared=false;
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
	/// SC block
	SC_mix_weight=0.5;
	SC_mix_weight_set=false;
	SC_criterion_a=0.001;
	SC_criterion_a_set=false;
	SC_criterion_b=0.001;
	SC_criterion_b_set=false;
	SC_orbital_count=0;
	/// Ph block
	Ph_skip_points_set=false;
	Ph_skip_points=2;
	Ph_calc_points_set=false;
	Ph_calc_points=5;
	/// PA block
	PA_input_file_set=false;
	PA_input_file=new char[HDF5_file_length+1];
	PA_E_min_set=false;
	PA_E_max_set=false;
	PA_dE_set=false;
	PA_E_pixel_set=false;
	PA_initial_state_set=false;
	PA_initial_state=new char[PA_state_length+1];
	PA_final_state_set=false;
	PA_final_state=new char[PA_state_length+1];
	PA_final_state_step=0.01;
	PA_final_state_step_set=false;
	PA_polarization_set=false;
	PA_polarization=new char[PA_state_length+1];
	PA_theta_set=false;
	PA_phi_set=false;
	PA_AO_file_set=false;
	PA_AO_file=new char[HDF5_file_length+1];
	PA_sigma_max=5.0;
	PA_decay_max=10;
	PA_ext_up=0;
	PA_ext_ri=0;
	PA_ext_dn=0;
	PA_ext_le=0;
	PA_ext_set=false;
	PA_weighting=false;
	PA_weighting_set=false;
	PA_weighting_axis=new double[3];
	PA_weighting_axis_set=false;
	PA_weighting_shape=new char[PA_state_length+1];
	PA_weighting_shape_set=false;
	PA_weighting_origin_set=false;
	PA_weighting_width_set=false;
	PA_use_angstrom=true;
	PA_use_angstrom_set=false;
	PA_output_data=new char[PA_state_length+1];
	strcpy(PA_output_data, "PAD");
	PA_output_data_set=false;
	PA_reflection=false;
	PA_reflection_set=false;
	PA_reflection_coef=0.0;
	PA_reflection_coef_set=false;
	PA_include_neg_depth=true;
	PA_include_neg_depth_set=false;
	PA_excitation_energy=0;
	PA_excitation_energy_set=false;
	PA_FPFS_energy_step=0.01;
	PA_FPFS_energy_step_set=false;
	PA_FPFS_range=3;
	PA_FPFS_kRange=2.0;
	PA_FPFS_kRange_set=false;
	PA_FPFS=false;
	PA_FPFS_file=new char[HDF5_file_length+1];
	PA_FPFS_file_set=false;
	PA_VPS_file=new char[HDF5_file_length+1];
	PA_VPS_file_set=false;
	PA_Lebedev_order_ave=86;
	PA_Lebedev_order_int=266;
	PA_lp_max=5;
	PA_zero_threshold=1e-8;
	PA_theta_points=128;
	PA_ignore_core=false;
	PA_ignore_core_set=false;
	PA_add_nonorth_term=false;
	PA_add_nonorth_term_set=false;
	PA_calc_all_nonloc=false;
	PA_calc_all_nonloc_set=false;
	PA_calc_all_loc=true;
	PA_calc_all_loc_set=false;
	PA_ignore_nonlocal=false;
	PA_ignore_nonlocal_set=false;
	PA_FPFS_Numerov=false;
	PA_FPFS_Numerov_set=false;
	PA_FPFS_margin=15.0;
	PA_FPFS_edge_smoothing=false;
	PA_FPFS_edge_smoothing_set=false;
	PA_FPFS_smoothing_E=0;
	PA_FPFS_smoothing_E_set=false;
	PA_FPFS_smoothing_k=0;
	PA_FPFS_smoothing_k_set=false;
	PA_FPFS_bulk_set=false;
	PA_FPFS_bulk_min_ang=0.0;
	PA_FPFS_bulk_max_ang=0.0;
	PA_FPFS_bulk_count=0;
	PA_FPFS_bulk_kz_steps=100;
	PA_FPFS_bulk_kappaz_steps=200;
	PA_interpolate_wfn=true;
	PA_interpolate_wfn_coef=0.5;
	PA_FPFS_bulk_tolerance=0.000;
	PA_FPFS_bulk_buffer_size=128;
	PA_FPFS_perturbation=false;
	PA_FPFS_perturbation_set=false;
	PA_calc_complex_dispersion=true;
	PA_calc_complex_dispersion_set=false;
	PA_FPFS_complex_dispersion_criterion=0.02;
	PA_FPFS_dispersion_group_criterion_left=0.01;
	PA_FPFS_dispersion_group_criterion_right=0.03;
	PA_FPFS_gap_coefficient=0.3;
	
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
	delete PA_input_file;
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
	delete[] sprintf_buffer;
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
				At_v_x[i]=-Z*1.0/(mu*x_coordinates[i]);
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
	delete[] input_line_c;
	return status;
}

void modify_potential(int Z, double mu){
	if(strcmp(At_potential, "Thomas-Fermi")==0){
		for(int i=1; i<x_count; i++){
			// exclude i=0 because At_v_x[0]=0
			// only for Thomas-Fermi potential
			if(At_v_x[i]>1.0/Z){
				At_v_x[i]*=-Z*1.0/(mu*x_coordinates[i]);
			}else{
				At_v_x[i]=-1.0/(mu*x_coordinates[i]);
			}
		}
	}
}

double load_potential_H5(Group atomG, double mu){
	int x_count2=r_att_int(atomG, "length");
	double x_coordinates2[x_count2];
	double v_x2[x_count2];
	r_data_1d(atomG, "Potential", x_count2, &v_x2[0]);
	r_att_1d(atomG, "x", x_count2, &x_coordinates2[0]);

	At_v_x[0]=v_x2[0];
	int i, j;
	for(i=1; i<x_count; i++){
		bool x_found=false;
		for(j=0; j<x_count2-1; j++){
			if(x_coordinates2[j]<x_coordinates[i] && x_coordinates[i]<x_coordinates2[j+1]){
				double dx1=x_coordinates[i]-x_coordinates2[j];
				double dx2=x_coordinates2[j+1]-x_coordinates[i];
				At_v_x[i]=(v_x2[j]*dx2+v_x2[j+1]*dx1)/(dx1+dx2);
				x_found=true;
				break;
			}else if(abs(x_coordinates2[j+1]-x_coordinates[i])<Data_read_error){
				At_v_x[i]=v_x2[j+1];
				x_found=true;
				break;
			}
		}
		if(!x_found){
			// if the potential is not in database, use V(r)=-1/r=-1/mu*x
			At_v_x[i]=-1.0/(mu*x_coordinates[i]);
		}
	}
	int mod_start_index=r_att_int(atomG, "mod_start_index");
	return x_coordinates2[mod_start_index];
}
