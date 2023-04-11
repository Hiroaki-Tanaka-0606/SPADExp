// calcPSF validation of self-consistent field calculation of an atom

#include <cstring>

#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"

int validation_scf_atom(){
	char* sprintf_buffer=new char[Log_length+1];
	int status=1;
	int i, j;
	int total_occupation=0;
	int SC_index=0;
	write_log((char*)"----Input value validation----");
	// Radial-grid
	setup_radial_grid();
	At_p_x=new double[x_count];
	At_p_diff_x=new double[x_count];
	At_v_x=new double[x_count];
	SC_sigma_x=new double[x_count];

	// properties from the Atomic-wfn block
	write_log((char*)"----Atomic-wfn block----");
	/// Z: must be Z_single
	if(!Z_single_set){
		write_log((char*)"Error: you must specify Z");
		status=0; goto FINALIZATION;
	}
	if(Z_min_set || Z_max_set){
		write_log((char*)"Error: you cannot specify multiple Z in SCF calculation");
		status=0; goto FINALIZATION;
	}

	sprintf(sprintf_buffer, "%32s = %d", "Atomic number Z", Z_single);
	write_log(sprintf_buffer);
	
	/// (Initial) potential: Thomas-Fermi
	if(At_potential_file_set){
		sprintf(sprintf_buffer, "%32s = %s", "Potential_file", At_potential_file);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: Potential_file is necessary for the initial potential");
		status=0; goto FINALIZATION;
	}
	strcpy(At_potential, "Thomas-Fermi");

	/// Solution (At_solution)
	if(strcmp(At_solution, "RK1")==0 || strcmp(At_solution, "RK4")==0 || strcmp(At_solution, "Numerov")==0){
		sprintf(sprintf_buffer, "%32s = %s", "Solution", At_solution);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: Solution should be 'RK1', 'RK4', or 'Numerov'");
		status=0; goto FINALIZATION;
	}

	/// Bisection_step
	if(Bisection_step<=0){
		write_log((char*)"Error: Bisection_step should be larger than 0");
		status=0; goto FINALIZATION;
	}
	sprintf(sprintf_buffer, "%32s = %-10.3e", "Bisection_step", Bisection_step);
	write_log(sprintf_buffer);

	/// Threshold
	sprintf(sprintf_buffer, "%32s = %-10.3e", "E_Threshold", At_E_threshold);
	write_log(sprintf_buffer);
	if(At_E_threshold<0){
		write_log((char*)"Error: E_Threshold should be positive");
		status=0; goto FINALIZATION;
	}
	
	/// Radius_factor
	sprintf(sprintf_buffer, "%32s = %-10.3e", "Radius_factor", At_radius_factor);
	write_log(sprintf_buffer);
	if(At_radius_factor<0){
		write_log((char*)"Error: Radius_factor should be positive");
		status=0; goto FINALIZATION;
	}

	// Occupation block
	write_log((char*)"----Occupation block----");
	if(!Oc_block_appeared){
		write_log((char*)"Error: Occupation block does not exist");
		status=0; goto FINALIZATION;
	}

	for(i=0; i<Occupation_count; i++){
		for(j=0; j<=i; j++){
			if(At_occupation[i][j]>0){
				sprintf(sprintf_buffer, "%26s%1d%4s%1d = %d", "Occupation at n=", (i+1), ", l=", j, At_occupation[i][j]);
				write_log(sprintf_buffer);
				total_occupation+=At_occupation[i][j];
				SC_orbital_count++;
			}
		}
	}
	if(total_occupation!=Z_single){
		write_log((char*)"Error: Sum of occupation is different from Z");
		status=0; goto FINALIZATION;
	}
	SC_n_list=new int[SC_orbital_count];
	SC_l_list=new int[SC_orbital_count];
	SC_eigen_list=new double[SC_orbital_count];
	SC_p_x=new double*[SC_orbital_count];
	for(i=0; i<Occupation_count; i++){
		for(j=0; j<=i; j++){
			if(At_occupation[i][j]>0){
				SC_n_list[SC_index]=i+1;
				SC_l_list[SC_index]=j;
				SC_eigen_list[SC_index]=0;
				SC_p_x[SC_index]=new double[x_count];
				SC_index++;
			}
		}
	}
		
	// properties from the SCF-atom block
	write_log((char*)"----SCF-atom block----");

	/// Mix_weight (SC_mix_weight)
	sprintf(sprintf_buffer, "%32s = %8.3f", "Mix_weight", SC_mix_weight);
	write_log(sprintf_buffer);
	if(SC_mix_weight<0 || SC_mix_weight>1){
		write_log((char*)"Error: Mix_weight should be between 0.0 and 1.0");
		status=0; goto FINALIZATION;
	}

	/// self-consistent criterion (SC_criterion_a and SC_criterion_b)
	sprintf(sprintf_buffer, "%32s = %10.3e, %10.3e", "SCF criterion (a, b)", SC_criterion_a, SC_criterion_b);
	write_log(sprintf_buffer);
	if(SC_criterion_a<0 || SC_criterion_b<0){
		write_log((char*)"Error: Criteria should be positive");
		status=0; goto FINALIZATION;
	}
	
	goto FINALIZATION;
 FINALIZATION:
	delete[] sprintf_buffer;
	return status;
}
