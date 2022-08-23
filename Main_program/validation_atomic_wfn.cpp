// calcPSF validation of atomic wavefunction calculation

#include <cstring>

#include "log.hpp"
#include "variables_ext.hpp"
#include "setup.hpp"

int validation_atomic_wfn(){
	char* sprintf_buffer=new char[Log_length+1];
	int status=1;
	write_log((char*)"----Input value validation----");
	// Radial-grid
	setup_radial_grid();
	At_p_x=new double[x_count];
	At_p_diff_x=new double[x_count];
	At_v_x=new double[x_count];

	// Atomic-wfn
	write_log((char*)"----Atomic-wfn block----");
	/// n
	if((n_min_set && !n_max_set) || (!n_min_set && n_max_set)){
		write_log((char*)"Error: when you specify either of n_min and n_max, you also need to specify the other");
		status=0; goto FINALIZATION;
	}
	if(n_min_set && n_single_set){
		write_log((char*)"Error: you cannot specify n in both min & max and single ways");
		status=0; goto FINALIZATION;
	}
	if(!n_min_set && !n_single_set){
		write_log((char*)"Error: you must specify n in either of min & max or single");
		status=0; goto FINALIZATION;
	}
	//// copy n_single to n_min and n_max if n_single_set
	if(n_single_set){
		n_min=n_single;
		n_max=n_single;
	}
	//// value validation: 1<=n_min<=n_max
	if(n_min>n_max){
		write_log((char*)"Error: n_min cannot be larger than n_max");
		status=0; goto FINALIZATION;
	}
	if(n_min<1){
		write_log((char*)"Error: n must be larger than 0");
		status=0; goto FINALIZATION;
	}
	if(n_min==n_max){
		sprintf(sprintf_buffer, "%32s = %3d", "Principal quantum number n", n_min);
		write_log(sprintf_buffer);
	}else{
		sprintf(sprintf_buffer, "%32s = %3d to %3d", "Principal quantum number n", n_min, n_max);
		write_log(sprintf_buffer);
	}
	
	/// l
	if((l_min_set && !l_max_set) || (!l_min_set && l_max_set)){
		write_log((char*)"Error: when you specify either of l_min and l_max, you also need to specify the other");
		status=0; goto FINALIZATION;
	}
	if(l_min_set && l_single_set){
		write_log((char*)"Error: you cannot specify l in both min & max and single ways");
		status=0; goto FINALIZATION;
	}
	if(!l_min_set && !l_single_set){
		write_log((char*)"Error: you must specify l in either of min & max or single");
		status=0; goto FINALIZATION;
	}
	//// copy l_single to l_min and l_max if l_single_set
	if(l_single_set){
		l_min=l_single;
		l_max=l_single;
	}
	//// value validation: 0<=l_min<=l_max
	if(l_min>l_max){
		write_log((char*)"Error: l_min cannot be larger than l_max");
		status=0; goto FINALIZATION;
	}
	if(l_min<0){
		write_log((char*)"Error: l must be positive");
		status=0; goto FINALIZATION;
	}
	if(l_min==l_max){
		sprintf(sprintf_buffer, "%32s = %3d", "Azimuthal quantum number l", l_min);
		write_log(sprintf_buffer);
	}else{
		sprintf(sprintf_buffer, "%32s = %3d to %3d", "Azimuthal quantum number l", l_min, l_max);
		write_log(sprintf_buffer);
	}

	/// Z
	if((Z_min_set && !Z_max_set) || (!Z_min_set && Z_max_set)){
		write_log((char*)"Error: when you specify either of Z_min and Z_max, you also need to specify the other");
		status=0; goto FINALIZATION;
	}
	if(Z_min_set && Z_single_set){
		write_log((char*)"Error: you cannot specify Z in both min & max and single ways");
		status=0; goto FINALIZATION;
	}
	if(!Z_min_set && !Z_single_set){
		write_log((char*)"Error: you must specify Z in either of min & max or single");
		status=0; goto FINALIZATION;
	}
	//// copy Z_single to Z_min and Z_max if Z_single_set
	if(Z_single_set){
		Z_min=Z_single;
		Z_max=Z_single;
	}
	//// value validation: 1<=Z_min<=Z_max
	if(Z_min>Z_max){
		write_log((char*)"Error: Z_min cannot be larger than Z_max");
		status=0; goto FINALIZATION;
	}
	if(Z_min<0){
		write_log((char*)"Error: Z must be larger than 0");
		status=0; goto FINALIZATION;
	}
	if(Z_min==Z_max){
		sprintf(sprintf_buffer, "%32s = %3d", "Atomic number Z", Z_min);
		write_log(sprintf_buffer);
	}else{
		sprintf(sprintf_buffer, "%32s = %3d to %3d", "Atomic number Z", Z_min, Z_max);
		write_log(sprintf_buffer);
	}

	/// Potential: H-like, Thomas-Fermi, file
	if(strcmp(At_potential, "H-like")==0 || strcmp(At_potential, "Thomas-Fermi")==0 || strcmp(At_potential, "file")==0){
		sprintf(sprintf_buffer, "%32s = %s", "Potential", At_potential);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: invalid potential");
		status=0; goto FINALIZATION;
	}
	if(strcmp(At_potential, "H-like")!=0){
		// Thomas-Fermi or file potentials need Potential_file
		if(At_potential_file_set){
			sprintf(sprintf_buffer, "%32s = %s", "Potential_file", At_potential_file);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Thomas-Fermi or file potentials need Potential_file");
			status=0; goto FINALIZATION;
		}
	}

	/// Solution (At_solution): RK1, Numerov
	if(strcmp(At_solution, "RK1")==0 /*|| strcmp(At_solution, "RK4")==0*/ || strcmp(At_solution, "Numerov")==0){
		sprintf(sprintf_buffer, "%32s = %s", "Solution", At_solution);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: Solution should be 'RK1' or 'Numerov'");
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
		
	goto FINALIZATION;
 FINALIZATION:
	delete sprintf_buffer;
	return status;
}
