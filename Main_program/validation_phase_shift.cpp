// calcPSF phase shift calculation

#include <cstring>
#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"

int validation_phase_shift(){
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

	/// Potential database (At_potential_file)
	if(At_potential_file_set){
		sprintf(sprintf_buffer, "%32s = %s", "Potential_file", At_potential_file);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Potential_file not found, V(r)=-1/r is used");
	}

	/// Solution (At_solution)
	if(strcmp(At_solution, "RK1")==0 || strcmp(At_solution, "RK4")==0 || strcmp(At_solution, "Numerov")==0){
		sprintf(sprintf_buffer, "%32s = %s", "Solution", At_solution);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: Solution should be 'RK1', 'RK4', or 'Numerov'");
		status=0; goto FINALIZATION;
	}

	// Excitation-energy block
	write_log((char*)"----Excitation-energy block----");
	if(!Ex_block_appeared){
		write_log((char*)"Error: Excitation-energy block is necessary");
		status=0; goto FINALIZATION;
	}
	for(i=0; i<Ex_energy_count; i++){
		sprintf(sprintf_buffer, "%28s[%2d] = %8.2f eV (%8.2f Eh)", "Excitation energy", i, Ex_energies[i]*Eh, Ex_energies[i]);
		write_log(sprintf_buffer);
	}

	// Orbital block
	write_log((char*)"----Orbital block----");
	if(!Or_block_appeared){
		write_log((char*)"Error: Orbital block is necessary");
		status=0; goto FINALIZATION;
	}
	for(i=0; i<Ph_orbital_count; i++){
		sprintf(sprintf_buffer, "%32s : l = %1d, EB = %8.2f eV (%8.2f Eh)", Ph_orbital_labels[i], Ph_l_list[i], Ph_binding_energies[i]*Eh, Ph_binding_energies[i]);
		write_log(sprintf_buffer);
	}

	// Phase-shift block
	write_log((char*)"----Phase-shift block----");
	sprintf(sprintf_buffer, "%32s = %d", "Skip_points", Ph_skip_points);
	write_log(sprintf_buffer);
	if(Ph_skip_points<0){
		write_log((char*)"Error: Skip_points should not be negative");
		status=0; goto FINALIZATION;
	}
	sprintf(sprintf_buffer, "%32s = %d", "Calc_points", Ph_calc_points);
	write_log(sprintf_buffer);
	if(Ph_calc_points<1){
		write_log((char*)"Error: Calc_points should be positive");
		status=0; goto FINALIZATION;
	}

		
	goto FINALIZATION;
 FINALIZATION:
	delete[] sprintf_buffer;
	return status;
	return 1;
}
