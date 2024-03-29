// calcPSF validation of Thomas-Fermi potential calculation

#include <cstring>

#include "log.hpp"
#include "variables_ext.hpp"
#include "setup.hpp"

int validation_Thomas_Fermi(){
	char* sprintf_buffer=new char[Log_length+1];
	int status=1;
	write_log((char*)"----Input value validation----");
	// Radial-grid
	setup_radial_grid();
	TF_phi=new double[x_count];
	TF_phi_diff=new double[x_count];

	// Thomas-Fermi
	write_log((char*)"----Thomas-Fermi block----");
	sprintf(sprintf_buffer, "%32s = %s", "Calculation_test", (TF_test?"true":"false"));
	write_log(sprintf_buffer);
	if(TF_test){
		// test mode
		// Initial_diff_(offset, delta, size) are necessary
		if(Initial_diff_offset_set && Initial_diff_delta_set && Initial_diff_size_set){
			sprintf(sprintf_buffer, "%32s = %-8.3f", "Initial_diff_offset", Initial_diff_offset);
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "%32s = %-8.3f", "Initial_diff_delta", Initial_diff_delta);
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "%32s = %d", "Initial_diff_size", Initial_diff_size);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: in the test mode, Initial_diff_offset, delta, and size are necessary");
			status=0; goto FINALIZATION;
		}
	}else{
		// actual mode
		// Initial_diff_(min, max) are necessary
		sprintf(sprintf_buffer, "%32s = %-8.3f", "Initial_diff_min", Initial_diff_min);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "%32s = %-8.3f", "Initial_diff_max", Initial_diff_max);
		write_log(sprintf_buffer);
		if(Initial_diff_min > Initial_diff_max){
			write_log((char*)"Error: Initial_diff_min should be smaller than Initial_diff_max");
			status=0; goto FINALIZATION;
		}
		sprintf(sprintf_buffer, "%32s = %-8.3e", "Threshold", TF_threshold);
		write_log(sprintf_buffer);
		if(TF_threshold<0){
			write_log((char*)"Error: Threshold should be positive");
			status=0; goto FINALIZATION;
		}
	}
	///solution: RK1 or RK4
	if(strcmp(TF_solution, "RK1")==0 || strcmp(TF_solution, "RK4")==0){
		sprintf(sprintf_buffer, "%32s = %s", "Solution", TF_solution);
		write_log(sprintf_buffer);
	}else{
		write_log((char*)"Error: Solution should be 'RK1' or 'RK4'");
		status=0; goto FINALIZATION;
	}

  write_log((char*)"");
	
	goto FINALIZATION;
 FINALIZATION:
	delete[] sprintf_buffer;
	return status;
}
