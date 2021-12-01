// calcPSF validation of photoemission structure factor calculation
#include <cstring>

#include "log.hpp"
#include "variables_ext.hpp"
#include "setup.hpp"

int validation_PSF(){
	char* sprintf_buffer=new char[Log_length+1];
	int status=1;
	write_log((char*)"----Input value validation----");

	// PSF
	write_log((char*)"----PSF block----");
	/// Input_file
	if(!PS_input_file_set){
		write_log((char*)"Error: Input_file is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %s", "Input_file", PS_input_file);
		write_log(sprintf_buffer);
	}
	/// E_min, E_max, and E_pixel
	if(!PS_E_min_set || !PS_E_max_set || !PS_E_pixel_set){
		write_log((char*)"Error: E_min, E_max, and E_pixel are necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %8.3f eV to %8.3f eV, step %8.3f eV", "Energy range", PS_E_min, PS_E_max, PS_E_pixel);
		write_log(sprintf_buffer);
	}
	/// dE
	if(!PS_dE_set){
		write_log((char*)"Error: dE is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %8.3f eV", "dE", PS_dE);
		write_log(sprintf_buffer);
	}
	/// Initial_state
	if(!PS_initial_state_set){
		write_log((char*)"Error: Initial_state is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PS_initial_state, "PAO")==0 || strcmp(PS_initial_state, "AO")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Initial_state", PS_initial_state);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Initial_state should be 'PAO' or 'AO'");
			status=0; goto FINALIZATION;
		}
	}
	/// Final_state
	if(!PS_final_state_set){
		write_log((char*)"Error: Final_state is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PS_final_state, "PW")==0 || strcmp(PS_final_state, "Calc")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Final_state", PS_final_state);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Final_state should be 'PW' or 'Calc'");
			status=0; goto FINALIZATION;
		}
	}
	/// Polarization
	if(!PS_polarization_set){
		write_log((char*)"Error: Polarization is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PS_polarization, "Linear")==0 || strcmp(PS_polarization, "RCircular")==0 || strcmp(PS_polarization, "LCircular")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Polarization", PS_polarization);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Polarization should be 'Linear', 'RCircular' or 'LCircular'");
			status=0; goto FINALIZATION;
		}
	}
	/// theta and phi
	if(!PS_theta_set || !PS_phi_set){
		write_log((char*)"Error: Theta and Phi are necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = (%8.3f deg, %8.3f deg)", "Polarization angle", PS_theta, PS_phi);
		write_log(sprintf_buffer);
	}
	/// AO_file
	if(!PS_AO_file_set){
		write_log((char*)"Error: Atomic_orbitals_file is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %s", "Atomic_orbitals_file", PS_AO_file);
		write_log(sprintf_buffer);
	}
		
	goto FINALIZATION;
 FINALIZATION:
	delete sprintf_buffer;
	return status;
}
