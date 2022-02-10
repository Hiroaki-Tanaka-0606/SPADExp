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
			if(strcmp(PS_final_state, "Calc")==0){
				// properties for final state calculations: see validation_phase_shift.cpp
				/// grid
				setup_radial_grid();
				At_p_x=new double[x_count];
				At_p_diff_x=new double[x_count];
				At_v_x=new double[x_count];
				write_log((char*)"----Radial-grid block ended----");
				/// Final_state_step (PS_final_state_step)
				sprintf(sprintf_buffer, "%32s = %8.3f", "Final_state_step", PS_final_state_step);
				write_log(sprintf_buffer);
				if(PS_final_state_step<0){
					write_log((char*)"Error: Final_state_step should be positive");
					status=0; goto FINALIZATION;
				}
				if(At_potential_file_set){
					sprintf(sprintf_buffer, "%32s = %s", "Potential_file (Atomic-wfn)", At_potential_file);
					write_log(sprintf_buffer);
				}else{
					write_log((char*)"Error: Potential_file is necessary");
					status=0; goto FINALIZATION;
				}
				/// Solution (At_solution)
				if(strcmp(At_solution, "RK1")==0 || strcmp(At_solution, "RK4")==0 || strcmp(At_solution, "Numerov")==0){
					sprintf(sprintf_buffer, "%32s = %s", "Solution (Atomic-wfn)", At_solution);
					write_log(sprintf_buffer);
				}else{
					write_log((char*)"Error: Solution should be 'RK1', 'RK4', or 'Numerov'");
					status=0; goto FINALIZATION;
				}
				/// Skip_ and Calc_ points
				sprintf(sprintf_buffer, "%32s = %d", "Skip_points (Phase-shift)", Ph_skip_points);
				write_log(sprintf_buffer);
				if(Ph_skip_points<0){
					write_log((char*)"Error: Skip_points should not be negative");
					status=0; goto FINALIZATION;
				}
				sprintf(sprintf_buffer, "%32s = %d", "Calc_points (Phase-shift)", Ph_calc_points);
				write_log(sprintf_buffer);
				if(Ph_calc_points<1){
					write_log((char*)"Error: Calc_points should be positive");
					status=0; goto FINALIZATION;
				}
	
			}
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
	/// Extend
	if(PS_ext_set){
		if(PS_ext_up<0 || PS_ext_ri<0 || PS_ext_dn<0 || PS_ext_le<0){
			write_log((char*)"Error: Extend should not be negative");
			status=0; goto FINALIZATION;
		}
		sprintf(sprintf_buffer, "%32s = %d %d %d %d", "Extend", PS_ext_up, PS_ext_ri, PS_ext_dn, PS_ext_le);
		write_log(sprintf_buffer);
	}
	/// Weighting
	sprintf(sprintf_buffer, "%32s = %s", "Weighting", PS_weighting ? "true" : "false");
	write_log(sprintf_buffer);
	if(PS_weighting){
		const char* unit=PS_use_angstrom ? "Ang" : "a.u.";
		/// Weighting_axis
		if(PS_weighting_axis_set){
			sprintf(sprintf_buffer, "%32s = (%8.2f, %8.2f, %8.2f) [%s]", "Weighting_axis", PS_weighting_axis[0], PS_weighting_axis[1], PS_weighting_axis[2], unit);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Weighting_axis is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_shape
		if(PS_weighting_shape_set){
			if(strcmp(PS_weighting_shape, "Rect")==0 || strcmp(PS_weighting_shape, "Exp")==0){
				sprintf(sprintf_buffer, "%32s = %s", "Weighting_shape", PS_weighting_shape);
				write_log(sprintf_buffer);
			}else{
				write_log((char*)"Error: Weighting_shape should be Rect or Exp");
				status=0; goto FINALIZATION;
			}
		}else{
			write_log((char*)"Error: Weighting_shape is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_origin
		if(PS_weighting_origin_set){
			if(PS_weighting_origin>=0){
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], parallel to the axis", "Distance of the origin", PS_weighting_origin, unit);
				write_log(sprintf_buffer);
			}else{
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], anti-parallel to the axis", "Distance of the origin", PS_weighting_origin, unit);
				write_log(sprintf_buffer);
			}
		}else{
			write_log((char*)"Error: Weighting_origin is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_width
		if(PS_weighting_width_set){
			if(PS_weighting_width>=0){
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], parallel to the axis", "Weighting_width", PS_weighting_width, unit);
				write_log(sprintf_buffer);
			}else{
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], anti-parallel to the axis", "Weighting_width", PS_weighting_width, unit);
				write_log(sprintf_buffer);
			}
		}else{
			write_log((char*)"Error: Weighting_width is necessary");
			status=0; goto FINALIZATION;
		}
	}
	/// Output_data
	sprintf(sprintf_buffer, "%32s = %s", "Output_data", PS_output_data);
	write_log(sprintf_buffer);
	if(strcmp(PS_output_data, "PSF")!=0 && strcmp(PS_output_data, "Band")!=0 /*&& strcmp(PS_output_data, "PSF_real")!=0 && strcmp(PS_output_data, "PSF_imag")!=0*/){
		write_log((char*)"Error: Invalid Output_data");
		status=0; goto FINALIZATION;
	}
	goto FINALIZATION;
 FINALIZATION:
	delete sprintf_buffer;
	return status;
}
