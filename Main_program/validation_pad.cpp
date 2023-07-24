// calcPSF validation of photoemission angular distribution calculation
#include <cstring>

#include "log.hpp"
#include "variables_ext.hpp"
#include "setup.hpp"

int validation_PAD(){
	char* sprintf_buffer=new char[Log_length+1];
	int status=1;
	write_log((char*)"----Input value validation----");

	// PAD
	write_log((char*)"----PAD block----");
	/// Input_file
	if(!PA_input_file_set){
		write_log((char*)"Error: Input_file is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %s", "Input_file", PA_input_file);
		write_log(sprintf_buffer);
	}
	/// E_min, E_max, and E_pixel
	if(!PA_E_min_set || !PA_E_max_set || !PA_E_pixel_set){
		write_log((char*)"Error: E_min, E_max, and E_pixel are necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %8.3f eV to %8.3f eV, step %8.3f eV", "Energy range", PA_E_min, PA_E_max, PA_E_pixel);
		write_log(sprintf_buffer);
	}
	/// dE
	if(!PA_dE_set){
		write_log((char*)"Error: dE is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %8.3f eV", "dE", PA_dE);
		write_log(sprintf_buffer);
	}
	/// Initial_state
	if(!PA_initial_state_set){
		write_log((char*)"Error: Initial_state is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PA_initial_state, "PAO")==0 || strcmp(PA_initial_state, "AO")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Initial_state", PA_initial_state);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Initial_state should be 'PAO' or 'AO'");
			status=0; goto FINALIZATION;
		}
	}
	sprintf(sprintf_buffer, "%32s = %s", "Interpolate_wfn", PA_interpolate_wfn?"true":"false");
	write_log(sprintf_buffer);
	if(PA_interpolate_wfn){
		sprintf(sprintf_buffer, "%32s = %.2f", "Interpolate_wfn_coef", PA_interpolate_wfn_coef);
		write_log(sprintf_buffer);
	}
	/// Final_state
	if(!PA_final_state_set){
		write_log((char*)"Error: Final_state is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PA_final_state, "PW")==0 || strcmp(PA_final_state, "Calc")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Final_state", PA_final_state);
			write_log(sprintf_buffer);
			if(strcmp(PA_final_state, "Calc")==0){
				// properties for final state calculations: see validation_phase_shift.cpp
				/// grid
				setup_radial_grid();
				At_p_x=new double[x_count];
				At_p_diff_x=new double[x_count];
				At_v_x=new double[x_count];
				write_log((char*)"----Radial-grid block ended----");
				/// Final_state_step (PA_final_state_step)
				sprintf(sprintf_buffer, "%32s = %8.3f", "Final_state_step", PA_final_state_step);
				write_log(sprintf_buffer);
				if(PA_final_state_step<0){
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
		}else if(strcmp(PA_final_state, "FP_PAO")==0 || strcmp(PA_final_state, "FP_AO")==0){
			PA_FPFS=true;
			sprintf(sprintf_buffer, "%32s = %s", "Final_state", PA_final_state);
			write_log(sprintf_buffer);
			// weighting axis and origin are used to determine empty atoms
			if(PA_weighting==false){
				write_log((char*)"Error: Weighting should be set to true in FP_PAO or FP_AO");
				//write_log((char*)"Error: Weighting should be set to true in FP_PAO");
				status=0; goto FINALIZATION;
			}
			sprintf(sprintf_buffer, "%32s = %s", "Ignore_nonlocal", PA_ignore_nonlocal?"true":"false");
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "%32s = %s", "FPFS_Numerov", PA_FPFS_Numerov?"true":"false");
			write_log(sprintf_buffer);
			// reflection cannot be used
			if(PA_reflection==true){
				write_log((char*)"Error: Reflection cannot be used in FP_PAO and FP_AO");
				//write_log((char*)"Error: Reflection cannot be used in FP_PAO");
				status=0; goto FINALIZATION;
			}
			/// Excitation energy
			if(PA_excitation_energy_set==false){
				write_log((char*)"Error: Excitation_energy is necessary");
				status=0; goto FINALIZATION;
			}
			sprintf(sprintf_buffer, "%32s = %.2f eV", "Excitation_energy", PA_excitation_energy);
			write_log(sprintf_buffer);
			/// VPS_file
			if(PA_VPS_file_set==false){
				write_log((char*)"Error: VPS_file is necessary");
				status=0; goto FINALIZATION;
			}
			sprintf(sprintf_buffer, "%32s = %s", "VPS_file", PA_VPS_file);
			write_log(sprintf_buffer);
			/// Nonloc_offset
			sprintf(sprintf_buffer, "%32s = %.2f Bohr", "Nonloc_offset", PA_FPFS_nonloc_offset);
			write_log(sprintf_buffer);
			/// FPFS_energy_step
			sprintf(sprintf_buffer, "%32s = %.3f eV", "FPFS_energy_step", PA_FPFS_energy_step);
			write_log(sprintf_buffer);
			if(PA_FPFS_Numerov){
				PA_FPFS_kRange=1.0;
				write_log((char*)"Note: FPFS_kRange is set to 1 because FPFS_Numerov is true");
			}
			/// FPFS_bulk
			if(PA_FPFS_bulk_set){
				sprintf(sprintf_buffer, "%32s = [%.3f, %.3f] (Ang)", "FPFS_bulk range", PA_FPFS_bulk_min_ang, PA_FPFS_bulk_max_ang);
				write_log(sprintf_buffer);
				sprintf(sprintf_buffer, "%32s = %d", "Number of bulk layers", PA_FPFS_bulk_count);
				write_log(sprintf_buffer);
				if(PA_FPFS_Numerov){
					write_log((char*)"Error: FPFS_Numerov = true and FPFS_bulk are incompatible");
					status=0; goto FINALIZATION;
				}
			}
			/// FPFS_kRange
			sprintf(sprintf_buffer, "%32s = %.3f", "FPFS_kRange", PA_FPFS_kRange);
			write_log(sprintf_buffer);
			/// FPFS_file
			if(PA_FPFS_file_set){
				sprintf(sprintf_buffer, "%32s = %s", "FPFS_file", PA_FPFS_file);
				write_log(sprintf_buffer);
			}
		}else{
			write_log((char*)"Error: Final_state should be 'PW', 'Calc', 'FP_PAO', or 'FP_AO'");
			//write_log((char*)"Error: Final_state should be 'PW', 'Calc', 'FP_PAO'");
			status=0; goto FINALIZATION;
		}
	}
	/// Polarization
	if(!PA_polarization_set){
		write_log((char*)"Error: Polarization is necessary");
		status=0; goto FINALIZATION;
	}else{
		if(strcmp(PA_polarization, "Linear")==0 || strcmp(PA_polarization, "RCircular")==0 || strcmp(PA_polarization, "LCircular")==0){
			sprintf(sprintf_buffer, "%32s = %s", "Polarization", PA_polarization);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Polarization should be 'Linear', 'RCircular' or 'LCircular'");
			status=0; goto FINALIZATION;
		}
	}
	/// theta and phi
	if(!PA_theta_set || !PA_phi_set){
		write_log((char*)"Error: Theta and Phi are necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = (%8.3f deg, %8.3f deg)", "Polarization angle", PA_theta, PA_phi);
		write_log(sprintf_buffer);
	}
	/// AO_file
	if(!PA_AO_file_set){
		write_log((char*)"Error: Atomic_orbitals_file is necessary");
		status=0; goto FINALIZATION;
	}else{
		sprintf(sprintf_buffer, "%32s = %s", "Atomic_orbitals_file", PA_AO_file);
		write_log(sprintf_buffer);
	}
	/// Extend
	if(PA_ext_set){
		if(PA_ext_up<0 || PA_ext_ri<0 || PA_ext_dn<0 || PA_ext_le<0){
			write_log((char*)"Error: Extend should not be negative");
			status=0; goto FINALIZATION;
		}
		sprintf(sprintf_buffer, "%32s = %d %d %d %d", "Extend", PA_ext_up, PA_ext_ri, PA_ext_dn, PA_ext_le);
		write_log(sprintf_buffer);
	}
	/// Weighting
	sprintf(sprintf_buffer, "%32s = %s", "Weighting", PA_weighting ? "true" : "false");
	write_log(sprintf_buffer);
	if(PA_weighting){
		const char* unit=PA_use_angstrom ? "Ang" : "a.u.";
		/// Weighting_axis
		if(PA_weighting_axis_set){
			sprintf(sprintf_buffer, "%32s = (%8.2f, %8.2f, %8.2f) [%s]", "Weighting_axis", PA_weighting_axis[0], PA_weighting_axis[1], PA_weighting_axis[2], unit);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Weighting_axis is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_shape
		if(PA_weighting_shape_set){
			if(strcmp(PA_weighting_shape, "Rect")==0 || strcmp(PA_weighting_shape, "Exp")==0 || strcmp(PA_weighting_shape, "Tri")==0 || strcmp(PA_weighting_shape, "Sqrt")==0){
				sprintf(sprintf_buffer, "%32s = %s", "Weighting_shape", PA_weighting_shape);
				write_log(sprintf_buffer);
			}else{
				write_log((char*)"Error: Weighting_shape should be Rect, Exp, Tri, or Sqrt");
				status=0; goto FINALIZATION;
			}
		}else{
			write_log((char*)"Error: Weighting_shape is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_origin
		if(PA_weighting_origin_set){
			if(PA_weighting_origin>=0){
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], parallel to the axis", "Distance of the origin", PA_weighting_origin, unit);
				write_log(sprintf_buffer);
			}else{
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], anti-parallel to the axis", "Distance of the origin", PA_weighting_origin, unit);
				write_log(sprintf_buffer);
			}
		}else{
			write_log((char*)"Error: Weighting_origin is necessary");
			status=0; goto FINALIZATION;
		}
		/// Weighting_width
		if(PA_weighting_width_set){
			if(PA_weighting_width>=0){
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], parallel to the axis", "Weighting_width", PA_weighting_width, unit);
				write_log(sprintf_buffer);
			}else{
				sprintf(sprintf_buffer, "%32s = %8.2f [%s], anti-parallel to the axis", "Weighting_width", PA_weighting_width, unit);
				write_log(sprintf_buffer);
			}
		}else{
			write_log((char*)"Error: Weighting_width is necessary");
			status=0; goto FINALIZATION;
		}
		/// Include_neg_depth
		sprintf(sprintf_buffer, "%32s = %s", "Include_neg_depth", PA_include_neg_depth ? "true": "false");
		write_log(sprintf_buffer);
	}
	/// Reflection
	if(PA_reflection){
		if(!PA_weighting){
			write_log((char*)"Error: Weighting should be true if Reflection is true");
			status=0; goto FINALIZATION;
		}
		if(PA_reflection_coef_set){
			sprintf(sprintf_buffer, "%32s = %8.2f", "Reflection coefficient", PA_reflection_coef);
			write_log(sprintf_buffer);
		}else{
			write_log((char*)"Error: Reflection_coef is necessary if Reflection is true");
			status=0; goto FINALIZATION;
		}
	}
	/// Ignore_core
	sprintf(sprintf_buffer, "%32s = %s", "Ignore_core", PA_ignore_core?"true":"false");
	write_log(sprintf_buffer);
	if(PA_ignore_core && !PA_VPS_file_set){
		write_log((char*)"Error: Ignoreing core needs the VPS file");
		status=0; goto FINALIZATION;
	}
	/// Add_nonorth_term
	sprintf(sprintf_buffer, "%32s = %s", "Add_nonorth_term", PA_add_nonorth_term?"true":"false");
	write_log(sprintf_buffer);
	/// Orth_correction
	sprintf(sprintf_buffer, "%32s = %s", "Orth_correction", PA_orth_correction?"true":"false");
	write_log(sprintf_buffer);
	if(PA_FPFS){
		/// Calc_all_nonloc
		sprintf(sprintf_buffer, "%32s = %s", "Calc_all_nonloc", PA_calc_all_nonloc?"true":"false");
		write_log(sprintf_buffer);
	}
	/// 
	/// Output_data
	sprintf(sprintf_buffer, "%32s = %s", "Output_data", PA_output_data);
	write_log(sprintf_buffer);
	if(strcmp(PA_output_data, "PAD")!=0 && strcmp(PA_output_data, "Band")!=0 /*&& strcmp(PA_output_data, "PAD_real")!=0 && strcmp(PA_output_data, "PAD_imag")!=0*/){
		write_log((char*)"Error: Invalid Output_data");
		status=0; goto FINALIZATION;
	}
	goto FINALIZATION;
 FINALIZATION:
	delete[] sprintf_buffer;
	return status;
}
