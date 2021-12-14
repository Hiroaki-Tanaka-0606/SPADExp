// calcPSF phase shift calculation

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include <H5Cpp.h>
#include <complex>

#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"
#include "calculation_atomic_wfn.hpp"
#include "HDF5_tools.hpp"

using namespace std;

double analyze_phase(double mu, double k_length, int lp, double mod_start_x);

void calc_phase_shift(){
	char* sprintf_buffer=new char[Log_length+1];
	int i,j,x;

	write_log((char*)"----Phase shift calculation----");
	int Z=Z_single;
	double mu=pow(3.0*M_PI/4.0, 2.0/3.0)/(2.0*pow(Z, 1.0/3.0));

	sprintf(sprintf_buffer, "TF scaling: mu=%10.5e", mu);
	write_log(sprintf_buffer);

	double mod_start_x=0;
	if(At_potential_file_set){
		write_log((char*)"Load potential from the database");
		H5File database(At_potential_file, H5F_ACC_RDONLY);

		char group_name[4];
		sprintf(group_name, "%03d", Z);
		Group atomG(database.openGroup(group_name));
		
		// prepare potential by interpolation
		mod_start_x=load_potential_H5(atomG, mu);

		sprintf(sprintf_buffer, "Modified potential starts from x = %8.3f", mod_start_x);
		write_log(sprintf_buffer);
	}else{
	  write_log((char*)"Potential is -1/r");
		for(i=0; i<x_count; i++){
			At_v_x[i]=-1.0/(mu*x_coordinates[i]);
		}
	}
	
	
	/*
		for(i=0; i<x_count; i++){
		sprintf(sprintf_buffer, "%8.3f %10.3e", x_coordinates[i], At_v_x[i]);
		write_log(sprintf_buffer);
		}
		cout << endl;
		for(i=0; i<x_count2; i++){
		sprintf(sprintf_buffer, "%8.3f %8.3f", x_coordinates2[i], v_x2[i]);
		write_log(sprintf_buffer);
		}
		return;*/

	fprintf(Output_file_obj, "# calcPSF\n");
	fprintf(Output_file_obj, "# Phase-shift calculation\n");
	fprintf(Output_file_obj, "# Z = %d\n", Z);
	fprintf(Output_file_obj, "\n");
	fprintf(Output_file_obj, "#    lp ");
	for(i=0; i<Ex_energy_count; i++){
		fprintf(Output_file_obj, "%6.1f ", Ex_energies[i]*Eh);
	}
	fprintf(Output_file_obj, " (eV)\n");
	
	double p_x_list[Ph_orbital_count][Ex_energy_count][2][x_count];
	
	for(i=0; i<Ph_orbital_count; i++){
		sprintf(sprintf_buffer, "%s orbital: l = %1d, EB = %8.2f", Ph_orbital_labels[i], Ph_l_list[i], Ph_binding_energies[i]);
		write_log(sprintf_buffer);
		int dl;
		for(dl=-1; dl<=1; dl+=2){
			int lp=dl+Ph_l_list[i];
			if(lp<0){
				continue;
			}						
			fprintf(Output_file_obj, "# %2s %2d ", Ph_orbital_labels[i], lp);
			for(j=0; j<Ex_energy_count; j++){
				double Ekin=Ex_energies[j]-Ph_binding_energies[i];
				double k_length=sqrt(2.0*Ekin);
				sprintf(sprintf_buffer, "Excitation energy = %8.2f, Kinetic energy = %8.2f, k = %8.3f", Ex_energies[j], Ekin, k_length);
				write_log(sprintf_buffer);
				if(Ekin<0){
					fprintf(Output_file_obj, "------ ");
					continue;
				}
				sprintf(sprintf_buffer, "Calculation for lp = %1d", lp);
				write_log(sprintf_buffer);
				Atomic_wfn_evolution(mu, lp, Ekin, false);

				/*
				for(x=0; x<x_count; x++){
					sprintf(sprintf_buffer, "%8.3f %8.3f", mu*x_coordinates[x], At_p_x[x]);
					write_log(sprintf_buffer);
					}return;*/
				double phase_shift=analyze_phase(mu, k_length, lp, mod_start_x);
				fprintf(Output_file_obj, "%6.4f ", phase_shift);
				sprintf(sprintf_buffer, "Averaged phase shift = %6.4f", phase_shift);
				write_log(sprintf_buffer);
				for(x=0; x<x_count; x++){
					p_x_list[i][j][(dl+1)/2][x]=At_p_x[x];
				}
			}
			fprintf(Output_file_obj, "\n");
		}
	}
	fprintf(Output_file_obj, "\n");

	// output of the wavefunctions
	for(j=0; j<Ex_energy_count; j++){
		fprintf(Output_file_obj, "# Excitation energy = %8.2f eV (%8.2f Eh)\n", Ex_energies[j]*Eh, Ex_energies[j]);
		fprintf(Output_file_obj, "# r ");
		for(i=0; i<Ph_orbital_count; i++){
			fprintf(Output_file_obj, "%s- %s+ ", Ph_orbital_labels[i], Ph_orbital_labels[i]);
		}
		fprintf(Output_file_obj, "\n");
		for(x=0; x<x_count; x++){
			fprintf(Output_file_obj, "%10.4f ", x_coordinates[x]*mu);
			for(i=0; i<Ph_orbital_count; i++){
				fprintf(Output_file_obj, "%10.3e %10.3e ", p_x_list[i][j][0][x], p_x_list[i][j][1][x]);
			}
			fprintf(Output_file_obj, "\n");
		}
		fprintf(Output_file_obj, "\n");
	}
	return;
}

double analyze_phase(double mu, double k_length, int lp, double mod_start_x){
	char* sprintf_buffer=new char[Log_length+1];
	bool increasing;
	bool before_mod=true;
	int local_top_count=0;
	double r;
	double n_calc;
	double delta;
	complex<double> phase_sum(0, 0);
	int calc_count=0;
	double amplitude_sum=0;
	int x;
	// p(r) ~ sin(phase)
	// phase = kr + 1/k log(2kr) - lp*pi/2 + delta
	// phase/2pi = {kr + 1/k log(2kr)}/2pi - lp/4 + delta/2pi
	// n_calc=n-delta/2pi
	// n=ceil(n_calc)
	// delta=2pi*(n-n_calc)
	for(x=0; x<x_count-1; x++){
		if(x_coordinates[x]<mod_start_x){
			continue;
		}
		if(before_mod){
			increasing=(At_p_x[x+1]>At_p_x[x]);
			before_mod=false;
			continue;
		}
		if(increasing){
			// search local maximum
			// phase = (2n+0.5)*pi
			// phase/2pi = n + 0.25
			if(At_p_x[x+1]<At_p_x[x]){
				// sprintf(sprintf_buffer, "Local maximum at i = %5d, x = %8.3f, P(r) = %10.3e", x, x_coordinates[x], At_p_x[x]);
				// write_log(sprintf_buffer);
				increasing=false;
				local_top_count++;
				if(local_top_count>Ph_skip_points){
					r=mu*x_coordinates[x];
					n_calc=(k_length*r + log(2.0*k_length*r)/k_length)/(2.0*M_PI)-lp/4.0-0.25;
					//n_calc=(k_length*r)/(2.0*M_PI)-lp/4.0-0.25;
					delta=2.0*M_PI*(ceil(n_calc)-n_calc);
					// sprintf(sprintf_buffer, "delta = %8.3f", delta);
					// write_log(sprintf_buffer);
					phase_sum+=complex<double>(cos(delta), sin(delta));
					amplitude_sum+=At_p_x[x];
					calc_count++;
				}
				if(local_top_count>=Ph_calc_points+Ph_skip_points){
					break;
				}
			}
		}else{
			// search local minimum
			// phase = (2n-0.5)*pi
			// phase/2pi = n - 0.25
			if(At_p_x[x+1]>At_p_x[x]){
				// sprintf(sprintf_buffer, "Local minimum at i = %5d, x = %8.3f, P(r) = %10.3e", x, x_coordinates[x], At_p_x[x]);
				// write_log(sprintf_buffer);
				increasing=true;
				local_top_count++;
				if(local_top_count>Ph_skip_points){
					r=mu*x_coordinates[x];
					n_calc=(k_length*r + log(2.0*k_length*r)/k_length)/(2.0*M_PI)-lp/4.0+0.25;
					//n_calc=(k_length*r)/(2.0*M_PI)-lp/4.0+0.25;
					delta=2.0*M_PI*(ceil(n_calc)-n_calc);
					// sprintf(sprintf_buffer, "delta = %8.3f", delta);
					// write_log(sprintf_buffer);
					phase_sum+=complex<double>(cos(delta), sin(delta));
					amplitude_sum+=(-At_p_x[x]);
					calc_count++;
				}
				if(local_top_count>=Ph_calc_points+Ph_skip_points){
					break;
				}
			}
		}
	}
	complex<double> phase_average;
	double amplitude_average=amplitude_sum/(calc_count*1.0);
	phase_average=phase_sum/(calc_count*1.0);
	double phase_output=arg(phase_average);
	if(phase_output<0){
		phase_output+=2*M_PI;
	}

	// normalize so that P(r)=1/k sin(kr+..)
	for(x=0; x<x_count; x++){
		At_p_x[x]/=(amplitude_average*k_length);		
	}
	return phase_output;
}
