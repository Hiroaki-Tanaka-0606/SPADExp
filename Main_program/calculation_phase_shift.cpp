// calcPSF phase shift calculation

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include <H5Cpp.h>

#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"
#include "calculation_atomic_wfn.hpp"
#include "HDF5_tools.hpp"

using namespace std;

void calc_phase_shift(){
	char* sprintf_buffer=new char[Log_length+1];
	int i,j,x;

	write_log((char*)"----Phase shift calculation----");
	int Z=Z_single;
	double mu=pow(3.0*M_PI/4.0, 2.0/3.0)/(2.0*pow(Z, 1.0/3.0));

	sprintf(sprintf_buffer, "TF scaling: mu=%10.5e", mu);
	write_log(sprintf_buffer);

	write_log((char*)"Load potential from the database");
	H5File database(At_potential_file, H5F_ACC_RDONLY);

	char group_name[4];
	sprintf(group_name, "%03d", Z);
	Group atomG(database.openGroup(group_name));

	// prepare potential by interpolation
	int x_count2=r_att_int(atomG, "length");
	double x_coordinates2[x_count2];
	double v_x2[x_count2];
	r_data_1d(atomG, "Potential", x_count2, &v_x2[0]);
	r_att_1d(atomG, "x", x_count2, &x_coordinates2[0]);

	int mod_start_index=r_att_int(atomG, "mod_start_index");
	double mod_start_x=x_coordinates2[mod_start_index];

	At_v_x[0]=v_x2[0];
	for(i=1; i<x_count; i++){
		for(j=0; j<x_count2-1; j++){
			if(x_coordinates2[j]<x_coordinates[i] && x_coordinates[i]<x_coordinates2[j+1]){
				double dx1=x_coordinates[i]-x_coordinates2[j];
				double dx2=x_coordinates2[j+1]-x_coordinates[i];
				At_v_x[i]=(v_x2[j]*dx2+v_x2[j+1]*dx1)/(dx1+dx2);
				break;
			}else if(abs(x_coordinates2[j+1]-x_coordinates[i])<Data_read_error){
				At_v_x[i]=v_x2[j+1];
				break;
			}
		}
		// if the potential is not in database, use V(r)=-1/r=-1/mu*x
		At_v_x[i]=-1.0/(mu*x_coordinates[i]);
	}
	
	/*
	for(i=0; i<x_count; i++){
		At_v_x[i]=0;
		}*/
	
	sprintf(sprintf_buffer, "Modified potential starts from i = %d, x = %8.3f", mod_start_index, mod_start_x);
	write_log(sprintf_buffer);
	
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

	for(i=0; i<Ph_orbital_count; i++){
		sprintf(sprintf_buffer, "%s orbital: l = %1d, EB = %8.2f", Ph_orbital_labels[i], Ph_l_list[i], Ph_binding_energies[i]);
		write_log(sprintf_buffer);
		for(j=0; j<Ex_energy_count; j++){
			double Ekin=Ex_energies[j]-Ph_binding_energies[i];
			double k_length=sqrt(2.0*Ekin);
			sprintf(sprintf_buffer, "Excitation energy = %8.2f, Kinetic energy = %8.2f, k = %8.3f", Ex_energies[j], Ekin, k_length);
			write_log(sprintf_buffer);
			if(Ekin<0){
				continue;
			}
			int dl;
			for(dl=-1; dl<1; dl+=2){
				int lp=dl+Ph_l_list[i];
				if(lp<0){
					continue;
				}						
				sprintf(sprintf_buffer, "Calculation for lp = %1d", lp);
				write_log(sprintf_buffer);
				Atomic_wfn_evolution(mu, lp, Ekin, false);
				
				
				for(x=0; x<x_count; x++){
					sprintf(sprintf_buffer, "%8.3f %8.3f", mu*x_coordinates[x], At_p_x[x]);
					write_log(sprintf_buffer);
				}return;
				
				bool increasing;
				int local_top_count=0;
				double r;
				double n_calc;
				double delta;
				// p(r) ~ sin(phase)
				// phase = kr - 1/k log(2kr) - lp*pi/2 + delta
				// phase/2pi = {kr - 1/k log(2kr)}/pi - lp/4 + delta/2pi
				for(x=mod_start_index; x<x_count-1; x++){
					if(x==mod_start_index){
						increasing=(At_p_x[x+1]>At_p_x[x]);
						continue;
					}
					if(increasing){
						// search local maximum
						// phase = (2n+0.5)*pi
						// phase/2pi = n + 0.25
						if(At_p_x[x+1]<At_p_x[x]){
							sprintf(sprintf_buffer, "Local maximum at i = %d, x = %8.3f", x, x_coordinates[x]);
							write_log(sprintf_buffer);
							increasing=false;
							local_top_count++;
							if(local_top_count>Ph_skip_points){
								r=mu*x_coordinates[x];
								n_calc=(k_length*r - log(2.0*k_length*r)/k_length)/(2.0*M_PI)-lp/4.0-0.25;
								//n_calc=(k_length*r)/(2.0*M_PI)-lp/4.0-0.25;
								delta=2.0*M_PI*(n_calc-floor(n_calc));
								sprintf(sprintf_buffer, "delta = %8.3f", delta);
								//write_log(sprintf_buffer);
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
							sprintf(sprintf_buffer, "Local minimum at i = %d, x = %8.3f", x, x_coordinates[x]);
							//write_log(sprintf_buffer);
							increasing=true;
							local_top_count++;
							if(local_top_count>Ph_skip_points){
								r=mu*x_coordinates[x];
								n_calc=(k_length*r - log(2.0*k_length*r)/k_length)/(2.0*M_PI)-lp/4.0+0.25;
								//n_calc=(k_length*r)/(2.0*M_PI)-lp/4.0+0.25;
								delta=2.0*M_PI*(n_calc-floor(n_calc));
								sprintf(sprintf_buffer, "delta = %8.3f", delta);
								//write_log(sprintf_buffer);
							}
							if(local_top_count>=Ph_calc_points+Ph_skip_points){
								break;
							}
						}
					}
				}
				
				return;
			}
		}
	}
	return;
}
