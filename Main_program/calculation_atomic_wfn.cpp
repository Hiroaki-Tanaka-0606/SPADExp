// calcPSF calculation of atomic wavefunction

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>

#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"

void Atomic_wfn_evolution(double mu, int l, double E, bool reverse);
double calc_atomic_wfn(double mu, int l, int nodes);
int count_nodes(int i_max);
int matching_radius_index(double mu, int l, double E);
double calc_second_deriv_coeff(double mu, int l, double E, int index);
double calc_atomic_wfn2(double mu, int l, int nodes, double eigen_guess);

using namespace std;

void sequence_atomic_wfn(){
	int Z_current;
	int n_current;
	int l_current;
	int i;
	int setup_result;
	char* sprintf_buffer=new char[Log_length+1];
	double mu;
	int nodes_min, nodes_max, nodes_current;
	double E_n_l;
	write_log((char*)"----Sequential calculation of atomic wavefunctions----");

	fprintf(Output_file_obj, "# calcPSF\n");
	fprintf(Output_file_obj, "# Sequential calculation of atomic wavefunctions\n");
	fprintf(Output_file_obj, "# Z ");
	for(n_current=n_min; n_current<=n_max; n_current++){
		for(l_current=l_min; l_current<=l_max, l_current<=n_current-1; l_current++){
			fprintf(Output_file_obj, "(%d,%d) ", n_current, l_current);
		}
	}
	fprintf(Output_file_obj, "\n");
				
	for(Z_current=Z_min; Z_current<=Z_max; Z_current++){
		sprintf(sprintf_buffer, "---- Z = %3d ----", Z_current);
		write_log(sprintf_buffer);
		mu=pow(3.0*M_PI/4.0, 2.0/3.0)/(2.0*pow(Z_current, 1.0/3.0));
		sprintf(sprintf_buffer, "TF scaling: mu = %10.5e", mu);
		write_log(sprintf_buffer);

		fprintf(Output_file_obj, "%3d ", Z_current);
		// setup the potential
		setup_result=setup_potential(Z_current, mu);
		if(setup_result!=1){
			return;
		}
		modify_potential(Z_current, mu);
		/*
		 for(i=0; i<x_count; i++){
			sprintf(sprintf_buffer, "v(%8.3f) = %10.5e", x_coordinates[i], At_v_x[i]);
		 	write_log(sprintf_buffer);
			}*/

		for(n_current=n_min; n_current<=n_max; n_current++){
			for(l_current=l_min; (l_current<=l_max && l_current<=n_current-1); l_current++){
				
				sprintf(sprintf_buffer, "---- n = %3d, l = %3d ( %3d nodes ) ----", n_current, l_current, n_current-l_current-1);
				write_log(sprintf_buffer);
				E_n_l=calc_atomic_wfn(mu, l_current, n_current-l_current-1);
				fprintf(Output_file_obj, "%12.5e ", E_n_l);

				/* debug: wfn */
				for(i=0; i<x_count; i++){
					sprintf(sprintf_buffer, "%10.5e %10.5e", x_coordinates[i], At_p_x[i]);
					write_log(sprintf_buffer);
				}
			}
		}
		fprintf(Output_file_obj, "\n");
	}
	delete sprintf_buffer;
}

double calc_atomic_wfn(double mu, int l, int nodes){
	char* sprintf_buffer=new char[Log_length+1];
  
	double current_step=Bisection_step;

	double E_max=0.0;
	double E_min=-current_step;
	double E_center;
	int matching_index;
	int num_nodes;
	write_log((char*)"Search adequate E_initial (trial)");
	int current_trials=1;
	/*
	// for debug
	for(E_min=-0.01; E_min>-1; E_min-=0.01){
		matching_index=matching_radius_index(mu, l, E_min);
		Atomic_wfn_evolution(mu, l, E_min, false);
		num_nodes=count_nodes(matching_index);
		sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_min, num_nodes);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "matching_index = %3d ( x = %8.3f )", matching_index, x_coordinates[matching_index]);
		write_log(sprintf_buffer);
		for(int i=0; i<x_count; i++){
			sprintf(sprintf_buffer, "%8.3f %12.5e", x_coordinates[i], At_p_x[i]);
			//write_log(sprintf_buffer);
		}
		//break;
	}
	return 0;*/

	E_max=Bisection_step*floor((l*(l+1)*1.0/(2*pow(mu*x_coordinates[x_count-1], 2))+At_v_x[x_count-1])/Bisection_step);
	E_min=E_max;
	
	// first calculation (outward)
	matching_index=matching_radius_index(mu, l, E_min);
  Atomic_wfn_evolution(mu, l, E_min, false);
	num_nodes=count_nodes(matching_index);
	sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_min, num_nodes);
	write_log(sprintf_buffer);
	// find the energy which evolves the wfn with nodes <= nodes_min
	while(num_nodes>nodes){
		E_max=E_min;
		E_min-=current_step;
		current_step*=2;
		current_trials++;
		matching_index=matching_radius_index(mu, l, E_min);
		Atomic_wfn_evolution(mu, l, E_min, false);
		num_nodes=count_nodes(matching_index);
		//sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_min, num_nodes);
		//write_log(sprintf_buffer);
		//sprintf(sprintf_buffer, "matching_index = %3d ( x = %8.3f )", matching_index, x_coordinates[matching_index]);
		//write_log(sprintf_buffer);
	}
	sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_min, num_nodes);
	write_log(sprintf_buffer);
	if(num_nodes==nodes){
		E_center=E_min;
	}
	write_log((char*)"Search adequate E_initial (bisection)");
	// find the energy which evolves the wfn with nodes = nodes_min (bisection method)
	while(num_nodes!=nodes || E_min/E_max>At_bisection_threshold){
		E_center=(E_max+E_min)/2.0;
		current_trials++;
		matching_index=matching_radius_index(mu, l, E_center);
		Atomic_wfn_evolution(mu, l, E_center, false);
		num_nodes=count_nodes(matching_index);
		//sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_center, num_nodes);
		//write_log(sprintf_buffer);
		//sprintf(sprintf_buffer, "matching_index = %3d ( x = %8.3f )", matching_index, x_coordinates[matching_index]);
		//write_log(sprintf_buffer);
		if(num_nodes<=nodes){
			E_min=E_center;
		}else{
			E_max=E_center;
		}
		
	}
	sprintf(sprintf_buffer, "Step %3d: E = %12.5e, nodes = %3d", current_trials, E_center, num_nodes);
	write_log(sprintf_buffer);
	
	return calc_atomic_wfn2(mu, l, nodes, E_center);
}

double calc_atomic_wfn2(double mu, int l, int nodes, double eigen_guess){
	char* sprintf_buffer=new char[Log_length+1];
	// integration to the outward
	double E_center=eigen_guess;
	int matching_index=matching_radius_index(mu, l, E_center);
  Atomic_wfn_evolution(mu, l, E_center, false);
	int num_nodes=count_nodes(matching_index);

	
	write_log((char*)"Search the eigenvalue by perturbation");
	int current_trials=1;

	double E_diff=At_E_threshold*2;
	double log_diff_outward, log_diff_inward, int_outward, int_inward;
	double interval_outward, interval_inward;
	double delta_E;
	int i;
	while((current_trials<At_min_iteration) || (E_diff>At_E_threshold && current_trials<At_max_iteration)){
		// at this time the outward evolution is in At_p_x
		interval_outward=x_coordinates[matching_index]-x_coordinates[matching_index-1];
		interval_inward=x_coordinates[matching_index+1]-x_coordinates[matching_index];

		// outward properties
		log_diff_outward=(At_p_x[matching_index]-At_p_x[matching_index-1])/(interval_outward*At_p_x[matching_index]);
		int_outward=0.0;
		for(i=0; i<matching_index; i++){
			int_outward+=pow(At_p_x[i], 2)*(x_coordinates[i+1]-x_coordinates[i]);
		}
		int_outward/=pow(At_p_x[matching_index], 2);

		// inward integration
		Atomic_wfn_evolution(mu, l, E_center, true);
		// inward properties
		log_diff_inward=(At_p_x[matching_index+1]-At_p_x[matching_index])/(interval_inward*At_p_x[matching_index]);
		int_inward=0.0;
		for(i=matching_index; i<x_count-1; i++){
			int_inward+=pow(At_p_x[i], 2)*(x_coordinates[i+1]-x_coordinates[i]);
		}
		int_inward/=pow(At_p_x[matching_index], 2);

		delta_E=(log_diff_outward-log_diff_inward)/(int_outward+int_inward);
		if(!isfinite(delta_E)){
			write_log((char*)"Warning: delta_E is not finite");
			delta_E=(log_diff_outward-log_diff_inward);
		}
		
		E_diff=fabs(delta_E/E_center);
		//sprintf(sprintf_buffer, "%12.5e %12.5e %12.5e %12.5e", log_diff_outward, log_diff_inward, int_outward, int_inward);
		//write_log(sprintf_buffer);
		//sprintf(sprintf_buffer, "Step %3d: E = %12.5e, dE = %12.5e, E_diff = %12.5e", current_trials, E_center, delta_E, E_diff);
		//write_log(sprintf_buffer);
		//sprintf(sprintf_buffer, "matching_index = %3d ( x = %8.3f )", matching_index, x_coordinates[matching_index]);
		//write_log(sprintf_buffer);

		// check if the next E_center has the same number of nodes
		matching_index=matching_radius_index(mu, l, E_center+delta_E);
		Atomic_wfn_evolution(mu, l, E_center+delta_E, false);
		num_nodes=count_nodes(matching_index);
		while(num_nodes!=nodes){
		  sprintf(sprintf_buffer, "Warning: E+dE (%12.5e) has %3d nodes, but the target is %3d", E_center+delta_E, num_nodes, nodes);
			write_log(sprintf_buffer);
			delta_E/=2.0;
			matching_index=matching_radius_index(mu, l, E_center+delta_E);
			Atomic_wfn_evolution(mu, l, E_center+delta_E, false);
			num_nodes=count_nodes(matching_index);
		}
		E_center+=delta_E;
		current_trials++;
	}
	sprintf(sprintf_buffer, "Step %3d: E = %12.5e, dE = %12.5e, E_diff = %12.5e", current_trials, E_center, delta_E, E_diff);
	write_log(sprintf_buffer);
	
	// obtain the solution
	/// copy the outward integration	
	double* At_p_outward=new double[x_count];
	for(i=0; i<x_count; i++){
		At_p_outward[i]=At_p_x[i];
	}
	/// inward integration
	matching_index=matching_radius_index(mu, l, E_center);
	Atomic_wfn_evolution(mu, l, E_center, true);	
	/// continuity condition at matching_index
	double coef=At_p_x[matching_index]/At_p_outward[matching_index];
	for(i=matching_index; i<x_count; i++){
		At_p_x[i]/=coef;
	}
	/// normalization w.r.t. r (not x)
	double p_int=0.0;
	for(i=0; i<matching_index; i++){
		p_int+=pow(At_p_outward[i], 2)*(x_coordinates[i+1]-x_coordinates[i]);
	}
	for(i=matching_index; i<x_count; i++){
		p_int+=pow(At_p_x[i], 2)*(x_coordinates[i+1]-x_coordinates[i]);
	}
	coef=sqrt(p_int*mu);
	for(i=0; i<matching_index; i++){
		At_p_x[i]=At_p_outward[i]/coef;
	}
	for(i=matching_index; i<x_count; i++){
		At_p_x[i]/=coef;
	}

	sprintf(sprintf_buffer, "E = %12.5e, matching_index = %3d (x = %8.3f)", E_center, matching_index, x_coordinates[matching_index]);
	write_log(sprintf_buffer);
	/*
	for(i=0; i<x_count; i++){
		sprintf(sprintf_buffer, "%8.3f %12.5e", x_coordinates[i], At_p_x[i]);
		write_log(sprintf_buffer);
		}*/
										 
	

  
	delete sprintf_buffer;
	return E_center;
}

int count_nodes(int i_max){
	int num_nodes=0;
	for(int i=2; i<=i_max; i++){
		if(At_p_x[i]*At_p_x[i-1]<0){
			num_nodes++;
		}
	}
	return num_nodes;
}
										 

int matching_radius_index(double mu, int l, double E){
	int i=x_count-1;
	while(calc_second_deriv_coeff(mu, l, E, i)>0 && i>=0){
		i--;
	}
	if(i>=x_count-2){
		return x_count-2;
	}else{
		return i+1;
	}
}

// Schroedinger differential equation of radial wave function
// (p, p_diff)'=f(p, p_diff)=(p_diff, (l(l+1)/x^2+2*mu^2(v(x)-E))*p)
// 
double calc_second_deriv_coeff(double mu, int l, double E, int index){
	return l*(l+1)*1.0/pow(x_coordinates[index], 2)+2*pow(mu, 2)*(At_v_x[index]-E);
	//return E;
}

void Atomic_wfn_evolution(double mu, int l, double E, bool reverse){
	// calculation range
	int last_grid=0;
	double last_x=Radial_grid_intervals[0]*Radial_grid_points[0];
	int last_x_index=Radial_grid_points[0];
	
	int matching_index=matching_radius_index(mu, l, E);
	double calc_last_x=x_coordinates[matching_index]*At_radius_factor;
	while(last_x<calc_last_x && last_grid < Radial_grid_count){
		last_grid++;
		last_x+=Radial_grid_intervals[last_grid]*Radial_grid_points[last_grid];
		last_x_index+=Radial_grid_points[last_grid];
	}	
	// initialize
	int i, j;
	for(i=0; i<x_count; i++){
		At_p_x[i]=0;
		At_p_diff_x[i]=0;
	}
	if(!reverse){
		// outward: 0 to inf
		At_p_x[0]=0;
		At_p_diff_x[0]=At_initial_diff;
	}else{
		// inward: inf to 0
		At_p_x[last_x_index]=0;
		At_p_diff_x[last_x_index]=At_initial_diff;
	}
	
	
	if(strcmp(At_solution, "RK1")==0){
		// RK1: Euler method
		if(!reverse){
			// outward: 0 to inf
			for(i=1; i<=last_x_index; i++){
				At_p_x[i]=At_p_x[i-1]+At_p_diff_x[i-1]*(x_coordinates[i]-x_coordinates[i-1]);
				if(i==1){
					// to avoid the division by zero
					At_p_diff_x[i]=At_p_diff_x[i-1];
				}else{
					At_p_diff_x[i]=At_p_diff_x[i-1]+At_p_x[i-1]*calc_second_deriv_coeff(mu, l, E, i-1)*(x_coordinates[i]-x_coordinates[i-1]);
				}
			}
		}else{
			// inward: inf to 0
			for(i=last_x_index; i>0; i--){
				At_p_x[i-1]=At_p_x[i]-At_p_diff_x[i]*(x_coordinates[i]-x_coordinates[i-1]);
				At_p_diff_x[i-1]=At_p_diff_x[i]-At_p_x[i]*calc_second_deriv_coeff(mu, l, E, i)*(x_coordinates[i]-x_coordinates[i-1]);		
			}
		}
	}else if(strcmp(At_solution, "Numerov")==0){
		// Numerov method
		double interval;
		int current_index;
		double x;
		double axi;
		double axim1;
		double axip1;
		if(!reverse){
			// outward: 0 to inf
			current_index=0;
			for(i=0; i<=last_grid; i++){
				interval=Radial_grid_intervals[i];
				for(j=0; j<Radial_grid_points[i]; j++){
					x=x_coordinates[current_index];
					if(j==0){
						At_p_x[current_index+1]=At_p_x[current_index]+At_p_diff_x[current_index]*interval;
						if(i==0){
							At_p_diff_x[current_index+1]=At_p_diff_x[current_index];
						}else{
							At_p_diff_x[current_index+1]=At_p_diff_x[current_index]+At_p_x[current_index]*calc_second_deriv_coeff(mu, l, E, current_index)*interval;
						}
						current_index++;
						continue;
					}
					if(j==1 && i==0){
						At_p_x[current_index+1]=At_p_x[current_index]+At_p_diff_x[current_index]*interval;
						At_p_diff_x[current_index+1]=At_p_diff_x[current_index]+At_p_x[current_index]*calc_second_deriv_coeff(mu, l, E, current_index)*interval;
						current_index++;
						continue;
					}
					axi=-calc_second_deriv_coeff(mu, l, E, current_index);
					axip1=-calc_second_deriv_coeff(mu, l, E, current_index+1);
					axim1=-calc_second_deriv_coeff(mu, l, E, current_index-1);

					At_p_x[current_index+1]=(2*(1.0-pow(interval, 2)*5.0/12.0*axi)*At_p_x[current_index]-(1.0+pow(interval, 2)/12.0*axim1)*At_p_x[current_index-1])/(1.0+pow(interval, 2)/12.0*axip1);
					At_p_diff_x[current_index+1]=(At_p_x[current_index+1]-At_p_x[current_index])/interval;
					
					current_index++;
				}
			}
		}else{
			// inward: inf to 0
			current_index=last_x_index;
			for(i=last_grid; i>=0; i--){
				interval=Radial_grid_intervals[i];
				for(j=0; j<Radial_grid_points[i]; j++){
					x=x_coordinates[current_index];
					if(j==0 || (i==0 && j==Radial_grid_points[i]-1)){
						At_p_x[current_index-1]=At_p_x[current_index]-At_p_diff_x[current_index]*interval;
						At_p_diff_x[current_index-1]=At_p_diff_x[current_index]-At_p_x[current_index]*calc_second_deriv_coeff(mu, l, E, current_index)*interval;
						current_index--;
						continue;
					}
					axi=-calc_second_deriv_coeff(mu, l, E, current_index);
					axip1=-calc_second_deriv_coeff(mu, l, E, current_index-1);
					axim1=-calc_second_deriv_coeff(mu, l, E, current_index+1);

					At_p_x[current_index-1]=(2*(1.0-pow(interval, 2)*5.0/12.0*axi)*At_p_x[current_index]-(1.0+pow(interval, 2)/12.0*axim1)*At_p_x[current_index+1])/(1.0+pow(interval, 2)/12.0*axip1);
					At_p_diff_x[current_index-1]=(At_p_x[current_index]-At_p_x[current_index-1])/interval;
					
					current_index--;
				}
			}
		}
	}	  
}
