// calcPSF SCF calculation of an atom

#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>

#include "variables_ext.hpp"
#include "log.hpp"
#include "setup.hpp"
#include "calculation_atomic_wfn.hpp"

using namespace std;

void calc_atom(double mu, double* v_x, bool first_calculation, bool use_guess);
void calc_potential(int Z, double mu, double* v_x_unmod, double* v_x);
void mix_potential(double* v_x, double* v_x_next);
double calc_criterion_a(double* v1_X, double* v2_x);
double calc_criterion_b(double* v1_x, double* v2_x, double mu);

void scf_calc_atom(){
	char* sprintf_buffer=new char[Log_length+1];
	int i,j;

	write_log((char*)"----SCF calculation of an atom----");
	int Z=Z_single;
	double mu=pow(3.0*M_PI/4.0, 2.0/3.0)/(2.0*pow(Z, 1.0/3.0));
	double* At_v_x_next=new double[x_count];
	// before modification (v_x=-1/r if (0>)v_x_unmod>-1/r)
	double* At_v_x_unmod=new double[x_count];
	double* At_v_x_next_unmod=new double[x_count];
		
	sprintf(sprintf_buffer, "TF scaling: mu = %10.5e", mu);
	write_log(sprintf_buffer);

	// first calculation: use TF potential
	int setup_result=setup_potential(Z, mu);
	if(setup_result!=1){
		return;
	}
	// copy the unmodified potential
	for(i=0; i<x_count; i++){
		At_v_x_unmod[i]=At_v_x[i];
	}
	modify_potential(Z, mu);
	/*
	for(i=0; i<x_count; i++){
		sprintf(sprintf_buffer, "%8.3f %8.3f", x_coordinates[i], At_v_x[i]);
		write_log(sprintf_buffer);
		}*/	
	calc_atom(mu, At_v_x_next, true, false);
	calc_potential(Z, mu, At_v_x_next_unmod, At_v_x_next);
	mix_potential(At_v_x_unmod, At_v_x_next_unmod);
	mix_potential(At_v_x, At_v_x_next);
	/*for(i=0; i<x_count; i++){
		sprintf(sprintf_buffer, "%8.3f %8.3f %8.3f %8.3f", x_coordinates[i], At_v_x[i], At_v_x_next_unmod[i], At_v_x_next[i]);
		write_log(sprintf_buffer);
		}*/
	double alpha=calc_criterion_a(At_v_x, At_v_x_next);
	double beta=calc_criterion_b(At_v_x_unmod, At_v_x_next_unmod, mu);
	sprintf(sprintf_buffer, "Criterion a = %8.3e, b = %8.3e", alpha, beta);
	write_log(sprintf_buffer);
	int trial_count=0;
	while(alpha>SC_criterion_a || beta>SC_criterion_b){
		trial_count++;
		sprintf(sprintf_buffer, "----Trial N=%d ----", trial_count);
		write_log(sprintf_buffer);
		calc_atom(mu, At_v_x_next, false, /*trial_count>10*/ false);
		// copy of unmodified potential
		// copy of modified potential is in calc_atom
		for(i=0; i<x_count; i++){
			At_v_x_unmod[i]=At_v_x_next_unmod[i];
		}
		calc_potential(Z, mu, At_v_x_next_unmod, At_v_x_next);
		mix_potential(At_v_x, At_v_x_next);
		mix_potential(At_v_x_unmod, At_v_x_next_unmod);
		
		alpha=calc_criterion_a(At_v_x, At_v_x_next);
		beta=calc_criterion_b(At_v_x_unmod, At_v_x_next_unmod, mu);

		sprintf(sprintf_buffer, "Criterion a = %8.3e, b = %8.3e", alpha, beta);
		write_log(sprintf_buffer);
		// break;
	}
	/*
	for(i=0; i<x_count; i+=4){
		sprintf(sprintf_buffer, "%8.3f %8.5f", x_coordinates[i], -At_v_x_next[i]*mu*x_coordinates[i]/Z);
		write_log(sprintf_buffer);
		}*/
	delete sprintf_buffer;
}

void calc_atom(double mu, double* v_x, bool first_calculation, bool use_guess){
	int i, j;
	char* sprintf_buffer=new char[Log_length+1];
	double* eigen_guess=new double[SC_orbital_count];
	// clear sigma_x
	for(i=0; i<x_count; i++){
		SC_sigma_x[i]=0;
	}
	if(!first_calculation){
		for(i=0; i<SC_orbital_count; i++){
			sprintf(sprintf_buffer, "Preparation for n=%d, l=%d", SC_n_list[i], SC_l_list[i]);
			write_log(sprintf_buffer);
			eigen_guess[i]=SC_eigen_list[i];
			sprintf(sprintf_buffer, "Eigenvalue in the last calculations: %8.3f", eigen_guess[i]);
			write_log(sprintf_buffer);
			for(j=0; j<x_count; j++){
				if(j!=0){
					eigen_guess[i]+=pow(SC_p_x[i][j], 2)*(v_x[j]-At_v_x[j])*mu*(x_coordinates[j]-x_coordinates[j-1]);
				}
			}
			sprintf(sprintf_buffer, "Guessed eigenvalue: %8.3f", eigen_guess[i]);
			write_log(sprintf_buffer);
		}
		for(j=0; j<x_count; j++){
			At_v_x[j]=v_x[j];
		}
	}

	for(i=0; i<SC_orbital_count; i++){
		int n=SC_n_list[i];
		int l=SC_l_list[i];
		sprintf(sprintf_buffer, "Calculation for n=%d, l=%d", n, l);
		write_log(sprintf_buffer);
		if(!use_guess){
			SC_eigen_list[i]=calc_atomic_wfn(mu, l, n-l-1);
		}else{
			SC_eigen_list[i]=calc_atomic_wfn2(mu, l, n-l-1, eigen_guess[i]);
			// SC_eigen_list[i]=calc_atomic_wfn(mu, l, n-l-1);
		}
		// copy the wfn
		for(j=0; j<x_count; j++){
			SC_p_x[i][j]=At_p_x[j];
			SC_sigma_x[j]+=At_occupation[n-1][l]*pow(At_p_x[j], 2);
		}
	}
}

void calc_potential(int Z, double mu, double* v_x_unmod, double* v_x){
	// V=-Z/r                     potential from nucleus
	//   +1/r int_0^r sigma(R) dR Hartree potential from the inside
	//   +int_r^inf sigma(R)/R dR Hartree potential from the outside
 	//   -3(3/8pi rho(r))^1/3     LDA exchange potential
	// r=mu x
	// rho(r)=sigma(r)/4pi r^2
	int i, j;
	for(i=0; i<x_count; i++){
		if(i==0){
			v_x_unmod[i]=0;
			v_x[i]=0;
			continue;
		}
		double r=mu*x_coordinates[i];
		double V1=-Z*1.0/r;
		double V2=0;
		for(j=0; j<i; j++){
			V2+=SC_sigma_x[j]*(x_coordinates[j+1]-x_coordinates[j]);
		}
		V2*=1.0/x_coordinates[i]; // mu/(mu*x)=x
		double V3=0;
		for(j=i; j<x_count-1; j++){
			V3+=SC_sigma_x[j]*(x_coordinates[j+1]-x_coordinates[j])/x_coordinates[j];
		}
		double V4=-3.0*pow(3.0/(8.0*M_PI)*SC_sigma_x[i]/(4.0*M_PI*pow(x_coordinates[i]*mu, 2.0)), 1.0/3.0);
		v_x_unmod[i]=V1+V2+V3+V4;
	}

	// modification
	for(i=1; i<x_count; i++){
		double v_x_threshold=-1.0/(mu*x_coordinates[i]);
		if(v_x_unmod[i]>v_x_threshold){
			v_x[i]=v_x_threshold;
		}else{
			v_x[i]=v_x_unmod[i];
		}
	}
}

void mix_potential(double* v_x, double* v_x_next){
	int i;
	for(i=0; i<x_count; i++){
		v_x_next[i]=(1.0-SC_mix_weight)*v_x_next[i]+SC_mix_weight*v_x[i];
	}
}


double calc_criterion_a(double* v1_x, double* v2_x){
	double alpha=-10.0;
	int i;
	for(i=1; i<x_count; i++){
		double diff=(v1_x[i]-v2_x[i])/v1_x[i];
		if(abs(diff)>alpha){
			alpha=abs(diff);
		}
		// cout << diff << endl;
	}
	return alpha;
}


double calc_criterion_b(double* v1_x, double* v2_x, double mu){
	double beta=0.0;
	int i;
	for(i=0; i<x_count; i++){
		double diff=mu*x_coordinates[i]*(v1_x[i]-v2_x[i]);
		if(abs(diff)>beta){
			beta=abs(diff);
		}
	}
	return beta;
}
