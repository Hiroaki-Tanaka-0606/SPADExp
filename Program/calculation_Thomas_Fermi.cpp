// calcPSF calculation of Thomas-Fermi potential calculation

#include <cmath>

#include "variables_ext.hpp"
#include "log.hpp"

void Thomas_Fermi_evolution();

void test_Thomas_Fermi(){
	// perform Thomas-Fermi differential equation calculation for some initial values
	write_log((char*)"----Test calculation----");
	int i, j;
	char* sprintf_buffer=new char[Log_length+1];
	double** TF_phis=new double*[Initial_diff_size+1];
	double Initial_diff;

	fprintf(Output_file_obj, "# calcPSF\n");
	fprintf(Output_file_obj, "# Thomas-Fermi calculation, test-mode\n");

	for(i=0; i<= Initial_diff_size; i++){
		Initial_diff=Initial_diff_offset+Initial_diff_delta*i;
		TF_phi[0]=1;
		TF_phi_diff[0]=Initial_diff;
		Thomas_Fermi_evolution();
		TF_phis[i]=new double[x_count];
		for(j=0; j<x_count; j++){
			TF_phis[i][j]=TF_phi[j];
		}
		
		sprintf(sprintf_buffer, "Step %3d: Initial_diff = %8.3f", (i+1), Initial_diff);
		fprintf(Output_file_obj, "# %s\n", sprintf_buffer);
		write_log(sprintf_buffer);
	}

	fprintf(Output_file_obj, "# x phi_i(x)\n");
	for(i=0; i<x_count; i++){
		fprintf(Output_file_obj, "%11.5f ", x_coordinates[i]);
		for(j=0; j<=Initial_diff_size; j++){
			fprintf(Output_file_obj, "%8.5f ", TF_phis[j][i]);
		}
		fprintf(Output_file_obj, "\n");
	}
		
	
	for(i=0; i<=Initial_diff_size; i++){
		delete TF_phis[i];
	}
	delete TF_phis;
	delete sprintf_buffer;
}

double trial_Thomas_Fermi(double Initial_diff){
	TF_phi[0]=1;
	TF_phi_diff[0]=Initial_diff;
	Thomas_Fermi_evolution();
	return TF_phi[x_count-1];
}

void calc_Thomas_Fermi(){
	// find the adequate Initial_diff using bisection method
	write_log((char*)"----Actual calculation----");
	fprintf(Output_file_obj, "# calcPSF\n");
	fprintf(Output_file_obj, "# Thomas-Fermi calculation, actual-mode\n");

	char* sprintf_buffer=new char[Log_length+1];

	double phi_edge_min, phi_edge_max, phi_edge_center;

	
	phi_edge_min=trial_Thomas_Fermi(Initial_diff_min);
	phi_edge_max=trial_Thomas_Fermi(Initial_diff_max);

	if(phi_edge_min*phi_edge_max>0){
		write_log((char*)"Error: cannot start the bisection method");
		delete sprintf_buffer;
		return;
	}
	double Initial_diff_center=(Initial_diff_min+Initial_diff_max)/2.0;
	phi_edge_center=trial_Thomas_Fermi(Initial_diff_center);
	int current_step=1;
	
	while(phi_edge_center<0 || phi_edge_center>TF_threshold){
		if(phi_edge_center*phi_edge_min<0){
			Initial_diff_max=Initial_diff_center;
		}else{
			Initial_diff_min=Initial_diff_center;
		}

	
		sprintf(sprintf_buffer, "Step %3d: Initial_diff = %11.8f, phi_edge = %11.5e", current_step, Initial_diff_center, phi_edge_center);
		write_log(sprintf_buffer);
		
		Initial_diff_center=(Initial_diff_min+Initial_diff_max)/2.0;
		phi_edge_center=trial_Thomas_Fermi(Initial_diff_center);
		current_step++;
	}

	write_log((char*)"Calculation finished");
	sprintf(sprintf_buffer, "Step %3d: Initial_diff = %11.8f, phi_edge = %11.5e", current_step, Initial_diff_center, phi_edge_center);
	write_log(sprintf_buffer);

	fprintf(Output_file_obj, "# Initial diff = %11.8f\n", Initial_diff_center);
	fprintf(Output_file_obj, "# x phi(x)\n");
	int i;
	for(i=0; i<x_count; i++){
		fprintf(Output_file_obj, "%11.5f %8.5e\n", x_coordinates[i], TF_phi[i]);
	}
	delete sprintf_buffer;
}

void Thomas_Fermi_evolution(){
	int i;
	for(i=1; i<x_count; i++){
		TF_phi[i]=TF_phi[i-1]+TF_phi_diff[i-1]*(x_coordinates[i]-x_coordinates[i-1]);
		if(i==1){
			// to avoid the division by zero
			TF_phi_diff[i]=TF_phi_diff[i-1];
		}else if(TF_phi[i-1]<0){
			TF_phi_diff[i]=0;
		}else{
			TF_phi_diff[i]=TF_phi_diff[i-1]+pow(TF_phi[i-1], 1.5)/sqrt(x_coordinates[i-1])*(x_coordinates[i]-x_coordinates[i-1]);
		}
	}
}
