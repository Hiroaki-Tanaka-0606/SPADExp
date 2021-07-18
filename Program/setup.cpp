// calcPSF setup (initialization and finalization)

#include "variables_ext.hpp"
#include <cstring>

void initialize(){
	// constants
	buffer_length=1024;
	Calculation_length=64;
	Log_file_length=1024;
	Log_length=256;
	Log_buffer=new char[Log_length+1];
	Solution_length=16;
	
	// variables
	/// blocks
	Ct_block_appeared=false;
	TF_block_appeared=false;
	/// Ct block
	Calculation=new char[Calculation_length+1];
	Calculation_set=false;
	Log_file=new char[Log_file_length+1];
	Log_file_set=false;
	Console_log=true;
	Console_log_set=false;
	Output_file=new char[Log_file_length+1];
	Output_file_set=false;
	/// TF block
  TF_test=false;
	TF_test_set=false;
	Initial_diff_offset_set=false;
	Initial_diff_delta_set=false;
	Initial_diff_size_set=false;
	Initial_diff_min=-1.49; // -1.52 for RK1
	Initial_diff_max=-1.48; // -1.51 for RK1
	Initial_diff_min_set=false;
	Initial_diff_max_set=false;
	TF_threshold=1e-5;
	TF_threshold_set=false;
	TF_solution=new char[Solution_length+1];
	strcpy(TF_solution, "RK4");
	TF_solution_set=false;
	
}

void finalize(){
	delete Calculation;
	if(Log_file_set && Log_file_obj!=NULL){
		fclose(Log_file_obj);
	}
	if(Output_file_set && Output_file_obj!=NULL){
		fclose(Output_file_obj);
	}
	if(Radial_grid_count>0){
		delete Radial_grid_intervals;
		delete Radial_grid_points;
	}
}

void initialize_radial_grid(){
	Radial_grid_count=11;
	Radial_grid_intervals=new double[11];
	Radial_grid_intervals[0]=0.0025;
	Radial_grid_intervals[1]=0.005;
	Radial_grid_intervals[2]=0.01;
	Radial_grid_intervals[3]=0.02;
	Radial_grid_intervals[4]=0.04;
	Radial_grid_intervals[5]=0.08;
	Radial_grid_intervals[6]=0.16;
	Radial_grid_intervals[7]=0.32;
	Radial_grid_intervals[8]=0.64;
	Radial_grid_intervals[9]=1.28;
	Radial_grid_intervals[10]=2.56;
	Radial_grid_points=new int[11];
	Radial_grid_points[0]=40;
	Radial_grid_points[1]=40;
	Radial_grid_points[2]=40;
	Radial_grid_points[3]=40;
	Radial_grid_points[4]=40;
	Radial_grid_points[5]=40;
	Radial_grid_points[6]=40;
	Radial_grid_points[7]=40;
	Radial_grid_points[8]=40;
	Radial_grid_points[9]=40;
	Radial_grid_points[10]=40;
}

void generate_radial_grid(){
	int i, j, index;
	int points_count=1;
	for(i=0; i<Radial_grid_count; i++){
		points_count+=Radial_grid_points[i];
	}
	x_coordinates=new double[points_count];
	x_coordinates[0]=0;
	index=1;
	double x_current=0;
	for(i=0; i<Radial_grid_count; i++){
		for(j=0; j<Radial_grid_points[i]; j++){
			x_coordinates[index]=x_coordinates[index-1]+Radial_grid_intervals[i];
			index++;
		}
	}
	TF_phi=new double[points_count];
	TF_phi_diff=new double[points_count];
	x_count=points_count;
}
