// calcPSF tools for OpenMX
// postproc_main

// include libraries from include path
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <string>
#include <ctime>
#include "H5Cpp.h"

// include libraries from current directory
#include "inputTools.hpp"
#include "HDF5_tools.hpp"

// namespace
using namespace std;
using namespace H5;


int main(int argc, const char** argv){
	cout << "Starting postproc tool for OpenMX..." << endl;
	// argv[1] should be an input file name
	if(argc<2){
		printf("Usage: %s input_file\n", argv[0]);
		return 0;
	}

	ifstream input(argv[1]);
	if(!input){
		printf("Error: cannot open %s\n", argv[0]);
		return 0;
	}
	load_file(&input);
	input.close();

	// buffer sizes
	int itemSize=32;
	int valSize=128;
	int pathSize=1024;
	
	// variables
	int i, j;
	bool curved;
	int dimension;
	double origin_frac[3];
	double x_frac[3];
	double y_frac[3];
	double origin_orth[3];
	double x_orth[3];
	double y_orth[3];
	double xMin, xMax;
	double yMin, yMax;
	int x_count;
	int y_count;
	double cell[3][3];
	double bandCell[3][3];

	char sysName[valSize+1];
	char currentDir[valSize+1];
	char vectorUnit[valSize+1];
	char atomUnit[valSize+1];

	int specNum;
	int atomNum;

	if(load_str("System.Name", sysName, valSize)==1){
		printf("System.Name is %s\n", sysName);
	}else{
		return 0;
	}

	if(load_str("System.CurrentDir", currentDir, valSize)==1){
		printf("System.CurrentDir is %s\n", currentDir);
	}else{
		strcpy(currentDir, "./");
	}

	if(load_int("Species.Number", &specNum)==1){
		printf("Species.Number is %d\n", specNum);
	}else{
		return 0;
	}

	int line_species_start=find_str("<Definition.of.Atomic.Species");
	int line_species_end=find_str("Definition.of.Atomic.Species>");
	char species[specNum][4][valSize];
	if(line_species_start>=0 && line_species_end>=0 &&
		 line_species_end-line_species_start-1==specNum){
		for(i=0; i<specNum; i++){
			char* range=get_line(line_species_start+1+i);
			if(sscanf(range, "%s %s %s", species[i][0], species[i][1], species[i][3])!=3){
				printf("Error in parsing Atomic.Species\n");
				return 0;
			}
			for(j=0; j<strlen(species[i][1]); j++){
				if(species[i][1][j]=='-'){
					species[i][1][j]='\0';
					strcpy(species[i][2], &species[i][1][j+1]);
				}
			}
			// printf("%s %s %s %s\n", species[i][0], species[i][1], species[i][2], species[i][3]);
		}
	}else{
		printf("Error in Definition.of.Atomic.Species");
		return 0;
	}

	if(load_int("Atoms.Number", &atomNum)==1){
		printf("Atoms.Number is %d\n", atomNum);
	}else{
		return 0;
	}

	int line_atom_start=find_str("<Atoms.SpeciesAndCoordinates");
	int line_atom_end=find_str("Atoms.SpeciesAndCoordinates>");
	char atoms[atomNum][valSize];
	double coordinates[atomNum][3];
	if(line_atom_start>=0 && line_atom_end>=0 &&
		 line_atom_end-line_atom_start-1==atomNum){
		for(i=0; i<atomNum; i++){
			char* range=get_line(line_atom_start+1+i);
			if(sscanf(range, "%*s %s %lf %lf %lf", atoms[i], &coordinates[i][0], &coordinates[i][1], &coordinates[i][2])!=4){
				printf("Error in parsing Atoms.SpeciesAndCoordinates\n");
				return 0;
			}
		}
	}else{
		printf("Error in Atoms.SpeciesAndCoordinates\n");
		return 0;
	}

	if(load_str("Atoms.SpeciesAndCoordinates.Unit", atomUnit, valSize)==1){
		printf("Atoms.UnitVectors.Unit is %s\n", atomUnit);
	}else{
		return 0;
	}
	
	
	
	if(load_str("Atoms.UnitVectors.Unit", vectorUnit, valSize)==1){
		printf("Atoms.UnitVectors.Unit is %s\n", vectorUnit);
	}else{
		return 0;
	}

	if(load_int("calcPSF.dimension", &dimension)==1){
		printf("calcPSF.dimension is %d\n", dimension);
		if(!(1<=dimension && dimension<=2)){
			printf("Error: invalid dimension\n");
			return 0;
		}	
	}else{
		return 0;
	}
	
	if(load_logical("calcPSF.curved", &curved)==1){
		printf("calcPSF.curved is %s\n", (curved ? "true": "false"));
	}else{
		return 0;
	}

	if(load_doublev("calcPSF.origin", 3, origin_frac)==1){
		// printf("Origin is (%.3f, %.3f, %.3f) in fractional coordinates\n", origin_frac[0], origin_frac[1], origin_frac[2]);
	}else{
		return 0;
	}

	int line_range_start=find_str("<calcPSF.range");
	int line_range_end=find_str("calcPSF.range>");

	if(line_range_start>=0 && line_range_end>=0){
		if(line_range_end-line_range_start-1==dimension){
			char* range=get_line(line_range_start+1);
			if(sscanf(range, "%lf %lf %lf %lf %lf %d", &x_frac[0], &x_frac[1], &x_frac[2], &xMin, &xMax, &x_count)!=6){
				printf("Error in parsing x range");
				return 0;
			}
			if(x_count<1){
				printf("Error: x_count must be positive\n");
				return 0;
			}
			range=get_line(line_range_start+2);
			if(dimension>1){
				if(sscanf(range, "%lf %lf %lf %lf %lf %d", &y_frac[0], &y_frac[1], &y_frac[2], &yMin, &yMax, &y_count)!=6){
					printf("Error in parsing y range");
					return 0;
				}
				if(y_count<1){
					printf("Error: y_count must be positive\n");
					return 0;
				}
			}
		}else{
			printf("Dimension error\n");
			return 0;
		}
	}else{
			return 0;
	}

	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", x_frac[0], x_frac[1], x_frac[2], xMin, xMax, x_count);
	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", y_frac[0], y_frac[1], y_frac[2], yMin, yMax, y_count);

	int total_count= (dimension==1) ? x_count : x_count*y_count;
	double kps[total_count][3];
	int line_kp_start=find_str("<MO.kpoint");
	int line_kp_end=find_str("MO.kpoint>");

	if(line_kp_start>=0 && line_kp_end>=0){
		if(line_kp_end-line_kp_start-1==total_count){
			for(i=0; i<total_count; i++){
				char* range=get_line(line_kp_start+1+i);
				if(sscanf(range, "%lf %lf %lf", &kps[i][0], &kps[i][1], &kps[i][2])!=3){
					printf("Error in pasing MO.kpoint\n");
					return 0;
				}
			}
		}else{
			printf("Error in MO.kpoint\n");
		}
	}else{
		return 0;
	}
	

	int line_bandCell_start=find_str("<Band.KPath.UnitCell");
	int line_bandCell_end=find_str("Band.KPath.UnitCell>");

	int line_cell_start=find_str("<Atoms.UnitVectors");
	int line_cell_end=find_str("Atoms.UnitVectors>");


	if(line_cell_end-line_cell_start-1!=3){
		printf("Cell error\n");
		return 0;
	}

	for(i=0; i<3; i++){
		char* cell_str=get_line(line_cell_start+1+i);
		if(sscanf(cell_str, "%lf %lf %lf", &cell[i][0], &cell[i][1], &cell[i][2])!=3){
			printf("Error in parsing cell\n");
			return 0;
		}
	}

	if(line_bandCell_start>=0 && line_bandCell_end>=0){
		if(line_bandCell_end-line_bandCell_start-1!=3){
			printf("BandCell error\n");
			return 0;
		}
		for(i=0; i<3; i++){
			char* cell_str=get_line(line_bandCell_start+1+i);
			if(sscanf(cell_str, "%lf %lf %lf", &bandCell[i][0], &bandCell[i][1], &bandCell[i][2])!=3){
				printf("Error in parsing bandCell\n");
				return 0;
			}
		}
	}else{
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				bandCell[i][j]=cell[i][j];
			}
		}
	}

	printf("Unit vectors (real space) for atoms\n");
	for(i=0; i<3; i++){
	 	printf("a_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), cell[i][0], cell[i][1], cell[i][2]);
	}
	
	printf("Unit vectors (real space) for bands\n");
	for(i=0; i<3; i++){
	 	printf("a_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), bandCell[i][0], bandCell[i][1], bandCell[i][2]);
	}


	cout << "Loading from the input file finished." << endl;
	cout << endl;

	// open the output HDF5 file
	char outFilePath[pathSize];
	sprintf(outFilePath, "%s%s.hdf5", currentDir, sysName);
	printf("Open %s\n", outFilePath);

	H5File output(outFilePath, H5F_ACC_TRUNC);

	// variable-length string
	StrType str_var(PredType::C_S1, H5T_VARIABLE);
	str_var.setCset(H5T_CSET_UTF8);

	// variables for data processing
	Attribute at;
	int data_i[1];
	string data_s[1];
	bool data_b[1];

	// att Datetime @root
	time_t datetime_now=time(NULL);
	struct tm *timeptr=localtime(&datetime_now);
	char time_str[valSize+1];
	strftime(time_str, valSize, "%Y-%m-%d %H:%M:%S", timeptr);
	w_att_str(output.openGroup("/"), "Datetime", time_str);

	// group /Input
	Group InputG(output.createGroup("/Input"));

	// group /Input/Atomic.Species
	Group SpecG(InputG.createGroup("Atomic.Species"));

	// atts @/Input/Atomic.Species
	w_att_int(SpecG, "Length", specNum);

	// data @/Input/Atomic.Species
	char spec_str[specNum][itemSize+1];
	for(i=0; i<specNum; i++){
		for(j=0; j<4; j++){
			if(strlen(species[i][j])>itemSize){
				printf("Too long item in Atomic.Species");
				return 0;
			}
		}
	}
	for(i=0; i<specNum; i++){
	  strcpy(spec_str[i], species[i][0]);
	}
	w_data_1c(SpecG, "Labels", specNum, itemSize+1, (char**)spec_str);

	for(i=0; i<specNum; i++){
		strcpy(spec_str[i], species[i][1]);
	}
	w_data_1c(SpecG, "PAOs", specNum, itemSize+1, (char**)spec_str);

	for(i=0; i<specNum; i++){
		strcpy(spec_str[i], species[i][2]);
	}
	w_data_1c(SpecG, "Orbits", specNum, itemSize+1, (char**)spec_str);
	
	for(i=0; i<specNum; i++){
		strcpy(spec_str[i], species[i][3]);
	}
	w_data_1c(SpecG, "Pseudopotentials", specNum, itemSize+1, (char**)spec_str);

	// group /Input/Atoms.SpeciesAndCoordinates
	Group AtomG(InputG.createGroup("Atoms.SpeciesAndCoordinates"));

	// atts @/Input/Atoms.~~
	w_att_int(AtomG, "Length", atomNum);
	w_att_str(AtomG, "Unit", atomUnit);

	// data @/Input/Atoms.~~
	char atom_str[atomNum][itemSize+1];
	for(i=0; i<atomNum; i++){
		if(strlen(atoms[i])>itemSize){
			printf("Too long item in Atoms.SpeciesAndCoordinates");
			return 0;
		}
	}
	for(i=0; i<atomNum; i++){
		strcpy(atom_str[i], atoms[i]);
	}
	w_data_1c(AtomG, "Labels", atomNum, itemSize+1, (char**) atom_str);
	w_data_2d(AtomG, "Coordinates", atomNum, 3, (double**) coordinates);
	
	
	// group /Input/UnitCells
	Group UnitG(InputG.createGroup("UnitCells"));

	// atts @/Input/UnitCells
	w_att_str(UnitG, "Unit", vectorUnit);
	w_att_2d(UnitG, "Atoms", 3, 3, (double**)cell);
	w_att_2d(UnitG, "Bands", 3, 3, (double**)bandCell);


	// group /Input/Kpath
	Group KpathG(InputG.createGroup("Kpath"));

	// atts @/Input/Kpath
  w_att_int(KpathG, "Dimension", dimension);
	w_att_bool(KpathG, "Curved", curved);
	w_att_1d(KpathG, "Origin", 3, origin_frac);
	w_att_1d(KpathG, "Xvector", 3, x_frac);
	w_att_1d(KpathG, "Yvector", 3, y_frac);
	w_att_int(KpathG, "Xcount", x_count);
	w_att_int(KpathG, "Ycount", y_count);
	w_att_int(KpathG, "Count", total_count);
	double xRange[2]={xMin, xMax};
	w_att_1d(KpathG, "Xrange", 2, xRange);
	double yRange[2]={yMin, yMax};
	w_att_1d(KpathG, "Yrange", 2, yRange);

	// data @/Input/Kpath
	w_data_2d(KpathG, "Coordinates", total_count, 3, (double**)kps);
	
	// DataSpace ds(rank, [x,y,z]);



	
	output.close();
	cout << "Finished writing data." << endl;
}
