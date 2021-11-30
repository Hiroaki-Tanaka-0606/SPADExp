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
#include <complex>
#include <H5Cpp.h>

// include libraries from current directory
#include "inputTools.hpp"
#include "HDF5_tools.hpp"

// namespace
using namespace std;
using namespace H5;

bool readLine(ifstream* fp, char* line_c, int bufSize);

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
	int bufSize=4096;

	// diff criterion
	double diff_threshold1=0.01; // kp coordinates
	
	// variables
	int i, j, k, l;
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
	char spinPol[valSize+1];
	int spinPol_i; // 0: off, 1: on, 2: nc (noncollinear)

	int specNum;
	int atomNum;

	int maxN;
	int minN;

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
	int numOrbits[specNum];
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
			numOrbits[i]=0;
			for(j=0; j<strlen(species[i][2]); j+=2){
				char orbitLabel=species[i][2][j];
				char numLabel[1]={species[i][2][j+1]};
				if(orbitLabel=='s'){
					numOrbits[i]+=atoi(numLabel);
				}else if(orbitLabel=='p'){
					numOrbits[i]+=3*atoi(numLabel);
				}else if(orbitLabel=='d'){
					numOrbits[i]+=5*atoi(numLabel);
				}else if(orbitLabel=='f'){
					numOrbits[i]+=7*atoi(numLabel);
				}else{
					printf("Error: invalid orbit label %c\n", orbitLabel);
					return 0;
				}
			}
			// printf("%s %s %s %s\n", species[i][0], species[i][1], species[i][2], species[i][3]);
			// cout << numOrbits[i] << endl;
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
	int atomOrbits[atomNum];
	int orbitIndices[atomNum];
	int totalOrbits=0;
	if(line_atom_start>=0 && line_atom_end>=0 &&
		 line_atom_end-line_atom_start-1==atomNum){
		for(i=0; i<atomNum; i++){
			char* range=get_line(line_atom_start+1+i);
			if(sscanf(range, "%*s %s %lf %lf %lf", atoms[i], &coordinates[i][0], &coordinates[i][1], &coordinates[i][2])!=4){
				printf("Error in parsing Atoms.SpeciesAndCoordinates\n");
				return 0;
			}
			bool specFound=false;
			for(j=0; j<specNum; j++){
				if(strcmp(atoms[i], species[j][0])==0){
					atomOrbits[i]=numOrbits[j];
					totalOrbits+=numOrbits[j];
					specFound=true;
					// cout << j << endl;
				}
			}
			if(!specFound){
				printf("Error: cannot find atom %s\n", atoms[i]);
				return 0;
			}

				
		}
	}else{
		printf("Error in Atoms.SpeciesAndCoordinates\n");
		return 0;
	}
	j=0;
	orbitIndices[0]=0;
	for(i=0; i<atomNum-1; i++){
		j+=atomOrbits[i];
		orbitIndices[i+1]=j;
	}
	// for(i=0; i<atomNum; i++){
	// 	cout << orbitIndices[i] << endl;
	// }
	

	if(load_str("Atoms.SpeciesAndCoordinates.Unit", atomUnit, valSize)==1){
		printf("Atoms.SpeciesAndCoordinates.Unit is %s\n", atomUnit);
	}else{
		return 0;
	}
	
	// cout << totalOrbits << endl;
	
	if(load_str("Atoms.UnitVectors.Unit", vectorUnit, valSize)==1){
		printf("Atoms.UnitVectors.Unit is %s\n", vectorUnit);
	}else{
		return 0;
	}

	if(load_str("scf.SpinPolarization", spinPol, valSize)==1){
		printf("scf.SpinPolarization is %s\n", spinPol);
		if(strncasecmp("off", spinPol, strlen(spinPol))==0 && strlen(spinPol)==3){
			spinPol_i=0;
		}else if(strncasecmp("on", spinPol, strlen(spinPol))==0 && strlen(spinPol)==2){
			spinPol_i=1;
		}else if(strncasecmp("nc", spinPol, strlen(spinPol))==0 && strlen(spinPol)==2){
			spinPol_i=2;
		}else{
			printf("Error: invalid scf.SpinPolarization\n");
			return 0;
		}
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

	if(load_int("calcPSF.minN", &minN)==1){
		printf("calcPSF.minN is %d\n", minN);
	}else{
		return 0;
	}

	if(load_int("calcPSF.maxN", &maxN)==1){
		printf("calcPSF.maxN is %d\n", maxN);
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
	w_data_2d(UnitG, "Atoms", 3, 3, (double**)cell);
	w_data_2d(UnitG, "Bands", 3, 3, (double**)bandCell);


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
	w_att_int(KpathG, "minN", minN);
	w_att_int(KpathG, "maxN", maxN);

	// data @/Input/Kpath
	w_data_2d(KpathG, "Coordinates", total_count, 3, (double**)kps);
	
	cout << "Finished writing /Input" << endl;
	cout << endl;
	
	// open the output of OpenMX
	char outOpenMXPath[pathSize];
	sprintf(outOpenMXPath, "%s%s.out", currentDir, sysName);
	printf("Open %s\n", outOpenMXPath);

	// data loaded from the output file
	int numBands= (spinPol_i==2) ? totalOrbits*2 : totalOrbits;
	int numBands_out= (spinPol_i==2) ? (maxN-minN+1)*2 : (maxN-minN+1);
	double Band[total_count][numBands_out];
	double BandUp[total_count][numBands_out];
	double BandDn[total_count][numBands_out];
	double eigen_lowerTop=0; // largest eigenvalue of band minN-1
	double eigen_upperBottom=0; // smallest eigenvalue of band maxN+1

	int digit=(spinPol_i==2)?4:2;
	double LCAO_all[total_count][numBands_out][totalOrbits][digit];
	double LCAOUp_all[total_count][numBands_out][totalOrbits][2];
	double LCAODn_all[total_count][numBands_out][totalOrbits][2];

	ifstream outOpenMX(outOpenMXPath);
	string line;
	char line_c[bufSize];
	int line_number=0;
	bool LCAO_found=false;
	while(getline(outOpenMX, line)){
		line_number++;
		if(line=="        Eigenvalues (Hartree) and LCAO coefficients        "){
			printf("The LCAO block starts from line %d\n", line_number+4);
			LCAO_found=true;
			break;
		}
	}
	if(!LCAO_found){
		printf("Error: cannot find the LCAO block\n");
		return 0;
	}

	for(i=0; i<3; i++){
		getline(outOpenMX, line);
		line_number++;
	}
	// here, the next line is the first line of the LCAO block

	int kp_index;
	int actual_line_length;
	char kpBuf[valSize];
	int indexBuf;
	bool kp_found;
	double EF_Eh; // in units of the Hartree energy (27.2 eV)
	for(kp_index=1; kp_index<=total_count; kp_index++){
		// cout << kp_index << endl << endl;
		// # of k-point = XX
		kp_found=false;
		while(readLine(&outOpenMX, line_c, bufSize)){
			line_number++;
			if(sscanf(line_c, "%s %*s %*s %*s %d", kpBuf, &indexBuf)==2 &&
				 strlen(kpBuf)==1 && kpBuf[0]=='#'){
				if(indexBuf!=kp_index){
					printf("Error: kp index mismatch\n");
					return 0;
				}
				printf("kp #%d starts\n", kp_index);
				kp_found=true;
				break;
			}
		}
		if(!kp_found){
			printf("Error: cannot find kp #%d\n", kp_index);
			return 0;
		}
		// kp coordinates check
		readLine(&outOpenMX, line_c, bufSize);
		 double kp_read[3];
		 if(sscanf(line_c, "%*s %lf %*s %lf %*s %lf", &kp_read[0], &kp_read[1], &kp_read[2])==3){
			 /* for(i=0; i<3; i++){
		 		if(abs(kp_read[i]-kps[kp_index-1][i])>diff_threshold1){
		 			printf("Warning: kp coordinate mismatch\n");
		 		}
				}*/
		 }else{
			printf("Error: cannot find kp coordinate\n");
			return 0;
		}
		readLine(&outOpenMX, line_c, bufSize);
		readLine(&outOpenMX, line_c, bufSize);
		// chemical potential (only when kp_index==1)
		if(kp_index==1){
			if(sscanf(line_c, "%*s %*s %*s %*s %lf", &EF_Eh)==1){
				printf("The Fermi energy (in units of the Hartree energy) is %.5f\n", EF_Eh);
			}else{
				printf("Error in parsing Chemical Potential\n");
				return 0;
			}
		}
		// skip 6 rows if spinPol==off or nc, else 7 (spinPol==on)
		// <=> The former: spinPol_i==0 or 2, the latter: spinPol_i==1
		int skipRows= (spinPol_i==1) ? 7 : 6;
		for(i=0; i<skipRows; i++){
			readLine(&outOpenMX, line_c, bufSize);
		}

		int band_remaining=numBands;
		int band_index=0;
		int lcao_remaining=numBands;
		int lcao_index=0;
		bool loadingUp=true;
		bool loadingUp_l=true;
		double eigen[4];
		double lcao[8];
		while(band_remaining>0){
			readLine(&outOpenMX, line_c, bufSize);
			// parse eigenvalues
			switch(spinPol_i){
			case 0:
				// spinPol off: 4 values, only (U)
				if(sscanf(line_c, "%lf %lf %lf %lf", &eigen[0], &eigen[1], &eigen[2], &eigen[3])==min(4, band_remaining)){
					for(i=0; (i<4 && band_remaining>0); i++){
						if(band_index>=minN && band_index<=maxN){
							Band[kp_index-1][band_index-minN]=eigen[i];
						}
						if(band_index==minN-1){
							if(kp_index==1 || eigen[i]>eigen_lowerTop){
								eigen_lowerTop=eigen[i];
							}
						}
						if(band_index==maxN+1){
							if(kp_index==1 || eigen[i]<eigen_upperBottom){
								eigen_upperBottom=eigen[i];
							}
						}
						band_index++;
						band_remaining--;
					}
				}else{
					printf("Error in parsing eigenvalues\n");
					return 0;
				}
				break;
			case 1:
				// sponPol on: 4 values, (U)s following (D)s
				if(sscanf(line_c, "%lf %lf %lf %lf", &eigen[0], &eigen[1], &eigen[2], &eigen[3])==min(4, band_remaining)){
					for(i=0; (i<4 && band_remaining>0); i++){
						if(band_index>=minN && band_index<=maxN){
							if(loadingUp){
								BandUp[kp_index-1][band_index-minN]=eigen[i];
							}else{
								BandDn[kp_index-1][band_index-minN]=eigen[i];
							}
						}
						if(band_index==minN-1){
							if(kp_index==1 || eigen[i]>eigen_lowerTop){
								eigen_lowerTop=eigen[i];
							}
						}
						if(band_index==maxN+1){
							if(kp_index==1 || eigen[i]<eigen_upperBottom){
								eigen_upperBottom=eigen[i];
							}
						}
						band_index++;
						band_remaining--;
					}
					if(band_remaining==0 && loadingUp){
						band_remaining=numBands;
						band_index=0;
						loadingUp=false;
					}
				}else{
					printf("Error in parsing eigenvalues\n");
					return 0;
				}
				break;
			case 2:
				// spinPol nc: 2 values
				if(sscanf(line_c, "%lf %lf", &eigen[0], &eigen[1])==min(2, band_remaining)){
					for(i=0; (i<2 && band_remaining>0); i++){
						if(band_index>=minN && band_index<=maxN){
							Band[kp_index-1][band_index-minN]=eigen[i];
						}
						if(band_index==minN-1){
							if(kp_index==1 || eigen[i]>eigen_lowerTop){
								eigen_lowerTop=eigen[i];
							}
						}
						if(band_index==maxN+1){
							if(kp_index==1 || eigen[i]<eigen_upperBottom){
								eigen_upperBottom=eigen[i];
							}
						}
						band_index++;
						band_remaining--;
					}
				}else{
					printf("Error in parsing eigenvalues\n");
					return 0;
				}
				break;
			}

			// skip 3 rows
			for(i=0; i<3; i++){
				readLine(&outOpenMX, line_c, bufSize);
			}

			// read LCAO coefficients
			for(i=0; i<totalOrbits; i++){
				readLine(&outOpenMX, line_c, bufSize);
				// cout << line_c << endl;
				int sscanf_result;
				bool flag=false;
				for(j=0; j<atomNum; j++){
					if(i==orbitIndices[j]){
						flag=true;
					}
				}
				if(flag){
					sscanf_result=sscanf(line_c, "%*d %*s %*d %*s %lf %lf %lf %lf %lf %lf %lf %lf",
															 &lcao[0], &lcao[1], &lcao[2], &lcao[3], &lcao[4], &lcao[5], &lcao[6], &lcao[7]);
				}else{
					sscanf_result=sscanf(line_c, "%*d %*s %lf %lf %lf %lf %lf %lf %lf %lf",
															 &lcao[0], &lcao[1], &lcao[2], &lcao[3], &lcao[4], &lcao[5], &lcao[6], &lcao[7]);
				}
				switch(spinPol_i){
				case 0:
					// spinPol off
					if(sscanf_result==2*min(4, lcao_remaining)){
						int lcao_remaining_tmp=lcao_remaining;
						for(j=0; (j<4 && lcao_remaining_tmp>0); j++){
							if(lcao_index+j>=minN && lcao_index+j<=maxN){
								LCAO_all[kp_index-1][lcao_index+j-minN][i][0]=lcao[2*j];
								LCAO_all[kp_index-1][lcao_index+j-minN][i][1]=lcao[2*j+1];
							}
							lcao_remaining_tmp--;
						}
					}else{
						printf("Error in parsing LCAO coefficients\n");
						return 0;
					}
					break;
				case 1:
					// spinPol on
					if(sscanf_result==2*min(4, lcao_remaining)){
						int lcao_remaining_tmp=lcao_remaining;
						for(j=0; (j<4 && lcao_remaining_tmp>0); j++){
							if(lcao_index+j>=minN && lcao_index+j<=maxN){
								if(loadingUp_l){
									LCAOUp_all[kp_index-1][lcao_index+j-minN][i][0]=lcao[2*j];
									LCAOUp_all[kp_index-1][lcao_index+j-minN][i][1]=lcao[2*j+1];
								}else{
									LCAODn_all[kp_index-1][lcao_index+j-minN][i][0]=lcao[2*j];
									LCAODn_all[kp_index-1][lcao_index+j-minN][i][1]=lcao[2*j+1];
								}
							}
							lcao_remaining_tmp--;
						}
					}else{
						printf("Error in parsing LCAO coefficients\n");
						// cout << flag << endl;
						return 0;
					}
					break;
				case 2:
					// spinPol nc
					if(sscanf_result==4*min(2, lcao_remaining)){
						int lcao_remaining_tmp=lcao_remaining;
						for(j=0; (j<2 && lcao_remaining_tmp>0); j++){
							if(lcao_index+j>=minN && lcao_index+j<=maxN){
								LCAO_all[kp_index-1][lcao_index+j-minN][i][0]=lcao[4*j];
								LCAO_all[kp_index-1][lcao_index+j-minN][i][1]=lcao[4*j+1];
								LCAO_all[kp_index-1][lcao_index+j-minN][i][2]=lcao[4*j+2];
								LCAO_all[kp_index-1][lcao_index+j-minN][i][3]=lcao[4*j+3];
							}
							lcao_remaining_tmp--;
						}
					}else{
						printf("Error in parsing LCAO coefficients\n");
						return 0;
					}
					break;
				}
			}

			// update lcao_index and lcao_remaining
			if(spinPol_i==2){
				if(lcao_remaining<2){
					lcao_index+=lcao_remaining;
					lcao_remaining=0;
				}else{
					lcao_index+=2;
					lcao_remaining-=2;
				}
			}else{
				if(lcao_remaining<4){
					lcao_index+=lcao_remaining;
					lcao_remaining=0;
				}else{
					lcao_index+=4;
					lcao_remaining-=4;
				}
				if(lcao_remaining==0 && spinPol_i==1 && loadingUp_l){
					loadingUp_l=false;
					lcao_remaining=numBands;
					lcao_index=0;
					// cout << "!!" << endl;
				}
			}
			// cout << lcao_remaining << endl;
			// skip 2 rows
			for(i=0; i<2; i++){
				readLine(&outOpenMX, line_c, bufSize);
			}
		}
			
			
		
	}

	// output lowerTop and upperBottom
	if(minN>0){
		printf("Largest eigenenergy in band minN-1 is %8.3f\n", (eigen_lowerTop-EF_Eh)*27.2114);
	}
	if(maxN<numBands-1){
		printf("Smallest eigenenergy in band maxN+1 is %8.3f\n", (eigen_upperBottom-EF_Eh)*27.2114);
	}

	
	// group /Output
	Group OutputG(output.createGroup("/Output"));


	
	// band data
	if(spinPol_i==1){
		w_data_2d(OutputG, "BandUp", total_count, numBands_out, (double**)BandUp);
		w_data_2d(OutputG, "BandDn", total_count, numBands_out, (double**)BandDn);
	}else{
		w_data_2d(OutputG, "Band", total_count, numBands_out, (double**)Band);
	}

	// LCAO data
	// For the conversion from px, py, pz to p(l=-1, 0, 1), see openmx/source/AngularF.c
	Group LCAOG(OutputG.createGroup("LCAO"));
	int specIndex=0;
	char groupName[itemSize];
	char datasetName[itemSize];
	int lcao_index=0;
	int twoLp1;
	int ip, jp, kp;
	for(i=0; i<atomNum; i++){
		for(j=0; j<specNum; j++){
			if(strcmp(atoms[i], species[j][0])==0){
				specIndex=j;
				break;
			}
		}
		sprintf(groupName, "%d_%s", (i+1), atoms[i]);
		Group atomG(LCAOG.createGroup(groupName));
		for(j=0; j<strlen(species[specIndex][2]); j+=2){
			char orbitLabel=species[specIndex][2][j];
			char numLabel[1]={species[specIndex][2][j+1]};
			for(k=0; k<atoi(numLabel); k++){
				// printf("%s %c %d\n", atoms[i], orbitLabel, k);
				if(orbitLabel=='s'){
					twoLp1=1;
				}else if(orbitLabel=='p'){
					twoLp1=3;
				}else if(orbitLabel=='d'){
					twoLp1=5;
				}else if(orbitLabel=='f'){
					twoLp1=7;
				}
				if(spinPol_i!=1){
					sprintf(datasetName, "%c%d", orbitLabel, k);
					int digit=(spinPol_i==2)?4:2;
					double LCAO_ext[total_count][numBands_out][twoLp1][digit];
					for(l=0; l<twoLp1; l++){
						for(ip=0; ip<total_count; ip++){
							for(jp=0; jp<numBands_out; jp++){
								for(kp=0; kp<digit; kp++){
									LCAO_ext[ip][jp][l][kp]=LCAO_all[ip][jp][lcao_index+l][kp];
								}
							}
						}
					}
					lcao_index+=l;
					w_data_4d(atomG, datasetName, total_count, numBands_out, twoLp1, digit, (double****)LCAO_ext);
				}else{
					// spinPol_i=1
					double LCAOUp_ext[total_count][numBands_out][twoLp1][2];
					double LCAODn_ext[total_count][numBands_out][twoLp1][2];
					for(l=0; l<twoLp1; l++){
						for(ip=0; ip<total_count; ip++){
							for(jp=0; jp<numBands_out; jp++){
								for(kp=0; kp<2; kp++){
									LCAOUp_ext[ip][jp][l][kp]=LCAOUp_all[ip][jp][lcao_index+l][kp];
									LCAODn_ext[ip][jp][l][kp]=LCAODn_all[ip][jp][lcao_index+l][kp];
								}
							}
						}
					}
					lcao_index+=l;
					sprintf(datasetName, "%c%dUp", orbitLabel, k);
					w_data_4d(atomG, datasetName, total_count, numBands_out, twoLp1, 2, (double****)LCAOUp_ext);
					sprintf(datasetName, "%c%dDn", orbitLabel, k);
					w_data_4d(atomG, datasetName, total_count, numBands_out, twoLp1, 2, (double****)LCAODn_ext);
					
				}
			}
		}
		
	}

	// atts @/Output
	w_att_double(OutputG, "EF_Eh", EF_Eh);
	w_att_str(OutputG, "Spin", spinPol);
	
	outOpenMX.close();
	output.close();
	cout << "Finished writing data." << endl;
}


	
bool readLine(ifstream* fp, char* line_c, int bufSize){
	string line_s;
	if(getline(*fp, line_s)){
		int actual_line_length=line_s.length();
		if(actual_line_length>bufSize){
			printf("Error: line longer than %d characters\n", bufSize);
			return false;
		}
		line_s.copy(line_c, actual_line_length);
		line_c[actual_line_length]='\0';
		return true;
	}else{
		return false;
	}
}
