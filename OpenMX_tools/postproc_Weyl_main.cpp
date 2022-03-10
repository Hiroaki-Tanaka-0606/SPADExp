// SPADExp tools for OpenMX
// postproc_Weyl_main

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
void outProd(double* a, double* b, double* c);
double inProd(double* a, double* b);

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
	int bufSize=262144;

	// diff criterion
	double diff_threshold1=0.01; // kp coordinates

	// Hartree energy (eV)
	double Eh_eV=27.2114;
	
	// variables
	int i, j, k, l;
	int dimension;
	double x_frac[3];
	double y_frac[3];
	double z_frac[3];
	double x_orth[3];
	double y_orth[3];
	double z_orth[3];
	double xMin, xMax;
	double yMin, yMax;
	double zMin, zMax;
	int x_count;
	int y_count;
	int z_count;
	double cell[3][3];
	double rec_cell[3][3];

	char sysName[valSize+1];
	char currentDir[valSize+1];
	int maxN;
	int minN;

	if(load_str("System.Name", sysName, valSize)==1){
		printf("System.Name is %s\n", sysName);
	}else{
		return 0;
	}

	if(load_str("System.CurrentDirectory", currentDir, valSize)==1){
		printf("System.CurrentDirectory is %s\n", currentDir);
	}else{
		strcpy(currentDir, "./");
	}


	if(load_int("calcPSF.Weyl.dimension", &dimension)==1){
		printf("calcPSF.Weyl.dimension is %d\n", dimension);
		if(!(1<=dimension && dimension<=3)){
			printf("Error: invalid dimension\n");
			return 0;
		}	
	}else{
		return 0;
	}
	

	int line_range_start=find_str("<calcPSF.Weyl.range");
	int line_range_end=find_str("calcPSF.Weyl.range>");

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
			range=get_line(line_range_start+3);
			if(dimension>2){
				if(sscanf(range, "%lf %lf %lf %lf %lf %d", &z_frac[0], &z_frac[1], &z_frac[2], &zMin, &zMax, &z_count)!=6){
					printf("Error in parsing z range");
					return 0;
				}
				if(y_count<1){
					printf("Error: z_count must be positive\n");
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
	int total_count=x_count;
	if(dimension>1){
		total_count*=y_count;
	}
	if(dimension>2){
		total_count*=z_count;
	}

	if(load_int("calcPSF.Weyl.minN", &minN)==1){
		printf("calcPSF.Weyl.minN is %d\n", minN);
	}else{
		return 0;
	}

	if(load_int("calcPSF.Weyl.maxN", &maxN)==1){
		printf("calcPSF.Weyl.maxN is %d\n", maxN);
	}else{
		return 0;
	}
  
	int line_bandCell_start=find_str("<Band.KPath.UnitCell");
	int line_bandCell_end=find_str("Band.KPath.UnitCell>");

	int line_cell_start=find_str("<Atoms.UnitVectors");
	int line_cell_end=find_str("Atoms.UnitVectors>");


	if(line_bandCell_start>=0 && line_bandCell_end>=0){
		line_cell_start=line_bandCell_start;
		line_cell_end=line_bandCell_end;
	}
	
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

	printf("Unit vectors (real space)\n");
	for(i=0; i<3; i++){
	 	printf("a_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), cell[i][0], cell[i][1], cell[i][2]);
	}	// calculate the reciprocal unit cell
	double vBuf[3];
	outProd(cell[1], cell[2], vBuf);
	double det=inProd(cell[0], vBuf);

	outProd(cell[1], cell[2], vBuf);
	for(i=0; i<3; i++){
		rec_cell[0][i]=2*M_PI*vBuf[i]/det;
	}
	outProd(cell[2], cell[0], vBuf);
	for(i=0; i<3; i++){
		rec_cell[1][i]=2*M_PI*vBuf[i]/det;
	}
	outProd(cell[0], cell[1], vBuf);
	for(i=0; i<3; i++){
		rec_cell[2][i]=2*M_PI*vBuf[i]/det;
	}
	
	printf("Unit vectors (reciprocal space)\n");
	for(i=0; i<3; i++){
	 	printf("b_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), rec_cell[i][0], rec_cell[i][1], rec_cell[i][2]);
	}
	printf("\n");

	for(i=0; i<3; i++){
		x_orth[i]=0;
		y_orth[i]=0;
		z_orth[i]=0;
		for(j=0; j<3; j++){
			x_orth[i]+=rec_cell[j][i]*x_frac[j];
			y_orth[i]+=rec_cell[j][i]*y_frac[j];
			z_orth[i]+=rec_cell[j][i]*z_frac[j];
		}
	}

	double r_X=sqrt(inProd(x_orth, x_orth));
	double r_Y=sqrt(inProd(y_orth, y_orth));
	double r_Z=sqrt(inProd(z_orth, z_orth));
	double dkx=0;
	double dky=0;
	double dkz=0;
	if(x_count>1){
		dkx=r_X*(xMax-xMin)/(x_count-1);
	}
	if(dimension>1 && y_count>1){
		dky=r_Y*(yMax-yMin)/(y_count-1);
	}
	if(dimension>2 && z_count>1){
		dkz=r_Z*(zMax-zMin)/(z_count-1);
	}
	
	cout << "Loading from the input file finished." << endl;
	cout << endl;

	// open the output HDF5 file
	char outFilePath[pathSize];
	sprintf(outFilePath, "%s%s.Band.hdf5", currentDir, sysName);
	printf("Open %s\n", outFilePath);

	H5File output(outFilePath, H5F_ACC_TRUNC);

	// att Datetime @root
	time_t datetime_now=time(NULL);
	struct tm *timeptr=localtime(&datetime_now);
	char time_str[valSize+1];
	strftime(time_str, valSize, "%Y-%m-%d %H:%M:%S", timeptr);
	Group rootG(output.openGroup("/"));
	w_att_str(rootG, "Datetime", time_str);
	w_att_int(rootG, "Dimension", dimension);
	w_att_1d(rootG, "Xvector", 3, x_frac);
	w_att_1d(rootG, "Yvector", 3, y_frac);
	w_att_1d(rootG, "Zvector", 3, z_frac);
	w_att_int(rootG, "Xcount", x_count);
	w_att_int(rootG, "Ycount", y_count);
	w_att_int(rootG, "Zcount", z_count);
	w_att_int(rootG, "Count", total_count);
	double xRange[2]={xMin, xMax};
	w_att_1d(rootG, "Xrange", 2, xRange);
	double yRange[2]={yMin, yMax};
	w_att_1d(rootG, "Yrange", 2, yRange);
	double zRange[2]={zMin, zMax};
	w_att_1d(rootG, "Zrange", 2, zRange);

	// minN and maxN: the index starts from zero
	w_att_int(rootG, "minN", minN);
	w_att_int(rootG, "maxN", maxN);

	if(dimension==1){
		double offset1[1]={r_X*xMin};
		double delta1[1]={dkx};
		int size1[1]={x_count};
		w_att_1d(rootG, "Offset", 1, offset1);
		w_att_1d(rootG, "Delta", 1, delta1);
		w_att_1i(rootG, "Size", 1, size1);
	}else if(dimension==2){
		double offset2[2]={r_X*xMin, r_Y*yMin};
		double delta2[2]={dkx, dky};
		int size2[2]={x_count, y_count};
		w_att_1d(rootG, "Offset", 2, offset2);
		w_att_1d(rootG, "Delta", 2, delta2);
		w_att_1i(rootG, "Size", 2, size2);
	}else if(dimension==3){
		double offset3[3]={r_X*xMin, r_Y*yMin, r_Z*zMin};
		double delta3[3]={dkx, dky, dkz};
		int size3[3]={x_count, y_count, z_count};
		w_att_1d(rootG, "Offset", 3, offset3);
		w_att_1d(rootG, "Delta", 3, delta3);
		w_att_1i(rootG, "Size", 3, size3);
	}

	// open the output of OpenMX
	char outBandPath[pathSize];
	sprintf(outBandPath, "%s%s.Band", currentDir, sysName);
	printf("Open %s\n", outBandPath);

	int numBands_out=maxN-minN+1;
	double* Band_alloc=new double[total_count*numBands_out];
	double** Band=new double*[total_count];
	for(i=0; i<total_count; i++){
		Band[i]=&Band_alloc[i*numBands_out];
	}
  
	ifstream outBand(outBandPath);
	string line;
	char line_c[bufSize+1];
	int sscanf_result;
	
	// row 1: numBands numSpin EF(Eh)
	int numBands;
	int numSpin;
	double EF_Eh;
	readLine(&outBand, line_c, bufSize);
	sscanf_result=sscanf(line_c, "%d %d %lf", &numBands, &numSpin, &EF_Eh);
	if(sscanf_result!=3){
		printf("Error in parsing row 1");
		return 0;
	}
	if(numSpin!=0){
		printf("Error: this program does not accept numSpin!=0 data");
		return 0;
	}
	
	// row 2: unit cell (real space), skip
	readLine(&outBand, line_c, bufSize);
	// row 3: number of Kpath
	readLine(&outBand, line_c, bufSize);
	int numKpath;
	sscanf_result=sscanf(line_c, "%d", &numKpath);
	if(sscanf_result!=1){
		printf("Error in parsing row 3");
		return 0;
	}
	// row 4~: Kpath (numKpath rows), skip
	for(i=0; i<numKpath; i++){
		readLine(&outBand, line_c, bufSize);
	}

	// band dispersion
	for(i=0; i<total_count; i++){
		// top row: numBands kx ky kz, skip
		readLine(&outBand, line_c, bufSize);
		// bottom row: eigenvalues (Eh)
		readLine(&outBand, line_c, bufSize);
		char sscanf_format1[bufSize+1];
		char sscanf_format2[bufSize+1];
		// printf("k=%d\n", i);
		for(j=0; j<bufSize; j++){
			sscanf_format1[j]='\0';
			sscanf_format2[j]='\0';
		}
		// prepare the sscanf format to load #minN-1 eigenvalue
		// the index starts from zero
		strcpy(sscanf_format1, "%lf");
		for(j=0; j<minN; j++){
			sprintf(sscanf_format2, "%%*lf %s", sscanf_format1);
			strcpy(sscanf_format1, sscanf_format2);
		}
		for(j=minN; j<=maxN; j++){
			double eigen;
			sscanf_result=sscanf(line_c, sscanf_format1, &eigen);
			if(sscanf_result!=1){
				printf("Error in parsing eigenvalue k=%d, j=%d\n", i, j);
				printf("%s\n", sscanf_format1);
				return 0;
			}
			double eigen_eV=(eigen-EF_Eh)*Eh_eV;
			Band[i][j-minN]=eigen_eV;
			// prepare the format for the next eigenvalue
			sprintf(sscanf_format2, "%%*lf %s", sscanf_format1);
			strcpy(sscanf_format1, sscanf_format2);
		}
		
	}
	/*
	for(i=0; i<total_count; i++){
		for(j=minN; j<=maxN; j++){
			printf("%11.5f ", Band[i][j-minN]);
		}
		printf("\n");
		}*/

	if(dimension==1){
		w_data_2d(rootG, "Band", x_count, numBands_out, (double**)&Band[0][0]);
	}else if(dimension==2){
		w_data_3d(rootG, "Band", x_count, y_count, numBands_out, (double***)&Band[0][0]);
	}else if(dimension==3){
		w_data_4d(rootG, "Band", x_count, y_count, z_count, numBands_out, (double****)&Band[0][0]);
	}

	double minEigen[numBands_out];
	double maxEigen[numBands_out];
	for(i=0; i<numBands_out; i++){
		for(j=0; j<total_count; j++){
			if(j==0){
				minEigen[i]=Band[j][i];
				maxEigen[i]=Band[j][i];
				continue;
			}
			if(Band[j][i]<minEigen[i]){
				minEigen[i]=Band[j][i];
			}
			if(Band[j][i]>maxEigen[i]){
				maxEigen[i]=Band[j][i];
			}
		}
		printf("Band #%4d: min=%11.5f max=%11.5f (eV)\n", i+minN, minEigen[i], maxEigen[i]);
	}
	
	double minDiff[numBands_out-1];
	int minDiff_index[numBands_out-1];
	for(i=0; i<numBands_out-1; i++){
		for(j=0; j<total_count; j++){
			double diff=Band[j][i+1]-Band[j][i];
			if(j==0){
				minDiff[i]=diff;
				minDiff_index[i]=j;
				continue;
			}
			if(diff<minDiff[i]){
				minDiff[i]=diff;
				minDiff_index[i]=j;
			}
		}
		printf("Minimum bandgap between #%4d and #%4d: %11.5f eV (%11.5f - %11.5f) at i=%6d ", i+minN, i+minN+1, minDiff[i], Band[minDiff_index[i]][i+1], Band[minDiff_index[i]][i], minDiff_index[i]);
		int minI=minDiff_index[i];
		int minIx, minIy, minIz;
		if(dimension==1){
			printf("(ix=%4d)\n", minI);
		}else if(dimension==2){
			minIx=minI/y_count;
			minIy=minI-minIx*y_count;
			printf("(ix=%4d, iy=%4d)\n", minIx, minIy);
		}else if(dimension==3){
			minIx=minI/(y_count*z_count);
			int minI2=minI-minIx*(y_count*z_count);
			minIy=minI2/z_count;
			minIz=minI2-minIy*z_count;
			printf("(ix=%4d, iy=%4d, iz=%4d)\n", minIx, minIy, minIz);
		}
		
	}

	cout << "Loading from the output file finished." << endl << endl;
	
	return 0;
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


void outProd(double* a, double* b, double* c){
	// c=a\times b
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

double inProd(double* a, double* b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}
