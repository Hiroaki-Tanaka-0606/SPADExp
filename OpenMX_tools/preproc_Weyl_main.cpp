// calcPSF tools for OpenMX
// preproc_Weyl_main

// include libraries from include path
#include <iostream>
#include <fstream>
#define _USE_MATH_DEFINES
#include <cmath>

// include libraries from current directory
#include "inputTools.hpp"

// namespace
using namespace std;

void outProd(double* a, double* b, double* c);
double inProd(double* a, double* b);
void calcKp(double p, double q, double r, double* X, double* Y, double* Z, double* output);
double orth_threshold=0.001;

int main(int argc, const char** argv){
	cout << "Starting preproc tool for OpenMX..." << endl;
	// argv[1] should be an input file name
	if(argc<3){
		printf("Usage: %s input_file output_file\n", argv[0]);
		return 0;
	}

	ifstream input(argv[1]);
	if(!input){
		printf("Error: cannot open %s\n", argv[0]);
		return 0;
	}
	load_file(&input);
	input.close();

	// variables
	int i, j;
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

	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", x_frac[0], x_frac[1], x_frac[2], xMin, xMax, x_count);
	// printf("(%.3f, %.3f, %.3f) [%.3f, %.3f] %d\n", y_frac[0], y_frac[1], y_frac[2], yMin, yMax, y_count);

	int line_bandCell_start=find_str("<Band.KPath.UnitCell");
	int line_bandCell_end=find_str("Band.KPath.UnitCell>");

	int line_cell_start=find_str("<Atoms.UnitVectors");
	int line_cell_end=find_str("Atoms.UnitVectors>");

	if(line_bandCell_start>=0 && line_bandCell_end>=0){
		line_cell_start=line_bandCell_start;
		line_cell_end=line_bandCell_end;
	}

	if(line_cell_end-line_cell_start-1!=3){
		printf("Cell error");
		return 0;
	}

	for(i=0; i<3; i++){
		char* cell_str=get_line(line_cell_start+1+i);
		if(sscanf(cell_str, "%lf %lf %lf", &cell[i][0], &cell[i][1], &cell[i][2])!=3){
			printf("Error in parsing cell");
			return 0;
		}
	}

	printf("Unit vectors (real space)\n");
	for(i=0; i<3; i++){
	 	printf("a_%d = (%8.3f, %8.3f, %8.3f)\n", (i+1), cell[i][0], cell[i][1], cell[i][2]);
	}

	// calculate the reciprocal unit cell
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

	printf("X = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", x_orth[0], x_orth[1], x_orth[2], r_X);
	printf("Y = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", y_orth[0], y_orth[1], y_orth[2], r_Y);
	printf("Z = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", z_orth[0], z_orth[1], z_orth[2], r_Z);

	double XY=inProd(x_orth, y_orth);
	double YZ=inProd(y_orth, z_orth);
	double ZX=inProd(z_orth, x_orth);
	
	if(abs(XY)>orth_threshold){
		printf("Error: X and Y are not orthogonal\n");
		return 0;
	}
	if(abs(YZ)>orth_threshold){
		printf("Error: Y and Z are not orthogonal\n");
		return 0;
	}
	if(abs(ZX)>orth_threshold){
		printf("Error: Z and X are not orthogonal\n");
		return 0;
	}

	printf("X, Y, and Z are orthogonal\n\n");

	// prepare the output file
	// neglect Band.dispersion, Band.Nkpath, and <Band.kpath> block
	ofstream output(argv[2]);
	copyInput_Weyl(&output);

	// Band.dispersion and Band.Nkpath
	output << endl;
	output << "#" << endl;
	output << "# Band output" << endl;
	output << "# written by preproc_Weyl.o" << endl;
	output << "#" << endl;
	output << endl;
	output << "Band.dispersion on" << endl;
	output << "Band.Nkpath ";
	if(dimension==1){
		output << 1 << endl;
	}else if(dimension==2){
		output << x_count << endl;
	}else if(dimension==3){
		output << (x_count*y_count) << endl;
	}
		
	double p, q;
	double dp, dq;
	double kp[3];
	double r_Kp;
	double pathStart[3];
	double pathEnd[3];
	char writeBuf[1024];
	if(x_count<2){
		dp=0;
	}else{
		dp=(xMax-xMin)/(x_count-1);
	}
	if(y_count<2){
		dq=0;
	}else{
		dq=(yMax-yMin)/(y_count-1);
	}
	output << "<Band.kpath" << endl;

	int current_path=0;
	if(dimension==3){
		for(i=0; i<x_count; i++){
			for(j=0; j<y_count; j++){
				p=xMin+dp*i;
				q=yMin+dq*j;
				// cout << p << " " << q << endl;
				calcKp(p, q, (double)zMin, x_frac, y_frac, z_frac, pathStart);
				calcKp(p, q, (double)zMax, x_frac, y_frac, z_frac, pathEnd);
				sprintf(writeBuf, "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f p%d p%d", z_count, pathStart[0], pathStart[1], pathStart[2], pathEnd[0], pathEnd[1], pathEnd[2], current_path*2, current_path*2+1);
				output << writeBuf << endl;
				current_path++;
			}
		}
	}else if(dimension==2){
		q=0;
		for(i=0; i<x_count; i++){
			p=xMin+dp*i;
			// cout << p << " " << q << endl;
			calcKp(p, (double)yMin, 0, x_frac, y_frac, z_frac, pathStart);
			calcKp(p, (double)yMax, 0, x_frac, y_frac, z_frac, pathEnd);
			sprintf(writeBuf, "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f p%d p%d", y_count, pathStart[0], pathStart[1], pathStart[2], pathEnd[0], pathEnd[1], pathEnd[2], current_path*2, current_path*2+1);
			output << writeBuf << endl;
			current_path++;
		}
	}else if(dimension==1){
		calcKp((double)xMin, 0, 0, x_frac, y_frac, z_frac, pathStart);
		calcKp((double)xMax, 0, 0, x_frac, y_frac, z_frac, pathEnd);
		sprintf(writeBuf, "%4d %11.5f %11.5f %11.5f %11.5f %11.5f %11.5f p%d p%d", x_count, pathStart[0], pathStart[1], pathStart[2], pathEnd[0], pathEnd[1], pathEnd[2], current_path*2, current_path*2+1);
		output << writeBuf << endl;
		current_path++;
	}
	output << "Band.kpath>" << endl;

	cout << "Preparation finished." << endl;
	
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

void calcKp(double p, double q, double r, double* X, double* Y, double* Z, double* output){
	// calculate pX+qY+rZ
	int i;
	for(i=0; i<3; i++){
		output[i]=p*X[i]+q*Y[i]+r*Z[i];
	}
}
