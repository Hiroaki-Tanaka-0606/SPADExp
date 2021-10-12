// calcPSF tools for OpenMX
// preproc_main

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
double calcKp(double p, double q, double r_X, double r_Y, double r_O,
						double* origin_frac, double* x_frac, double* y_frac, bool curved, double* kp);
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
	double rec_cell[3][3];

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
		origin_orth[i]=0;
		x_orth[i]=0;
		y_orth[i]=0;
		for(j=0; j<3; j++){
			origin_orth[i]+=rec_cell[j][i]*origin_frac[j];
			x_orth[i]+=rec_cell[j][i]*x_frac[j];
			y_orth[i]+=rec_cell[j][i]*y_frac[j];
		}
	}

	double r_O=sqrt(inProd(origin_orth, origin_orth));
	double r_X=sqrt(inProd(x_orth, x_orth));
	double r_Y=sqrt(inProd(y_orth, y_orth));

	printf("O = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", origin_orth[0], origin_orth[1], origin_orth[2], r_O);
	printf("X = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", x_orth[0], x_orth[1], x_orth[2], r_X);
	printf("Y = (%8.3f, %8.3f, %8.3f) length = %8.3f\n", y_orth[0], y_orth[1], y_orth[2], r_Y);
	double OX=inProd(origin_orth, x_orth);
	double XY=inProd(x_orth, y_orth);
	double YO=inProd(y_orth, origin_orth);
	if(abs(OX)>orth_threshold){
		printf("Error: O and X are not orthogonal\n");
		return 0;
	}
	if(abs(XY)>orth_threshold){
		printf("Error: X and Y are not orthogonal\n");
		return 0;
	}
	if(abs(YO)>orth_threshold){
		printf("Error: Y and O are not orthogonal\n");
		return 0;
	}

	printf("O, X, and Y are orthogonal\n\n");

	// prepare the output file
	// neglect MO.fileout, MO.Nkpoint, and <MO.kpoint> block
	ofstream output(argv[2]);
	copyInput(&output);

	// MO.fileout and MO.Nkpoint
	output << endl;
	output << "#" << endl;
	output << "# MO output" << endl;
	output << "# written by preproc.o" << endl;
	output << "#" << endl;
	output << endl;
	output << "MO.fileout on" << endl;
	output << "MO.Nkpoint ";
	if(dimension==1){
		output << x_count << endl;
	}else{
		output << (x_count*y_count) << endl;
	}

	// K=pX+qY+rO
	// r=1 if not curved, else so that |K| = r_O
	double p, q;
	double dp, dq;
	double kp[3];
	double r_Kp;
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
	output << "<MO.kpoint" << endl;
	cout << "kpoints in fractional coordinates" << endl;
	if(dimension==2){
		for(j=0; j<y_count; j++){
			for(i=0; i<x_count; i++){
				p=xMin+dp*i;
				q=yMin+dq*j;
				// cout << p << " " << q << endl;
				r_Kp=calcKp(p, q, r_X, r_Y, r_O, origin_frac, x_frac, y_frac, curved, kp);
				sprintf(writeBuf, "%11.5f %11.5f %11.5f # r = %11.5f", kp[0], kp[1], kp[2], r_Kp);
				output << writeBuf << endl;					
			}
		}
	}else{
		q=0;
		for(i=0; i<x_count; i++){
			p=xMin+dp*i;
			// cout << p << " " << q << endl;
			r_Kp=calcKp(p, q, r_X, r_Y, r_O, origin_frac, x_frac, y_frac, curved, kp);
			sprintf(writeBuf, "%11.5f %11.5f %11.5f", kp[0], kp[1], kp[2]);
			output << writeBuf << endl;
			printf("%s # r = %11.5f\n", writeBuf, r_Kp);
		}
	}
	output << "MO.kpoint>" << endl;

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

double calcKp(double p, double q, double r_X, double r_Y, double r_O,
						double* origin_frac, double* x_frac, double* y_frac, bool curved, double* kp){
	double r;
	if(!curved){
		r=1.0;
	}else{
		// (p*r_X)^2+(q*r_Y)^2+(r*r_O)^2=r_O^2, 0<r
		r=sqrt(pow(r_O, 2)-pow(p*r_X, 2)-pow(q*r_Y, 2))/r_O;
	}
	double r_Kp=sqrt(pow(p*r_X, 2)+pow(q*r_Y, 2)+pow(r*r_O, 2));
	int i;
	for(i=0; i<3; i++){
		kp[i]=p*x_frac[i]+q*y_frac[i]+r*origin_frac[i];
	}
	return r_Kp;
	
}
