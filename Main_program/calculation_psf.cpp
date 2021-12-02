// calcPSF photoemission structure factor calculation

#include <iostream>
#include <complex>
#include <chrono>
#include <string>
#include <omp.h>
#include <H5Cpp.h>
#include <algorithm>
#include <cstring>
#include "log.hpp"
#include "variables_ext.hpp"
#include "HDF5_tools.hpp"
#include "convert_LCAO.hpp"

using namespace std;
using namespace H5;
// test program of BLAS and OMP parallelization
extern "C"{
	double ddot_(
							 int* N,
							 double* CX,
							 int* INCX,
							 double* CY,
							 int* INCY);
	// in case of complex<double>, the result is the first argument
	void zdotu_(
												 complex<double>* result,
												 int* N,
												 complex<double>* CX,
												 int* INCX,
												 complex<double>* CY,
												 int* INCY);
}

complex<double> zdot(int* N, complex<double>* X, complex<double>* Y){
	int INC=1;
  complex<double> result;
	zdotu_(&result, N, &X[0], &INC, &Y[0], &INC);
	return result;
}

complex<double> zdot2(int* N, complex<double>* X, complex<double>* Y){
	complex<double> ret=complex<double>(0, 0);
	int i;
	for(i=0; i<*N; i++){
		ret+=X[i]*Y[i];
	}
	return ret;
}

double ddot(int* N, double* X, double* Y){
	int INC=1;
	return ddot_(N, &X[0], &INC, &Y[0], &INC);
}

void calculate_PSF(){
	// 1. load input
	write_log((char*)"----Load input data----");
	char* sprintf_buffer=new char[Log_length+1];
	string fileName(PS_input_file);
	H5File input(fileName, H5F_ACC_RDONLY);
	int size1, size2, size3, size4, length;
	int i, j;
	
	Group InputG(input.openGroup("/Input"));
	/// 1-1. Atomic.Species
	Group SpecG(InputG.openGroup("Atomic.Species"));
	//// Length (attribute)
	int atom_spec_length=r_att_int(SpecG, "Length");
	//// Labels
	s_data_1c(SpecG, "Labels", &size1, &length);
	char atom_spec_label[size1][length];
	r_data_1c(SpecG, "Labels", size1, length, (char**) atom_spec_label);
	//// Orbitals
	s_data_1c(SpecG, "Orbits", &size1, &length);
	char atom_spec_orbits[size1][length];
	r_data_1c(SpecG, "Orbits", size1, length, (char**) atom_spec_orbits);
	//// PAOs
	s_data_1c(SpecG, "PAOs", &size1, &length);
	char atom_spec_PAO[size1][length];
	r_data_1c(SpecG, "PAOs", size1, length, (char**) atom_spec_PAO);
	//// Pseudopotentials
	s_data_1c(SpecG, "Pseudopotentials", &size1, &length);
	char atom_spec_PS[size1][length];
	r_data_1c(SpecG, "Pseudopotentials", size1, length, (char**) atom_spec_PS);

	write_log((char*)"----Atomic species----");
	for(i=0; i<atom_spec_length; i++){
		sprintf(sprintf_buffer, "%32s : %s&%s, %s", atom_spec_label[i], atom_spec_PAO[i], atom_spec_PS[i], atom_spec_orbits[i]);
		write_log(sprintf_buffer);
	}

	/// 1-2. Atoms.SpeciesAndCoordinates
	Group AtomG(InputG.openGroup("Atoms.SpeciesAndCoordinates"));
	//// Length (attribute)
	int atom_length=r_att_int(AtomG, "Length");
	//// Unit (attribute): ang, au, or frac
	string atom_unit=r_att_str(AtomG, "Unit");
	char atom_unit_c[atom_unit.length()+1];
	atom_unit.copy(atom_unit_c, atom_unit.length());
	atom_unit_c[atom_unit.length()]='\0';
	//// Labels
	s_data_1c(AtomG, "Labels", &size1, &length);
	char atom_labels[size1][length];
	r_data_1c(AtomG, "Labels", size1, length, (char**) atom_labels);
	//// coordinates (in unit of atom_unit)
	s_data_2d(AtomG, "Coordinates", &size1, &size2);
	double atom_coordinates[size1][size2];
	r_data_2d(AtomG, "Coordinates", size1, size2, (double**) atom_coordinates);

	write_log((char*)"----Atoms----");
	for(i=0; i<atom_length; i++){
		sprintf(sprintf_buffer, "%32s : (%8.3f, %8.3f, %8.3f)", atom_labels[i], atom_coordinates[i][0], atom_coordinates[i][1], atom_coordinates[i][2]);
		write_log(sprintf_buffer);
	}
	sprintf(sprintf_buffer, "in unit of %s", atom_unit_c);
	write_log(sprintf_buffer);

	/// 1-3. UnitCells
	Group UnitG(InputG.openGroup("UnitCells"));
	//// Unit (attribute): ang or au
	string unit=r_att_str(UnitG, "Unit");
	char unit_c[atom_unit.length()+1];
	unit.copy(unit_c, unit.length());
	unit[unit.length()]='\0';
	//// Bands (in unit of unit)
	s_data_2d(UnitG, "Bands", &size1, &size2);
	double band_cell[size1][size2];
	r_data_2d(UnitG, "Bands", size1, size2, (double**) band_cell);
	//// Atoms (in unit of unit)
	s_data_2d(UnitG, "Atoms", &size1, &size2);
	double atom_cell[size1][size2];
	r_data_2d(UnitG, "Atoms", size1, size2, (double**) atom_cell);

	write_log((char*)"----Unit cell for atoms----");
	for(i=0; i<3; i++){
		sprintf(sprintf_buffer, "%31s%1d = (%8.3f, %8.3f, %8.3f)", "a", (i+1), atom_cell[i][0], atom_cell[i][1], atom_cell[i][2]);
		write_log(sprintf_buffer);
	}
	sprintf(sprintf_buffer, "in unit of %s", unit_c);
	write_log(sprintf_buffer);
	write_log((char*)"----Unit cell for bands----");
	for(i=0; i<3; i++){
		sprintf(sprintf_buffer, "%31s%1d = (%8.3f, %8.3f, %8.3f)", "a", (i+1), band_cell[i][0], band_cell[i][1], band_cell[i][2]);
		write_log(sprintf_buffer);
	}
	sprintf(sprintf_buffer, "in unit of %s", unit_c);
	write_log(sprintf_buffer);

	/// 1-4. Kpath
	Group KpathG(InputG.openGroup("Kpath"));
	//// Dimension (attribute): 1 or 2
	int dimension=r_att_int(KpathG, "Dimension");
	//// Curved
	bool curved=r_att_bool(KpathG, "Curved");
	//// Origin, Xvector, Yvector in fractional unit of reciprocal unit vectors
	//// Xcount, Ycount, Xrange, Yrange
	double origin_frac[3];
	double x_frac[3];
	double y_frac[3];
	int x_count, y_count;
	double x_range[2];
	double y_range[2];
	r_att_1d(KpathG, "Origin", 3, (double*) origin_frac);
	r_att_1d(KpathG, "Xvector", 3, (double*) x_frac);
	x_count=r_att_int(KpathG, "Xcount");
	r_att_1d(KpathG, "Xrange", 2, (double*) x_range);
	if(dimension==2){
		r_att_1d(KpathG, "Yvector", 3, (double*) y_frac);
		y_count=r_att_int(KpathG, "Ycount");
		r_att_1d(KpathG, "Yrange", 2, (double*) y_range);
	}
	//// Count
	int total_count=r_att_int(KpathG, "Count");
	//// min and max band indicies in the input
	int min_N=r_att_int(KpathG, "minN");
	int max_N=r_att_int(KpathG, "maxN");
	//// k points
	double k_points[total_count][3];
	r_data_2d(KpathG, "Coordinates", total_count, 3, (double**)k_points);

	write_log((char*)"----k points----");
	sprintf(sprintf_buffer, "%32s = (%8.3f, %8.3f, %8.3f)", "Origin", origin_frac[0], origin_frac[1], origin_frac[2]);
	write_log(sprintf_buffer);
	sprintf(sprintf_buffer, "%32s = (%8.3f, %8.3f, %8.3f)", "X", x_frac[0], x_frac[1], x_frac[2]);
	write_log(sprintf_buffer);
	sprintf(sprintf_buffer, "%32s = %8.3f to %8.3f, %d points", "X range", x_range[0], x_range[1], x_count);
	write_log(sprintf_buffer);
	if(dimension==2){
		sprintf(sprintf_buffer, "%32s = (%8.3f, %8.3f, %8.3f)", "Y", y_frac[0], y_frac[1], y_frac[2]);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "%32s = %8.3f to %8.3f, %d points", "Y range", y_range[0], y_range[1], y_count);
		write_log(sprintf_buffer);
	}
	sprintf(sprintf_buffer, "%32s = %d", "Total number of points", total_count);
	write_log(sprintf_buffer);

	/// 1-5. output
	Group OutputG=input.openGroup("/Output");
	//// EF (attribute) in unit of Eh
	double EF_Eh=r_att_double(OutputG, "EF_Eh");
	//// Spin (attribute)
	string spin=r_att_str(OutputG, "Spin");
	transform(spin.begin(), spin.end(), spin.begin(), ::tolower);
	//// spin index: 0 (off), 1(on), 2(nc)
	int spin_i;
	if(spin==string("off")){
		spin_i=0;
	}else if(spin==string("on")){
		spin_i=1;
	}else if(spin==string("nc")){
		spin_i=2;
	}else{
		write_log((char*)"Error in spin");
		return;
	}
	//// Band dispersion
	double** band;
	double** band_up;
	double** band_dn;
	write_log((char*)"----Loading band dispersion----");
	// size of the first index is total_count
	int numBands; // size of the second index
	if(spin_i==0 || spin_i==2){
		// off or nc: use band
		s_data_2d(OutputG, "Band", &size1, &size2);
		numBands=size2;
		double* band_alloc=new double[size1*size2];
		r_data_2d(OutputG, "Band", size1, size2, (double**) band_alloc);
		band=new double*[size1];
		for(i=0; i<size1; i++){
			band[i]=&band_alloc[0]+size2*i;
		}
	}else{
		// on: use band_up and band_dn
		s_data_2d(OutputG, "BandUp", &size1, &size2);
		numBands=size2;
		double* band_up_alloc=new double[size1*size2];
		double* band_dn_alloc=new double[size1*size2];
		r_data_2d(OutputG, "BandUp", size1, size2, (double**) band_up_alloc);
		r_data_2d(OutputG, "BandDn", size1, size2, (double**) band_dn_alloc);
		band_up=new double*[size1];
		band_dn=new double*[size1];
		for(i=0; i<size1; i++){
			band_up[i]=&band_up_alloc[0]+size2*i;
			band_dn[i]=&band_dn_alloc[0]+size2*i;
		}
	}
	//// LCAO
	write_log((char*)"----Loading LCAO coefficients----");
	chrono::system_clock::time_point start=chrono::system_clock::now();
	
	Group LCAOG(OutputG.openGroup("LCAO"));

	int ia; // for atom
	int is; // for species
	int io; // for orbit letter index
	int im; // for multiplicity of the orbit label
	int is_search; // for searching species
	int ip; // for orbit index in LCAO[ia], s0, s1, p0, p1, ... or s0Up, s0Dn, s1Up, s1Dn, ...
	char groupName[32];
	char dsName[32];

	complex<double>****** LCAO=new complex<double>*****[atom_length];
	int numOrbits[atom_length];
	int* lList[atom_length];
	for(ia=0; ia<atom_length; ia++){
		for(is_search=0; is_search<atom_spec_length; is_search++){
			if(strcmp(atom_labels[ia], atom_spec_label[is_search])==0){
				is=is_search;
				break;
			}	
		}
		sprintf(groupName, "%d_%s", (ia+1), atom_labels[ia]);
		Group atomG(LCAOG.openGroup(groupName));
		int orbitLength=strlen(atom_spec_orbits[is]);
		int numOrbits_calc=0;
		for(io=0; io<orbitLength; io+=2){
			char numLabel[2]={atom_spec_orbits[is][io+1], '\0'};
			numOrbits_calc+=atoi(numLabel);
		}
		if(spin_i==0 || spin_i==2){
			numOrbits[ia]=numOrbits_calc;
		}else{
			numOrbits[ia]=2*numOrbits_calc;
		}
		LCAO[ia]=new complex<double>****[numOrbits[ia]];
		lList[ia]=new int[numOrbits[ia]];
		ip=0;
		for(io=0; io<orbitLength; io+=2){
			char orbitLabel=atom_spec_orbits[is][io];
			char numLabel[2]={atom_spec_orbits[is][io+1], '\0'};
			int numLabel_i=atoi(numLabel);
			int twoLp1;
			if(orbitLabel=='s'){
				twoLp1=1;
			}else if(orbitLabel=='p'){
				twoLp1=3;
			}else if(orbitLabel=='d'){
				twoLp1=5;
			}else if(orbitLabel=='f'){
				twoLp1=7;
			}else{
				write_log((char*)"Orbit label error");
				return;
			}
			for(im=0; im<numLabel_i; im++){
				// LCAO data: [total_count][numBands][twoLp1][digit]
				// digit=2 [re, im] in on and off, 4 [re_up, im_up, re_dn, im_dn] in nc
				// LCAO after conversion is complex<double>[total_count][numBands][twoLp1][digit/2]
				if(spin_i==0 || spin_i==2){
					sprintf(dsName, "%c%d", orbitLabel, im);
					// cout << dsName << endl;
					s_data_4d(atomG, dsName, &size1, &size2, &size3, &size4);
					if(size1!=total_count || size2!=numBands || size3!=twoLp1){
						write_log((char*)"LCAO size mismatch");
						return;
					}
					double* LCAO_raw=new double[size1*size2*size3*size4];
					r_data_4d(atomG, dsName, size1, size2, size3, size4, (double****) LCAO_raw);
					LCAO[ia][ip]=convert_LCAO(size1, size2, size3, size4, LCAO_raw);
					delete LCAO_raw;
					/*
					if(ia==0){
						for(int iii=0; iii<size3; iii++){
							for(int jjj=0; jjj<size4/2; jjj++){
								cout << LCAO[ia][ip][0][0][iii][jjj];
							}
							cout << endl;
						}
						}*/
					lList[ia][ip]=(twoLp1-1)/2;
					ip++;
				}else{
					sprintf(dsName, "%c%dUp", orbitLabel, im);
					s_data_4d(atomG, dsName, &size1, &size2, &size3, &size4);
					if(size1!=total_count || size2!=numBands || size3!=twoLp1){
						write_log((char*)"LCAO size mismatch");
						return;
					}
					double* LCAO_rawUp=new double[size1*size2*size3*size4];
					r_data_4d(atomG, dsName, size1, size2, size3, size4, (double****) LCAO_rawUp);
					LCAO[ia][ip]=convert_LCAO(size1, size2, size3, size4, LCAO_rawUp);
					/*
					if(ia==0){
						for(int iii=0; iii<size3; iii++){
							for(int jjj=0; jjj<size4/2; jjj++){
								cout << LCAO[ia][ip][0][0][iii][jjj];
							}
							cout << endl;
						}
						}*/
					lList[ia][ip]=(twoLp1-1)/2;
					ip++;
					delete LCAO_rawUp;
					
					sprintf(dsName, "%c%dDn", orbitLabel, im);
					s_data_4d(atomG, dsName, &size1, &size2, &size3, &size4);
					if(size1!=total_count || size2!=numBands || size3!=twoLp1){
						write_log((char*)"LCAO size mismatch");
						return;
					}
				  double* LCAO_rawDn=new double[size1*size2*size3*size4];
					r_data_4d(atomG, dsName, size1, size2, size3, size4, (double****) LCAO_rawDn);
					LCAO[ia][ip]=convert_LCAO(size1, size2, size3, size4, LCAO_rawDn);
					delete LCAO_rawDn;
					lList[ia][ip]=(twoLp1-1)/2;
					ip++;
							
				}
			}
		}
	}
	chrono::system_clock::time_point end=chrono::system_clock::now();
	double duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("LCAO load time: %.3f [ms]\n", duration);	
	return;
}







void PSF_test(){
	int N=128;
	complex<double> matrix[N][N];
	complex<double> vec[N];
	int i,j;
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			matrix[i][j]=complex<double>(i+2*j, 0);
		}
		vec[i]=complex<double>(i*i, 0);
	}
	/*
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%4.1f ", matrix[i][j].real());
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++){
		printf("%4.1f ", vec[i].real());
	}
	printf("\n");*/

	
	complex<double> dot_product;
	int tMax=10000;
	int t;
	chrono::system_clock::time_point start=chrono::system_clock::now();
	for(t=0; t<tMax; t++){
		for(i=0; i<N; i++){
			dot_product=zdot(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product.real());
		}
	}
	chrono::system_clock::time_point end=chrono::system_clock::now();
	double duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Time1: %.3f [ms]\n", duration);

	complex<double> dot_product2;
	start=chrono::system_clock::now();
#pragma omp parallel shared(N, matrix, vec) private(dot_product2)
#pragma omp for
	for(i=0; i<N; i++){
		for(t=0; t<tMax; t++){
			dot_product2=zdot(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product2.real());
		}
	}
	end=chrono::system_clock::now();
	duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	printf("Time2: %.3f [ms]\n", duration);

	return;

	/*
	start=clock();
	for(t=0; t<tMax; t++){
		for(i=0; i<N; i++){
			dot_product=zdot2(&N, matrix[i], vec);
			// printf("M[%d] * V = %4.1f\n", i, dot_product.real());
		}
	}
  end=clock();
	duration=(double)(end-start)/CLOCKS_PER_SEC*1000;
	printf("Time3: %.3f [ms]\n", duration);
	*/
	/*
	double dmatrix[N][N];
	double dvec[N];
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			dmatrix[i][j]=i+2*j;
		}
		dvec[i]=i*i;
	}
	for(i=0; i<N; i++){
		for(j=0; j<N; j++){
			printf("%4.1f ", dmatrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
	for(i=0; i<N; i++){
		printf("%4.1f ", dvec[i]);
	}
	printf("\n");
	
	double ddot_product;
	for(i=0; i<N; i++){
		ddot_product=ddot(&N, dmatrix[i], dvec);
		printf("M[%d] * V = %4.1f\n", i, ddot_product);
		}*/
	
}
