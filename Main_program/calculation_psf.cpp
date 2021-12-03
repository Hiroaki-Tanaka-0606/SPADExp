// calcPSF photoemission structure factor calculation

#include <iostream>
#include <complex>
#include <chrono>
#include <string>
#include <omp.h>
#include <H5Cpp.h>
#include <algorithm>
#include <cstring>
#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include "log.hpp"
#include "variables_ext.hpp"
#include "HDF5_tools.hpp"
#include "physical_tools.hpp"

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
	string file_name(PS_input_file);
	H5File input(file_name, H5F_ACC_RDONLY);
	int size1, size2, size3, size4, length;
	int i, j, k;

	int item_size=32;
	int val_size=128;
	
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
	//// analysis of orbitals
	int num_orbits[atom_spec_length];
	int* l_list[atom_spec_length];
	char** orbital_list[atom_spec_length];
	for(i=0; i<atom_spec_length; i++){
		int orbit_length=strlen(atom_spec_orbits[i]);
		int num_orbits_calc=0;
		for(j=0; j<orbit_length; j+=2){
			char num_label[2]={atom_spec_orbits[i][j+1], '\0'};
			num_orbits_calc+=atoi(num_label);
		}
		num_orbits[i]=num_orbits_calc;
		l_list[i]=new int[num_orbits_calc];
		orbital_list[i]=new char*[num_orbits_calc];
		int io=0;
		for(j=0; j<orbit_length; j+=2){
			char l_char=atom_spec_orbits[i][j];
			int l_int=0;
			if(l_char=='s'){
				l_int=0;
			}else if(l_char=='p'){
				l_int=1;
			}else if(l_char=='d'){
				l_int=2;
			}else if(l_char=='f'){
				l_int=3;
			}else{
				return;
			}
			char num_label[2]={atom_spec_orbits[i][j+1], '\0'};
			int num_label_i=atoi(num_label);
			for(int p=0; p<num_label_i; p++){
				l_list[i][io]=l_int;
				orbital_list[i][io]=new char[item_size];
				sprintf(orbital_list[i][io], "%c%d", l_char, p);
				io++;
			}						
		}
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
	transform(atom_unit.begin(), atom_unit.end(), atom_unit.begin(), ::tolower);
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
	transform(unit.begin(), unit.end(), unit.begin(), ::tolower);
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
	double** band; // [total_count][num_bands]
	double** band_up;
	double** band_dn;
	write_log((char*)"----Loading band dispersion----");
	// size of the first index is total_count
	int num_bands; // size of the second index
	if(spin_i==0 || spin_i==2){
		// off or nc: use band
		s_data_2d(OutputG, "Band", &size1, &size2);
		num_bands=size2;
		double* band_alloc=new double[size1*size2];
		r_data_2d(OutputG, "Band", size1, size2, (double**) band_alloc);
		band=new double*[size1];
		for(i=0; i<size1; i++){
			band[i]=&band_alloc[0]+size2*i;
		}
	}else{
		// on: use band_up and band_dn
		s_data_2d(OutputG, "BandUp", &size1, &size2);
		num_bands=size2;
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
	int io; // for orbital
	int is_search; // for searching species
	char group_name[item_size];
	char ds_name[item_size];

	complex<double>****** LCAO=new complex<double>*****[atom_length]; // [atom_length][num_orbits2][[[[data]]]]
	int atom_spec_index[atom_length];
	for(ia=0; ia<atom_length; ia++){
		for(is_search=0; is_search<atom_spec_length; is_search++){
			if(strcmp(atom_labels[ia], atom_spec_label[is_search])==0){
				is=is_search;
				atom_spec_index[ia]=is_search;
				break;
			}	
		}
		sprintf(group_name, "%d_%s", (ia+1), atom_labels[ia]);
		Group atomG(LCAOG.openGroup(group_name));
		int num_orbits2=num_orbits[is];
		if(spin_i==1){
			num_orbits2*=2;
		}
		LCAO[ia]=new complex<double>****[num_orbits2];
		for(io=0; io<num_orbits[is]; io++){
			int twoLp1=l_list[is][io]*2+1;
			// LCAO data: [total_count][numBands][twoLp1][digit]
			// digit=2 [re, im] in on and off, 4 [re_up, im_up, re_dn, im_dn] in nc
			// LCAO after conversion is complex<double>[total_count][numBands][twoLp1][digit/2]
			if(spin_i==0 || spin_i==2){
				sprintf(ds_name, "%s", orbital_list[is][io]);
				// cout << ds_name << endl;
				s_data_4d(atomG, ds_name, &size1, &size2, &size3, &size4);
				if(size1!=total_count || size2!=num_bands || size3!=twoLp1){
					write_log((char*)"LCAO size mismatch");
					return;
				}
				double* LCAO_raw=new double[size1*size2*size3*size4];
				r_data_4d(atomG, ds_name, size1, size2, size3, size4, (double****) LCAO_raw);
				LCAO[ia][io]=convert_LCAO(size1, size2, size3, size4, LCAO_raw);
				delete LCAO_raw;
				/*
					if(ia==0){
					for(int iii=0; iii<size3; iii++){
					for(int jjj=0; jjj<size4/2; jjj++){
					cout << LCAO[ia][io][0][0][iii][jjj];
					}
					cout << endl;
					}
					}*/
				}else{
					sprintf(ds_name, "%sUp", orbital_list[is][io]);
					s_data_4d(atomG, ds_name, &size1, &size2, &size3, &size4);
					if(size1!=total_count || size2!=num_bands || size3!=twoLp1){
						write_log((char*)"LCAO size mismatch");
						return;
					}
					double* LCAO_rawUp=new double[size1*size2*size3*size4];
					r_data_4d(atomG, ds_name, size1, size2, size3, size4, (double****) LCAO_rawUp);
					LCAO[ia][io*2]=convert_LCAO(size1, size2, size3, size4, LCAO_rawUp);
					/*
					if(ia==0){
						for(int iii=0; iii<size3; iii++){
							for(int jjj=0; jjj<size4/2; jjj++){
								cout << LCAO[ia][io][0][0][iii][jjj];
							}
							cout << endl;
						}
						}*/
					delete LCAO_rawUp;
					
					sprintf(ds_name, "%sDn", orbital_list[is][io]);
					s_data_4d(atomG, ds_name, &size1, &size2, &size3, &size4);
					if(size1!=total_count || size2!=num_bands || size3!=twoLp1){
						write_log((char*)"LCAO size mismatch");
						return;
					}
				  double* LCAO_rawDn=new double[size1*size2*size3*size4];
					r_data_4d(atomG, ds_name, size1, size2, size3, size4, (double****) LCAO_rawDn);
					LCAO[ia][io*2+1]=convert_LCAO(size1, size2, size3, size4, LCAO_rawDn);
					delete LCAO_rawDn;							
			}
		}
	}
	chrono::system_clock::time_point end=chrono::system_clock::now();
	double duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	sprintf(sprintf_buffer, "LCAO load time: %.3f [ms]", duration);
	write_log(sprintf_buffer);

	/*
	for(i=0; i<atom_length; i++){
		for(j=0; j<num_orbits[i]; j++){
			cout << orbital_list[i][j] << " ";
		}
		cout << endl;
		}*/
	
	
	// 2. load initial state
	write_log((char*)"----Load initial states----");

	double*** wfn_phi_rdr=new double**[atom_spec_length]; // [is][io][r]
	double** wfn_r=new double*[atom_spec_length]; // [is][r]
	char AO_group_name[item_size*2];
	int wfn_length[atom_spec_length];

	string AO_file_name(PS_AO_file);
	H5File AO_file(AO_file_name, H5F_ACC_RDONLY);
	
	for(is=0; is<atom_spec_length; is++){
		wfn_phi_rdr[is]=new double*[num_orbits[is]];
		sprintf(AO_group_name, "%s&%s", atom_spec_PAO[is], atom_spec_PS[is]);
		sprintf(sprintf_buffer, "----Open %s----", AO_group_name);
		write_log(sprintf_buffer);
		Group AOG(AO_file.openGroup(AO_group_name));
		wfn_length[is]=r_att_int(AOG, "length");
		wfn_r[is]=new double[wfn_length[is]];
		r_att_1d(AOG, "r", wfn_length[is], &wfn_r[is][0]);
		for(io=0; io<num_orbits[is]; io++){
			s_data_2d(AOG, orbital_list[is][io], &size1, &size2);
			if(size2!=wfn_length[is]){
				cout << "Size error" << endl;
				return;
			}
			double wfn_both[2][wfn_length[is]];
			r_data_2d(AOG, orbital_list[is][io], 2, wfn_length[is], (double**) wfn_both);
			wfn_phi_rdr[is][io]=new double[wfn_length[is]];
			for(j=0; j<wfn_length[is]-1; j++){
				if(strcmp(PS_initial_state, "PAO")==0){
					wfn_phi_rdr[is][io][j]=wfn_both[0][j]*wfn_r[is][j]*(wfn_r[is][j+1]-wfn_r[is][j]);
				}else{
					wfn_phi_rdr[is][io][j]=wfn_both[1][j]*wfn_r[is][j]*(wfn_r[is][j+1]-wfn_r[is][j]);
				}
			}
			wfn_phi_rdr[is][io][wfn_length[is]-1]=0;
		} 
	}

	// 3. calculate final states
	if(strcmp(PS_final_state, "Calc")==0){
		write_log((char*)"Caculation of the final state is not yet implemented");
		return;
	}
	
	// 4. psf calculation
	write_log((char*)"----PSF Calculation----");
  start=chrono::system_clock::now();
	int num_points_E=ceil((PS_E_max-PS_E_min)/PS_E_pixel+1);
	sprintf(sprintf_buffer, "%32s = %d", "Number of points along E", num_points_E);
	write_log(sprintf_buffer);

	/// 4-1. calculate atom position and k point in au (or au^-1) orthonormal basis
	double rec_cell[3][3];
	double op[3]; // for outer product
	double det; // for determinant
	/// reciprocal unit cell
	outer_product(band_cell[1], band_cell[2], op);
	det=inner_product(band_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[0][i]=2*M_PI*op[i]/det;
	}

	outer_product(band_cell[2], band_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[1][i]=2*M_PI*op[i]/det;
	}

	outer_product(band_cell[0], band_cell[1], op);
	for(i=0; i<3; i++){
		rec_cell[2][i]=2*M_PI*op[i]/det;
	}
	/// resclale to au or au^-1
	//// unit is already transformed to lowercase
	if(unit==string("ang")){
		for(i=0; i<3; i++){
			for(j=0; j<3; j++){
				rec_cell[i][j]*=au_ang;
				atom_cell[i][j]/=au_ang;
			}
		}
	}
	/// x and y vectors
	double x_vector[3]={0, 0, 0};
	double y_vector[3]={0, 0, 0};
	double x_length;
	double y_length;
	double dx_length;
	double dy_length;
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			x_vector[j]+=x_frac[i]*rec_cell[i][j];
			if(dimension==2){
				y_vector[j]+=y_frac[i]*rec_cell[i][j];
			}
		}
	}
	x_length=sqrt(inner_product(x_vector, x_vector));
	dx_length=x_length/x_count;
	if(dimension==2){
		y_length=sqrt(inner_product(y_vector, y_vector));
		dy_length=y_length/y_count;
	}

	/// atom positions
	//// atom_unit is already transformed to lowercase
	if(atom_unit==string("au")){
		// nothing to do
	}else if(atom_unit==string("ang")){
		for(i=0; i<atom_length; i++){
			for(j=0; j<3; j++){
				atom_coordinates[i][j]/=au_ang;
			}
		}
	}else if(atom_unit==string("frac")){
		for(i=0; i<atom_length; i++){
			double atom_coordinates_frac[3];
			for(j=0; j<3; j++){
				atom_coordinates_frac[j]=atom_coordinates[i][j];
				atom_coordinates[i][j]=0;
			}
			for(j=0; j<3; j++){
				for(k=0; k<3; k++){
					atom_coordinates[i][k]+=atom_coordinates_frac[j]*atom_cell[j][k];
				}
			}
		}
	}
	/*
	for(i=0; i<atom_length; i++){
		sprintf(sprintf_buffer, "%8.3f %8.3f %8.3f", atom_coordinates[i][0], atom_coordinates[i][1], atom_coordinates[i][2]);
		write_log(sprintf_buffer);
		}*/
	/// k points
	for(i=0; i<total_count; i++){
		double k_points_frac[3];
		for(j=0; j<3; j++){
			k_points_frac[j]=k_points[i][j];
			k_points[i][j]=0;
		}
		for(j=0; j<3; j++){
			for(k=0; k<3; k++){
				k_points[i][k]+=k_points_frac[j]*rec_cell[j][k];
			}
		}
	}

	/// 4-2. calculation of operator coefficients
	complex<double> Y_coeff[3];
	operator_coefficient(PS_polarization, PS_theta, PS_phi, Y_coeff);
	for(i=0; i<3; i++){
		sprintf(sprintf_buffer, "%29s[%1d] = %8.3f + %8.3f*i", "Y_coeff", i, Y_coeff[i].real(), Y_coeff[i].imag());
		write_log(sprintf_buffer);
	}
	

	/// 4-3. tail profile
	int tail_index=floor(PS_dE*PS_sigma_max/PS_E_pixel);
	double tail[tail_index+1];
	for(i=0; i<=tail_index; i++){
		tail[i]=1.0/(PS_dE*sqrt(2*M_PI))*exp(-1.0/2.0*(i*PS_E_pixel/PS_dE)*(i*PS_E_pixel/PS_dE));
	}
  
	/// 4-4. calculation for each k point
	double dispersion2[total_count][num_points_E]={0};
	int ik; // for k point
	int ib; // for band index
	int sp; // for spin: 0=Up, 1=Dn
	int sp_max=(spin_i==0 || spin_i==2) ? 1 : 2;
	int il; // for l
	int ir; // for r
	complex<double> m1jlp[5]={complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1), complex<double>(1, 0)};
	complex<double> Ylm_k[5][9];
	double Gaunt_arr[4][7][5][9]={0}; // [l][m+l][lp][mp+lp];
	int ip, jp;
	for(i=0; i<4; i++){
		for(j=0; j<7; j++){
			for(ip=0; ip<5; ip++){
				for(jp=0; jp<9; jp++){
					Gaunt_arr[i][j][ip][jp]=Gaunt(ip, jp-ip, i, j-i);
				}
			}
		}
	}
  end=chrono::system_clock::now();
  duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	sprintf(sprintf_buffer, "PSF preparation time: %.3f [ms]", duration);
	write_log(sprintf_buffer);
	//#pragma omp parallel firstprivate(Y_coeff, sp_max, num_bands, EF_Eh, Eh, PS_E_min, PS_E_pixel, num_points_E, tail_index,  m1jlp, Gaunt_arr, atom_length, atom_spec_index, num_orbits, spin_i,  atom_coordinates) private(ib, sp, Ylm_k, ia, is, io, il, ir, j)
	//#pragma omp for
	for(ik=0; ik<total_count; ik++){
		// cout << ik << endl;
		spherical_harmonics(k_points[ik], &Ylm_k[0][0]);
		double k_length=sqrt(inner_product(k_points[ik], k_points[ik]));
		/*
		for(int ii=0; ii<5; ii++){
			for(int jj=0; jj<9; jj++){
				char buffer[1024];
				sprintf(buffer, "%8.3f %8.3f", Ylm_k[ii][jj].real(), Ylm_k[ii][jj].imag());
				write_log(buffer);
			}
			write_log((char*)"");
			}*/
		for(sp=0; sp<sp_max; sp++){
			for(ib=0; ib<num_bands; ib++){
				double eigen;
				if(spin_i==0 || spin_i==2){
					eigen=(band[ik][ib]-EF_Eh)*Eh;
				}else if(sp==0){
					eigen=(band_up[ik][ib]-EF_Eh)*Eh;
				}else{
					eigen=(band_dn[ik][ib]-EF_Eh)*Eh;
				}
				
				int eigen_index=round((eigen-PS_E_min)/PS_E_pixel);
				if(eigen_index-tail_index>=num_points_E){
					break;
				}
				if(eigen_index+tail_index<0){
					continue;
				}
				complex<double> PSF_1(0, 0);
				complex<double> PSF_2(0, 0);
				for(ia=0; ia<atom_length; ia++){
					double kt=inner_product(k_points[ik], atom_coordinates[ia]);
					complex<double> atom_phase(cos(kt), -sin(kt));
					// cout << kt << endl;
					is=atom_spec_index[ia];
					int num_orbits2=num_orbits[is];
					if(spin_i==1){
						num_orbits2*=2;
					}
					double final_states[5][wfn_length[is]];
					if(strcmp(PS_final_state, "PW")==0){
						for(il=0; il<5; il++){
							for(ir=0; ir<wfn_length[is]; ir++){
								final_states[il][ir]=wfn_r[is][ir]*sp_bessel(il, wfn_r[is][ir]*k_length);
							}
						}
					}
					for(io=0; io<num_orbits2; io++){
						int io2=io; // for initial state, l_list
						if(spin_i==1){
							// io=0, 2, ... -> Up (sp=0), io=1, 3, ... -> Dn (sp=1)
							if(sp==0 && io%2==1){
								continue;
							}
							if(sp==1 && io%2==0){
								continue;
							}
							io2=io/2;
						}
						complex<double>** LCAO_use=LCAO[ia][io][ik][ib]; // [twoLp1][digit]
						// l: azimuthal quantum number of initial state
						// lp: of final state, lp=l+dl
						// m: magnetic quantum number of initial state
						// mpl: m+l (0, ..., 2l)
						// jp1: j+1
						// mp: of final state, mp=m+j
						// mpplp: mp+lp=m+j+lp
						int l=l_list[is][io2];
						int dl; // -1 or +1
						for(dl=-1; dl<=1; dl+=2){
							int lp=l+dl;
							if(lp<0){
								continue;
							}
							double radial_part=ddot(&wfn_length[is], &final_states[lp][0], &wfn_phi_rdr[is][io2][0]);
							int mpl;
							complex<double> coeff2_1(0, 0);
							complex<double> coeff2_2(0, 0);
							for(mpl=0; mpl<=2*l; mpl++){
								int m=mpl-l;
								int jp1St=max(-1, -(m+lp))+1; // include
								int jp1En=min(1, lp-m)+2; // not include
								int jp1;
								complex<double> coeff1(0, 0);
								for(jp1=jp1St; jp1<jp1En; jp1++){
									int mpplp=jp1-1+m+lp;
									coeff1+=Ylm_k[lp][mpplp]*Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1];
								}
								// cout << coeff1 << endl;
								coeff2_1+=LCAO_use[mpl][0]*coeff1;
								if(spin_i==2){
									coeff2_2+=LCAO_use[mpl][1]*coeff1;
								}
							}
							// cout << endl;
							PSF_1+=m1jlp[lp]*radial_part*atom_phase*coeff2_1;
							if(spin_i==2){
								PSF_2+=m1jlp[lp]*radial_part*atom_phase*coeff2_2;
							}
						}
					}
				}
				// cout << PSF_1 << endl;
				double PSF_norm=norm(PSF_1);
				if(spin_i==2){
					PSF_norm+=norm(PSF_2);
				}
				for(j=-tail_index; j<=tail_index; j++){
					if(eigen_index+j>=0 && eigen_index+j<num_points_E){
						dispersion2[ik][eigen_index+j]+=tail[abs(j)]*PSF_norm;
					}
				}
				
			} // end of for(ib)
		} // end of for(sp)
	} // end of for(ik) and omp
	
  end=chrono::system_clock::now();
  duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	sprintf(sprintf_buffer, "PSF calculation time: %.3f [ms]", duration);
	write_log(sprintf_buffer);

	double disp_max=0.0;
	for(i=0; i<total_count; i++){
		for(j=0; j<num_points_E; j++){
			if(disp_max<dispersion2[i][j]){
				disp_max=dispersion2[i][j];
			}
		}
	}

	sprintf(sprintf_buffer, "max intentity: %8.3e", disp_max);
	write_log(sprintf_buffer);

	/// 4-5 redimension (if dimension==2)
	double*** dispersion3;
	if(dimension==2){
		double* dispersion3_alloc=new double[total_count*num_points_E];
		dispersion3=new double**[x_count];
		for(i=0; i<x_count; i++){
			dispersion3[i]=new double*[y_count];
			for(j=0; j<y_count; j++){
				dispersion3[i][j]=&dispersion3_alloc[i*y_count*num_points_E+j*num_points_E];
			}
		}
		for(i=0; i<x_count; i++){
			for(j=0; j<y_count; j++){
				for(k=0; k<num_points_E; k++){
					dispersion3[i][j][k]=dispersion2[i+j*x_count][k];
				}
			}
		}
	}

	/// 5. output
	write_log((char*)"----Output----");
	H5File output(Output_file, H5F_ACC_TRUNC);

	// att Datetime @root
	time_t datetime_now=time(NULL);
	struct tm *timeptr=localtime(&datetime_now);
	char time_str[val_size+1];
	strftime(time_str, val_size, "%Y-%m-%d %H:%M:%S", timeptr);
	w_att_str(output.openGroup("/"), "Datetime", time_str);

	if(dimension==1){
		w_data_2d(output.openGroup("/"), "Dispersion", total_count, num_points_E, (double**) dispersion2);
	}else{
		w_data_3d(output.openGroup("/"), "Dispersion", x_count, y_count, num_points_E, (double***) &dispersion3[0][0][0]);
	}

	return;
}
