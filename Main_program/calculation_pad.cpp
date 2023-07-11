// SPADExp photoemission angular distribution calculation

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
#include <filesystem>
#include "log.hpp"
#include "variables_ext.hpp"
#include "HDF5_tools.hpp"
#include "physical_tools.hpp"
#include "calculation_phase_shift.hpp"
#include "calculation_atomic_wfn.hpp"
#include "setup.hpp"
#include "sphere_lebedev_rule.hpp"
#include "allocation_tools.hpp"

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
}

double ddot(int* N, double* X, double* Y){
	int INC=1;
	return ddot_(N, &X[0], &INC, &Y[0], &INC);
}


void calculate_PAD(){
	// 1. load input
	write_log((char*)"----Load input data----");
	char* sprintf_buffer=new char[Log_length+1];
	string file_name(PA_input_file);
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
	// cout << atom_unit_c << endl;
	// cout << atom_unit.length() << endl;
	// cout << atom_unit << endl;
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
	double kx_frac[3];
	double ky_frac[3];
	int kx_count, ky_count;
	double kx_range[2];
	double ky_range[2];
	r_att_1d(KpathG, "Origin", 3, (double*) origin_frac);
	r_att_1d(KpathG, "Xvector", 3, (double*) kx_frac);
	kx_count=r_att_int(KpathG, "Xcount");
	r_att_1d(KpathG, "Xrange", 2, (double*) kx_range);
	if(dimension==2){
		r_att_1d(KpathG, "Yvector", 3, (double*) ky_frac);
		ky_count=r_att_int(KpathG, "Ycount");
		r_att_1d(KpathG, "Yrange", 2, (double*) ky_range);
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
	sprintf(sprintf_buffer, "%32s = (%8.3f, %8.3f, %8.3f)", "kx", kx_frac[0], kx_frac[1], kx_frac[2]);
	write_log(sprintf_buffer);
	sprintf(sprintf_buffer, "%32s = %8.3f to %8.3f, %d points", "kx range", kx_range[0], kx_range[1], kx_count);
	write_log(sprintf_buffer);
	if(dimension==2){
		sprintf(sprintf_buffer, "%32s = (%8.3f, %8.3f, %8.3f)", "ky", ky_frac[0], ky_frac[1], ky_frac[2]);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "%32s = %8.3f to %8.3f, %d points", "ky range", ky_range[0], ky_range[1], ky_count);
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
			// LCAO[ia][io] data: [total_count][numBands][twoLp1][digit]
			// digit=2 [re, im] in on and off, 4 [re_up, im_up, re_dn, im_dn] in nc
			// LCAO[ia][io] after conversion is complex<double>[total_count][numBands][twoLp1][digit/2]
			// io=0, 2, ... -> Up (sp=0), io=1, 3, ... -> Dn (sp=1)
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
				delete[] LCAO_raw;
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
				delete[] LCAO_rawUp;
					
				sprintf(ds_name, "%sDn", orbital_list[is][io]);
				s_data_4d(atomG, ds_name, &size1, &size2, &size3, &size4);
				if(size1!=total_count || size2!=num_bands || size3!=twoLp1){
					write_log((char*)"LCAO size mismatch");
					return;
				}
				double* LCAO_rawDn=new double[size1*size2*size3*size4];
				r_data_4d(atomG, ds_name, size1, size2, size3, size4, (double****) LCAO_rawDn);
				LCAO[ia][io*2+1]=convert_LCAO(size1, size2, size3, size4, LCAO_rawDn);
				delete[] LCAO_rawDn;							
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
	/// interpolate wave functions	
	Group PotentialG(OutputG.openGroup("Potential"));
	int VKS_count_for_wfn[3];
	r_att_1i(PotentialG, "Count", 3, &VKS_count_for_wfn[0]);
	double wfn_interpolate_dz=sqrt(inner_product(atom_cell[0], atom_cell[0]))/VKS_count_for_wfn[0]*PA_interpolate_wfn_coef;
  
	double*** wfn_phi_rdr=new double**[atom_spec_length]; // [is][io][r]
	double*** wfn_phi_PAO=new double**[atom_spec_length]; // [is][io][r] for FPFS
	double*** wfn_phi_AO=new double**[atom_spec_length];  // [is][io][r] for FPFS
	double*** wfn_phi_dr=new double**[atom_spec_length];  // [is][io][r] for nonorth_term
	double*** wfn_phi=new double**[atom_spec_length];     // [is][io][r] for orth_correction
	double** wfn_r=new double*[atom_spec_length]; // [is][r] original values
	double** wfn_r_reduced=new double*[atom_spec_length];
	double** wfn_r_use=new double*[atom_spec_length];
	char AO_group_name[item_size*2];
	int wfn_length[atom_spec_length];
	int wfn_length_reduced[atom_spec_length];
	int wfn_length_use[atom_spec_length];
	int Z[atom_spec_length];
	bool empty_atoms[atom_spec_length];

	string AO_file_name(PA_AO_file);
	H5File AO_file(AO_file_name, H5F_ACC_RDONLY);
	
	for(is=0; is<atom_spec_length; is++){
		wfn_phi_rdr[is]=new double*[num_orbits[is]];
		wfn_phi_dr[is]=new double*[num_orbits[is]];
		wfn_phi_PAO[is]=new double*[num_orbits[is]];
		wfn_phi_AO[is]=new double*[num_orbits[is]];
		if(PA_orth_correction){
			wfn_phi[is]=new double*[num_orbits[is]];
		}
		sprintf(AO_group_name, "%s&%s", atom_spec_PAO[is], atom_spec_PS[is]);
		sprintf(sprintf_buffer, "----Open %s----", AO_group_name);
		write_log(sprintf_buffer);
		bool empty_atom=false;
		if(strcmp(atom_spec_PS[is], "E")==0){
			write_log((char*)"This orbital is for an empty atom");
			empty_atom=true;
			int is_for_read=-1;
			for(int is2=0; is2<atom_spec_length; is2++){
				if(is2==is){
					continue;
				}
				if(strcmp(atom_spec_PAO[is], atom_spec_PAO[is2])==0){
					is_for_read=is2;
					break;
				}
			}
			if(is_for_read==-1){
				write_log((char*)"Orbital not found from occupied atoms");
				return;
			}
			sprintf(AO_group_name, "%s&%s", atom_spec_PAO[is_for_read], atom_spec_PS[is_for_read]);
		}
		empty_atoms[is]=empty_atom;
		// printf("%s\n", AO_group_name);
		
		Group AOG(AO_file.openGroup(AO_group_name));
		// read wfn_r
		wfn_length[is]=r_att_int(AOG, "length");
		Z[is]=r_att_int(AOG, "Z");
		wfn_r[is]=new double[wfn_length[is]];
		r_att_1d(AOG, "r", wfn_length[is], &wfn_r[is][0]);
		// generate wfn_r_reduced
		if(PA_interpolate_wfn){
			wfn_length_reduced[is]=floor(wfn_r[is][wfn_length[is]-1]/wfn_interpolate_dz)+1;
			wfn_r_reduced[is]=new double[wfn_length_reduced[is]];
			for(int ir=0; ir<wfn_length_reduced[is]; ir++){
				wfn_r_reduced[is][ir]=(ir+0.5)*wfn_interpolate_dz;
			}
		}

		wfn_length_use[is]=PA_interpolate_wfn?wfn_length_reduced[is]:wfn_length[is];
		wfn_r_use[is]=new double[wfn_length_use[is]];
		for(j=0; j<wfn_length_use[is]; j++){
			wfn_r_use[is][j]=PA_interpolate_wfn?wfn_r_reduced[is][j]:wfn_r[is][j];
		}
		// load wfn_phi
		for(io=0; io<num_orbits[is]; io++){
			s_data_2d(AOG, orbital_list[is][io], &size1, &size2);
			if(size2!=wfn_length[is]){
				cout << "Size error" << endl;
				return;
			}
			double wfn_both[2][wfn_length[is]];
			r_data_2d(AOG, orbital_list[is][io], 2, wfn_length[is], (double**) wfn_both);
			if(PA_interpolate_wfn){
				wfn_phi_PAO[is][io]=interpolate_wfn(wfn_length[is], &wfn_both[0][0], &wfn_r[is][0], wfn_length_reduced[is], wfn_interpolate_dz);
				wfn_phi_AO[is][io]=interpolate_wfn(wfn_length[is], &wfn_both[1][0], &wfn_r[is][0], wfn_length_reduced[is], wfn_interpolate_dz);
			}else{
				wfn_phi_PAO[is][io]=new double[wfn_length_use[is]];
				wfn_phi_AO[is][io]=new double[wfn_length_use[is]];
				for(j=0; j<wfn_length_use[is]; j++){
					wfn_phi_PAO[is][io][j]=wfn_both[0][j];
					wfn_phi_AO[is][io][j]=wfn_both[1][j];
				}
			}
			wfn_phi_rdr[is][io]=new double[wfn_length_use[is]];
			wfn_phi_dr[is][io]=new double[wfn_length_use[is]];
			if(PA_orth_correction){
				wfn_phi[is][io]=new double[wfn_length_use[is]];
			}
			for(j=0; j<wfn_length_use[is]-1; j++){
				if(empty_atom==true || strcmp(PA_initial_state, "PAO")==0){
					wfn_phi_rdr[is][io][j]=wfn_phi_PAO[is][io][j]*wfn_r_use[is][j]*(wfn_r_use[is][j+1]-wfn_r_use[is][j]);
					wfn_phi_dr[is][io][j]=wfn_phi_PAO[is][io][j]*(wfn_r_use[is][j+1]-wfn_r_use[is][j]);
					if(PA_orth_correction){
						wfn_phi[is][io][j]=wfn_phi_PAO[is][io][j];
					}
				}else{
					wfn_phi_rdr[is][io][j]=wfn_phi_AO[is][io][j]*wfn_r_use[is][j]*(wfn_r_use[is][j+1]-wfn_r_use[is][j]);
					wfn_phi_dr[is][io][j]=wfn_phi_AO[is][io][j]*(wfn_r_use[is][j+1]-wfn_r_use[is][j]);
					if(PA_orth_correction){
						wfn_phi[is][io][j]=wfn_phi_AO[is][io][j];
					}
				}
			}
			wfn_phi_rdr[is][io][wfn_length_use[is]-1]=0;
			wfn_phi_dr[is][io][wfn_length_use[is]-1]=0;
			/*for(j=0; j<wfn_length[is]; j++){
				printf("%f %f\n", wfn_r[is][j], wfn_phi_rdr[is][io][j]);
				}*/
		} 
	}

	// 3. pad calculation
	write_log((char*)"----PAD Calculation----");
  start=chrono::system_clock::now();
	int num_points_E=ceil((PA_E_max-PA_E_min)/PA_E_pixel+1);
	sprintf(sprintf_buffer, "%32s = %d", "Number of points along E", num_points_E);
	write_log(sprintf_buffer);

	/// 3-1. calculate atom position and k point in au (or au^-1) orthonormal basis
	double rec_cell[3][3];
	double op[3]; // for outer product
	double det; // for determinant
	/// reciprocal unit cell
	outer_product(band_cell[1], band_cell[2], op);
	det=inner_product(band_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[0][i]=2.0*M_PI*op[i]/det;
	}

	outer_product(band_cell[2], band_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[1][i]=2.0*M_PI*op[i]/det;
	}

	outer_product(band_cell[0], band_cell[1], op);
	for(i=0; i<3; i++){
		rec_cell[2][i]=2.0*M_PI*op[i]/det;
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
	double kx_vector[3]={0, 0, 0};
	double ky_vector[3]={0, 0, 0};
	double kx_length;
	double ky_length;
	double dkx_length;
	double dky_length;
	double kx_distance[3];
	double ky_distance[3];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			kx_vector[j]+=kx_frac[i]*rec_cell[i][j];
			if(dimension==2){
				ky_vector[j]+=ky_frac[i]*rec_cell[i][j];
			}
		}
	}
	for(i=0; i<3; i++){
		kx_distance[i]=kx_vector[i]*(kx_range[1]-kx_range[0]);
		if(dimension==2){
			ky_distance[i]=ky_vector[i]*(ky_range[1]-ky_range[0]);
		}
	}
	
	kx_length=sqrt(inner_product(kx_vector, kx_vector));
	dkx_length=kx_length*(kx_range[1]-kx_range[0])/(kx_count-1);
	if(dimension==2){
		ky_length=sqrt(inner_product(ky_vector, ky_vector));
		dky_length=ky_length*(ky_range[1]-ky_range[0])/(ky_count-1);
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

	/// weighting
	/// Rect:
	//// if width>0 [origin, origin+width]*1
	//// if width<0 [origin+width, origin]*1
	/// Exp:
	//// if width>0 [origin, origin+width*decay_max]*exp(-(x-origin)/(|width|*2))
	//// if width<0 [origin+width*decay_max, origin]*exp((x-origin)/(|width|*2))
	//// both exp part are the same
	/// Tri:
	//// if width>0 [origin, origin+width]*(x-(origin+width))*(-1/width)
	//// if width<0 [origin-|width|, origin]*(x-(origin-|width|))*(1/|width|)
	/// When Include_neg_depth = true, z satisfying (width>0 && origin>z) || (width<0 && origin<z) becomes 1 in Exp, Tri, and Sqrt

	/// for FPFS, assuming width<0, the bottom edge of the calculation range is range_min-margin
	
	double atom_weighting[atom_length];
	bool atom_weighting_flag[atom_length];
	bool atom_above_surface[atom_length];
	double FPFS_z_bottom=0.0;
	int FPFS_z_start=0;

	// axis_au is also used in reflection correction
	double axis_au[3];
	if(PA_weighting==true){
		double origin_au=PA_weighting_origin;
		double width_au=PA_weighting_width;
		if(PA_use_angstrom){
			origin_au/=au_ang;
			width_au/=au_ang;
			for(i=0; i<3; i++){
				axis_au[i]=PA_weighting_axis[i]/au_ang;
			}
		}else{
			for(i=0; i<3; i++){
				axis_au[i]=PA_weighting_axis[i];
			}
		}
		// normalize axis
		double axis_length=sqrt(inner_product(axis_au, axis_au));
		for(i=0; i<3; i++){
			axis_au[i]/=axis_length;
		}
		double range_min, range_max;
		if(strcmp(PA_weighting_shape, "Rect")==0){
			// Rect
			if(width_au>=0){
				range_min=origin_au;
				range_max=origin_au+width_au;
			}else{
				range_min=origin_au+width_au;
				range_max=origin_au;
			}
		}else if(strcmp(PA_weighting_shape, "Exp")==0){
			// Exp
			if(width_au>=0){
				range_min=origin_au;
				range_max=origin_au+width_au*PA_decay_max;
			}else{
				range_min=origin_au+width_au*PA_decay_max;
				range_max=origin_au;
			}
		}else{
			// Tri or Sqrt
			if(width_au>=0){
				range_min=origin_au;
				range_max=origin_au+width_au;
			}else{
				range_min=origin_au+width_au;
				range_max=origin_au;
			}
		}
		sprintf(sprintf_buffer, "%32s = [%8.2f, %8.2f] [a.u.]", "Weighting range", range_min, range_max);
		write_log(sprintf_buffer);
		for(i=0; i<atom_length; i++){
			double signed_length=inner_product(axis_au, atom_coordinates[i]);
			if(range_min<signed_length && signed_length<range_max){
				atom_weighting_flag[i]=true;
				if(strcmp(PA_weighting_shape, "Rect")==0){
					atom_weighting[i]=1;
				}else if(strcmp(PA_weighting_shape, "Exp")==0){
					atom_weighting[i]=exp(-(signed_length-origin_au)/(width_au*2));
				}else if(strcmp(PA_weighting_shape, "Tri")==0){
					atom_weighting[i]=(signed_length-(origin_au+width_au))*(-1.0/width_au);
				}else{ // Sqrt
					atom_weighting[i]=sqrt((signed_length-(origin_au+width_au))*(-1.0/width_au));
				}
			}else{
				atom_weighting_flag[i]=false;
				atom_weighting[i]=0;
			}
		}
		if(PA_include_neg_depth==true){
			for(i=0; i<atom_length; i++){
				double signed_length=inner_product(axis_au, atom_coordinates[i]);
				if((width_au>0 && signed_length<origin_au) || (width_au<0 && signed_length>origin_au)){
					atom_weighting_flag[i]=true;
					atom_weighting[i]=1.0;
				}
			}
		}
		for(i=0; i<atom_length; i++){
			double signed_length=inner_product(axis_au, atom_coordinates[i]);
			atom_above_surface[i]=(width_au>0 && signed_length<origin_au) || (width_au<0 && signed_length>origin_au);
		}
		write_log((char*)"----Weighting----");
		for(i=0; i<atom_length; i++){
			if(atom_weighting_flag[i]==true){
				sprintf(sprintf_buffer, "%32s (%8.3f, %8.3f, %8.3f) [a.u.]: %10.3e%s", atom_labels[i], atom_coordinates[i][0], atom_coordinates[i][1], atom_coordinates[i][2], atom_weighting[i], atom_above_surface[i]?", above surface":"");
				write_log(sprintf_buffer);
			}else{
				sprintf(sprintf_buffer, "%32s (%8.3f, %8.3f, %8.3f) [a.u.]: out of range%s", atom_labels[i], atom_coordinates[i][0], atom_coordinates[i][1], atom_coordinates[i][2], atom_above_surface[i]?", above surface":"");
				write_log(sprintf_buffer);
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
	/// extend k space (or path)
	double** k_points_ext;
	int* k_index_reduced;
	int total_count_ext=total_count;
	int kx_count_ext=kx_count;
	int ky_count_ext=ky_count;
	if(PA_ext_set){
		if(curved){
			write_log((char*)"Error: Extend cannot be used in the curved dispersion");
			return;
		}
		if(dimension==1){
			// E-k: use ext_ri and ext_le
			total_count_ext=(total_count-1)*(1+PA_ext_ri+PA_ext_le)+1;
			double* k_points_ext_alloc=new double[total_count_ext*3];
			k_points_ext=new double*[total_count_ext];
			k_index_reduced=new int[total_count_ext];
			for(i=-PA_ext_le; i<=PA_ext_ri; i++){
				for(j=0; j<total_count; j++){
					int index_ext=(i+PA_ext_le)*(total_count-1)+j;
					k_index_reduced[index_ext]=j;
					k_points_ext[index_ext]=&k_points_ext_alloc[3*index_ext];
					for(k=0; k<3; k++){
						k_points_ext[index_ext][k]=k_points[j][k]+i*kx_distance[k];
					}
				}
			}/*
				 for(i=0; i<total_count_ext; i++){
				 cout << k_index_reduced[i] << endl;
				 }*/
		}else{
			// E-k-k: use ext_up, ri, dn, and le
			kx_count_ext=((kx_count-1)*(1+PA_ext_ri+PA_ext_le)+1);
			ky_count_ext=((ky_count-1)*(1+PA_ext_up+PA_ext_dn)+1);
			total_count_ext=kx_count_ext*ky_count_ext;
			double* k_points_ext_alloc=new double[total_count_ext*3];
			k_points_ext=new double*[total_count_ext];
			k_index_reduced=new int[total_count_ext];
			int ix, jx, iy, jy;
			for(iy=-PA_ext_dn; iy<=PA_ext_up; iy++){
				for(jy=0; jy<ky_count; jy++){
					for(ix=-PA_ext_le; ix<=PA_ext_ri; ix++){
						for(jx=0; jx<kx_count; jx++){
							int index_original=jy*kx_count+jx;
							int index_ext=((iy+PA_ext_dn)*(ky_count-1)+jy)*kx_count_ext+(ix+PA_ext_le)*(kx_count-1)+jx;
							k_index_reduced[index_ext]=index_original;
							k_points_ext[index_ext]=&k_points_ext_alloc[3*index_ext];
							for(k=0; k<3; k++){
								k_points_ext[index_ext][k]=k_points[index_original][k]+iy*ky_distance[k]+ix*kx_distance[k];
							}
						}
					}
				}
			}
		}
	}

	/// 3-2. calculation of operator coefficients
	complex<double> Y_coeff[3];
	operator_coefficient(PA_polarization, PA_theta, PA_phi, Y_coeff);
	write_log((char*)"----Operator coefficients----");
	for(i=0; i<3; i++){
		sprintf(sprintf_buffer, "%29s[%1d] = %8.3f + %8.3f*i", "Y_coeff", i, Y_coeff[i].real(), Y_coeff[i].imag());
		write_log(sprintf_buffer);
	}
	complex<double> e_vec[3];
	electric_field_vector(PA_polarization, PA_theta, PA_phi, e_vec);
	write_log((char*)"----Electric field vector----");
	for(i=0; i<3; i++){
		sprintf(sprintf_buffer, "%29s[%1d] = %8.3f + %8.3f*i", "e_vec", i, e_vec[i].real(), e_vec[i].imag());
		write_log(sprintf_buffer);
	}
	
	/// 3-3. tail profile
	int tail_index=floor(PA_dE*PA_sigma_max/PA_E_pixel);
	double tail[tail_index+1];
	for(i=0; i<=tail_index; i++){
		tail[i]=1.0/(PA_dE*sqrt(2.0*M_PI))*exp(-1.0/2.0*(i*PA_E_pixel/PA_dE)*(i*PA_E_pixel/PA_dE));
		// printf("%d %f\n", i, tail[i]);
	}


	// 3-4. final states
	/// for Calc
	int k_index_min=-1;
	int k_index_max=-1;
	int final_states_count;
	double**** final_states_calc; // [k_index][atom_spec][l][r] for Calc
	double*** final_states_phase; // [k_index][atom_spec][l] for Calc
	/// for FPFS
	int** final_states_EScale;                // [total_count_ext][FPIndex]=round(E/FPFS_energy_step)
	int** final_states_spin;                  // [total_count_ext][FPIndex]=up(0) or down(1), always zero for spin_i==0 (off)
	double*** final_states_k;                 // [total_count_ext][FPIndex]=k vector (au^-1)
	double** final_states_zgels_norm;         // [total_count_ext][FPIndex]=zgels norm (in FPFS_bulk)
	int* final_states_FP_size;                // [total_count_ext]=FPIndex_size
	int E_min_scale;                          // round(E_min_tail/PA_FPFS_energy_step);
	int E_max_scale;                          // round(E_max_tail/PA_FPFS_energy_step);
	//// for local part
	complex<double>**** final_states_FP_loc;  // [total_count_ext][FPIndex][ig][iz] for FP_PAO and FP_AO
	int** final_states_FP_g_size;             // [total_count_ext][FPIndex]=ig size
	int**** final_states_FP_g;                // [total_count_ext][FPIndex][ig][2]=[n1,n2]
	double**** final_states_FP_g_vec;         // [total_count_ext][FPIndex][ig][3]=g vector
	complex<double>*** final_states_FP_loc_edge; // [total_count_ext][FPIndex][ig]
	//// for nonlocal part
	complex<double>***** final_states_FP_nonloc; // [total_count_ext][FPIndex][ia][l][mpl] @rc
	double***** final_states_FP_nonloc_r;     // [EScale][sp][ia][l][ir] based on wfn_r, satisfying the value at the end is 1.0
	double*** final_states_FP_norm1;          // [total_count_ext][FPIndex][ia]
  double*** final_states_FP_norm2;          // [total_count_ext][FPIndex][ia]
	//// for bulk
	int** final_states_FP_g_bulk;             // [ig'][3]
	double** final_states_FP_g_vec_bulk;      // [ig'][3]
	int final_states_FP_g_size_bulk;          // ig' size
	complex<double>**** final_states_FP_bulk; // [total_count_ext][FPIndex][in][ig']
	//complex<double>***** final_states_FP_bulk_z; // [total_count_ext][FPIndex][in][ig][iz]
	double*** final_states_FP_bulk_kz;        // [total_count_ext][FPIndex][in]
	double*** final_states_FP_bulk_kappaz;    // [total_count_ext][FPIndex][in]
	int** final_states_FP_bulk_count;         // [total_count_ext][FPIndex] = in size
	complex<double>*** final_states_FP_bulk_coefs; // [total_count_ext][FPIndex][in] = linear combination coefficients
	double*** FP_bulk_dispersion_up; // [total_count_ext][ikz][ig']
	double*** FP_bulk_dispersion_dn; // [total_count_ext][ikz][ig']
	//complex<double>*** FP_bulk_dispersion_complex; // [ikz][ikappaz][ig']
	double* FP_bulk_dispersion_kz=new double[PA_FPFS_bulk_kz_steps]; // [ikz]
	double* FP_bulk_dispersion_kappaz; // [ikappaz]
	double* FP_bulk_dispersion_mkappaz;
	int FP_bulk_kappaz_count=0;
	int FP_bulk_kappaz_border=0;
	complex<double>*** FP_bulk_dispersion_c_up; // [total_count_ext][ikappaz, ib]
	complex<double>*** FP_bulk_dispersion_c_dn;
	complex<double>*** FP_bulk_dispersion_c_BZ_up;
	complex<double>*** FP_bulk_dispersion_c_BZ_dn;
	complex<double>*** FP_bulk_dispersion_mc_up;
	complex<double>*** FP_bulk_dispersion_mc_dn;
	complex<double>*** FP_bulk_dispersion_mc_BZ_up;
	complex<double>*** FP_bulk_dispersion_mc_BZ_dn;
	int** FP_bulk_dispersion_c_count_up; // [total_count_ext][ikappaz]=ib size
	int** FP_bulk_dispersion_c_count_dn;
	int** FP_bulk_dispersion_c_BZ_count_up;
	int** FP_bulk_dispersion_c_BZ_count_dn;
	int** FP_bulk_dispersion_mc_count_up;
	int** FP_bulk_dispersion_mc_count_dn;
	int** FP_bulk_dispersion_mc_BZ_count_up;
	int** FP_bulk_dispersion_mc_BZ_count_dn;
	int*** FP_bulk_connection_c_up; // [total_count_ext][ikappaz, ib]=connected ib at ikappaz+1
	int*** FP_bulk_connection_c_dn;
	int*** FP_bulk_connection_c_BZ_up;
	int*** FP_bulk_connection_c_BZ_dn;
	int*** FP_bulk_connection_mc_up;
	int*** FP_bulk_connection_mc_dn;
	int*** FP_bulk_connection_mc_BZ_up;
	int*** FP_bulk_connection_mc_BZ_dn;
	
	//// VPS-related variables
	int VPS_l_length[atom_spec_length];       // [is]=il value
	int VPS_r_length[atom_spec_length];       // [is]=ir value
	int VPS_j_length[atom_spec_length];       // [is]=2 if j_dependent else 1
	int* VPS_l[atom_spec_length];             // [is][il] = orbital angular momenta
	double* VPS_r[atom_spec_length];          // [is][ir] = r value
	double** VPS_E[atom_spec_length];         // [is][il][ij] = projector energies
	double* VPS_loc[atom_spec_length];        // [is][ir] = local potentials (-> not used in calculations)
	double*** VPS_nonloc[atom_spec_length];   // [is][il][ij][ir] = nonlocal potentials
	double VPS_cutoff[atom_spec_length];      // [is]=max of radial cutoff
	int vps_cutoff_index[atom_spec_length];   // [is]=radial cutoff index based on VPS_r
	int wfn_cutoff_index[atom_spec_length];   // [is]=radial cutoff index based on wfn_r
	double* VPS_E_ave[atom_spec_length];      // [is][il] = projector energies averaged w.r.t. ij
	double** VPS_nonloc_ave[atom_spec_length];// [is][il][ir] = averaged nonlocal potentials, ir corresponds to VPS_r
	double*** VKS0;                           // [ix][iy][iz] = Kohn-Sham potential for up spin (wo/ nonlocal)
	double*** VKS1;                           // [ix][iy][iz] = Kohn-Sham potential for dn spin (wo/ nonlocal)
	double* VKS0_r[atom_length];              // [ia][ir] (average) = radial Kohn-Sham potential for up spin ir ~ VPS_r
	double* dVKS0_r[atom_length];             // [ia][ir] (stdev)
	double* VKS1_r[atom_length];              // [ia][ir] (average) = radial Kohn-Sham potential for dn spin
	double* dVKS1_r[atom_length];             // [ia][ir] (stdev)
	int VKS_count[3];                         // = Size of VKS0 and VKS1
	// Fourier components
	complex<double>** Vgg0;                   // [ig][iz] = Fourier-expanded potential for up spin
	complex<double>** Vgg1;                   // [ig][iz] = Fourier-expanded potential for dn spin
	double** Vgg0_abs;                        // [ig][iz] = |Vgg0[ig][iz]|
	double** Vgg1_abs;
	int** Vgg_list;                           // [ig][2] = linear combination coefficients for g vector
	double** Vgg_vector;                      // [ig][3] = g
	int Vgg_count=0;                          // ig size
	double final_states_dz;                   // = atom_cell[0] length / VKS_count[0]
	/// for bulk calculations
	double FPFS_bulk_min=PA_FPFS_bulk_min_ang/au_ang;
	double FPFS_bulk_max=PA_FPFS_bulk_max_ang/au_ang;
	int FPFS_bulk_count=PA_FPFS_bulk_count;
	double FPFS_bulk_height=(FPFS_bulk_max-FPFS_bulk_min)/(FPFS_bulk_count*1.0);
	double FPFS_gz_length=2.0*M_PI/FPFS_bulk_height;
	complex<double>* Vgg0_bulk;
	double* Vgg0_abs_bulk;
	complex<double>* Vgg1_bulk;
	double* Vgg1_abs_bulk;
	int** Vgg_list_bulk;
	double** Vgg_vector_bulk;
	int Vgg_count_bulk=0;
	int z_count_bulk;
	double dz_bulk;
	complex<double>** Vgg0_average; // [ig][iz']
	double** Vgg0_average_abs;
	double** Vgg0_stdev;
	complex<double>** Vgg1_average;
  double** Vgg1_stdev;
	double** Vgg1_average_abs;
	
	
	// Determining the g range for VPS
	// k=sqrt(2*hn)
	// g_max=2*k*FPFS_kRange, 2 is necessary because g=g1-g2
	// the search area is (n1, n2) \in k/min(b1, b2)*FPFS_kRange*FPFS_range
	
	// Determing the g range for wave function
	// |k+G|/2 < |(k//, kz)|
	// the search area is (n1, n2) \in k/min(b1, b2)*FPFS_kRange*FPFS_range

	/// for orthogonality correction
	double** Self_radial_int[atom_spec_length]; // [is][io1][io2] = \int r P[io1]P[io2] dr
	
	if(strcmp(PA_final_state, "Calc")==0){
		write_log((char*)"----Calculate final states----");
		for(i=0; i<total_count_ext; i++){
			double* k_point;
			if(PA_ext_set){
				k_point=k_points_ext[i];
			}else{
				k_point=k_points[i];
			}
			double k_length=sqrt(inner_product(k_point, k_point));
			int k_index=round(k_length/PA_final_state_step);
			if(k_index_min<0 || k_index_min>k_index){
				k_index_min=k_index;
			}
			if(k_index_max<0 || k_index_max<k_index){
				k_index_max=k_index;
			}
		}
		sprintf(sprintf_buffer, "%32s = [%d, %d] * %8.3f", "Rounded k_length range", k_index_min, k_index_max, PA_final_state_step);
		write_log(sprintf_buffer);
		final_states_count=k_index_max-k_index_min+1;

		final_states_calc=new double***[final_states_count];
		for(i=0; i<final_states_count; i++){
			final_states_calc[i]=new double**[atom_spec_length];
			for(j=0; j<atom_spec_length; j++){
				double* final_states_alloc=new double[5*wfn_length_use[j]];
				final_states_calc[i][j]=new double*[5];
				for(k=0; k<5; k++){
					final_states_calc[i][j][k]=&final_states_alloc[k*wfn_length_use[j]];
					// cout << final_states_calc[i][j][k] << endl;
				}
			}
		}
		
		double* final_states_phase_alloc=new double[final_states_count*atom_spec_length*5];
		final_states_phase=new double**[final_states_count];
		for(i=0; i<final_states_count; i++){
			final_states_phase[i]=new double*[atom_spec_length];
			for(j=0; j<atom_spec_length; j++){
				final_states_phase[i][j]=&final_states_phase_alloc[i*atom_spec_length*5+j*5];
			}
		}

		H5File database(At_potential_file, H5F_ACC_RDONLY);
		char group_name[4];
		
		for(j=0; j<atom_spec_length; j++){
			double mu=pow(3.0*M_PI/4.0, 2.0/3.0)/(2.0*pow(Z[j], 1.0/3.0));
			// cout << x_count << " " << x_coordinates[x_count-1] << endl;
			sprintf(sprintf_buffer, "Z = %3d, mu = %8.3f", Z[j], mu);
			write_log(sprintf_buffer);
			// open group
			sprintf(group_name, "%03d", Z[j]);
			Group atomG(database.openGroup(group_name));
			// prepare potential by interpolation
			double mod_start_x=load_potential_H5(atomG, mu);
			
			for(i=0; i<final_states_count; i++){
				double k_length=(k_index_min+i)*PA_final_state_step;
				double Ekin=k_length*k_length*0.5;
				int l;
				for(l=0; l<5; l++){
					Atomic_wfn_evolution(mu, l, Ekin, false);
					final_states_phase[i][j][l]=analyze_phase(mu, k_length, l, mod_start_x);
					// interpolate final states so that they correspond to wfn_r_use[j]
					int x1, x2;
					for(x1=0; x1<wfn_length_use[j]; x1++){
						bool r_found=false;
						for(x2=0; x2<x_count-1; x2++){
							if(mu*x_coordinates[x2] < wfn_r_use[j][x1] && wfn_r_use[j][x1] < mu*x_coordinates[x2+1]){
								double dr1=wfn_r_use[j][x1]-mu*x_coordinates[x2];
								double dr2=mu*x_coordinates[x2+1]-wfn_r_use[j][x1];
								final_states_calc[i][j][l][x1]=(At_p_x[x2]*dr2+At_p_x[x2+1]*dr1)/(dr1+dr2);
								r_found=true;
								break;
							}
						}
						if(!r_found){
							final_states_calc[i][j][l][x1]=0;
							// cout << "Not found" << endl;
						}
					}
				}
			}
		}

		/*		
					for(k=0; k<5; k++){
					for(i=0; i<final_states_count; i++){
					for(j=0; j<atom_spec_length; j++){
					cout << final_states_phase[i][j][k] << " ";
					}		
					}
					cout << endl;
					}return;
					for(i=0; i<wfn_length[0]; i++){
					cout << wfn_r[0][i] << " " << final_states_calc[0][0][0][i] << endl;
					}*/
				
		// return;
	}else if(PA_FPFS){
		write_log((char*)"----Calculate First-principles final states----");
		if(spin_i==2){
			write_log((char*)"Error: this software can not perform noncollinear FPFS calculations.");
			return;
		}
		// load pseudopotentials
		write_log((char*)"----Load Pseudopotentials----");
		H5File VPS_db(PA_VPS_file, H5F_ACC_RDONLY);
		for(int is=0; is<atom_spec_length; is++){
			Group VPSG(VPS_db.openGroup(atom_spec_PS[is]));
			sprintf(sprintf_buffer, "Group: %s", atom_spec_PS[is]);
			write_log(sprintf_buffer);
			VPS_l_length[is]=r_att_int(VPSG, "num_vps_proj");
			VPS_r_length[is]=r_att_int(VPSG, "num_grid");
			VPS_j_length[is]=(r_att_bool(VPSG, "j_dependent")) ? 2:1;
			VPS_cutoff[is]=r_att_double(VPSG, "max_cutoff");
			VPS_l[is]=new int[VPS_l_length[is]];
			VPS_r[is]=new double[VPS_r_length[is]];
			VPS_loc[is]=new double[VPS_r_length[is]];
			double* buffer=new double[VPS_l_length[is]*VPS_j_length[is]]; // for VPS_E
			VPS_E[is]=new double*[VPS_l_length[is]];
			VPS_E_ave[is]=new double[VPS_l_length[is]];
			double* buffer2=new double[VPS_l_length[is]*VPS_j_length[is]*VPS_r_length[is]]; // for VPS_nonloc
			VPS_nonloc[is]=new double**[VPS_l_length[is]];
			VPS_nonloc_ave[is]=new double*[VPS_l_length[is]];
			for(int il=0; il<VPS_l_length[is]; il++){
				VPS_E[is][il]=&(buffer[il*VPS_j_length[is]]);
				VPS_nonloc[is][il]=new double*[VPS_j_length[is]];
				for(int ij=0; ij<VPS_j_length[is]; ij++){
					VPS_nonloc[is][il][ij]=&(buffer2[il*VPS_j_length[is]*VPS_r_length[is]+ij*VPS_r_length[is]]);
				}
			}
			// printf("l %d j %d r %d\n", VPS_l_length[is], VPS_j_length[is], VPS_r_length[is]);
			r_att_1d(VPSG, "r", VPS_r_length[is], &VPS_r[is][0]);
			r_data_1i(VPSG, "orbital_angular_momenta", VPS_l_length[is], &VPS_l[is][0]);
			r_data_2d(VPSG, "project_energies", VPS_l_length[is], VPS_j_length[is], (double**)&VPS_E[is][0][0]);
			r_data_1d(VPSG, "local_potential", VPS_r_length[is], &VPS_loc[is][0]);
			if(!empty_atoms[is]){
				r_data_3d(VPSG, "nonlocal_potentials", VPS_l_length[is], VPS_j_length[is], VPS_r_length[is], (double***)&VPS_nonloc[is][0][0][0]);
			}
			/*for(int ir=0; ir<VPS_r_length[is]; ir++){
				printf("%7.4f %8.4f\n", VPS_r[is][ir], VPS_loc[is][ir]);
				}*/
			if(empty_atoms[is]){
				vps_cutoff_index[is]=0;
				wfn_cutoff_index[is]=0;
			}else{
				for(int ir=0; ir<VPS_r_length[is]; ir++){
					if(VPS_r[is][ir]>VPS_cutoff[is]){
						vps_cutoff_index[is]=ir;
						break;
					}
				}
				for(int ir=0; ir<wfn_length_use[is]; ir++){
					if(wfn_r_use[is][ir]>VPS_cutoff[is]){
						wfn_cutoff_index[is]=ir;
						break;
					}
				}
			}
			sprintf(sprintf_buffer, "VPS cutoff #%-4d = %.3f Bohr ~ %.3f Bohr", vps_cutoff_index[is], VPS_r[is][vps_cutoff_index[is]], VPS_cutoff[is]);
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "Wfn cutoff #%-4d = %.3f Bohr ~ %.3f Bohr", wfn_cutoff_index[is], wfn_r_use[is][wfn_cutoff_index[is]], VPS_cutoff[is]);
			write_log(sprintf_buffer);
			if(empty_atoms[is]){
				continue;
			}
			for(int il=0; il<VPS_l_length[is]; il++){
				VPS_nonloc_ave[is][il]=new double[vps_cutoff_index[is]];
			}
			// averaging
			for(int il=0; il<VPS_l_length[is]; il++){
				if(VPS_j_length[is]==1){
					// l-dependent
					VPS_E_ave[is][il]=VPS_E[is][il][0];
					for(int ir=0; ir<vps_cutoff_index[is]; ir++){
						VPS_nonloc_ave[is][il][ir]=VPS_nonloc[is][il][0][ir];
					}
				}else{
					// j-dependent
					double ld=VPS_l[is][il]*1.0;
					VPS_E_ave[is][il]=((ld+1.0)*VPS_E[is][il][0]+ld*VPS_E[is][il][1])/(2.0*ld+1.0);
					// printf("%6.2f %6.2f %6.2f\n", VPS_E[is][il][0], VPS_E[is][il][1], VPS_E_ave[is][il]);
					for(int ir=0; ir<vps_cutoff_index[is]; ir++){
						VPS_nonloc_ave[is][il][ir]=((ld+1.0)*VPS_nonloc[is][il][0][ir]+ld*VPS_nonloc[is][il][1][ir])/(2.0*ld+1.0);
					}
				}
			}
		}
		write_log((char*)"----Load Kohn-Sham potentials----");
		// load
		// the ignore_VKS flag is for debug
		bool ignore_VKS=false;
		
		Group PotentialG(OutputG.openGroup("Potential"));
		r_att_1i(PotentialG, "Count", 3, &VKS_count[0]);
		VKS0=alloc_dcube(VKS_count[0], VKS_count[1], VKS_count[2]);
		if(!ignore_VKS){
		  r_data_3d(PotentialG, "V0", VKS_count[0], VKS_count[1], VKS_count[2], (double***)&VKS0[0][0][0]);
		}
		if(spin_i==1 || spin_i==2){
			VKS1=alloc_dcube(VKS_count[0], VKS_count[1], VKS_count[2]);
			if(!ignore_VKS){
			  r_data_3d(PotentialG, "V1", VKS_count[0], VKS_count[1], VKS_count[2], (double***)&VKS1[0][0][0]);
			}
		}
		// caculate average for the nonlocal areas
		double Lebedev_r[3][PA_Lebedev_order_ave];
		double Lebedev_w[PA_Lebedev_order_ave];
		ld_by_order(PA_Lebedev_order_ave, Lebedev_r[0], Lebedev_r[1], Lebedev_r[2], Lebedev_w);
		/*
			for(int i=0; i<PA_Lebedev_order_ave; i++){
			printf("(%7.4f, %7.4f, %7.4f) w=%7.4f, r=%7.4f\n", Lebedev_r[0][i], Lebedev_r[1][i], Lebedev_r[2][i], Lebedev_w[i], sqrt(Lebedev_r[0][i]*Lebedev_r[0][i]+Lebedev_r[1][i]*Lebedev_r[1][i]+Lebedev_r[2][i]*Lebedev_r[2][i]));
			}*/
		double r_ref[3];
		for(int ia=0; ia<atom_length; ia++){
			int is=atom_spec_index[ia];
			VKS0_r[ia]=new double[VPS_r_length[is]];
			dVKS0_r[ia]=new double[VPS_r_length[is]];
			if(spin_i>0){
				VKS1_r[ia]=new double[VPS_r_length[is]];
				dVKS1_r[ia]=new double[VPS_r_length[is]];
			}
			for(int ir=0; ir<VPS_r_length[is]; ir++){
				VKS0_r[ia][ir]=0.0;
				dVKS0_r[ia][ir]=0.0;
				if(spin_i>0){
					VKS1_r[ia][ir]=0.0;
					dVKS1_r[ia][ir]=0.0;
				}
				double ri=VPS_r[is][ir];
				for(int iL=0; iL<PA_Lebedev_order_ave; iL++){
					// set the coordinate
					for(int ix=0; ix<3; ix++){
						r_ref[ix]=atom_coordinates[ia][ix]+Lebedev_r[ix][iL]*ri;
					}
					double VKS_ref=interpolate_potential(r_ref, VKS_count, &VKS0[0][0][0], &atom_cell[0][0]);
					VKS0_r[ia][ir]+=VKS_ref;
					dVKS0_r[ia][ir]+=VKS_ref*VKS_ref;
					if(spin_i>0){
						VKS_ref=interpolate_potential(r_ref, VKS_count, &VKS1[0][0][0], &atom_cell[0][0]);
						VKS1_r[ia][ir]+=VKS_ref;
						dVKS1_r[ia][ir]+=VKS_ref*VKS_ref;
					}
				}
				VKS0_r[ia][ir]/=PA_Lebedev_order_ave;
				dVKS0_r[ia][ir]/=PA_Lebedev_order_ave;
				dVKS0_r[ia][ir]-=VKS0_r[ia][ir]*VKS0_r[ia][ir];
				if(spin_i>0){
					VKS1_r[ia][ir]/=PA_Lebedev_order_ave;
					dVKS1_r[ia][ir]/=PA_Lebedev_order_ave;
					dVKS1_r[ia][ir]-=VKS1_r[ia][ir]*VKS1_r[ia][ir];
				}
				/*
					if(is==0){
					printf("%10.6f %10.6f %10.6f\n", VPS_r[is][ir], VKS_r[ia][ir], dVKS_r[ia][ir]);
					}*/
			}
			/*
				if(is==0){
				printf("\n");
				}*/
		}

		// Fourier expansion of VKS
		write_log((char*)"----Fourier expansion of the Kohn-Sham potential----");
		double khn_approx=sqrt(2.0*PA_excitation_energy/Eh);
		double b1_length=sqrt(inner_product(rec_cell[1], rec_cell[1]));
		double b2_length=sqrt(inner_product(rec_cell[2], rec_cell[2]));
		double b_length=min(b1_length, b2_length);
		double gMax=2*khn_approx*PA_FPFS_kRange;
		int n_range=ceil(2*khn_approx/b_length*PA_FPFS_range*PA_FPFS_kRange);
		int n_range_z;
		sprintf(sprintf_buffer, "khn = %.3f", khn_approx);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "b   = %.3f", b_length);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "gMax= %.3f", gMax);
		write_log(sprintf_buffer);
		sprintf(sprintf_buffer, "n   = %d", n_range);
		write_log(sprintf_buffer);
		bool g_included[n_range*2+1][n_range*2+1];
		double g12[3];
		double g12_length;
		for(int n1=-n_range; n1<=n_range; n1++){
			int n1p=n1+n_range;
			for(int n2=-n_range; n2<=n_range; n2++){
				int n2p=n2+n_range;
				for(int ix=0; ix<3; ix++){
					g12[ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
				}
				g12_length=sqrt(inner_product(g12, g12));
				if(g12_length<gMax){
					g_included[n1p][n2p]=true;
					Vgg_count++;
				}else{
					g_included[n1p][n2p]=false;
				}
			}
		}
		sprintf(sprintf_buffer, "Vgg size = %d", Vgg_count);
		write_log(sprintf_buffer);
		Vgg0=alloc_zmatrix(Vgg_count, VKS_count[0]);
		Vgg0_abs=alloc_dmatrix(Vgg_count, VKS_count[0]);
		if(spin_i==1 || spin_i==2){
			Vgg1=alloc_zmatrix(Vgg_count, VKS_count[0]);
			Vgg1_abs=alloc_dmatrix(Vgg_count, VKS_count[0]);
		}
		Vgg_list=alloc_imatrix(Vgg_count, 2);
		Vgg_vector=alloc_dmatrix(Vgg_count, 3);
		
		int ig_count=0;
		for(int n1=-n_range; n1<=n_range; n1++){
			int n1p=n1+n_range;
			for(int n2=-n_range; n2<=n_range; n2++){
				int n2p=n2+n_range;
				if(g_included[n1p][n2p]==true){
					Vgg_list[ig_count][0]=n1;
					Vgg_list[ig_count][1]=n2;
					for(int ix=0; ix<3; ix++){
						Vgg_vector[ig_count][ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
					}
					//printf("%8.3f %8.3f %8.3f\n", Vgg_vector[ig_count][0], Vgg_vector[ig_count][1], Vgg_vector[ig_count][2]);
					ig_count++;
				}
			}
		}

#pragma omp parallel firstprivate(atom_cell)
#pragma omp for
		for(int ig=0; ig<Vgg_count; ig++){
		  Fourier_expansion_z(&VKS0[0][0][0], VKS_count, Vgg_vector[ig], &atom_cell[0][0], Vgg0[ig], Vgg0_abs[ig]);
			if(spin_i==1 || spin_i==2){
				Fourier_expansion_z(&VKS1[0][0][0], VKS_count, Vgg_vector[ig], &atom_cell[0][0], Vgg1[ig], Vgg1_abs[ig]);
			}
		}
		/*
			for(int ig=0; ig<Vgg_count; ig++){
			printf("(n1, n2)=(%d, %d)\n", Vgg_list[ig][0], Vgg_list[ig][1]);
			for(int iz=0; iz<VKS_count[0]; iz++){
			printf("Vg = (%f, %f)\n", Vgg0[ig][iz].real(), Vgg0[ig][iz].imag());
			}
			}*/
		final_states_dz=sqrt(inner_product(atom_cell[0], atom_cell[0]))/VKS_count[0];
		sprintf(sprintf_buffer, "dz = %.3f", final_states_dz);
		write_log(sprintf_buffer);
		if(PA_FPFS_bulk_set){
			FPFS_z_start=floor(FPFS_bulk_max/final_states_dz);
			sprintf(sprintf_buffer, "z_start = %d", FPFS_z_start);
			write_log(sprintf_buffer);
		}

		if(PA_FPFS_bulk_set){
			write_log((char*)"----Kohn-Sham potential calculations for bulk states----");
			z_count_bulk=round(FPFS_bulk_height/final_states_dz);
			dz_bulk=FPFS_bulk_height/(z_count_bulk*1.0);
			sprintf(sprintf_buffer, "z_count for bulk = %d", z_count_bulk);
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "dz for bulk = %.3f", dz_bulk);
			write_log(sprintf_buffer);
			sprintf(sprintf_buffer, "2pi/c for bulk = %.3f", FPFS_gz_length);
			write_log(sprintf_buffer);

			Vgg0_average=alloc_zmatrix(Vgg_count, z_count_bulk);
			Vgg0_stdev=alloc_dmatrix(Vgg_count, z_count_bulk);
			Vgg0_average_abs=alloc_dmatrix(Vgg_count, z_count_bulk);
			if(spin_i==1 || spin_i==2){
				Vgg1_average=alloc_zmatrix(Vgg_count, z_count_bulk);
				Vgg1_stdev=alloc_dmatrix(Vgg_count, z_count_bulk);
				Vgg1_average_abs=alloc_dmatrix(Vgg_count, z_count_bulk);
			}

			for(int ig=0; ig<Vgg_count; ig++){
				calc_bulk_potential(Vgg0[ig], final_states_dz, VKS_count[0], FPFS_bulk_min, FPFS_bulk_height, dz_bulk, z_count_bulk, FPFS_bulk_count, Vgg0_average[ig], Vgg0_stdev[ig]);
				/*
				for(int iz=0; iz<z_count_bulk; iz++){
					printf("(%8.4f, %8.4f) +- %8.4f\n", Vgg0_average[ig][iz].real(), Vgg0_average[ig][iz].imag(), Vgg0_stdev[ig][iz]);
				}
				printf("\n");*/
				if(spin_i==1 || spin_i==2){
					calc_bulk_potential(Vgg1[ig], final_states_dz, VKS_count[0], FPFS_bulk_min, FPFS_bulk_height, dz_bulk, z_count_bulk, FPFS_bulk_count, Vgg1_average[ig], Vgg1_stdev[ig]);
				}
				for(int iz=0; iz<z_count_bulk; iz++){
					Vgg0_average_abs[ig][iz]=abs(Vgg0_average[ig][iz]);
					if(spin_i==1 || spin_i==2){
						Vgg1_average_abs[ig][iz]=abs(Vgg1_average[ig][iz]);
					}
				}
			}

			// Fourier expansion along z
			n_range_z=ceil(2*khn_approx/FPFS_gz_length*PA_FPFS_kRange);
			double g123[3];
			double g123_length;
			
			for(int n1=-n_range; n1<=n_range; n1++){
				int n1p=n1+n_range;
				for(int n2=-n_range; n2<=n_range; n2++){
					int n2p=n2+n_range;
					for(int n3=-n_range_z; n3<=n_range_z; n3++){
						for(int ix=0; ix<3; ix++){
							g123[ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
						}
						g123[2]+=n3*FPFS_gz_length;
						g123_length=sqrt(inner_product(g123, g123));
						if(g123_length<gMax){
							Vgg_count_bulk++;
						}
					}
				}
			}
			sprintf(sprintf_buffer, "Vgg_bulk size = %d", Vgg_count_bulk);
			write_log(sprintf_buffer);
			Vgg0_bulk=new complex<double>[Vgg_count_bulk];
			Vgg0_abs_bulk=new double[Vgg_count_bulk];
			if(spin_i==1 || spin_i==2){
				Vgg1_bulk=new complex<double>[Vgg_count_bulk];
				Vgg1_abs_bulk=new double[Vgg_count_bulk];
			}
			Vgg_list_bulk=alloc_imatrix(Vgg_count_bulk, 3);
			Vgg_vector_bulk=alloc_dmatrix(Vgg_count_bulk, 3);

			ig_count=0;
			for(int n1=-n_range; n1<=n_range; n1++){
				for(int n2=-n_range; n2<=n_range; n2++){
					for(int n3=-n_range_z; n3<=n_range_z; n3++){
						for(int ix=0; ix<3; ix++){
							g123[ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
						}
						g123[2]+=n3*FPFS_gz_length;
						g123_length=sqrt(inner_product(g123, g123));
						if(g123_length<gMax){
							Vgg_list_bulk[ig_count][0]=n1;
							Vgg_list_bulk[ig_count][1]=n2;
							Vgg_list_bulk[ig_count][2]=n3;
							for(int ix=0; ix<3; ix++){
								Vgg_vector_bulk[ig_count][ix]=g123[ix];
							}
							// printf("%8.3f %8.3f %8.3f\n", Vgg_vector_bulk[ig_count][0], Vgg_vector_bulk[ig_count][1], Vgg_vector_bulk[ig_count][2]);
							ig_count++;
						}
					}
				}
			}

			for(int igb=0; igb<Vgg_count_bulk; igb++){
				int ig_found=-1;
				for(int ig=0; ig<Vgg_count; ig++){
					if(Vgg_list_bulk[igb][0]==Vgg_list[ig][0] && Vgg_list_bulk[igb][1]==Vgg_list[ig][1]){
						ig_found=ig;
						break;
					}
				}
				if(ig_found<0){
					write_log((char*)"Error: Vgg not found");
					return;
				}
				Vgg0_bulk[igb]=Fourier_expansion_1D(Vgg0_average[ig_found], Vgg_vector_bulk[igb][2], dz_bulk, z_count_bulk);
				Vgg0_abs_bulk[igb]=abs(Vgg0_bulk[igb]);
				// printf("%8.4f %8.4f %8.4f %8.4f\n", Vgg_vector_bulk[igb][0], Vgg_vector_bulk[igb][1], Vgg_vector_bulk[igb][2], Vgg0_abs_bulk[igb]);
				if(spin_i==1 || spin_i==2){
					Vgg1_bulk[igb]=Fourier_expansion_1D(Vgg1_average[ig_found], Vgg_vector_bulk[igb][2], dz_bulk, z_count_bulk);
					Vgg1_abs_bulk[igb]=abs(Vgg1_bulk[igb]);
				}
			}
			/*
			for(int igb=0; igb<Vgg_count_bulk; igb++){
				int igb_found=-1;
				for(int igb2=0; igb2<Vgg_count_bulk; igb2++){
					if(Vgg_list_bulk[igb][0]==-Vgg_list_bulk[igb2][0] &&
						 Vgg_list_bulk[igb][1]==-Vgg_list_bulk[igb2][1] &&
						 Vgg_list_bulk[igb][2]==-Vgg_list_bulk[igb2][2]){
						igb_found=igb2;
						break;
					}
				}
				complex<double> Vg1g2=Vgg0_bulk[igb];
				complex<double> Vg2g1=Vgg0_bulk[igb_found];
				double re_ave=(Vg1g2.real()+Vg2g1.real())*0.5;
				double im_ave=(Vg1g2.imag()-Vg2g1.imag())*0.5;
				Vgg0_bulk[igb]=complex<double>(re_ave, im_ave);
				Vgg0_bulk[igb_found]=complex<double>(re_ave, -im_ave);
				}*/

		} // end of if(PA_FPFS_bulk_set)
		
		write_log((char*)"----Energy scale calculations----");
		sprintf(sprintf_buffer, "Fermi level = %.2f Eh = %.3f eV", EF_Eh, EF_Eh*Eh);
		write_log(sprintf_buffer);
		double E_min_tail=PA_E_min-PA_E_pixel*(tail_index+0.5);
		double E_max_tail=PA_E_max+PA_E_pixel*(tail_index+0.5);
		E_min_scale=round(E_min_tail/PA_FPFS_energy_step);
		E_max_scale=round(E_max_tail/PA_FPFS_energy_step);
		// printf("Scale: [%d, %d]\n", E_min_scale, E_max_scale);
		int scale_width=E_max_scale-E_min_scale+1;
		final_states_EScale=new int*[total_count_ext];
		final_states_spin=new int*[total_count_ext];
		final_states_FP_size=new int[total_count_ext];
		final_states_k=new double**[total_count_ext];
		
		final_states_FP_loc=new complex<double>***[total_count_ext];
		final_states_FP_g_size=new int*[total_count_ext];
		final_states_FP_g=new int***[total_count_ext];
		final_states_FP_g_vec=new double***[total_count_ext];
		if(!PA_FPFS_bulk_set){
			final_states_FP_loc_edge=new complex<double>**[total_count_ext];
		}

		final_states_FP_nonloc=new complex<double>****[total_count_ext];
		final_states_FP_norm1=new double**[total_count_ext];
		final_states_FP_norm2=new double**[total_count_ext];

		if(PA_FPFS_bulk_set){
			final_states_FP_bulk=new complex<double>***[total_count_ext];
			final_states_FP_bulk_count=new int*[total_count_ext];
			//final_states_FP_bulk_z=new complex<double>****[total_count_ext];
			final_states_FP_bulk_kz=new double**[total_count_ext];
			final_states_FP_bulk_kappaz=new double**[total_count_ext];
			final_states_FP_bulk_coefs=new complex<double>**[total_count_ext];
			if(PA_FPFS_file_set){
				FP_bulk_dispersion_c_up=new complex<double>**[total_count_ext];
				FP_bulk_dispersion_c_BZ_up=new complex<double>**[total_count_ext];
				FP_bulk_dispersion_mc_up=new complex<double>**[total_count_ext];
				FP_bulk_dispersion_mc_BZ_up=new complex<double>**[total_count_ext];
				FP_bulk_dispersion_c_count_up=new int*[total_count_ext];
				FP_bulk_dispersion_c_BZ_count_up=new int*[total_count_ext];
				FP_bulk_dispersion_mc_count_up=new int*[total_count_ext];
				FP_bulk_dispersion_mc_BZ_count_up=new int*[total_count_ext];
				FP_bulk_connection_c_up=new int**[total_count_ext];
				FP_bulk_connection_c_BZ_up=new int**[total_count_ext];
				FP_bulk_connection_mc_up=new int**[total_count_ext];
				FP_bulk_connection_mc_BZ_up=new int**[total_count_ext];
				if(spin_i>0){
					FP_bulk_dispersion_c_dn=new complex<double>**[total_count_ext];
					FP_bulk_dispersion_c_BZ_dn=new complex<double>**[total_count_ext];
					FP_bulk_dispersion_mc_dn=new complex<double>**[total_count_ext];
					FP_bulk_dispersion_mc_BZ_dn=new complex<double>**[total_count_ext];
					FP_bulk_dispersion_c_count_dn=new int*[total_count_ext];
					FP_bulk_dispersion_c_BZ_count_dn=new int*[total_count_ext];
					FP_bulk_dispersion_mc_count_dn=new int*[total_count_ext];
					FP_bulk_dispersion_mc_BZ_count_dn=new int*[total_count_ext];
					FP_bulk_connection_c_dn=new int**[total_count_ext];
					FP_bulk_connection_c_BZ_dn=new int**[total_count_ext];
					FP_bulk_connection_mc_dn=new int**[total_count_ext];
					FP_bulk_connection_mc_BZ_dn=new int**[total_count_ext];
				}
			}
		}
		if(!PA_FPFS_Numerov){
			final_states_zgels_norm=new double*[total_count_ext];
		}
		// final_states_FP_bulk_dispersion=new double*[total_count_ext];
		int digit=(spin_i==2)?2:1;
		
		final_states_FP_nonloc_r=new double****[scale_width];
		int org_indices_spin_count=spin_i==0?1:2;
		int org_indices[scale_width*atom_length*org_indices_spin_count][3];
		int org_indices_count=0;
		for(int ie=0; ie<scale_width; ie++){
			final_states_FP_nonloc_r[ie]=new double***[org_indices_spin_count];
			for(int sp=0; sp<org_indices_spin_count; sp++){
				final_states_FP_nonloc_r[ie][sp]=new double**[atom_length];
				for(int ia=0; ia<atom_length; ia++){
					int is=atom_spec_index[ia];
					if(empty_atoms[is]){
						continue;
					}
					if(!PA_calc_all_nonloc && !atom_weighting_flag[ia]){
						continue;
					}		
					final_states_FP_nonloc_r[ie][sp][ia]=alloc_dmatrix(5, wfn_cutoff_index[is]);
					org_indices[org_indices_count][0]=ie;
					org_indices[org_indices_count][1]=sp;
					org_indices[org_indices_count][2]=ia;
					org_indices_count++;
				}
			}
		}

		// if PA_FPFS_Numerov=false, prepare the matrix buffer
		complex<double>*** left_matrix_buffer;
		complex<double>*** right_matrix_buffer;
	  if(!PA_FPFS_Numerov){
			int FP_g_count=0;
			double g_test[3];
			for(int n1=-n_range; n1<=n_range; n1++){
				for(int n2=-n_range; n2<=n_range; n2++){
					g_test[2]=0.0;
					for(int p=0; p<=1; p++){
						g_test[p]=rec_cell[1][p]*n1+rec_cell[2][p]*n2;
					}
					double g_length=sqrt(inner_product(g_test, g_test));
					if(g_length<khn_approx*PA_FPFS_kRange){
						FP_g_count++;
					}
				}
			}
			int z_count_new=VKS_count[0]-FPFS_z_start;
			int eq_dim=z_count_new*FP_g_count;
			int nrhs=2*FP_g_count;

			int num_threads=omp_get_max_threads();
			left_matrix_buffer=new complex<double>**[num_threads];
			right_matrix_buffer=new complex<double>**[num_threads];
			for(int it=0; it<num_threads; it++){
				left_matrix_buffer[it]=&alloc_zmatrix(eq_dim)[0];
				right_matrix_buffer[it]=&alloc_zmatrix(nrhs, eq_dim)[0];
			}
			
			sprintf(sprintf_buffer, "FP_g_count: %d", FP_g_count);
			write_log(sprintf_buffer);
		}

		complex<double>** Vgg0_matrix_bulk;
		complex<double>** Vgg1_matrix_bulk;
		complex<double>*** bulk_matrix_buffer;
		complex<double>*** bulk_VR_buffer;
		double kappaz_border;
		int kappaz_border_index;
		if(PA_FPFS_bulk_set){
			int FP_g_count_bulk=0;
			double g123[3];
			double g123_length;
			for(int n1=-n_range; n1<=n_range; n1++){
				int n1p=n1+n_range;
				for(int n2=-n_range; n2<=n_range; n2++){
					int n2p=n2+n_range;
					for(int n3=-n_range_z; n3<=n_range_z; n3++){
						for(int ix=0; ix<3; ix++){
							g123[ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
						}
						g123[2]+=n3*FPFS_gz_length;
						g123_length=sqrt(inner_product(g123, g123));
						if(g123_length<khn_approx*PA_FPFS_kRange){
							FP_g_count_bulk++;
						}
					}
				}
			}
			sprintf(sprintf_buffer, "FP_g_count_bulk: %d", FP_g_count_bulk);
			write_log(sprintf_buffer);
			final_states_FP_g_size_bulk=FP_g_count_bulk;
			final_states_FP_g_bulk=alloc_imatrix(FP_g_count_bulk, 3);
			final_states_FP_g_vec_bulk=alloc_dmatrix(FP_g_count_bulk, 3);
			FP_bulk_dispersion_up=alloc_dcube(total_count_ext, PA_FPFS_bulk_kz_steps, FP_g_count_bulk);
			if(spin_i==1 || spin_i==2){
				FP_bulk_dispersion_dn=alloc_dcube(total_count_ext, PA_FPFS_bulk_kz_steps, FP_g_count_bulk);
			}
			double dkz=FPFS_gz_length/(PA_FPFS_bulk_kz_steps*1.0);
			for(int ikz=0; ikz<PA_FPFS_bulk_kz_steps; ikz++){
				FP_bulk_dispersion_kz[ikz]=dkz*(ikz-PA_FPFS_bulk_kz_steps*0.5);
				//FP_bulk_dispersion_kz[ikz]=FPFS_gz_length*((ikz*1.0)/(PA_FPFS_bulk_kz_steps*1.0)-0.5);
			}

			// kappaz
			kappaz_border=PA_FPFS_gap_upper_limit/khn_approx;
			double dkappaz1=kappaz_border/(PA_FPFS_bulk_kappaz_steps_left*1.0);
			int kappaz_count1=PA_FPFS_bulk_kappaz_steps_left;
			double dkappaz2=FPFS_gz_length/(PA_FPFS_bulk_kappaz_steps_right*1.0);
			int kappaz_count2=ceil((khn_approx*PA_FPFS_kRange-kappaz_border)/dkappaz2);
			
			FP_bulk_kappaz_count=kappaz_count1+kappaz_count2;
			kappaz_border_index=kappaz_count1;
			FP_bulk_kappaz_border=kappaz_border_index;
			sprintf(sprintf_buffer, "kappaz count: %d", FP_bulk_kappaz_count);
			write_log(sprintf_buffer);
			FP_bulk_dispersion_kappaz=new double[FP_bulk_kappaz_count];
			FP_bulk_dispersion_mkappaz=new double[kappaz_count1];
			for(int ikz=0; ikz<kappaz_count1; ikz++){
				FP_bulk_dispersion_kappaz[ikz]=dkappaz1*(ikz+1.0);
				FP_bulk_dispersion_mkappaz[ikz]=-dkappaz1*(ikz+1.0);
			}
			for(int ikz=0; ikz<kappaz_count2; ikz++){
				FP_bulk_dispersion_kappaz[ikz+kappaz_count1]=kappaz_border+dkappaz2*(ikz+1.0);
			}

			//for(int ikz=0; ikz<FP_bulk_kappaz_count; ikz++){
			//	printf("%8.4f\n", FP_bulk_dispersion_kappaz[ikz]);
			//}
			
			int ig_count=0;
			for(int n1=-n_range; n1<=n_range; n1++){
				for(int n2=-n_range; n2<=n_range; n2++){
					for(int n3=-n_range_z; n3<=n_range_z; n3++){
						for(int ix=0; ix<3; ix++){
							g123[ix]=n1*rec_cell[1][ix]+n2*rec_cell[2][ix];
						}
						g123[2]+=n3*FPFS_gz_length;
						g123_length=sqrt(inner_product(g123, g123));
						if(g123_length<khn_approx*PA_FPFS_kRange){
							final_states_FP_g_bulk[ig_count][0]=n1;
							final_states_FP_g_bulk[ig_count][1]=n2;
							final_states_FP_g_bulk[ig_count][2]=n3;
							for(int ix=0; ix<3; ix++){
								final_states_FP_g_vec_bulk[ig_count][ix]=g123[ix];
							}
							ig_count++;
						}
					}
				}
			}			  
			Vgg0_matrix_bulk=alloc_zmatrix(FP_g_count_bulk);
			Vgg1_matrix_bulk=alloc_zmatrix(FP_g_count_bulk);
			for(int igb1=0; igb1<FP_g_count_bulk; igb1++){
				for(int igb2=0; igb2<FP_g_count_bulk; igb2++){
					int n1_12=final_states_FP_g_bulk[igb1][0]-final_states_FP_g_bulk[igb2][0];
					int n2_12=final_states_FP_g_bulk[igb1][1]-final_states_FP_g_bulk[igb2][1];
					int n3_12=final_states_FP_g_bulk[igb1][2]-final_states_FP_g_bulk[igb2][2];
					bool Vgg_found=false;
					for(int igb=0; igb<Vgg_count_bulk; igb++){
						if(Vgg_list_bulk[igb][0]==n1_12 && Vgg_list_bulk[igb][1]==n2_12 && Vgg_list_bulk[igb][2]==n3_12){
							Vgg0_matrix_bulk[igb1][igb2]=Vgg0_bulk[igb];
							if(spin_i==1 || spin_i==2){
								Vgg1_matrix_bulk[igb1][igb2]=Vgg1_bulk[igb];
							}
							Vgg_found=true;
							break;
						}
					}
					if(Vgg_found==false){
						write_log((char*)"Error: Vgg not found");
						return;
					}
				}
			}
			
			int num_threads=omp_get_max_threads();
			bulk_matrix_buffer=new complex<double>**[num_threads];
			bulk_VR_buffer=new complex<double>**[num_threads];
			for(int it=0; it<num_threads; it++){
				bulk_matrix_buffer[it]=alloc_zmatrix(FP_g_count_bulk);
				bulk_VR_buffer[it]=alloc_zmatrix(FP_g_count_bulk);
			}
		} // end of if(PA_FPFS_bulk_set)

		// Note: n_range are determined in the Fourier expansion of VKS
#pragma omp parallel firstprivate(scale_width, PA_ext_set, spin_i, num_bands, EF_Eh, Eh, PA_FPFS_energy_step, E_min_scale, E_max_scale, PA_excitation_energy, atom_above_surface, atom_length, n_range, digit, PA_FPFS_kRange, khn_approx, FPFS_z_start, n_range_z, final_states_FP_g_size_bulk) private(j)
#pragma omp for
		for(i=0; i<total_count_ext; i++){
			char* sprintf_buffer2=new char[Log_length+1];
			complex<double>** left_matrix;
			complex<double>** right_matrix;
			complex<double>** bulk_matrix;
			complex<double>** bulk_VR;
			if(!PA_FPFS_Numerov){
				int threadId=omp_get_thread_num();
				left_matrix=left_matrix_buffer[threadId];
				right_matrix=right_matrix_buffer[threadId];
			}
			if(PA_FPFS_bulk_set){
				int threadId=omp_get_thread_num();
				bulk_matrix=bulk_matrix_buffer[threadId];
				bulk_VR=bulk_VR_buffer[threadId];
			}
		  
			int count, count_up, count_dn;
			bool band_exists[scale_width];
			bool band_exists_up[scale_width];
			bool band_exists_dn[scale_width];
			double* k_point;
			int k_index=i;
			if(PA_ext_set){
				k_point=k_points_ext[i];
				k_index=k_index_reduced[k_index];
			}else{
				k_point=k_points[i];
			}
			int sp;
			int sp_max=(spin_i==0 || spin_i==2)? 1 : 2;
			for(sp=0; sp<sp_max;sp++){
				// initialization
				for(j=0; j<scale_width; j++){
					band_exists[j]=false;
					band_exists_up[j]=false;
					band_exists_dn[j]=false;
				}
			}
			count=0;
			count_up=0;
			count_dn=0;
			// fill
			for(j=0; j<num_bands; j++){
				double eigen;
				if(spin_i==0 || spin_i==2){
					eigen=(band[k_index][j]-EF_Eh)*Eh;
				}else if(sp==0){
					eigen=(band_up[k_index][j]-EF_Eh)*Eh;
				}else{
					eigen=(band_dn[k_index][j]-EF_Eh)*Eh;
				}
				int eigen_scale=round(eigen/PA_FPFS_energy_step);
				// printf("Eigen scale: %d\n", eigen_scale);
				int eigen_scale_offset=eigen_scale-E_min_scale;
				if(E_min_scale<= eigen_scale && eigen_scale<=E_max_scale){
					if(spin_i==0 || spin_i==2){
						if(band_exists[eigen_scale_offset]==false){
							count++;
							band_exists[eigen_scale_offset]=true;
						}
					}else if(sp==0){
						if(band_exists_up[eigen_scale_offset]==false){
							count_up++;
							band_exists_up[eigen_scale_offset]=true;
						}
					}else{
						if(band_exists_dn[eigen_scale_offset]==false){
							band_exists_dn[eigen_scale_offset]=true;
							count_dn++;
						}
					}
				}
			}
			// allocate matrices
			int FPIndex_size=0;
			if(spin_i==0){
				FPIndex_size=count;
			}else if(spin_i==1){
				FPIndex_size=count_up+count_dn;
			}else{
				FPIndex_size=count*2;
			}


			
			final_states_FP_size[i]=FPIndex_size;
			final_states_EScale[i]=new int[FPIndex_size];
			final_states_spin[i]=new int[FPIndex_size];
			final_states_k[i]=new double*[FPIndex_size];
			final_states_FP_nonloc[i]=new complex<double>***[FPIndex_size];
			final_states_FP_norm1[i]=new double*[FPIndex_size];
			final_states_FP_norm2[i]=new double*[FPIndex_size];
			
			final_states_FP_loc[i]=new complex<double>**[FPIndex_size];
			final_states_FP_g[i]=new int**[FPIndex_size];
			final_states_FP_g_vec[i]=new double**[FPIndex_size];
			final_states_FP_g_size[i]=new int[FPIndex_size];
			if(!PA_FPFS_bulk_set){
				final_states_FP_loc_edge[i]=new complex<double>*[FPIndex_size];
			}

			if(PA_FPFS_bulk_set){
				final_states_FP_bulk[i]=new complex<double>**[FPIndex_size];
				final_states_FP_bulk_count[i]=new int[FPIndex_size];
				//final_states_FP_bulk_z[i]=new complex<double>***[FPIndex_size];
				final_states_FP_bulk_kz[i]=new double*[FPIndex_size];
				final_states_FP_bulk_kappaz[i]=new double*[FPIndex_size];
				final_states_FP_bulk_coefs[i]=new complex<double>*[FPIndex_size];
			}
			if(!PA_FPFS_Numerov){
				final_states_zgels_norm[i]=new double[FPIndex_size];
			}

			for(j=0; j<FPIndex_size; j++){
				final_states_k[i][j]=new double[3];
				final_states_FP_nonloc[i][j]=alloc_zcube(atom_length, 5, 9);
				final_states_FP_norm1[i][j]=new double[atom_length];
				final_states_FP_norm2[i][j]=new double[atom_length];
				if(!PA_FPFS_Numerov){
					final_states_zgels_norm[i][j]=0.0;
				}
				for(int ia=0; ia<atom_length; ia++){
					int is=atom_spec_index[ia];
					if(empty_atoms[is]){
						continue;
					}
				}
			}
			int FPIndex=0;
			// fill
			for(sp=0; sp<sp_max;sp++){
				for(j=0; j<num_bands; j++){
					double eigen;
					if(spin_i==0 || spin_i==2){
						eigen=(band[k_index][j]-EF_Eh)*Eh;
					}else if(sp==0){
						eigen=(band_up[k_index][j]-EF_Eh)*Eh;
					}else{
						eigen=(band_dn[k_index][j]-EF_Eh)*Eh;
					}
					int eigen_scale=round(eigen/PA_FPFS_energy_step);
					if(E_min_scale<= eigen_scale && eigen_scale<=E_max_scale){
						bool needAddition=false;
						if(FPIndex==0){
							needAddition=true;
						}else{
							if((spin_i==0 || spin_i==2) && final_states_EScale[i][FPIndex-1]!=eigen_scale){
								needAddition=true;
							}else if(spin_i==1 && (final_states_EScale[i][FPIndex-1]!=eigen_scale || final_states_spin[i][FPIndex-1]!=sp)){
								needAddition=true;
							}
						}
						if(needAddition){
							if(spin_i==0){
								final_states_EScale[i][FPIndex]=eigen_scale;
								final_states_spin[i][FPIndex]=0;
								FPIndex++;
							}else if(spin_i==1){
								final_states_EScale[i][FPIndex]=eigen_scale;
								final_states_spin[i][FPIndex]=sp;
								FPIndex++;
							}else{
								final_states_EScale[i][FPIndex]=eigen_scale;
								final_states_spin[i][FPIndex]=0;
								final_states_EScale[i][FPIndex+1]=eigen_scale;
								final_states_spin[i][FPIndex+1]=1;
								FPIndex+=2;
							}
						}
					}
				}
			}
			 
			sprintf(sprintf_buffer2, "k=%4d, count=%3d", i, FPIndex_size);
			write_log(sprintf_buffer2);
			// bulk band calculations
			complex<double>** dispersion_c_up;
		  complex<double>** dispersion_c_dn;
			int** connection_c_up;
			int** connection_c_dn;
			int* dispersion_c_count_up;
			int* dispersion_c_count_dn;
			
			complex<double>** dispersion_mc_up;
		  complex<double>** dispersion_mc_dn;
			int** connection_mc_up;
			int** connection_mc_dn;
			int* dispersion_mc_count_up;
			int* dispersion_mc_count_dn;
			
			complex<double>** dispersion_c_BZ_up;
		  complex<double>** dispersion_c_BZ_dn;
			int** connection_c_BZ_up;
			int** connection_c_BZ_dn;
			int* dispersion_c_BZ_count_up;
			int* dispersion_c_BZ_count_dn;
			
			complex<double>** dispersion_mc_BZ_up;
		  complex<double>** dispersion_mc_BZ_dn;
			int** connection_mc_BZ_up;
			int** connection_mc_BZ_dn;
			int* dispersion_mc_BZ_count_up;
			int* dispersion_mc_BZ_count_dn;
			if(PA_FPFS_bulk_set && FPIndex_size>0){
				dispersion_c_count_up=new int[FP_bulk_kappaz_count];
				dispersion_c_count_dn=new int[FP_bulk_kappaz_count];
				dispersion_c_BZ_count_up=new int[kappaz_border_index];
				dispersion_c_BZ_count_dn=new int[kappaz_border_index];
				
				dispersion_mc_count_up=new int[FP_bulk_kappaz_count];
				dispersion_mc_count_dn=new int[FP_bulk_kappaz_count];
				dispersion_mc_BZ_count_up=new int[FP_bulk_kappaz_count];
				dispersion_mc_BZ_count_dn=new int[FP_bulk_kappaz_count];
				for(sp=0; sp<sp_max; sp++){
					complex<double>** Vgg_use;
					double** dispersion_use;
					complex<double>*** dispersion_c_pointer;
					int* dispersion_c_count_pointer;
					int*** connection_c_pointer;
					if(sp==0){
						Vgg_use=Vgg0_matrix_bulk;
						dispersion_use=FP_bulk_dispersion_up[i];
						dispersion_c_pointer=&dispersion_c_up;
						dispersion_c_count_pointer=dispersion_c_count_up;
						connection_c_pointer=&connection_c_up;
					}else{
						Vgg_use=Vgg1_matrix_bulk;
						dispersion_use=FP_bulk_dispersion_dn[i];
						dispersion_c_pointer=&dispersion_c_dn;
						dispersion_c_count_pointer=dispersion_c_count_dn;
						connection_c_pointer=&connection_c_dn;
					}
					calc_bulk_dispersion(k_point, PA_FPFS_bulk_kz_steps, FP_bulk_dispersion_kz, final_states_FP_g_size_bulk,
															 final_states_FP_g_vec_bulk, Vgg_use, dispersion_use, bulk_matrix);
					/*
					for(int ig=0; ig<final_states_FP_g_size_bulk; ig++){
						for(int ikz=0; ikz<PA_FPFS_bulk_kz_steps; ikz++){
							printf("%10.6f %10.6f\n", FP_bulk_dispersion_kz[ikz], FP_bulk_dispersion_up[i][ikz][ig]);
						}
						printf("\n");
						}*/

				  
					calc_bulk_dispersion_complex(k_point, 0.0, FP_bulk_kappaz_count, kappaz_border_index, FP_bulk_dispersion_kappaz, 
																			 final_states_FP_g_size_bulk, final_states_FP_g_bulk, final_states_FP_g_vec_bulk, Vgg_use,
																			 dispersion_c_pointer, dispersion_c_count_pointer, connection_c_pointer, bulk_matrix);
					/*
					for(int ikappaz=0; ikappaz<FP_bulk_kappaz_count-1; ikappaz++){
						for(int ib=0; ib<dispersion_c_count_pointer[ikappaz]; ib++){
							printf("%10.6f %10.6f %10.6f\n", FP_bulk_dispersion_kappaz[ikappaz], (*dispersion_c_pointer)[ikappaz][ib].real(), (*dispersion_c_pointer)[ikappaz][ib].imag());
							if((*connection_c_pointer)[ikappaz][ib]>=0){
								printf("%10.6f %10.6f %10.6f\n", FP_bulk_dispersion_kappaz[ikappaz+1], (*dispersion_c_pointer)[ikappaz+1][(*connection_c_pointer)[ikappaz][ib]].real(), (*dispersion_c_pointer)[ikappaz+1][(*connection_c_pointer)[ikappaz][ib]].imag());
							}
							printf("\n");
						}
						printf("\n");
						}*/
						// BZ, positive kappaz
					if(sp==0){
						Vgg_use=Vgg0_matrix_bulk;
						dispersion_c_pointer=&dispersion_c_BZ_up;
						dispersion_c_count_pointer=dispersion_c_BZ_count_up;
						connection_c_pointer=&connection_c_BZ_up;
					}else{
						Vgg_use=Vgg1_matrix_bulk;
						dispersion_c_pointer=&dispersion_c_BZ_dn;
						dispersion_c_count_pointer=dispersion_c_BZ_count_dn;
						connection_c_pointer=&connection_c_BZ_dn;
					}
					calc_bulk_dispersion_complex(k_point, FPFS_gz_length*0.5, kappaz_border_index, kappaz_border_index, FP_bulk_dispersion_kappaz, 
																			 final_states_FP_g_size_bulk, final_states_FP_g_bulk, final_states_FP_g_vec_bulk, Vgg_use,
																			 dispersion_c_pointer, dispersion_c_count_pointer, connection_c_pointer, bulk_matrix);

					// zero, negative kappaz
					if(sp==0){
						Vgg_use=Vgg0_matrix_bulk;
						dispersion_c_pointer=&dispersion_mc_up;
						dispersion_c_count_pointer=dispersion_mc_count_up;
						connection_c_pointer=&connection_mc_up;
					}else{
						Vgg_use=Vgg1_matrix_bulk;
						dispersion_c_pointer=&dispersion_mc_dn;
						dispersion_c_count_pointer=dispersion_mc_count_dn;
						connection_c_pointer=&connection_mc_dn;
					}
					calc_bulk_dispersion_complex(k_point, 0, kappaz_border_index, kappaz_border_index, FP_bulk_dispersion_mkappaz, 
																			 final_states_FP_g_size_bulk, final_states_FP_g_bulk, final_states_FP_g_vec_bulk, Vgg_use,
																			 dispersion_c_pointer, dispersion_c_count_pointer, connection_c_pointer, bulk_matrix);

					// BZ, negative kappaz
					if(sp==0){
						Vgg_use=Vgg0_matrix_bulk;
						dispersion_c_pointer=&dispersion_mc_BZ_up;
						dispersion_c_count_pointer=dispersion_mc_BZ_count_up;
						connection_c_pointer=&connection_mc_BZ_up;
					}else{
						Vgg_use=Vgg1_matrix_bulk;
						dispersion_c_pointer=&dispersion_mc_BZ_dn;
						dispersion_c_count_pointer=dispersion_mc_BZ_count_dn;
						connection_c_pointer=&connection_mc_BZ_dn;
					}
					calc_bulk_dispersion_complex(k_point, FPFS_gz_length*0.5, kappaz_border_index, kappaz_border_index, FP_bulk_dispersion_mkappaz, 
																			 final_states_FP_g_size_bulk, final_states_FP_g_bulk, final_states_FP_g_vec_bulk, Vgg_use,
																			 dispersion_c_pointer, dispersion_c_count_pointer, connection_c_pointer, bulk_matrix);
					//write_log((char*)"Bulk complex band calculation finished");

					/*
					for(int ikappaz=0; ikappaz<kappaz_border_index-1; ikappaz++){
						for(int ib=0; ib<dispersion_c_count_pointer[ikappaz]; ib++){
							printf("%10.6f %10.6f %10.6f\n", FP_bulk_dispersion_kappaz[ikappaz], (*dispersion_c_pointer)[ikappaz][ib].real(), (*dispersion_c_pointer)[ikappaz][ib].imag());
							if((*connection_c_pointer)[ikappaz][ib]>=0){
								printf("%10.6f %10.6f %10.6f\n", FP_bulk_dispersion_kappaz[ikappaz+1], (*dispersion_c_pointer)[ikappaz+1][(*connection_c_pointer)[ikappaz][ib]].real(), (*dispersion_c_pointer)[ikappaz+1][(*connection_c_pointer)[ikappaz][ib]].imag());
							}
							printf("\n");
						}
						printf("\n");
						}*/
				} // end of for(sp)
			} // end of if(PA_FPFS_bulk && FPIndex_size>0)
			/*
				for(j=0; j<FPIndex_size; j++){
				printf("Scale=%4d, Energy=%7.3f eV, spin=%d\n",
				final_states_EScale[i][j], final_states_EScale[i][j]*PA_FPFS_energy_step, final_states_spin[i][j]);
				}*/			
			// calculate the final states from the KS potential
			for(j=0; j<FPIndex_size; j++){
				// printf("k=%4d, FPIndex=%3d, local part\n", i, j);
				int eigen_scale=final_states_EScale[i][j];
				double kinetic_energy_Eh=(eigen_scale*PA_FPFS_energy_step+PA_excitation_energy)/Eh+EF_Eh;
				int eigen_spin=final_states_spin[i][j];
				// printf("Scale=%d, Ekin=%.2f, Spin=%d\n", eigen_scale, kinetic_energy_Eh, eigen_spin);
				double k_length_au=sqrt(2*kinetic_energy_Eh);
				double k_au[3]={k_point[0], k_point[1], 0.0}; // contains only in-plane components
				double kz_square=k_length_au*k_length_au-inner_product(k_au, k_au);
				if(kz_square<0){
					write_log((char*)"Warning: negative kz^2");
					final_states_FP_g_size[i][j]=0;
					continue;
				}
				double kz=sqrt(kz_square);
				final_states_k[i][j][2]=kz;
				for(int j1=0; j1<2; j1++){
					final_states_k[i][j][j1]=k_au[j1];
				}

				
				int FP_g_count=0;
				double g_test[3];
				double kpg_test[3];
				for(int n1=-n_range; n1<=n_range; n1++){
				  for(int n2=-n_range; n2<=n_range; n2++){
				    g_test[2]=0.0;
				    kpg_test[2]=0.0;
				    for(int p=0; p<=1; p++){
				      g_test[p]=rec_cell[1][p]*n1+rec_cell[2][p]*n2;
				      kpg_test[p]=g_test[p]+k_au[p];
				    }
				    double kpg_length=sqrt(inner_product(kpg_test, kpg_test));
				    double g_length=sqrt(inner_product(g_test, g_test));
				    if((PA_FPFS_Numerov && kpg_length<k_length_au) || (!PA_FPFS_Numerov && g_length<khn_approx*PA_FPFS_kRange)){
				      FP_g_count++;
				    }
				  }
				}

				
				int FP_bulk_count;
				// bulk eigenstate calculations
				if(PA_FPFS_bulk_set){
					complex<double>** Vgg_use;
					double** dispersion_use;
					complex<double>** dispersion_c_use;
					complex<double>** dispersion_c_BZ_use;
					int** connection_c_use;
					int** connection_c_BZ_use;
					int* dispersion_c_count_use;
					int* dispersion_c_BZ_count_use;
					complex<double>** dispersion_mc_use;
					complex<double>** dispersion_mc_BZ_use;
					int** connection_mc_use;
					int** connection_mc_BZ_use;
					int* dispersion_mc_count_use;
					int* dispersion_mc_BZ_count_use;
					if(final_states_spin[i][j]==0){
						Vgg_use=Vgg0_matrix_bulk;
						dispersion_use=FP_bulk_dispersion_up[i];
						dispersion_c_use=dispersion_c_up;
						dispersion_c_count_use=dispersion_c_count_up;
						connection_c_use=connection_c_up;
						dispersion_c_BZ_use=dispersion_c_BZ_up;
						dispersion_c_BZ_count_use=dispersion_c_BZ_count_up;
						connection_c_BZ_use=connection_c_BZ_up;
						dispersion_mc_use=dispersion_mc_up;
						dispersion_mc_count_use=dispersion_mc_count_up;
						connection_mc_use=connection_mc_up;
						dispersion_mc_BZ_use=dispersion_mc_BZ_up;
						dispersion_mc_BZ_count_use=dispersion_mc_BZ_count_up;
						connection_mc_BZ_use=connection_mc_BZ_up;
					}else{
						Vgg_use=Vgg1_matrix_bulk;
						dispersion_use=FP_bulk_dispersion_dn[i];
						dispersion_c_use=dispersion_c_dn;
						dispersion_c_count_use=dispersion_c_count_dn;
						connection_c_use=connection_c_dn;
						dispersion_c_BZ_use=dispersion_c_BZ_dn;
						dispersion_c_BZ_count_use=dispersion_c_BZ_count_dn;
						connection_c_BZ_use=connection_c_BZ_up;
						dispersion_mc_use=dispersion_mc_dn;
						dispersion_mc_count_use=dispersion_mc_count_dn;
						connection_mc_use=connection_mc_dn;
						dispersion_mc_BZ_use=dispersion_mc_BZ_dn;
						dispersion_mc_BZ_count_use=dispersion_mc_BZ_count_dn;
						connection_mc_BZ_use=connection_mc_BZ_up;
					}
					FP_bulk_count=solve_final_states_bulk(kinetic_energy_Eh, k_au, FPFS_gz_length, final_states_FP_g_size_bulk, final_states_FP_g_bulk,
																								final_states_FP_g_vec_bulk, Vgg_use, PA_FPFS_bulk_kz_steps, FP_bulk_dispersion_kz,
																								FP_bulk_kappaz_count, kappaz_border_index, FP_g_count,
																								FP_bulk_dispersion_kappaz, FP_bulk_dispersion_mkappaz,
																								dispersion_use, dispersion_c_use, dispersion_c_count_use, connection_c_use,
																								dispersion_c_BZ_use, dispersion_c_BZ_count_use, connection_c_BZ_use,
																								dispersion_mc_use, dispersion_mc_count_use, connection_mc_use,
																								dispersion_mc_BZ_use, dispersion_mc_BZ_count_use, connection_mc_BZ_use,
																								&final_states_FP_bulk[i][j], &final_states_FP_bulk_kz[i][j], &final_states_FP_bulk_kappaz[i][j],
																								bulk_matrix, bulk_VR);
					final_states_FP_bulk_count[i][j]=FP_bulk_count;
					final_states_FP_bulk_coefs[i][j]=new complex<double>[FP_bulk_count];
					// continue;
				} // end of if(PA_FPFS_bulk_set)
				
				// printf("k=(%.2f, %.2f, %.2f)\n", final_states_k[i][j][0], final_states_k[i][j][1], final_states_k[i][j][2]);
				
				// solve the Schroedinger equation using the Fourier expansion in the xy plane
				// composite list
				final_states_FP_g_size[i][j]=FP_g_count;
				final_states_FP_g[i][j]=alloc_imatrix(FP_g_count, 2);
				final_states_FP_g_vec[i][j]=alloc_dmatrix(FP_g_count, 3);
				if(!PA_FPFS_bulk_set){
					final_states_FP_loc_edge[i][j]=new complex<double>[FP_g_count];
					for(int ig=0; ig<FP_g_count; ig++){
						final_states_FP_loc_edge[i][j][ig]=complex<double>(0.0, 0.0);
					}
				}
				complex<double>*** final_states_FP_bulk_z;
				if(PA_FPFS_bulk_set){
					final_states_FP_bulk_z=alloc_zcube(FP_bulk_count, FP_g_count, VKS_count[0]);
				}

				int ig_count=0;
				for(int n1=-n_range; n1<=n_range; n1++){
				  for(int n2=-n_range; n2<=n_range; n2++){
				    g_test[2]=0.0;
				    kpg_test[2]=0.0;
				    for(int p=0; p<=1; p++){
				      g_test[p]=rec_cell[1][p]*n1+rec_cell[2][p]*n2;
				      kpg_test[p]=g_test[p]+k_au[p];
				    }
				    double kpg_length=sqrt(inner_product(kpg_test, kpg_test));
				    double g_length=sqrt(inner_product(g_test, g_test));
				    if((PA_FPFS_Numerov && kpg_length<k_length_au) || (!PA_FPFS_Numerov && g_length<khn_approx*PA_FPFS_kRange)){
				      final_states_FP_g[i][j][ig_count][0]=n1;
				      final_states_FP_g[i][j][ig_count][1]=n2;
				      for(int p=0; p<2; p++){
								final_states_FP_g_vec[i][j][ig_count][p]=g_test[p];
				      }
				      final_states_FP_g_vec[i][j][ig_count][2]=0.0;
				      ig_count++;
				    }
				  }
				}
				
				if(PA_FPFS_bulk_set){
					// fill final_states_FP_bulk_z
					int bulk_start=ceil(FPFS_bulk_min/final_states_dz);
					for(int in=0; in<FP_bulk_count; in++){
						bool FP_g_error=false;
						for(int igb=0; igb<final_states_FP_g_size_bulk; igb++){
							int ig_found=-1;
							for(int ig=0; ig<FP_g_count; ig++){
								if(final_states_FP_g_bulk[igb][0]==final_states_FP_g[i][j][ig][0] &&
									 final_states_FP_g_bulk[igb][1]==final_states_FP_g[i][j][ig][1]){
									ig_found=ig;
									break;
								}
							}
							if(ig_found<0){
								write_log((char*)"Error: FP_g not found");
								FP_g_error=true;
								break;
							}
							double kpgz=final_states_FP_bulk_kz[i][j][in]+final_states_FP_g_vec_bulk[igb][2];
							for(int iz=bulk_start; iz<=FPFS_z_start; iz++){
								double z_bulk=final_states_dz*iz-FPFS_bulk_max;
								complex<double> phase(cos(kpgz*z_bulk), sin(kpgz*z_bulk));
								double tail=exp(final_states_FP_bulk_kappaz[i][j][in]*z_bulk);
								final_states_FP_bulk_z[in][ig_found][iz]+=phase*final_states_FP_bulk[i][j][in][igb]*tail;
							}
						}
						if(FP_g_error){
							continue;
						}
						// printf("kz: %f\n", final_states_FP_bulk_kz[i][j][in]);
					}
					/*
					for(int in=0; in<FP_bulk_count; in++){
						for(int ig=0; ig<FP_g_count; ig++){
							printf("%d %d\n\n", final_states_FP_g[i][j][ig][0], final_states_FP_g[i][j][ig][1]);
							for(int iz=0; iz<VKS_count[0]; iz++){
								printf("%4d %8.4f %8.4f\n", iz, final_states_FP_bulk_z[in][ig][iz].real(), final_states_FP_bulk_z[in][ig][iz].imag());
							}
							printf("\n");
						}
						}*/
				}
				/*
					for(int ig=0; ig<final_states_FP_g_size[i][j]; ig++){
					printf("(n1, n2) = (%2d, %2d), g = (%6.3f, %6.3f, %6.3f)\n",
					final_states_FP_g[i][j][ig][0],
					final_states_FP_g[i][j][ig][1],
					final_states_FP_g_vec[i][j][ig][0],
					final_states_FP_g_vec[i][j][ig][1],
					final_states_FP_g_vec[i][j][ig][2]);
					}*/
				
				// prapare the Vgg matrix
				complex<double>*** Vgg_matrix=alloc_zpmatrix(FP_g_count);
				int V00_index=-1;
				bool Vgg_error=false;
				for(int ig1=0; ig1<FP_g_count; ig1++){
					if(final_states_FP_g[i][j][ig1][0]==0 && final_states_FP_g[i][j][ig1][1]==0){
						V00_index=ig1;
					}
					for(int ig2=0; ig2<FP_g_count; ig2++){
						int n1_12=final_states_FP_g[i][j][ig1][0]-final_states_FP_g[i][j][ig2][0];
						int n2_12=final_states_FP_g[i][j][ig1][1]-final_states_FP_g[i][j][ig2][1];
						bool Vgg_found=false;
						for(int ig=0; ig<Vgg_count; ig++){
							if(n1_12==Vgg_list[ig][0] && n2_12==Vgg_list[ig][1]){
								if(final_states_spin[i][j]==0){
									Vgg_matrix[ig1][ig2]=&Vgg0[ig][0];
								}else{
									Vgg_matrix[ig1][ig2]=&Vgg1[ig][0];
								}
								Vgg_found=true;
								/*
									for(int iz=0; iz<VKS_count[0]; iz++){
									printf("%f %f\n", Vgg_matrix[ig1][ig2][iz].real(), Vgg_matrix[ig1][ig2][iz].imag());
									}*/
								break;
							}
						}
						if(Vgg_found==false){
							write_log((char*)"Error: Vgg not found");
							Vgg_error=true;
						}
					}
				}
				if(Vgg_error==true){
					continue;
				}

				// prepare the solution buffer
				final_states_FP_loc[i][j]=alloc_zmatrix(FP_g_count, VKS_count[0]);
				if(PA_FPFS_Numerov){
					if(FP_g_count<1){
						write_log((char*)"Warning: No g vector is found");
					}else{
						// solve the differential equation by the Numerov method
						solve_final_state_Numerov(kinetic_energy_Eh, k_au, kz, FP_g_count, VKS_count[0],
																			final_states_dz, FPFS_z_start, V00_index, Vgg_matrix, final_states_FP_g_vec[i][j],
																			final_states_FP_loc[i][j]);
					}
				}else{
					if(!PA_FPFS_bulk_set){
						final_states_zgels_norm[i][j]
							=solve_final_state_Matrix(kinetic_energy_Eh, k_au, kz, FP_g_count, VKS_count[0],
																				final_states_dz, FPFS_z_start, V00_index, Vgg_matrix, final_states_FP_g_vec[i][j],
																				final_states_FP_loc[i][j], left_matrix, right_matrix, final_states_FP_loc_edge[i][j]);
					}else{
						final_states_zgels_norm[i][j]
							=solve_final_state_from_bulk(kinetic_energy_Eh, k_au, kz, FP_g_count, VKS_count[0], FP_bulk_count,
																					 final_states_dz, FPFS_z_start, V00_index, Vgg_matrix, final_states_FP_g_vec[i][j],
																					 final_states_FP_bulk_z, final_states_FP_loc[i][j], left_matrix, right_matrix,
																					 final_states_FP_bulk_coefs[i][j]);
						//solve_final_state_from_bulk_perturbation(kinetic_energy_Eh, k_au, kz, FP_g_count, VKS_count[0], FP_bulk_count,
						//																			 final_states_dz, FPFS_z_start, V00_index, Vgg_matrix, final_states_FP_g_vec[i][j],
						//																				 final_states_FP_bulk_z, final_states_FP_loc[i][j],
						//																				 final_states_FP_bulk_coefs[i][j]);
						delete_zcube(final_states_FP_bulk_z);
					}
				}
				if(FP_g_count>0){
					delete_zpmatrix(Vgg_matrix);
				}
				// printf("k=%4d, FPIndex=%3d, nonlocal part\n", i, j);
			} // for(j<FPIndex_size)
			if(PA_FPFS_bulk_set){
				if(FPIndex_size>0){
					if(!PA_FPFS_file_set){
						delete_zmatrix(dispersion_c_up);
						delete_imatrix(connection_c_up);
						delete_zmatrix(dispersion_c_BZ_up);
						delete_imatrix(connection_c_BZ_up);
						delete_zmatrix(dispersion_mc_up);
						delete_imatrix(connection_mc_up);
						delete_zmatrix(dispersion_mc_BZ_up);
						delete_imatrix(connection_mc_BZ_up);
						if(spin_i>0){
							delete_zmatrix(dispersion_c_dn);
							delete_imatrix(connection_c_dn);
							delete_zmatrix(dispersion_c_BZ_dn);
							delete_imatrix(connection_c_BZ_dn);
							delete_zmatrix(dispersion_mc_dn);
							delete_imatrix(connection_mc_dn);
							delete_zmatrix(dispersion_mc_BZ_dn);
							delete_imatrix(connection_mc_BZ_dn);
						}
						delete[] dispersion_c_count_up;
						delete[] dispersion_c_count_dn;
						delete[] dispersion_c_BZ_count_up;
						delete[] dispersion_c_BZ_count_dn;
						delete[] dispersion_mc_count_up;
						delete[] dispersion_mc_count_dn;
						delete[] dispersion_mc_BZ_count_up;
						delete[] dispersion_mc_BZ_count_dn;
					}else{
						FP_bulk_dispersion_c_up[i]    =dispersion_c_up;
						FP_bulk_dispersion_c_BZ_up[i] =dispersion_c_BZ_up;
						FP_bulk_dispersion_mc_up[i]   =dispersion_mc_up;
						FP_bulk_dispersion_mc_BZ_up[i]=dispersion_mc_BZ_up;
						FP_bulk_dispersion_c_count_up[i]    =dispersion_c_count_up;
						FP_bulk_dispersion_c_BZ_count_up[i] =dispersion_c_BZ_count_up;
						FP_bulk_dispersion_mc_count_up[i]   =dispersion_mc_count_up;
						FP_bulk_dispersion_mc_BZ_count_up[i]=dispersion_mc_BZ_count_up;
						FP_bulk_connection_c_up[i]    =connection_c_up;
						FP_bulk_connection_c_BZ_up[i] =connection_c_BZ_up;
						FP_bulk_connection_mc_up[i]   =connection_mc_up;
						FP_bulk_connection_mc_BZ_up[i]=connection_mc_BZ_up;
						if(spin_i>0){
							FP_bulk_dispersion_c_dn[i]    =dispersion_c_dn;
							FP_bulk_dispersion_c_BZ_dn[i] =dispersion_c_BZ_dn;
							FP_bulk_dispersion_mc_dn[i]   =dispersion_mc_dn;
							FP_bulk_dispersion_mc_BZ_dn[i]=dispersion_mc_BZ_dn;
							FP_bulk_dispersion_c_count_dn[i]    =dispersion_c_count_dn;
							FP_bulk_dispersion_c_BZ_count_dn[i] =dispersion_c_BZ_count_dn;
							FP_bulk_dispersion_mc_count_dn[i]   =dispersion_mc_count_dn;
							FP_bulk_dispersion_mc_BZ_count_dn[i]=dispersion_mc_BZ_count_dn;
							FP_bulk_connection_c_dn[i]    =connection_c_dn;
							FP_bulk_connection_c_BZ_dn[i] =connection_c_BZ_dn;
							FP_bulk_connection_mc_dn[i]   =connection_mc_dn;
							FP_bulk_connection_mc_BZ_dn[i]=connection_mc_BZ_dn;
						}
					}
				}
			}
			delete[] sprintf_buffer2;
		} // for(total_count_ext)

		// for debug
		/*
		for(int ik=0; ik<total_count_ext; ik++){
			for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
				for(int in=0; in<final_states_FP_bulk_count[ik][ifp]; in++){
					printf("%4d %8.4f\n", final_states_EScale[ik][ifp], final_states_FP_bulk_kz[ik][ifp][in]);
				}
			}
			printf("\n");
			}*/
		
		if(!PA_FPFS_Numerov){
			int num_threads=omp_get_max_threads();
			for(int it=0; it<num_threads; it++){
				delete_zmatrix(left_matrix_buffer[it]);
				delete_zmatrix(right_matrix_buffer[it]);
			}
		}
		if(PA_FPFS_bulk_set){
			int num_threads=omp_get_max_threads();
			for(int it=0; it<num_threads; it++){
				delete_zmatrix(bulk_matrix_buffer[it]);
				delete_zmatrix(bulk_VR_buffer[it]);
			}
		}

		if(!PA_ignore_core && !PA_ignore_nonlocal){
			write_log((char*)"---- Connection calculations between the core and valence ----");
#pragma omp parallel firstprivate(atom_length, PA_lp_max, PA_theta_points, PA_ext_set) private(j)
#pragma omp for
			for(i=0; i<total_count_ext; i++){
				double* k_point;
				int k_index=i;
				if(PA_ext_set){
					k_point=k_points_ext[i];
					k_index=k_index_reduced[k_index];
				}else{
					k_point=k_points[i];
				}
				double k_au[3]={k_point[0], k_point[1], 0.0};
				for(j=0; j<final_states_FP_size[i]; j++){
					// obtain the connection condition at rc
					complex<double> Ylm_k[6][11];
					complex<double> Ylm_kp[6][11];
					for(int ia=0; ia<atom_length; ia++){
						int is=atom_spec_index[ia];
						if(empty_atoms[is]){
							continue;
						}
						if(!PA_calc_all_nonloc && !atom_weighting_flag[ia]){
							continue;
						}		
						double tau_z=atom_coordinates[ia][2];
						double g_vec[3];
						double gp_vec[3];
						double kpg_vec[3];
						double kpgp_vec[3];
						final_states_FP_norm1[i][j][ia]=0.0;
						final_states_FP_norm2[i][j][ia]=0.0;
						complex<double> p1jlp[12]={
							complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1),
							complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1),
							complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1)
						};
						complex<double> m1jlp[12]={
							complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1),
							complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1),
							complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1)
						};
						int ir=wfn_cutoff_index[is];
						for(int l=0; l<5; l++){
							for(int m=-l; m<=l; m++){
								// by theta integral
								int mpl=m+l;
								final_states_FP_nonloc[i][j][ia][l][mpl]=complex<double>(0,0);
								for(int lp=0; lp<=PA_lp_max; lp++){
									if(abs(m)>lp){
										continue;
									}
									int mplp=m+lp;
									for(int ig=0; ig<final_states_FP_g_size[i][j]; ig++){
										complex<double> theta_integral(0,0);
										for(int ix=0; ix<3; ix++){
											g_vec[ix]=final_states_FP_g_vec[i][j][ig][ix];
											kpg_vec[ix]=g_vec[ix]+k_au[ix];
										}
										spherical_harmonics(kpg_vec, &Ylm_k[0][0]);
										double kpg_length=sqrt(inner_product(kpg_vec, kpg_vec));
										double kt=inner_product(kpg_vec, atom_coordinates[ia]);
										complex<double> atom_phase(cos(kt), sin(kt));
										for(int it=0; it<PA_theta_points; it++){
											double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
											double z_value=wfn_r_use[is][ir]*cos(theta)+tau_z;
											theta_integral+=sin(theta)*interpolate_fgz(z_value, final_states_FP_loc[i][j][ig], final_states_dz, VKS_count[0])
												*spherical_harmonic_theta(lp, m, theta)
												*spherical_harmonic_theta(l, m, theta)
												*M_PI/(PA_theta_points*1.0);
										}
										final_states_FP_nonloc[i][j][ia][l][mpl]+=4*M_PI*wfn_r_use[is][ir]*p1jlp[lp]*atom_phase
											*sp_bessel(lp, kpg_length*wfn_r_use[is][ir])*conj(Ylm_k[lp][mplp])*theta_integral;
									}
								}
								final_states_FP_norm1[i][j][ia]+=norm(final_states_FP_nonloc[i][j][ia][l][mpl])/wfn_r_use[is][ir]/wfn_r_use[is][ir];
								// by Lebedev integral
								/*
									double Lebedev_r[3][PA_Lebedev_order_int];
									double Lebedev_w[PA_Lebedev_order_int];
									ld_by_order(PA_Lebedev_order_int, Lebedev_r[0], Lebedev_r[1], Lebedev_r[2], Lebedev_w);
									complex<double> Ylm_Le[6][11];
									complex<double> nonloc_lebedev(0,0);
									double g_vec[3];
									double kpg_vec[3];
									for(int ile=0; ile<PA_Lebedev_order_int; ile++){
									double r_le[3];
									for(int ix=0; ix<3; ix++){
									r_le[ix]=Lebedev_r[ix][ile]*wfn_r[is][ir];
									}
									spherical_harmonics(r_le, &Ylm_Le[0][0]);
									for(int ix=0; ix<3; ix++){
									r_le[ix]+=atom_coordinates[ia][ix];
									}
								
									for(int ig=0; ig<final_states_FP_g_size[i][j]; ig++){
									for(int ix=0; ix<3; ix++){
									g_vec[ix]=final_states_FP_g_vec[i][j][ig][ix];
									kpg_vec[ix]=g_vec[ix]+k_au[ix];
									}
									double kpgt=inner_product(kpg_vec, r_le);
									nonloc_lebedev+=4*M_PI*complex<double>(cos(kpgt), sin(kpgt))*interpolate_fgz(r_le[2], final_states_FP_loc[i][j][ig], final_states_dz, VKS_count[0])
									*conj(Ylm_Le[l][mpl])*wfn_r[is][ir]*Lebedev_w[ile];
									}
									}*/
								// printf("Norm1 theta integral [%2d][%2d]= (%8.4f, %8.4f)\n", l, m, final_states_FP_nonloc[i][j][ia][l][mpl].real(), final_states_FP_nonloc[i][j][ia][l][mpl].imag());
								// printf("Norm1 Lebed integral [%2d][%2d]= (%8.4f, %8.4f)\n\n", l, m, nonloc_lebedev.real(), nonloc_lebedev.imag());
							
							} // for(m)
						} // for(l)
						// debug
						// printf("Norm1[%d][%d][%d]= %8.4f\n", i, j, ia, final_states_FP_norm1[i][j][ia]);
					
						//norm2
						// by theta integral
						complex<double> norm2_temporary(0,0);
						for(int l=0; l<=PA_lp_max; l++){
							for(int lp=0; lp<=PA_lp_max; lp++){
								int l_min=min(l, lp);
								for(int ig=0; ig<final_states_FP_g_size[i][j]; ig++){
									for(int ix=0; ix<3; ix++){
										g_vec[ix]=final_states_FP_g_vec[i][j][ig][ix];
										kpg_vec[ix]=g_vec[ix]+k_au[ix];
									}
									double kpg_length=sqrt(inner_product(kpg_vec, kpg_vec));
									spherical_harmonics(kpg_vec, &Ylm_k[0][0]);
									for(int igp=0; igp<final_states_FP_g_size[i][j]; igp++){
										for(int ix=0; ix<3; ix++){
											gp_vec[ix]=final_states_FP_g_vec[i][j][igp][ix];
											kpgp_vec[ix]=gp_vec[ix]+k_au[ix];
										}
										double kpgp_length=sqrt(inner_product(kpgp_vec, kpgp_vec));
										double ggt=inner_product(g_vec, atom_coordinates[ia])-inner_product(gp_vec, atom_coordinates[ia]);
										complex<double> g_phase(cos(ggt), sin(ggt));
										spherical_harmonics(kpgp_vec, &Ylm_kp[0][0]);
										for(int m=-l_min; m<=l_min; m++){
											int mpl=m+l;
											int mplp=m+lp;
											complex<double> theta_integral(0,0);
											for(int it=0; it<PA_theta_points; it++){
												double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
												double z_value=wfn_r_use[is][ir]*cos(theta)+tau_z;
												theta_integral+=
													sin(theta)*conj(interpolate_fgz(z_value, final_states_FP_loc[i][j][igp], final_states_dz, VKS_count[0]))
													*interpolate_fgz(z_value, final_states_FP_loc[i][j][ig], final_states_dz, VKS_count[0])
													*spherical_harmonic_theta(lp, m, theta)
													*spherical_harmonic_theta(l, m, theta)
													*M_PI/(PA_theta_points*1.0);
											}
											norm2_temporary+=16.0*M_PI*M_PI*p1jlp[l]*m1jlp[lp]*Ylm_kp[lp][mplp]*conj(Ylm_k[l][mpl])
												*sp_bessel(lp, kpgp_length*wfn_r_use[is][ir])*sp_bessel(l, kpg_length*wfn_r_use[is][ir])*theta_integral*g_phase;
										}
									}										
								}
							}
						}
						final_states_FP_norm2[i][j][ia]=abs(norm2_temporary);
						// by lebedev integral
						/*
							double Lebedev_r[3][PA_Lebedev_order_int];
							double Lebedev_w[PA_Lebedev_order_int];
							ld_by_order(PA_Lebedev_order_int, Lebedev_r[0], Lebedev_r[1], Lebedev_r[2], Lebedev_w);
							complex<double> nonloc_lebedev(0,0);
							for(int ile=0; ile<PA_Lebedev_order_int; ile++){
							double r_le[3];
							for(int ix=0; ix<3; ix++){
							r_le[ix]=Lebedev_r[ix][ile]*wfn_r[is][ir];
							r_le[ix]+=atom_coordinates[ia][ix];
							}
								
							for(int ig=0; ig<final_states_FP_g_size[i][j]; ig++){
							for(int ix=0; ix<3; ix++){
							g_vec[ix]=final_states_FP_g_vec[i][j][ig][ix];
							kpg_vec[ix]=g_vec[ix]+k_au[ix];
							}
							double kpgt=inner_product(kpg_vec, r_le);
							for(int igp=0; igp<final_states_FP_g_size[i][j]; igp++){
							for(int ix=0; ix<3; ix++){
							gp_vec[ix]=final_states_FP_g_vec[i][j][igp][ix];
							kpgp_vec[ix]=gp_vec[ix]+k_au[ix];
							}
							double kpgpt=inner_product(kpgp_vec, r_le);
							nonloc_lebedev+=4*M_PI*complex<double>(cos(kpgt), sin(kpgt))*complex<double>(cos(kpgpt), -sin(kpgpt))
							*interpolate_fgz(r_le[2], final_states_FP_loc[i][j][ig], final_states_dz, VKS_count[0])
							*conj(interpolate_fgz(r_le[2], final_states_FP_loc[i][j][igp], final_states_dz, VKS_count[0]))*Lebedev_w[ile];
							}
							}
							}*/
						//printf("Norm2 theta integral [%d][%d][%d]=(%8.4f,   0.0000)\n", i, j, ia, final_states_FP_norm2[i][j][ia]);
						//printf("Norm2 Lebed integral [%d][%d][%d]=(%8.4f, %8.4f)\n", i, j, ia, nonloc_lebedev.real(), nonloc_lebedev.imag());
					
					
						// debug
						// printf("Norm2[%d][%d][%d]= %8.4f\n", i, j, ia, final_states_FP_norm2[i][j][ia]);
					} //for (ia=0; ia<atom_length; ia++)
				} // for(j=0; j<FPIndex_size; j++)
			} // omp for(i=0; i<total_count_ext; i++)
		}

		// nonlocal core part
		if(!PA_ignore_core && !PA_ignore_nonlocal){
#pragma omp parallel firstprivate(E_min_scale, org_indices_count, EF_Eh, Eh, PA_excitation_energy)
#pragma omp for
			for(int iepa=0; iepa<org_indices_count; iepa++){
				char* sprintf_buffer2=new char[Log_length+1];
				int ie=org_indices[iepa][0];
				int sp=org_indices[iepa][1];
				int ia=org_indices[iepa][2];
				int is=atom_spec_index[ia];
				int eigen_scale=ie+E_min_scale;
				sprintf(sprintf_buffer2, "Index %6d/%6d, EScale: %4d, Spin: %1d, Atom: %3d", iepa+1, org_indices_count, eigen_scale, sp, ia);
				write_log(sprintf_buffer2);
				double kinetic_energy_Eh=(eigen_scale*PA_FPFS_energy_step+PA_excitation_energy)/Eh+EF_Eh;
				// obtain the nonlocal radial function
				double wfn_buffer[vps_cutoff_index[is]];
				for(int l=0; l<5; l++){
					if(sp==0){
						solve_nonlocal_wfn(kinetic_energy_Eh, l, vps_cutoff_index[is], VPS_r[is], VKS0_r[ia], VPS_l_length[is], VPS_l[is], VPS_E_ave[is], VPS_nonloc_ave[is], wfn_buffer);
					}else{
						solve_nonlocal_wfn(kinetic_energy_Eh, l, vps_cutoff_index[is], VPS_r[is], VKS1_r[ia], VPS_l_length[is], VPS_l[is], VPS_E_ave[is], VPS_nonloc_ave[is], wfn_buffer);
					}
					// conversion for wfn_r_use
					for(int ir=0; ir<wfn_cutoff_index[is]; ir++){
						if(wfn_r_use[is][ir]<VPS_r[is][0]){
							final_states_FP_nonloc_r[ie][sp][ia][l][ir]=wfn_buffer[0];
						}
						if(wfn_r_use[is][ir]>VPS_r[is][vps_cutoff_index[is]-1]){
							final_states_FP_nonloc_r[ie][sp][ia][l][ir]=wfn_buffer[vps_cutoff_index[is]-1];
						}	 
						for(int irp=0; irp<vps_cutoff_index[is]-1; irp++){
							if(VPS_r[is][irp] < wfn_r_use[is][ir] && wfn_r_use[is][ir] < VPS_r[is][irp+1]){
								final_states_FP_nonloc_r[ie][sp][ia][l][ir]=(wfn_buffer[irp]*(VPS_r[is][irp+1]-wfn_r_use[is][ir])+
																														 wfn_buffer[irp+1]*(wfn_r_use[is][ir]-VPS_r[is][irp]))/(VPS_r[is][irp+1]-VPS_r[is][irp]);
							}
						}
					}
				}
				delete[] sprintf_buffer2;
			} // for(iepa)
		} // if(!PA_ignore_core)
	} // if(FPFS)

	/// Load VPS file for Ignore_core
	if(PA_ignore_core && !PA_FPFS){
		write_log((char*)"----Load the pseudopotential files to determine the core----");
		H5File VPS_db(PA_VPS_file, H5F_ACC_RDONLY);
		for(int is=0; is<atom_spec_length; is++){
			Group VPSG(VPS_db.openGroup(atom_spec_PS[is]));
			VPS_cutoff[is]=r_att_double(VPSG, "max_cutoff");
			sprintf(sprintf_buffer, "%2s: %4.2f Bohr", atom_spec_label[is], VPS_cutoff[is]);
			write_log(sprintf_buffer);
		}
	}

	/// radial integration for orthogonality correction
	if(PA_orth_correction){
		write_log((char*)"----Radial integration for orthgonality correction----");
		for(int is=0; is<atom_spec_length; is++){
			Self_radial_int[is]=alloc_dmatrix(num_orbits[is]);
			for(int io1=0; io1<num_orbits[is]; io1++){
				for(int io2=0; io2<num_orbits[is]; io2++){
					Self_radial_int[is][io1][io2]=ddot(&wfn_length_use[is], &wfn_phi_rdr[is][io1][0], &wfn_phi[is][io2][0]);
					// printf("%5.2f ", Self_radial_int[is][io1][io2]);
				}
				// printf("\n");
			}
			// printf("\n");
		}
	}

	/// 3-5. calculation for each k point
	double dispersion2[total_count_ext][num_points_E];
	for(i=0; i<total_count_ext; i++){
		for(j=0; j<num_points_E; j++){
			dispersion2[i][j]=0.0;
		}
	}
	int ik; // for k point
	int ib; // for band index
	int sp; // for spin: 0=Up, 1=Dn
	int sp_max=(spin_i==0 || spin_i==2) ? 1 : 2;
	int il; // for l
	int ir; // for r
	complex<double> m1jlp[12]={
		complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1),
		complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1),
		complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1)
	};
	complex<double> Ylm_k[6][11];
	complex<double> Ylm_k2[6][11];
	double Gaunt_arr[4][7][5][9]; // [l][m+l][lp][mp+lp];
	int ip, jp;
	for(i=0; i<4; i++){
		for(j=0; j<7; j++){
			for(ip=0; ip<5; ip++){
				for(jp=0; jp<9; jp++){
					Gaunt_arr[i][j][ip][jp]=Gaunt(ip, jp-ip, i, j-i);
					// printf("%e ", Gaunt_arr[i][j][ip][jp]);
				}
			}
		}
	}
	
  end=chrono::system_clock::now();
  duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	sprintf(sprintf_buffer, "PAD preparation time: %.3f [ms]", duration);
	write_log(sprintf_buffer);
#pragma omp parallel firstprivate(Y_coeff, sp_max, num_bands, EF_Eh, Eh, PA_E_min, PA_E_pixel, num_points_E, tail_index,  m1jlp, Gaunt_arr, atom_length, atom_spec_index, num_orbits, spin_i,  atom_coordinates, PA_ext_set, PA_final_state_step, k_index_min, atom_weighting_flag, atom_weighting, PA_weighting, axis_au, PA_reflection, PA_reflection_coef, PA_FPFS, PA_FPFS_energy_step, E_min_scale) private(ib, sp, Ylm_k, Ylm_k2, ia, is, io, il, ir, j)
#pragma omp for
	for(ik=0; ik<total_count_ext; ik++){
		char* sprintf_buffer2=new char[Log_length+1];
		// cout << ik << endl;
		int ik_reduced=ik;
		double* k_point;
		double k_point2[3];
		if(PA_ext_set){
			ik_reduced=k_index_reduced[ik];
			k_point=k_points_ext[ik];
		}else{
			k_point=k_points[ik];
		}
		double k_length=sqrt(inner_product(k_point, k_point));
		spherical_harmonics(k_point, &Ylm_k[0][0]);
		
		if(PA_reflection){
			double k_inner=inner_product(k_point, axis_au);
			int ix;
			for(ix=0; ix<3; ix++){
				k_point2[ix]=k_point[ix]-axis_au[ix]*2*k_inner;
			}
			spherical_harmonics(k_point2, &Ylm_k2[0][0]);
		}
		// printf("%lf %lf %lf | %lf %lf %lf\n", k_point[0], k_point[1], k_point[2], k_point2[0], k_point2[1], k_point2[2]);
		/*
			for(int ii=0; ii<5; ii++){
			for(int jj=0; jj<9; jj++){
			char buffer[1024];
			sprintf(buffer, "%8.3f %8.3f", Ylm_k[ii][jj].real(), Ylm_k[ii][jj].imag());
			write_log(buffer);
			}
			write_log((char*)"");
			}*/

		// FPFS radial integration preparation
		// (FPFS_re, FPFS_im)=fgz(r*cos(theta)+tau_z)
		/*
		complex<double>***** FPFS_rt;
		if(PA_FPFS){
			
		  FPFS_rt=new complex<double>****[final_states_FP_size[ik]]; //[FPIndex][ig][ia][it][ir]
		  for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
				double kz_value=final_states_k[ik][ifp][2];
		    FPFS_rt[ifp]=new complex<double>***[final_states_FP_g_size[ik][ifp]];
		    for(int ig=0; ig<final_states_FP_g_size[ik][ifp]; ig++){
		      FPFS_rt[ifp][ig]=new complex<double>**[atom_length];
		      for(int ia=0; ia<atom_length; ia++){
						if(atom_weighting_flag[ia]==false){
							continue;
						}
						int is=atom_spec_index[ia];
						double tau_z=atom_coordinates[ia][2];
						FPFS_rt[ifp][ig][ia]=new complex<double>*[PA_theta_points];
						// printf("tau_z %7.3f r_max %7.3f\n", tau_z, wfn_r[is][wfn_length[is]-1]);
						for(int it=0; it<PA_theta_points; it++){
							double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
							FPFS_rt[ifp][ig][ia][it]=new complex<double>[wfn_length[is]];
							for(int ir=0; ir<wfn_length[is]; ir++){
								double z_value=wfn_r[is][ir]*cos(theta)+tau_z;
								FPFS_rt[ifp][ig][ia][it][ir]=interpolate_fgz(z_value, final_states_FP_loc[ik][ifp][ig], final_states_dz, VKS_count[0]);
							}
						}
		      }
		    }
				}
		}*/
		for(sp=0; sp<sp_max; sp++){
			for(ib=0; ib<num_bands; ib++){
				double eigen;
				if(spin_i==0 || spin_i==2){
					eigen=(band[ik_reduced][ib]-EF_Eh)*Eh;
				}else if(sp==0){
					eigen=(band_up[ik_reduced][ib]-EF_Eh)*Eh;
				}else{
					eigen=(band_dn[ik_reduced][ib]-EF_Eh)*Eh;
				}
				
				int eigen_index=round((eigen-PA_E_min)/PA_E_pixel);
				if(eigen_index-tail_index>=num_points_E){
					break;
				}
				if(eigen_index+tail_index<0){
					continue;
				}
				int eigen_scale=round(eigen/PA_FPFS_energy_step);
				// printf("Eigen scale: %d\n", eigen_scale);
				// find appropriate FPIndex
				int FPIndex_1=-1; // for all types of spin_i
				int FPIndex_2=-1; // for spin_i==2
				int FPIndex_size;
				if(PA_FPFS){
					FPIndex_size=final_states_FP_size[ik];
					for(j=0; j<FPIndex_size; j++){
						// cout << final_states_EScale[ik][j] << endl;
						if(final_states_EScale[ik][j]==eigen_scale){
							if(spin_i==0){
								FPIndex_1=j;
								break;
							}else if(spin_i==1){
								if(final_states_spin[ik][j]==sp){
									FPIndex_1=j;
									break;
								}
							}else if(spin_i==2){
								if(final_states_spin[ik][j]==0){
									FPIndex_1=j;
								}else{
									FPIndex_2=j;
								}
								if(FPIndex_1>=0 && FPIndex_2>=0){
									break;
								}
							}
						}
					}
					// printf("FPIndex = %d, %d\n", FPIndex_1, FPIndex_2);
					if(FPIndex_1<0 || (spin_i==2 && FPIndex_2<0)){
						write_log((char*)"Error: FP not found");
						continue;
					}
				}
				// PAD matrix element for up and down spins
				complex<double> PAD_1(0, 0);
				complex<double> PAD_2(0, 0);
				if(strcmp(PA_output_data, "Band")==0){
					PAD_1=complex<double>(1, 0);
					PAD_2=complex<double>(1, 0);
				}else{
					for(ia=0; ia<atom_length; ia++){
						// Norm of initial and final states for orthgonality correction
						complex<double> Norm_FI_1(0, 0);
						complex<double> Norm_FI_2(0, 0);
						if(PA_weighting==true && atom_weighting_flag[ia]==false){
							continue;
						}
						if(PA_FPFS==false){
							double kt=inner_product(k_point, atom_coordinates[ia]);
							complex<double> atom_phase(cos(kt), -sin(kt));
							if(PA_weighting==true){
								atom_phase*=atom_weighting[ia];
							}
							// cout << kt << endl;
							is=atom_spec_index[ia];
							int num_orbits2=num_orbits[is];
							if(spin_i==1){
								num_orbits2*=2;
							}
							double final_states[5][wfn_length_use[is]];
							for(il=0; il<5; il++){
								for(ir=0; ir<wfn_length_use[is]; ir++){
									final_states[il][ir]=0.0;
								}
							}
							int k_index;
							if(strcmp(PA_final_state, "PW")==0 && k_length>1e-5){
								for(il=0; il<5; il++){
									for(ir=0; ir<wfn_length_use[is]; ir++){
										if(PA_ignore_core && !empty_atoms[is] && wfn_r_use[is][ir]<VPS_cutoff[is]){
											final_states[il][ir]=0.0;
										}else{
											final_states[il][ir]=wfn_r_use[is][ir]*sp_bessel(il, wfn_r_use[is][ir]*k_length);
										}
									}
								}
							}else if(strcmp(PA_final_state, "Calc")==0){
								k_index=round(k_length/PA_final_state_step);
								for(il=0; il<5; il++){
									for(ir=0; ir<wfn_length_use[is]; ir++){
										if(PA_ignore_core && !empty_atoms[is] && wfn_r_use[is][ir]<VPS_cutoff[is]){
											final_states[il][ir]=0.0;
										}else{
											final_states[il][ir]=final_states_calc[k_index-k_index_min][is][il][ir];
										}
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
								complex<double>** LCAO_use=LCAO[ia][io][ik_reduced][ib]; // [twoLp1][digit]
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
									double radial_part=ddot(&wfn_length_use[is], &final_states[lp][0], &wfn_phi_rdr[is][io2][0]);
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
											if(PA_reflection==false){
												coeff1+=Ylm_k[lp][mpplp]*Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1];
											}else{
												coeff1+=(Ylm_k[lp][mpplp]+PA_reflection_coef*Ylm_k2[lp][mpplp])*Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1];
											}
											// cout << Ylm_k[lp][mpplp] << " ";
										}
										// cout << coeff1 << endl;
										coeff2_1+=LCAO_use[mpl][0]*coeff1;
										if(spin_i==2){
											coeff2_2+=LCAO_use[mpl][1]*coeff1;
										}
									}
									// cout << endl;
									complex<double> atom_phase2;
									if(strcmp(PA_final_state, "Calc")==0){
										atom_phase2=atom_phase*complex<double>(cos(final_states_phase[k_index-k_index_min][is][lp]), sin(final_states_phase[k_index-k_index_min][is][lp]));
									}else{
										atom_phase2=atom_phase;
									}
									// cout << atom_phase2 << endl;
									PAD_1+=m1jlp[lp]*radial_part*atom_phase2*coeff2_1;
									if(spin_i==2){
										PAD_2+=m1jlp[lp]*radial_part*atom_phase2*coeff2_2;
									}
								} // end of for(dl)
								if(PA_add_nonorth_term || PA_orth_correction){
									double e_vec_re[3];
									double e_vec_im[3];
									for(int ix=0; ix<3; ix++){
										e_vec_re[ix]=e_vec[ix].real();
										e_vec_im[ix]=e_vec[ix].imag();
									}
									double et_real=inner_product(e_vec_re, atom_coordinates[ia]);
									double et_imag=inner_product(e_vec_im, atom_coordinates[ia]);
									complex<double> et(et_real, et_imag);
									double radial_part=ddot(&wfn_length_use[is], &final_states[l][0], &wfn_phi_dr[is][io2][0]);
									complex<double> coeff2_1(0, 0);
									complex<double> coeff2_2(0, 0);
									int mpl;
									for(mpl=0; mpl<=2*l; mpl++){
										int m=mpl-l;
										complex<double> coeff1(0, 0);
										if(PA_reflection==false){
											coeff1+=Ylm_k[l][mpl];
										}else{
											coeff1+=(Ylm_k[l][mpl]+PA_reflection_coef*Ylm_k2[l][mpl]);
										}
										coeff2_1+=LCAO_use[mpl][0]*coeff1;
										if(spin_i==2){
											coeff2_2+=LCAO_use[mpl][1]*coeff1;
										}
									}
									// cout << endl;
									complex<double> atom_phase2;
									if(strcmp(PA_final_state, "Calc")==0){
										atom_phase2=atom_phase*complex<double>(cos(final_states_phase[k_index-k_index_min][is][l]), sin(final_states_phase[k_index-k_index_min][is][l]));
									}else{
										atom_phase2=atom_phase;
									}
									// cout << atom_phase2 << endl;
									if(PA_add_nonorth_term){
										PAD_1+=m1jlp[l]*radial_part*atom_phase2*coeff2_1*et;
										if(spin_i==2){
											PAD_2+=m1jlp[l]*radial_part*atom_phase2*coeff2_2*et;
										}
									}
									if(PA_orth_correction){
										Norm_FI_1+=m1jlp[l]*radial_part*atom_phase2*coeff2_1;
										if(spin_i==2){
											Norm_FI_2+=m1jlp[l]*radial_part*atom_phase2*coeff2_2;
										}
									}
								}
							}
						}else{
							// FPFS==true
							double g_vec[3];
							double kpg_vec[3]; // the z component is zero
							is=atom_spec_index[ia];
							double* final_states_re=new double[wfn_length_use[is]];
							double* final_states_im=new double[wfn_length_use[is]];
							double* sp_bessel_lp=new double[wfn_length_use[is]];
							if(final_states_FP_g_size[ik][FPIndex_1]==0){
								continue;
							}
							complex<double>** FPFS_rt=alloc_zmatrix(PA_theta_points, wfn_length_use[is]);
							for(int ig=0; ig<final_states_FP_g_size[ik][FPIndex_1]; ig++){
								// prepare FPFS_rt
								double kz_value=final_states_k[ik][FPIndex_1][2];
								double tau_z=atom_coordinates[ia][2];
								for(int it=0; it<PA_theta_points; it++){
									double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
									for(int ir=0; ir<wfn_length_use[is]; ir++){
										double z_value=wfn_r_use[is][ir]*cos(theta)+tau_z;
										FPFS_rt[it][ir]=interpolate_fgz(z_value, final_states_FP_loc[ik][FPIndex_1][ig], final_states_dz, VKS_count[0]);
									}
								}
				  
								// printf("ik %3d ib %3d ia %3d ig %3d\n", ik, ib, ia, ig);
								for(int ix=0; ix<3; ix++){
									g_vec[ix]=final_states_FP_g_vec[ik][FPIndex_1][ig][ix];
									kpg_vec[ix]=g_vec[ix]+k_point[ix];
								}
								// the z component is discarded
								g_vec[2]=0.0;
								kpg_vec[2]=0.0;
								// printf("kpg = %7.4f %7.4f %7.4f\n", kpg_vec[0], kpg_vec[1], kpg_vec[2]);
								spherical_harmonics(kpg_vec, &Ylm_k[0][0]);
								/*
									for(int il=0; il<6; il++){
									for(int im=0; im<=2*il; im++){
									printf("(%7.4f, %7.4f) ", Ylm_k[il][im].real(), Ylm_k[il][im].imag());
									}
									printf("\n");
									}
									printf("\n");*/
								double kpgt=inner_product(kpg_vec, atom_coordinates[ia]);
								complex<double> atom_phase(cos(kpgt), -sin(kpgt));
								atom_phase/=sqrt(2.0*M_PI);
								double kpg_length=sqrt(inner_product(kpg_vec, kpg_vec));
								if(PA_weighting==true){
									atom_phase*=atom_weighting[ia];
								}
								
								int num_orbits2=num_orbits[is];
								if(spin_i==1){
									num_orbits2*=2;
								}
								for(io=0; io<num_orbits2; io++){
									int io2=io; // for initial state radial wfn, l_list
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
									complex<double>** LCAO_use=LCAO[ia][io][ik_reduced][ib]; // [twoLp1][digit]
									// local (valence) part
									// l: azimuthal quantum number of initial state
									// lp: of final state
									// m: magnetic quantum number of initial state
									// mpl: m+l (0, ..., 2l)
									// jp1: j+1
									// mp: of final state, mp=m+j
									// mpplp: mp+lp=m+j+lp
									int l=l_list[is][io2];
									double* wfn_init=new double[wfn_length_use[is]];
									for(int ir=0; ir<wfn_length_use[is]; ir++){
										if(strcmp(PA_initial_state, "PAO")==0){
											wfn_init[ir]=wfn_phi_PAO[is][io2][ir];
										}else{
											if(ir<wfn_cutoff_index[is]){
												wfn_init[ir]=wfn_phi_AO[is][io2][ir];
											}else{
												wfn_init[ir]=wfn_phi_PAO[is][io2][ir];
											}
										}
										
									}
									for(int lp=0; lp<=PA_lp_max; lp++){
										// prepare the spherical Bessel function
										for(int ir=0; ir<wfn_length_use[is]; ir++){
											//if(PA_ignore_core && wfn_r[is][ir]<VPS_cutoff[is]){
											if(!empty_atoms[is] && ir<wfn_cutoff_index[is] && !PA_ignore_nonlocal){
												sp_bessel_lp[ir]=0.0;
											}else{
												sp_bessel_lp[ir]=sp_bessel(lp, kpg_length*wfn_r_use[is][ir]);
											}
											//printf("%d %12.7f %12.7f %12.7f\n", lp, wfn_r[is][ir], sp_bessel_lp[ir], kpg_length*wfn_r[is][ir]);
										}
										int mpl;
										for(mpl=0; mpl<=2*l; mpl++){
											// cout << mpl << endl;
											int m=mpl-l;
											int jp1St=max(-1, -(m+lp))+1; // include
											int jp1En=min(1, lp-m)+2; // not include
											int jp1;
											complex<double> coeff11(0, 0);
											//complex<double> coeff12(0, 0);
											//complex<double> coeff21(0, 0);
											//complex<double> coeff22(0, 0);
											for(jp1=jp1St; jp1<jp1En; jp1++){
												int mpplp=jp1-1+m+lp;
												// printf("l %2d, m %2d, lp %2d, mp %2d, j %2d\n", l, m, lp, mpplp-lp, jp1-1);
												if(spin_i==0 || spin_i==1){
													// digit=0 |1> -> |1> or |0> -> |0>
													// load final state
													complex<double> integral_rt(0, 0);
													for(int it=0; it<PA_theta_points; it++){
														double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
														// printf("%6.4f ", theta);
														for(ir=0; ir<wfn_length_use[is]-1; ir++){
															final_states_re[ir]=(wfn_r_use[is][ir+1]-wfn_r_use[is][ir])*wfn_r_use[is][ir]*wfn_r_use[is][ir]
																*sp_bessel_lp[ir]*FPFS_rt[it][ir].real();
															final_states_im[ir]=(wfn_r_use[is][ir+1]-wfn_r_use[is][ir])*wfn_r_use[is][ir]*wfn_r_use[is][ir]
																*sp_bessel_lp[ir]*FPFS_rt[it][ir].imag();
															//printf("%7.4f (%7.4f, %7.4f) (%7.4f, %7.4f) %7.4f\n", wfn_r[is][ir], final_states_re[ir], final_states_im[ir], FPFS_rt[FPIndex_1][ig][ia][it][ir].real(), FPFS_rt[FPIndex_1][ig][ia][it][ir].imag(), sp_bessel_lp[ir]);
														}
														final_states_re[wfn_length_use[is]-1]=0.0;
														final_states_im[wfn_length_use[is]-1]=0.0;
														double radial_part_re=ddot(&wfn_length_use[is], &final_states_re[0], &wfn_init[0]);
														double radial_part_im=ddot(&wfn_length_use[is], &final_states_im[0], &wfn_init[0]);
														complex<double> radial_part(radial_part_re, -radial_part_im);
														integral_rt+=radial_part*sin(theta)
															*spherical_harmonic_theta(lp, mpplp-lp, theta)
															*spherical_harmonic_theta(1, jp1-1, theta)
															*spherical_harmonic_theta(l, m, theta)*M_PI/(PA_theta_points*1.0);
													}
													//if(final_states_FP_g[ik][FPIndex_1][ig][0]==0 && final_states_FP_g[ik][FPIndex_1][ig][1]==0){
													//printf("lp %2d mp %2d l %2d m %2d real %7.4f imag %7.4f\n", lp, mpplp-lp, l, m, integral_rt.real(), integral_rt.imag());
													//printf("lp %2d mp %2d l %2d m %2d added real %7.4f imag %7.4f\n", lp, mpplp-lp, l, m, (Ylm_k[lp][mpplp]*integral_rt*Y_coeff[jp1]).real(), (Ylm_k[lp][mpplp]*integral_rt*Y_coeff[jp1]).imag());
													//printf("Ylm_k real %7.4f imag %7.4f\n", Ylm_k[lp][mpplp].real(), Ylm_k[lp][mpplp].imag());
													//}
													coeff11+=Ylm_k[lp][mpplp]*integral_rt*Y_coeff[jp1];
													// printf("\n");
												}else{
													// spin_i=2 (nc): not yet implemented
												}
											}
											coeff11*=atom_phase*m1jlp[lp];
											// cout << coeff1 << endl;
											if(spin_i==0 || spin_i==1){
												PAD_1+=LCAO_use[mpl][0]*coeff11;
											}else{
												//PAD_1+=LCAO_use[mpl][0]*coeff11+LCAO_use[mpl][1]*coeff21;
												//PAD_2+=LCAO_use[mpl][0]*coeff21+LCAO_use[mpl][1]*coeff22;
											}
											if(PA_add_nonorth_term || PA_orth_correction){
												if(lp<abs(m)){
													continue;
												}
												double e_vec_re[3];
												double e_vec_im[3];
												for(int ix=0; ix<3; ix++){
													e_vec_re[ix]=e_vec[ix].real();
													e_vec_im[ix]=e_vec[ix].imag();
												}
												double et_real=inner_product(e_vec_re, atom_coordinates[ia]);
												double et_imag=inner_product(e_vec_im, atom_coordinates[ia]);
												complex<double> et(et_real, et_imag);
												
												complex<double> integral_rt(0, 0);
												for(int it=0; it<PA_theta_points; it++){
													double theta=(it*1.0+0.5)/(PA_theta_points*1.0)*M_PI;
													// printf("%6.4f ", theta);
													for(ir=0; ir<wfn_length_use[is]-1; ir++){
														final_states_re[ir]=(wfn_r_use[is][ir+1]-wfn_r_use[is][ir])*wfn_r_use[is][ir]
															*sp_bessel_lp[ir]*FPFS_rt[it][ir].real();
														final_states_im[ir]=(wfn_r_use[is][ir+1]-wfn_r_use[is][ir])*wfn_r_use[is][ir]
															*sp_bessel_lp[ir]*FPFS_rt[it][ir].imag();
													}
													final_states_re[wfn_length_use[is]-1]=0.0;
													final_states_im[wfn_length_use[is]-1]=0.0;
													double radial_part_re=ddot(&wfn_length_use[is], &final_states_re[0], &wfn_init[0]);
													double radial_part_im=ddot(&wfn_length_use[is], &final_states_im[0], &wfn_init[0]);
													complex<double> radial_part(radial_part_re, -radial_part_im);
													integral_rt+=radial_part*sin(theta)
														*spherical_harmonic_theta(lp, m, theta)
														*spherical_harmonic_theta(l, m, theta)*M_PI/(PA_theta_points*1.0);
												}
												if(PA_add_nonorth_term){
													// Note: sqrt(2*pi) is due to atom_phase division by sqrt(2*pi) above
													if(spin_i==0 || spin_i==1){
														PAD_1+=LCAO_use[mpl][0]*atom_phase*sqrt(2.0*M_PI)*m1jlp[lp]*Ylm_k[lp][mpl]*integral_rt*et;
													}
												}
												if(PA_orth_correction){
													if(spin_i==0 || spin_i==1){
														Norm_FI_1+=LCAO_use[mpl][0]*atom_phase*sqrt(2.0*M_PI)*m1jlp[lp]*Ylm_k[lp][mpl]*integral_rt;
													}
												}
											}// end of Add_nonorth_term
										}// end of for(mpl)
									}// end of for(lp)
									delete[] wfn_init;
								}// end of for(io)
							}//end of for(ig)
							delete[] final_states_re;
							delete[] final_states_im;
							delete[] sp_bessel_lp;
							delete_zmatrix(FPFS_rt);
							
							//continue;
							// nonlocal (core) part
							// l: azimuthal quantum number of initial state
							// lp: of final state, lp=l+dl
							// m: magnetic quantum number of initial state
							// mpl: m+l (0, ..., 2l)
							// jp1: j+1
							// mp: of final state, mp=m+j
							// mpplp: mp+lp=m+j+lp
							if(!PA_ignore_core && !empty_atoms[is] && !PA_ignore_nonlocal){
								int ie=eigen_scale-E_min_scale;
								int num_orbits2=num_orbits[is];
								if(spin_i==1){
									num_orbits2*=2;
								}
								for(io=0; io<num_orbits2; io++){
									int io2=io; // for initial state radial wfn, l_list
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
									complex<double>** LCAO_use=LCAO[ia][io][ik_reduced][ib]; // [twoLp1][digit]
									int l=l_list[is][io2];
									int dl; // -1 or +1
									double atom_weight=1.0;
									if(PA_weighting){
										atom_weight=atom_weighting[ia];
									}
									for(dl=-1; dl<=1; dl+=2){
										int lp=l+dl;
										if(lp<0){
											continue;
										}
										double radial_part=ddot(&wfn_cutoff_index[is], &final_states_FP_nonloc_r[ie][sp][ia][lp][0], &wfn_phi_rdr[is][io2][0]);
										//printf("Radial: %f\n", radial_part);
										int mpl;
										complex<double> coeff2_1(0, 0);
										for(mpl=0; mpl<=2*l; mpl++){
											int m=mpl-l;
											int jp1St=max(-1, -(m+lp))+1; // include
											int jp1En=min(1, lp-m)+2; // not include
											int jp1;
											complex<double> coeff1(0, 0);
											for(jp1=jp1St; jp1<jp1En; jp1++){
												int mpplp=jp1-1+m+lp;
												// printf("l m lp mp: %2d %2d %2d %2d\n", l, m, lp, mpplp-lp);
												// printf("Edge: %f %f\n", final_states_FP_nonloc[ik][FPIndex_1][ia][lp][mpplp].real(), final_states_FP_nonloc[ik][FPIndex_1][ia][lp][mpplp].imag());
												coeff1+=Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1]*conj(final_states_FP_nonloc[ik][FPIndex_1][ia][lp][mpplp]);
												// cout << Ylm_k[lp][mpplp] << " ";
											}
											// cout << coeff1 << endl;
											coeff2_1+=LCAO_use[mpl][0]*coeff1;
										}
										// cout << endl;
										PAD_1+=radial_part*coeff2_1*atom_weight/(4.0*M_PI);
									} // end of for(dl)
									if(PA_add_nonorth_term || PA_orth_correction){
										double e_vec_re[3];
										double e_vec_im[3];
										for(int ix=0; ix<3; ix++){
											e_vec_re[ix]=e_vec[ix].real();
											e_vec_im[ix]=e_vec[ix].imag();
										}
										double et_real=inner_product(e_vec_re, atom_coordinates[ia]);
										double et_imag=inner_product(e_vec_im, atom_coordinates[ia]);
										complex<double> et(et_real, et_imag);
										double radial_part=ddot(&wfn_cutoff_index[is], &final_states_FP_nonloc_r[ie][sp][ia][l][0], &wfn_phi_dr[is][io2][0]);
										complex<double> coeff2_1(0, 0);
										int mpl;
										for(mpl=0; mpl<=2*l; mpl++){
											int m=mpl-l;
											coeff2_1+=LCAO_use[mpl][0]*conj(final_states_FP_nonloc[ik][FPIndex_1][ia][l][mpl]);
										}
										// cout << endl;
										if(PA_add_nonorth_term){
											PAD_1+=radial_part*coeff2_1*et*atom_weight/(4.0*M_PI);
										}
										if(PA_orth_correction){
											Norm_FI_1+=radial_part*coeff2_1*atom_weight/(4.0*M_PI);
										}
									} // end of add_nonorth_term
								} // end of for(io)
							} // end of PA_ignore_core
							// printf("ig end\n");
						} // end of if(FPFS==false)
						
						// Orthogonality correction
						if(PA_orth_correction){
							double sum_cia_1=0.0;
							double sum_cia_2=0.0;
							int num_orbits2=num_orbits[is];
							if(spin_i==1){
								num_orbits2*=2;
							}
							for(io=0; io<num_orbits2; io++){
								int io2=io; // for initial state radial wfn, l_list
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
								complex<double>** LCAO_use=LCAO[ia][io][ik_reduced][ib]; // [twoLp1][digit]
								int l=l_list[is][io2];
								for(int mpl=0; mpl<=2*l; mpl++){
									sum_cia_1+=norm(LCAO_use[mpl][0]);
									if(spin_i==2){
										sum_cia_2+=norm(LCAO_use[mpl][1]);
									}
								}
							}
							complex<double> t1=conj(Norm_FI_1)/sum_cia_1;
							// printf("t1: (%9.2e, %9.2e)\n", t1.real(), t1.imag()); 
							complex<double> t2(0,0);
							if(spin_i==2){
								t2=conj(Norm_FI_2)/sum_cia_2;
								// printf("t2: (%9.2e, %9.2e)\n", t2.real(), t2.imag());
							}
							// if(isfinite(t1.real()) && isfinite(t1.imag())){
							if(abs(sum_cia_1)>PA_zero_threshold){
								// io is for ket
								// iop is for bra
							  //printf("t1: (%9.2e, %9.2e)\n", t1.real(), t1.imag()); 
								complex<double> Self_matrix_element1(0,0);
								complex<double> Self_matrix_element2(0,0);
								for(io=0; io<num_orbits2; io++){
									int io2=io; // for initial state radial wfn, l_list
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
									complex<double>** LCAO_use=LCAO[ia][io][ik_reduced][ib]; // [twoLp1][digit]
									int l=l_list[is][io2];
									for(int mpl=0; mpl<=2*l; mpl++){
										int m=mpl-l;
										for(int iop=0; iop<num_orbits2; iop++){
											int iop2=iop;
											if(spin_i==1){
												if(sp==0 && iop%2==1){
													continue;
												}
												if(sp==1 && iop%2==0){
													continue;
												}
												iop2=iop/2;
											}
											int lp=l_list[is][iop2];
											complex<double>** LCAO_p_use=LCAO[ia][iop][ik_reduced][ib];
											if(lp==l+1 || lp==l-1){
												int jp1St=max(-1, -(m+lp))+1; // include
												int jp1En=min(1, lp-m)+2; // not include
												int jp1;
												for(jp1=jp1St; jp1<jp1En; jp1++){
													int mpplp=jp1-1+m+lp;
													Self_matrix_element1+=Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1]*conj(LCAO_p_use[mpplp][0])*LCAO_use[mpl][0]*Self_radial_int[is][iop2][io2];
													if(spin_i==2){
														Self_matrix_element2+=Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1]*conj(LCAO_p_use[mpplp][1])*LCAO_use[mpl][1]*Self_radial_int[is][iop2][io2];
													}
												}
											}
										} // end of for(iop)
									} // end of for(m)
								} // end of for(io)
								PAD_1-=conj(t1)*Self_matrix_element1;
								if(spin_i==2){
									PAD_2-=conj(t2)*Self_matrix_element2;
								}
							}
						} // end of orthogonality correction
						
					} // end of for(ia)
				} // end of if(output_data=="Band")
					// cout << PAD_1 << endl;
				double PAD_norm;
				if(strcmp(PA_output_data, "PAD")==0 || strcmp(PA_output_data, "Band")==0){
					PAD_norm=norm(PAD_1);
					if(spin_i==2){
						PAD_norm+=norm(PAD_2);
					}
				}
				
				// cout << "!" << PAD_norm << endl;
				for(j=-tail_index; j<=tail_index; j++){
					if(eigen_index+j>=0 && eigen_index+j<num_points_E){
						dispersion2[ik][eigen_index+j]+=tail[abs(j)]*PAD_norm;
					}
				}
				
			} // end of for(ib)
		} // end of for(sp)
		if(PA_FPFS){
			/*
			for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
				for(int ig=0; ig<final_states_FP_g_size[ik][ifp]; ig++){
					for(int ia=0; ia<atom_length; ia++){
						if(atom_weighting_flag[ia]==false){
							continue;
						}
						for(int it=0; it<PA_theta_points; it++){
							delete[] FPFS_rt[ifp][ig][ia][it];
						}
						delete[] FPFS_rt[ifp][ig][ia];
					}
					delete[] FPFS_rt[ifp][ig];
				}
				delete[] FPFS_rt[ifp];
			}
			delete[] FPFS_rt;*/
			sprintf(sprintf_buffer2, "k = %4d finished", ik);
			write_log(sprintf_buffer2);
		}
		delete[] sprintf_buffer2;
	} // end of for(ik) and omp
	
  end=chrono::system_clock::now();
  duration=chrono::duration_cast<chrono::milliseconds>(end-start).count();
	sprintf(sprintf_buffer, "PAD calculation time: %.3f [ms]", duration);
	write_log(sprintf_buffer);

	double disp_max=0.0;
	for(i=0; i<total_count_ext; i++){
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
		double* dispersion3_alloc=new double[total_count_ext*num_points_E];
		dispersion3=new double**[kx_count_ext];
		for(i=0; i<kx_count_ext; i++){
			dispersion3[i]=new double*[ky_count_ext];
			for(j=0; j<ky_count_ext; j++){
				dispersion3[i][j]=&dispersion3_alloc[i*ky_count_ext*num_points_E+j*num_points_E];
			}
		}
		for(i=0; i<kx_count_ext; i++){
			for(j=0; j<ky_count_ext; j++){
				for(k=0; k<num_points_E; k++){
					dispersion3[i][j][k]=dispersion2[i+j*kx_count_ext][k];
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
	Group rootG(output.openGroup("/"));
	w_att_str(rootG, "Datetime", time_str);
	w_att_int(rootG, "Dimension", dimension);

	if(dimension==1){
		w_data_2d(rootG, "Dispersion", total_count_ext, num_points_E, (double**) dispersion2);
	}else{
		w_data_3d(rootG, "Dispersion", kx_count_ext, ky_count_ext, num_points_E, (double***) &dispersion3[0][0][0]);
	}

	if(dimension==1){
		double kx_offset=kx_length*kx_range[0];
		double E_offset=PA_E_min;
		if(PA_ext_set){
			kx_offset=kx_length*(kx_range[0]-(kx_range[1]-kx_range[0])*PA_ext_le);
		}
		double offset2[2]={kx_offset, E_offset};
		w_att_1d(rootG, "Offset", 2, &offset2[0]);

		double delta2[2]={dkx_length, PA_E_pixel};
		w_att_1d(rootG, "Delta", 2, &delta2[0]);
		
		int Size2[2]={total_count_ext, num_points_E};
		w_att_1i(rootG, "Size", 2, &Size2[0]);
		
		w_att_1d(rootG, "Xvector", 3, &kx_vector[0]);
	}else{
		double kx_offset=kx_length*kx_range[0];
		double ky_offset=ky_length*ky_range[0];
		double E_offset=PA_E_min;
		if(PA_ext_set){
			kx_offset=kx_length*(kx_range[0]-(kx_range[1]-kx_range[0])*PA_ext_le);
			ky_offset=ky_length*(ky_range[0]-(ky_range[1]-ky_range[0])*PA_ext_dn);
		}
		double offset3[3]={kx_offset, ky_offset, E_offset};
		w_att_1d(rootG, "Offset", 3, &offset3[0]);

		double delta3[3]={dkx_length, dky_length, PA_E_pixel};
		w_att_1d(rootG, "Delta", 3, &delta3[0]);
		
		int Size3[3]={kx_count_ext, ky_count_ext, num_points_E};
		w_att_1i(rootG, "Size", 3, &Size3[0]);
		
		w_att_1d(rootG, "Xvector", 3, &kx_vector[0]);
		w_att_1d(rootG, "Yvector", 3, &ky_vector[0]);
	}

	w_att_double(rootG, "dE", PA_dE);
	w_att_str(rootG, "Initial_state", string(PA_initial_state));
	w_att_str(rootG, "Final_state", string(PA_final_state));
	if(strcmp(PA_final_state, "Calc")==0){
		w_att_double(rootG, "Final_state_step", PA_final_state_step);
	}
	w_att_str(rootG, "Polarization", string(PA_polarization));
	w_att_double(rootG, "Theta", PA_theta);
	w_att_double(rootG, "Phi", PA_phi);

	w_att_bool(rootG, "Weighting", PA_weighting);

	s_data_1c(AtomG, "Labels", &size1, &length);
	Group atomG_out(rootG.createGroup("Atoms"));
	w_data_1c(atomG_out, "Labels", size1, length, (char**) atom_labels);

	s_data_1c(SpecG, "Labels", &size1, &length);
	w_data_1c(atomG_out, "Species", size1, length, (char**) atom_spec_label);

	w_data_2d(atomG_out, "Coordinates", atom_length, 3, (double**) atom_coordinates);
	if(PA_weighting==true){
		w_data_1d(atomG_out, "Weighting", atom_length, (double*) atom_weighting);
	}

	w_data_2d(atomG_out, "UnitCell", 3, 3, (double**) atom_cell);

	double**** FP_loc_edge_export_re; // [sp][ig][ie][ik]
  double**** FP_loc_edge_export_im; // [sp][ig][ie][ik]
	int** FP_bulk_count_export_up; // [ie][ik];
	int** FP_bulk_count_export_dn;
	double** FP_zgels_norm_export_up; // [ie][ik];
	double** FP_zgels_norm_export_dn;
	
	int EScale_count=E_max_scale-E_min_scale+1;
	int g_count=0;
	if(PA_FPFS){
		Group FPFSG(rootG.createGroup("FPFS"));
		for(int ik=0; ik<total_count_ext; ik++){
			for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
				if(g_count<final_states_FP_g_size[ik][ifp]){
					g_count=final_states_FP_g_size[ik][ifp];
				}
			}
		}
		FP_loc_edge_export_re=new double***[sp_max];
		FP_loc_edge_export_im=new double***[sp_max];
		if(!PA_FPFS_bulk_set && !PA_FPFS_Numerov){
			for(int sp=0; sp<sp_max; sp++){
				double* buffer_re=new double[g_count*EScale_count*total_count_ext];
				FP_loc_edge_export_re[sp]=new double**[g_count];
				double* buffer_im=new double[g_count*EScale_count*total_count_ext];
				FP_loc_edge_export_im[sp]=new double**[g_count];
				for(int ig=0; ig<g_count; ig++){
					FP_loc_edge_export_re[sp][ig]=new double*[EScale_count];
					FP_loc_edge_export_im[sp][ig]=new double*[EScale_count];
					for(int ie=0; ie<EScale_count; ie++){
						FP_loc_edge_export_re[sp][ig][ie]=&buffer_re[ig*EScale_count*total_count_ext+ie*total_count_ext];
						FP_loc_edge_export_im[sp][ig][ie]=&buffer_im[ig*EScale_count*total_count_ext+ie*total_count_ext];
						for(int ik=0; ik<total_count_ext; ik++){
							FP_loc_edge_export_re[sp][ig][ie][ik]=0.0;
							FP_loc_edge_export_im[sp][ig][ie][ik]=0.0;
						}
					}
				}
			}
			for(int ik=0; ik<total_count_ext; ik++){
				for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
					int ie=final_states_EScale[ik][ifp]-E_min_scale;
					int sp=final_states_spin[ik][ifp];
					for(int ig=0; ig<g_count; ig++){
						FP_loc_edge_export_re[sp][ig][ie][ik]=final_states_FP_loc_edge[ik][ifp][ig].real();
						FP_loc_edge_export_im[sp][ig][ie][ik]=final_states_FP_loc_edge[ik][ifp][ig].imag();
					}
				}
			}
		}
		if(PA_FPFS_bulk_set){
			FP_bulk_count_export_up=alloc_imatrix(EScale_count, total_count_ext);
			if(spin_i>0){
				FP_bulk_count_export_dn=alloc_imatrix(EScale_count, total_count_ext);
			}
			for(int ik=0; ik<total_count_ext; ik++){
				for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
					int ie=final_states_EScale[ik][ifp]-E_min_scale;
					int sp=final_states_spin[ik][ifp];
					if(sp==0){
						FP_bulk_count_export_up[ie][ik]=final_states_FP_bulk_count[ik][ifp];
					}else{
						FP_bulk_count_export_dn[ie][ik]=final_states_FP_bulk_count[ik][ifp];
					}
				}
			}
			w_data_2i(FPFSG, "FP_bulk_count_up", EScale_count, total_count_ext, (int**)&FP_bulk_count_export_up[0][0]);
			if(spin_i>0){
				w_data_2i(FPFSG, "FP_bulk_count_dn", EScale_count, total_count_ext, (int**)&FP_bulk_count_export_dn[0][0]);
			}
		}
		if(!PA_FPFS_Numerov){
			FP_zgels_norm_export_up=alloc_dmatrix(EScale_count, total_count_ext);
			if(spin_i>0){
				FP_zgels_norm_export_dn=alloc_dmatrix(EScale_count, total_count_ext);
			}
			for(int ik=0; ik<total_count_ext; ik++){
				for(int ifp=0; ifp<final_states_FP_size[ik]; ifp++){
					int ie=final_states_EScale[ik][ifp]-E_min_scale;
					int sp=final_states_spin[ik][ifp];
					if(sp==0){
						FP_zgels_norm_export_up[ie][ik]=final_states_zgels_norm[ik][ifp];
					}else{
						FP_zgels_norm_export_dn[ie][ik]=final_states_zgels_norm[ik][ifp];
					}
				}
			}
			w_data_2d(FPFSG, "FP_zgels_norm_up", EScale_count, total_count_ext, (double**)&FP_zgels_norm_export_up[0][0]);
			if(spin_i>0){
				w_data_2d(FPFSG, "FP_zgels_norm_dn", EScale_count, total_count_ext, (double**)&FP_zgels_norm_export_dn[0][0]);
			}
		}
	  
		// FPFS part
		w_att_int(FPFSG, "E_min_scale", E_min_scale);
		w_att_int(FPFSG, "E_max_scale", E_max_scale);
		w_att_double(FPFSG, "FPFS_energy_step", PA_FPFS_energy_step);
		if(!PA_FPFS_bulk_set && !PA_FPFS_Numerov){
			for(int sp=0; sp<sp_max; sp++){
				if(spin_i==1){
					sprintf(group_name, "Local_edge_%s_real", sp==0?"Up":"Dn");
				}else{
					sprintf(group_name, "Local_edge_real");
				}
				w_data_3d(FPFSG, group_name, g_count, EScale_count, total_count_ext, (double***)&FP_loc_edge_export_re[sp][0][0][0]);
				if(spin_i==1){
					sprintf(group_name, "Local_edge_%s_imag", sp==0?"Up":"Dn");
				}else{
					sprintf(group_name, "Local_edge_imag");
				}
				w_data_3d(FPFSG, group_name, g_count, EScale_count, total_count_ext, (double***)&FP_loc_edge_export_im[sp][0][0][0]);
			}
		}
	}

	if(PA_FPFS && PA_FPFS_file_set){
		H5File output2(PA_FPFS_file, H5F_ACC_TRUNC);	// att Datetime @root
		Group rootG(output2.openGroup("/"));
		time_t datetime_now=time(NULL);
		struct tm *timeptr=localtime(&datetime_now);
		char time_str[val_size+1];
		strftime(time_str, val_size, "%Y-%m-%d %H:%M:%S", timeptr);
		w_att_str(rootG, "Datetime", time_str);
		w_att_int(rootG, "Dimension", dimension);

		if(dimension==1){
			w_data_2d(rootG, "Dispersion", total_count_ext, num_points_E, (double**) dispersion2);
		}else{
			w_data_3d(rootG, "Dispersion", kx_count_ext, ky_count_ext, num_points_E, (double***) &dispersion3[0][0][0]);
		}

		if(dimension==1){
			double kx_offset=kx_length*kx_range[0];
			double E_offset=PA_E_min;
			if(PA_ext_set){
				kx_offset=kx_length*(kx_range[0]-(kx_range[1]-kx_range[0])*PA_ext_le);
			}
			double offset2[2]={kx_offset, E_offset};
			w_att_1d(rootG, "Offset", 2, &offset2[0]);

			double delta2[2]={dkx_length, PA_E_pixel};
			w_att_1d(rootG, "Delta", 2, &delta2[0]);
		
			int Size2[2]={total_count_ext, num_points_E};
			w_att_1i(rootG, "Size", 2, &Size2[0]);
		
			w_att_1d(rootG, "Xvector", 3, &kx_vector[0]);
		}else{
			double kx_offset=kx_length*kx_range[0];
			double ky_offset=ky_length*ky_range[0];
			double E_offset=PA_E_min;
			if(PA_ext_set){
				kx_offset=kx_length*(kx_range[0]-(kx_range[1]-kx_range[0])*PA_ext_le);
				ky_offset=ky_length*(ky_range[0]-(ky_range[1]-ky_range[0])*PA_ext_dn);
			}
			double offset3[3]={kx_offset, ky_offset, E_offset};
			w_att_1d(rootG, "Offset", 3, &offset3[0]);

			double delta3[3]={dkx_length, dky_length, PA_E_pixel};
			w_att_1d(rootG, "Delta", 3, &delta3[0]);
		
			int Size3[3]={kx_count_ext, ky_count_ext, num_points_E};
			w_att_1i(rootG, "Size", 3, &Size3[0]);
		
			w_att_1d(rootG, "Xvector", 3, &kx_vector[0]);
			w_att_1d(rootG, "Yvector", 3, &ky_vector[0]);
		}

		w_att_double(rootG, "dE", PA_dE);
		w_att_str(rootG, "Initial_state", string(PA_initial_state));
		w_att_str(rootG, "Final_state", string(PA_final_state));
		if(strcmp(PA_final_state, "Calc")==0){
			w_att_double(rootG, "Final_state_step", PA_final_state_step);
		}
		w_att_str(rootG, "Polarization", string(PA_polarization));
		w_att_double(rootG, "Theta", PA_theta);
		w_att_double(rootG, "Phi", PA_phi);

		w_att_bool(rootG, "Weighting", PA_weighting);

		s_data_1c(AtomG, "Labels", &size1, &length);
		Group atomG_out(rootG.createGroup("Atoms"));
		w_data_1c(atomG_out, "Labels", size1, length, (char**) atom_labels);

		s_data_1c(SpecG, "Labels", &size1, &length);
		w_data_1c(atomG_out, "Species", size1, length, (char**) atom_spec_label);

		w_data_2d(atomG_out, "Coordinates", atom_length, 3, (double**) atom_coordinates);
		if(PA_weighting==true){
			w_data_1d(atomG_out, "Weighting", atom_length, (double*) atom_weighting);
		}

		w_data_2d(atomG_out, "UnitCell", 3, 3, (double**) atom_cell);

		// FPFS part
		Group FPFSG(rootG.createGroup("FPFS"));
		w_att_int(FPFSG, "Spin_i", spin_i);
		w_att_double(FPFSG, "Excitation_energy", PA_excitation_energy);
		w_att_double(FPFSG, "FPFS_energy_step", PA_FPFS_energy_step);
		if(PA_FPFS_bulk_set){
			w_att_double(FPFSG, "bulk_min", FPFS_bulk_min);
			w_att_double(FPFSG, "bulk_max", FPFS_bulk_max);
			w_att_int(FPFSG, "bulk_count", FPFS_bulk_count);
			w_att_int(FPFSG, "g_bulk_count", final_states_FP_g_size_bulk);
			w_data_2d(FPFSG, "g_vector_bulk", final_states_FP_g_size_bulk, 3, (double**)&final_states_FP_g_vec_bulk[0][0]);
		}
		int atom_props[atom_length];
		for(i=0; i<atom_length; i++){
			if(atom_above_surface[i]){
				atom_props[i]=1;
			}else{
				int is=atom_spec_index[i];
				if(empty_atoms[is]){
					atom_props[i]=-1;
				}else{
					atom_props[i]=0;
				}
			}
		}
		w_att_1i(FPFSG, "Atom_props", atom_length, &atom_props[0]);
		// orbitals
		//for(is=0; is<atom_spec_length; is++){
			//sprintf(group_name, "%d_%s", is+1, atom_spec_label[is]);
			//Group FPFSG_is(FPFSG.createGroup(group_name));
			// w_att_1i(FPFSG_is, "l_list", num_orbits[is], &l_list[is][0]);
			// w_att_1d(FPFSG_is, "r", wfn_length[is], &wfn_r[is][0]);
			/*
				for(io=0; io<num_orbits[is]; io++){
				sprintf(group_name, "orbital_%d", io);
				w_data_1d(FPFSG_is, group_name, wfn_length[is], &wfn_phi[is][io][0]);
				}*/
		//}

		// potentials
		Group PotG(FPFSG.createGroup("Potential"));
		w_data_3d(PotG, "V0_Kohn_Sham", VKS_count[0], VKS_count[1], VKS_count[2], (double***)&VKS0[0][0][0]);
		if(spin_i>0){
			w_data_3d(PotG, "V1_Kohn_Sham", VKS_count[0], VKS_count[1], VKS_count[2], (double***)&VKS1[0][0][0]);
		}
		w_data_2d(PotG, "Vgg_vector", Vgg_count, 3, (double**)&Vgg_vector[0][0]);
		w_data_2d(PotG, "Vgg0_abs", Vgg_count, VKS_count[0], (double**)&Vgg0_abs[0][0]);
		if(spin_i>0){
			w_data_2d(PotG, "Vgg1_abs", Vgg_count, VKS_count[0], (double**)&Vgg1_abs[0][0]);
		}
		if(!PA_ignore_core && !PA_ignore_nonlocal){
			for(ia=0; ia<atom_length; ia++){
				is=atom_spec_index[ia];
				sprintf(group_name, "%d_%s", ia+1, atom_spec_label[is]);
				Group PotG_ia(PotG.createGroup(group_name));
				w_data_1d(PotG_ia, "r", VPS_r_length[is], &VPS_r[is][0]);
				w_data_1d(PotG_ia, "V0", VPS_r_length[is], &VKS0_r[ia][0]);
				w_data_1d(PotG_ia, "V0_stdev", VPS_r_length[is], &dVKS0_r[ia][0]);
				if(spin_i>0){
					w_data_1d(PotG_ia, "V1", VPS_r_length[is], &VKS1_r[ia][0]);
					w_data_1d(PotG_ia, "V1_stdev", VPS_r_length[is], &dVKS1_r[ia][0]);
				}
			}
		}
		if(PA_FPFS_bulk_set){
			w_data_2d(PotG, "Vgg0_bulk_z_abs", Vgg_count, z_count_bulk, (double**)&Vgg0_average_abs[0][0]);
			w_data_2d(PotG, "Vgg0_bulk_z_stdev", Vgg_count, z_count_bulk, (double**)&Vgg0_stdev[0][0]);
			if(spin_i>0){
				w_data_2d(PotG, "Vgg1_bulk_z_abs", Vgg_count, z_count_bulk, (double**)&Vgg1_average_abs[0][0]);
				w_data_2d(PotG, "Vgg1_bulk_z_stdev", Vgg_count, z_count_bulk, (double**)&Vgg1_stdev[0][0]);
			}

			w_data_2d(PotG, "Vgg_vector_bulk", Vgg_count_bulk, 3, (double**)&Vgg_vector_bulk[0][0]);
			w_data_1d(PotG, "Vgg0_bulk_abs", Vgg_count_bulk, &Vgg0_abs_bulk[0]);
			if(spin_i>0){
				w_data_1d(PotG, "Vgg1_bulk_abs", Vgg_count_bulk, &Vgg1_abs_bulk[0]);
			}
		}
		// bulk band dispersion
		if(PA_FPFS_bulk_set){
			w_data_3d(FPFSG, "bulk_dispersion_up", total_count_ext, PA_FPFS_bulk_kz_steps, final_states_FP_g_size_bulk, (double***)&FP_bulk_dispersion_up[0][0][0]);
			if(spin_i>0){
				w_data_3d(FPFSG, "bulk_dispersion_dn", total_count_ext, PA_FPFS_bulk_kz_steps, final_states_FP_g_size_bulk, (double***)&FP_bulk_dispersion_dn[0][0][0]);
			}
			w_data_1d(FPFSG, "bulk_dispersion_kz", PA_FPFS_bulk_kz_steps, FP_bulk_dispersion_kz);
			w_data_1d(FPFSG, "bulk_dispersion_kappaz", FP_bulk_kappaz_count, FP_bulk_dispersion_kappaz);
			w_data_1d(FPFSG, "bulk_dispersion_mkappaz", FP_bulk_kappaz_border, FP_bulk_dispersion_mkappaz);
		}

		for(i=0; i<total_count_ext; i++){
			if(dimension==2){
				int kx=i%kx_count_ext;
				int ky=(i-kx)/kx_count_ext;
				sprintf(group_name, "%d_%d", kx, ky);
			}else{
				sprintf(group_name, "%d", i);
			}
			Group FPFSG_k(FPFSG.createGroup(group_name));
			w_att_int(FPFSG_k, "Count", final_states_FP_size[i]);
			w_att_1i(FPFSG_k, "Spin", final_states_FP_size[i], &final_states_spin[i][0]);
			w_data_1i(FPFSG_k, "Energy_scale", final_states_FP_size[i], &final_states_EScale[i][0]);
			if(PA_FPFS_bulk_set && final_states_FP_size[i]>0){
				w_data_1i(FPFSG_k, "Dispersion_c_count_up",     FP_bulk_kappaz_count,  &FP_bulk_dispersion_c_count_up[i][0]);
				w_data_1i(FPFSG_k, "Dispersion_c_BZ_count_up",  FP_bulk_kappaz_border, &FP_bulk_dispersion_c_BZ_count_up[i][0]);
				w_data_1i(FPFSG_k, "Dispersion_mc_count_up",    FP_bulk_kappaz_border, &FP_bulk_dispersion_mc_count_up[i][0]);
				w_data_1i(FPFSG_k, "Dispersion_mc_BZ_count_up", FP_bulk_kappaz_border, &FP_bulk_dispersion_mc_BZ_count_up[i][0]);
				w_data_2i(FPFSG_k, "Connection_c_up",     FP_bulk_kappaz_count,  FP_bulk_dispersion_c_count_up[i][0],     (int**)&FP_bulk_connection_c_up[i][0][0]);
				w_data_2i(FPFSG_k, "Connection_c_BZ_up",  FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_up[i][0],  (int**)&FP_bulk_connection_c_BZ_up[i][0][0]);
				w_data_2i(FPFSG_k, "Connection_mc_up",    FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_up[i][0],    (int**)&FP_bulk_connection_mc_up[i][0][0]);
				w_data_2i(FPFSG_k, "Connection_mc_BZ_up", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_up[i][0], (int**)&FP_bulk_connection_mc_BZ_up[i][0][0]);
				double*** dispersion_c_up;
				dispersion_c_up=alloc_dcube(FP_bulk_kappaz_count, FP_bulk_dispersion_c_count_up[i][0], 2);
				for(int ikappaz=0; ikappaz<FP_bulk_kappaz_count; ikappaz++){
					for(int ib=0; ib<FP_bulk_dispersion_c_count_up[i][ikappaz]; ib++){
						dispersion_c_up[ikappaz][ib][0]=FP_bulk_dispersion_c_up[i][ikappaz][ib].real();
						dispersion_c_up[ikappaz][ib][1]=FP_bulk_dispersion_c_up[i][ikappaz][ib].imag();
					}
				}
				w_data_3d(FPFSG_k, "Dispersion_c_up", FP_bulk_kappaz_count, FP_bulk_dispersion_c_count_up[i][0], 2, (double***)&dispersion_c_up[0][0][0]);
				delete_dcube(dispersion_c_up);
				
				dispersion_c_up=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_up[i][0], 2);
				for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
					for(int ib=0; ib<FP_bulk_dispersion_c_BZ_count_up[i][ikappaz]; ib++){
						dispersion_c_up[ikappaz][ib][0]=FP_bulk_dispersion_c_BZ_up[i][ikappaz][ib].real();
						dispersion_c_up[ikappaz][ib][1]=FP_bulk_dispersion_c_BZ_up[i][ikappaz][ib].imag();
					}
				}
				w_data_3d(FPFSG_k, "Dispersion_c_BZ_up", FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_up[i][0], 2, (double***)&dispersion_c_up[0][0][0]);
				delete_dcube(dispersion_c_up);

				dispersion_c_up=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_up[i][0], 2);
				for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
					for(int ib=0; ib<FP_bulk_dispersion_mc_count_up[i][ikappaz]; ib++){
						dispersion_c_up[ikappaz][ib][0]=FP_bulk_dispersion_mc_up[i][ikappaz][ib].real();
						dispersion_c_up[ikappaz][ib][1]=FP_bulk_dispersion_mc_up[i][ikappaz][ib].imag();
					}
				}
				w_data_3d(FPFSG_k, "Dispersion_mc_up", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_up[i][0], 2, (double***)&dispersion_c_up[0][0][0]);
				delete_dcube(dispersion_c_up);
				
				dispersion_c_up=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_up[i][0], 2);
				for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
					for(int ib=0; ib<FP_bulk_dispersion_mc_BZ_count_up[i][ikappaz]; ib++){
						dispersion_c_up[ikappaz][ib][0]=FP_bulk_dispersion_mc_BZ_up[i][ikappaz][ib].real();
						dispersion_c_up[ikappaz][ib][1]=FP_bulk_dispersion_mc_BZ_up[i][ikappaz][ib].imag();
					}
				}
				w_data_3d(FPFSG_k, "Dispersion_mc_BZ_up", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_up[i][0], 2, (double***)&dispersion_c_up[0][0][0]);
				delete_dcube(dispersion_c_up);
				if(spin_i>0){
					w_data_1i(FPFSG_k, "Dispersion_c_count_dn",     FP_bulk_kappaz_count,  &FP_bulk_dispersion_c_count_dn[i][0]);
					w_data_1i(FPFSG_k, "Dispersion_c_BZ_count_dn",  FP_bulk_kappaz_border, &FP_bulk_dispersion_c_BZ_count_dn[i][0]);
					w_data_1i(FPFSG_k, "Dispersion_mc_count_dn",    FP_bulk_kappaz_border, &FP_bulk_dispersion_mc_count_dn[i][0]);
					w_data_1i(FPFSG_k, "Dispersion_mc_BZ_count_dn", FP_bulk_kappaz_border, &FP_bulk_dispersion_mc_BZ_count_dn[i][0]);
					w_data_2i(FPFSG_k, "Connection_c_dn",     FP_bulk_kappaz_count,  FP_bulk_dispersion_c_count_dn[i][0],     (int**)&FP_bulk_connection_c_dn[i][0][0]);
					w_data_2i(FPFSG_k, "Connection_c_BZ_dn",  FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_dn[i][0],  (int**)&FP_bulk_connection_c_BZ_dn[i][0][0]);
					w_data_2i(FPFSG_k, "Connection_mc_dn",    FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_dn[i][0],    (int**)&FP_bulk_connection_mc_dn[i][0][0]);
					w_data_2i(FPFSG_k, "Connection_mc_BZ_dn", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_dn[i][0], (int**)&FP_bulk_connection_mc_BZ_dn[i][0][0]);
					double*** dispersion_c_dn;
					dispersion_c_dn=alloc_dcube(FP_bulk_kappaz_count, FP_bulk_dispersion_c_count_dn[i][0], 2);
					for(int ikappaz=0; ikappaz<FP_bulk_kappaz_count; ikappaz++){
						for(int ib=0; ib<FP_bulk_dispersion_c_count_dn[i][ikappaz]; ib++){
							dispersion_c_dn[ikappaz][ib][0]=FP_bulk_dispersion_c_dn[i][ikappaz][ib].real();
							dispersion_c_dn[ikappaz][ib][1]=FP_bulk_dispersion_c_dn[i][ikappaz][ib].imag();
						}
					}
					w_data_3d(FPFSG_k, "Dispersion_c_dn", FP_bulk_kappaz_count, FP_bulk_dispersion_c_count_dn[i][0], 2, (double***)&dispersion_c_dn[0][0][0]);
					delete_dcube(dispersion_c_dn);
				
					dispersion_c_dn=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_dn[i][0], 2);
					for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
						for(int ib=0; ib<FP_bulk_dispersion_c_BZ_count_dn[i][ikappaz]; ib++){
							dispersion_c_dn[ikappaz][ib][0]=FP_bulk_dispersion_c_BZ_dn[i][ikappaz][ib].real();
							dispersion_c_dn[ikappaz][ib][1]=FP_bulk_dispersion_c_BZ_dn[i][ikappaz][ib].imag();
						}
					}
					w_data_3d(FPFSG_k, "Dispersion_c_BZ_dn", FP_bulk_kappaz_border, FP_bulk_dispersion_c_BZ_count_dn[i][0], 2, (double***)&dispersion_c_dn[0][0][0]);
					delete_dcube(dispersion_c_dn);

					dispersion_c_dn=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_dn[i][0], 2);
					for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
						for(int ib=0; ib<FP_bulk_dispersion_mc_count_dn[i][ikappaz]; ib++){
							dispersion_c_dn[ikappaz][ib][0]=FP_bulk_dispersion_mc_dn[i][ikappaz][ib].real();
							dispersion_c_dn[ikappaz][ib][1]=FP_bulk_dispersion_mc_dn[i][ikappaz][ib].imag();
						}
					}
					w_data_3d(FPFSG_k, "Dispersion_mc_dn", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_count_dn[i][0], 2, (double***)&dispersion_c_dn[0][0][0]);
					delete_dcube(dispersion_c_dn);
				
					dispersion_c_dn=alloc_dcube(FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_dn[i][0], 2);
					for(int ikappaz=0; ikappaz<FP_bulk_kappaz_border; ikappaz++){
						for(int ib=0; ib<FP_bulk_dispersion_mc_BZ_count_dn[i][ikappaz]; ib++){
							dispersion_c_dn[ikappaz][ib][0]=FP_bulk_dispersion_mc_BZ_dn[i][ikappaz][ib].real();
							dispersion_c_dn[ikappaz][ib][1]=FP_bulk_dispersion_mc_BZ_dn[i][ikappaz][ib].imag();
						}
					}
					w_data_3d(FPFSG_k, "Dispersion_mc_BZ_dn", FP_bulk_kappaz_border, FP_bulk_dispersion_mc_BZ_count_dn[i][0], 2, (double***)&dispersion_c_dn[0][0][0]);
					delete_dcube(dispersion_c_dn);
				}
			}
			for(j=0; j<final_states_FP_size[i]; j++){
				sprintf(group_name, "%d", j);
				Group FPFSG_ke(FPFSG_k.createGroup(group_name));
				w_att_1d(FPFSG_ke, "k", 3, &final_states_k[i][j][0]);
				int g_size=final_states_FP_g_size[i][j];
				w_att_int(FPFSG_ke, "g_size", g_size);
				if(g_size>0){
					w_data_2d(FPFSG_ke, "g_vector", g_size, 3, (double**)&final_states_FP_g_vec[i][j][0][0]);
					int z_size=VKS_count[0];
					double final_states_FP_loc_export[g_size][z_size][2];
					for(int ig=0; ig<g_size; ig++){
						for(int iz=0; iz<z_size; iz++){
							final_states_FP_loc_export[ig][iz][0]=final_states_FP_loc[i][j][ig][iz].real();
							final_states_FP_loc_export[ig][iz][1]=final_states_FP_loc[i][j][ig][iz].imag();
						}
					}
					w_data_3d(FPFSG_ke, "FP_local", g_size, z_size, 2, (double***)&final_states_FP_loc_export[0][0][0]);
				}
				if(PA_FPFS_bulk_set){
					int FP_size_bulk=final_states_FP_bulk_count[i][j];
					int g_size_bulk=final_states_FP_g_size_bulk;
					w_att_int(FPFSG_ke, "bulk_count", FP_size_bulk);
					if(FP_size_bulk>0){
						double final_states_FP_bulk_export[FP_size_bulk][g_size_bulk][2];
						for(int in=0; in<FP_size_bulk; in++){
							for(int ig=0; ig<g_size_bulk; ig++){
								final_states_FP_bulk_export[in][ig][0]=final_states_FP_bulk[i][j][in][ig].real();
								final_states_FP_bulk_export[in][ig][1]=final_states_FP_bulk[i][j][in][ig].imag();
							}
						}
						w_data_3d(FPFSG_ke, "FP_bulk", FP_size_bulk, g_size_bulk, 2, (double***)&final_states_FP_bulk_export[0][0][0]);
						double final_states_FP_bulk_coefs_export[FP_size_bulk];
						for(int in=0; in<FP_size_bulk; in++){
							final_states_FP_bulk_coefs_export[in]=abs(final_states_FP_bulk_coefs[i][j][in]);
						}
						w_data_1d(FPFSG_ke, "FP_bulk_coefs_abs", FP_size_bulk, &final_states_FP_bulk_coefs_export[0]);
						w_data_1d(FPFSG_ke, "FP_bulk_kz", FP_size_bulk, &final_states_FP_bulk_kz[i][j][0]);
						w_data_1d(FPFSG_ke, "FP_bulk_kappaz", FP_size_bulk, &final_states_FP_bulk_kappaz[i][j][0]);
					}
				}
				if(!PA_ignore_core && !PA_ignore_nonlocal){
					double*** FP_nonloc_export=alloc_dcube(5, 9, 2);
					for(int ia=0; ia<atom_length; ia++){
						int is=atom_spec_index[ia];
						if(empty_atoms[is]){
							continue;
						}
						if(!PA_calc_all_nonloc && !atom_weighting_flag[ia]){
							continue;
						}
						for(int il=0; il<5; il++){
							for(int mpl=0; mpl<=il*2; mpl++){
								FP_nonloc_export[il][mpl][0]=final_states_FP_nonloc[i][j][ia][il][mpl].real();
								FP_nonloc_export[il][mpl][1]=final_states_FP_nonloc[i][j][ia][il][mpl].imag();
							}
						}
						sprintf(group_name, "Nonloc_%d_%s", ia+1, atom_spec_label[is]);
						w_data_3d(FPFSG_ke, group_name, 5, 9, 2, (double***)&FP_nonloc_export[0][0][0]);
					}
					delete_dcube(FP_nonloc_export);
				}
				/*
					for(ia=0; ia<atom_length; ia++){
					int is=atom_spec_index[ia];
					sprintf(group_name, "%d_%s", (ia+1), atom_labels[ia]);
					Group FPFSG_kea(FPFSG_ke.createGroup(group_name));
					int digit=(spin_i==2)?2:1;
					
					for(io=0; io<num_orbits[is]; io++){
					int l=l_list[is][io];
					sprintf(group_name, "orbital_%d", io);
					double LCAO_component[2*l+1][digit*2];
					for(int mpl=0; mpl<2*l+1; mpl++){
					LCAO_component[mpl][0]=final_states_LCAO[i][j][ia][io][mpl][0].real();
					LCAO_component[mpl][1]=final_states_LCAO[i][j][ia][io][mpl][0].imag();
					if(digit==2){
					LCAO_component[mpl][2]=final_states_LCAO[i][j][ia][io][mpl][1].real();
					LCAO_component[mpl][3]=final_states_LCAO[i][j][ia][io][mpl][1].imag();
					}
					}
					w_data_2d(FPFSG_kea, group_name, 2*l+1, digit*2, (double**)&LCAO_component[0][0]);
					}
					}*/
			}
		}
		if(!PA_FPFS_bulk_set && !PA_FPFS_Numerov){
			for(int sp=0; sp<sp_max; sp++){
				if(spin_i==1){
					sprintf(group_name, "Local_edge_%s_real", sp==0?"Up":"Dn");
				}else{
					sprintf(group_name, "Local_edge_real");
				}
				w_data_3d(FPFSG, group_name, g_count, EScale_count, total_count_ext, (double***)&FP_loc_edge_export_re[sp][0][0][0]);
				
				if(spin_i==1){
					sprintf(group_name, "Local_edge_%s_imag", sp==0?"Up":"Dn");
				}else{
					sprintf(group_name, "Local_edge_imag");
				}
				w_data_3d(FPFSG, group_name, g_count, EScale_count, total_count_ext, (double***)&FP_loc_edge_export_im[sp][0][0][0]);
			}
		}
		
		if(!PA_ignore_core && !PA_ignore_nonlocal){
			Group NonlocG(FPFSG.createGroup("Nonlocal"));
			for(int ie=0; ie<EScale_count; ie++){
				sprintf(group_name, "%d", ie);
				Group NonlocG_e(NonlocG.createGroup(group_name));
				for(int sp=0; sp<sp_max; sp++){
					sprintf(group_name, "%s", sp==0 ? "Up" : "Dn");
					Group NonlocG_es(NonlocG_e.createGroup(group_name));
					for(int ia=0; ia<atom_length; ia++){
						int is=atom_spec_index[ia];
						if(empty_atoms[is]){
							continue;
						}
						if(!PA_calc_all_nonloc && !atom_weighting_flag[ia]){
							continue;
						}		
						sprintf(group_name, "%d_%s", ia+1, atom_spec_label[is]);
						w_data_2d(NonlocG_es, group_name, 5, wfn_cutoff_index[is], (double**)&final_states_FP_nonloc_r[ie][sp][ia][0][0]);
					}
				}
			}
			for(int is=0; is<atom_spec_length; is++){
				sprintf(group_name, "%s_r", atom_spec_label[is]);
				w_data_1d(NonlocG, group_name, wfn_cutoff_index[is], &wfn_r_use[is][0]);
			}
		}
	}
	
	return;
}
