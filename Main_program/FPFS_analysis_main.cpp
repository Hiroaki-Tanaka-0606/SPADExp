// First-principles final state analysis

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
#include "HDF5_tools.hpp"
#include "physical_tools.hpp"
#include "variables.hpp"

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

int main(int argc, char** argv){
	if(argc<2){
		printf("Usage: %s .hdf5\n", argv[0]);
		return -1;
	}
	// arguments
	char* file_name=argv[1];

	// constants
	int FPFS_range=10;
	int item_size=64;

	// buffers
	char group_name[item_size];

	printf("----Load %s----\n", file_name);

	H5File input(file_name, H5F_ACC_RDONLY);
	Group rootG(input.openGroup("/"));
	
	int dimension=r_att_int(rootG, "Dimension");
	printf("%32s = %d\n", "Dimension", dimension);

	int size[3]={0, 0, 0};
	r_att_1i(rootG, "Size", dimension+1, &size[0]);
	printf("%32s = [%d, %d, %d]\n", "Size", size[0], size[1], size[2]);
	
	int total_count=size[0];
	if(dimension==2){
		total_count*=size[1];
	}
	printf("%32s = %d\n", "Total count", total_count);

	Group AtomG(rootG.openGroup("Atoms"));
	int atom_length;
	int label_length;
	s_data_1c(AtomG, "Labels", &atom_length, &label_length);
	char atom_labels[atom_length][label_length];
	r_data_1c(AtomG, "Labels", atom_length, label_length, (char**)atom_labels);

	int atom_spec_length;
	s_data_1c(AtomG, "Species", &atom_spec_length, &label_length);
	char atom_spec_label[atom_spec_length][label_length];
	r_data_1c(AtomG, "Species", atom_spec_length, label_length, (char**)atom_spec_label);

	int atom_spec_index[atom_length];
	for(int ia=0; ia<atom_length; ia++){
		for(int is_search=0; is_search<atom_spec_length; is_search++){
			if(strcmp(atom_labels[ia], atom_spec_label[is_search])==0){
				atom_spec_index[ia]=is_search;
				break;
			}	
		}
	}
	
	double atom_coordinates[atom_length][3];
	r_data_2d(AtomG, "Coordinates", atom_length, 3, (double**)atom_coordinates);

	Group FPFSG(rootG.openGroup("FPFS"));
	int atom_props[atom_length];
	r_att_1i(FPFSG, "Atom_props", atom_length, (int*)atom_props);
	
	for(int i=0; i<atom_length; i++){
		printf("%32s (%9.4f, %9.4f, %9.4f) [a.u.]: %s\n", atom_labels[i], atom_coordinates[i][0], atom_coordinates[i][1], atom_coordinates[i][2], atom_props[i]==-1 ? "Below the slab": (atom_props[i]==0 ? "In the slab": "Above the slab"));
	}
	
	double unit_cell[3][3];
	r_data_2d(AtomG, "UnitCell", 3, 3, (double**)unit_cell);
	for(int i=0; i<3; i++){
		printf("%31s%1d = (%9.4f, %9.4f, %9.4f) [a.u.]\n", "a", (i+1), unit_cell[i][0], unit_cell[i][1], unit_cell[i][2]);
	}

	double rec_cell[3][3];
	double acrossb[3];
	outer_product(unit_cell[1], unit_cell[2], acrossb);
	double abc=inner_product(acrossb, unit_cell[0]);
	for(int i=0; i<3; i++){
		rec_cell[0][i]=2*M_PI*acrossb[i]/abc;
	}
	outer_product(unit_cell[2], unit_cell[0], acrossb);
	for(int i=0; i<3; i++){
		rec_cell[1][i]=2*M_PI*acrossb[i]/abc;
	}
	outer_product(unit_cell[0], unit_cell[1], acrossb);
	for(int i=0; i<3; i++){
		rec_cell[2][i]=2*M_PI*acrossb[i]/abc;
	}

	for(int i=0; i<3; i++){
		printf("%31s%1d = (%9.4f, %9.4f, %9.4f) [/a.u.]\n", "b", (i+1), rec_cell[i][0], rec_cell[i][1], rec_cell[i][2]);
	}

	int spin_i=r_att_int(FPFSG, "Spin_i");
	printf("%32s = %d\n", "Spin degree of freedom", spin_i);
	int digit=spin_i==2?2:1;
	
	// orbitals
	int wfn_length[atom_spec_length];
	double** wfn_phi[atom_spec_length]; // [is][io][r]
	double* wfn_r[atom_spec_length];
	int* l_list[atom_spec_length];
	int num_orbits[atom_spec_length];
	for(int is=0; is<atom_spec_length; is++){
		sprintf(group_name, "%d_%s", is+1, atom_spec_label[is]);
		Group FPFSG_is(FPFSG.openGroup(group_name));
		// r
		s_att_1d(FPFSG_is, "r", &wfn_length[is]);
		wfn_r[is]=new double[wfn_length[is]];
		r_att_1d(FPFSG_is, "r", wfn_length[is], &wfn_r[is][0]);
		// l_list
		s_att_1i(FPFSG_is, "l_list", &num_orbits[is]);
		l_list[is]=new int[num_orbits[is]];
		r_att_1i(FPFSG_is, "l_list", num_orbits[is], &l_list[is][0]);
		// wfn_phi
		wfn_phi[is]=new double*[num_orbits[is]];
		for(int io=0; io<num_orbits[is]; io++){
			sprintf(group_name, "orbital_%d", io);
			wfn_phi[is][io]=new double[wfn_length[is]];
			r_data_1d(FPFSG_is, group_name, wfn_length[is], &wfn_phi[is][io][0]);
		}
	}
	
	complex<double> m1jlp[5]={complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1), complex<double>(1, 0)};
	complex<double> Ylm_k[5][9];
	for(int ik=0; ik<total_count; ik++){
		printf("k index = %d\n", ik);
		if(dimension==1){
			sprintf(group_name, "%d", ik);
		}else{
			int ikx=ik%size[0];
			int iky=(ik-ikx)/size[0];
			sprintf(group_name, "%d_%d", ikx, iky);
		}
		Group FPFSG_k(FPFSG.openGroup(group_name));
		int count=r_att_int(FPFSG_k, "Count");
		printf("  Number of final states = %d\n", count);
		int spin[count];
		r_att_1i(FPFSG_k, "Spin", count, &spin[0]);
		for(int j=0; j<count; j++){
			printf("  state index = %d\n", j);
			sprintf(group_name, "%d", j);
			Group FPFSG_ke(FPFSG_k.openGroup(group_name));
			double k[3];
			r_att_1d(FPFSG_ke, "k", 3, &k[0]);
			printf("    k vector = [%7.3f, %7.3f, %7.3f]\n", k[0], k[1], k[2]);
			printf("    spin = %s\n", spin[j]==0 ? "Up": "Dn");
			printf("    t kx ky kz  fourier components\n");
			double k_t[3]={k[0], k[1], k[2]};			
			for(double t=1.1; t>=-1.1; t-=0.1){
				k_t[0]=t*k[0];
				spherical_harmonics(k_t, &Ylm_k[0][0]);
				complex<double> final_state_fourier[digit];
				for(int id=0; id<digit; id++){
					final_state_fourier[id]=complex<double>(0, 0);
				}
				double k_t_length=sqrt(inner_product(k_t, k_t));
				// radial part calculations
				double** final_state_integral=new double*[atom_spec_length]; // [is][io]
				for(int is=0; is<atom_spec_length; is++){
					double* wfn_spBessel_rdr=new double[wfn_length[is]];
					final_state_integral[is]=new double[num_orbits[is]];
					for(int io=0; io<num_orbits[is]; io++){
						int final_state_l=l_list[is][io];
						for(int ir=0; ir<wfn_length[is]-1; ir++){
							wfn_spBessel_rdr[ir]=wfn_r[is][ir]*sp_bessel(final_state_l, wfn_r[is][ir]*k_t_length)*(wfn_r[is][ir+1]-wfn_r[is][ir]);
						}
						wfn_spBessel_rdr[wfn_length[is]-1]=0;
						final_state_integral[is][io]=ddot(&wfn_length[is], &wfn_spBessel_rdr[0], &wfn_phi[is][io][0]);
						// printf("is=%d, io=%d, mat=%e\n", is, io, final_state_integral[is][io]);
					}
				}
				// multiply LCAO
				for(int ia=0; ia<atom_length; ia++){
					if(atom_props[ia]!=1){
						continue;
					}
					int is=atom_spec_index[ia];
					double kt=inner_product(k_t, atom_coordinates[ia]);
					complex<double> atom_phase(cos(kt), -sin(kt));
					
					sprintf(group_name, "%d_%s", (ia+1), atom_labels[ia]);
					Group FPFSG_kea(FPFSG_ke.openGroup(group_name));
					for(int io=0; io<num_orbits[is]; io++){
						int l=l_list[is][io];
						complex<double> coeff=4*M_PI*atom_phase*m1jlp[l];
						double LCAO_component[2*l+1][digit*2];
						sprintf(group_name, "orbital_%d", io);
						r_data_2d(FPFSG_kea, group_name, 2*l+1, digit*2, (double**)&LCAO_component[0][0]);
						for(int mpl=0; mpl<2*l+1; mpl++){
							for(int id=0; id<digit; id++){
								complex<double> LCAO_element(LCAO_component[mpl][id*2], LCAO_component[mpl][id*2+1]);
								final_state_fourier[id]+=coeff*Ylm_k[l][mpl]*LCAO_element*final_state_integral[is][io];
							}
						}
					}
				}
				// export
				printf("%5.2f %7.3f %7.3f %7.3f  ", t, k_t[0], k_t[1], k_t[2]);
				for(int id=0; id<digit; id++){
					printf("%7.3f %7.3f  ", final_state_fourier[id].real(), final_state_fourier[id].imag());
				}
				printf("\n");
			}
		}
	}
	/*
	
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
				delete LCAO_raw;
				/*
					if(ia==0){
					for(int iii=0; iii<size3; iii++){
					for(int jjj=0; jjj<size4/2; jjj++){
					cout << LCAO[ia][io][0][0][iii][jjj];
					}
					cout << endl;
					}
					}
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
						}
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
		}
	
	
	// 2. load initial state
	write_log((char*)"----Load initial states----");

	double*** wfn_phi_rdr=new double**[atom_spec_length]; // [is][io][r]
	double*** wfn_phi=new double**[atom_spec_length]; //[is][io][r] for FPFS
	double** wfn_r=new double*[atom_spec_length]; // [is][r]
	char AO_group_name[item_size*2];
	int wfn_length[atom_spec_length];
	int Z[atom_spec_length];
	bool empty_atoms[atom_spec_length];

	string AO_file_name(PA_AO_file);
	H5File AO_file(AO_file_name, H5F_ACC_RDONLY);
	
	for(is=0; is<atom_spec_length; is++){
		wfn_phi_rdr[is]=new double*[num_orbits[is]];
		if(PA_FPFS){
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
		wfn_length[is]=r_att_int(AOG, "length");
		Z[is]=r_att_int(AOG, "Z");
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
			if(PA_FPFS){
				wfn_phi[is][io]=new double[wfn_length[is]];
			}
			for(j=0; j<wfn_length[is]-1; j++){
				if(empty_atom==true || strcmp(PA_initial_state, "PAO")==0){
					wfn_phi_rdr[is][io][j]=wfn_both[0][j]*wfn_r[is][j]*(wfn_r[is][j+1]-wfn_r[is][j]);
				}else{
					wfn_phi_rdr[is][io][j]=wfn_both[1][j]*wfn_r[is][j]*(wfn_r[is][j+1]-wfn_r[is][j]);
				}
				if(PA_FPFS){
					if(empty_atom==true || strcmp(PA_final_state, "FP_PAO")==0){
						wfn_phi[is][io][j]=wfn_both[0][j];
					}else{
						wfn_phi[is][io][j]=wfn_both[1][j];
					}
				}
			}
			wfn_phi_rdr[is][io][wfn_length[is]-1]=0;
			wfn_phi[is][io][wfn_length[is]-1]=0;
			/*for(j=0; j<wfn_length[is]; j++){
				printf("%f %f\n", wfn_r[is][j], wfn_phi_rdr[is][io][j]);
				}
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
	
	double atom_weighting[atom_length];
	bool atom_weighting_flag[atom_length];
	bool atom_above_surface[atom_length];

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
		}
	
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
				}
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
	

	/// 3-3. tail profile
	int tail_index=floor(PA_dE*PA_sigma_max/PA_E_pixel);
	double tail[tail_index+1];
	for(i=0; i<=tail_index; i++){
		tail[i]=1.0/(PA_dE*sqrt(2*M_PI))*exp(-1.0/2.0*(i*PA_E_pixel/PA_dE)*(i*PA_E_pixel/PA_dE));
		// printf("%d %f\n", i, tail[i]);
	}


	// 3-4. final states
	int k_index_min=-1;
	int k_index_max=-1;
	int final_states_count;
	double**** final_states_calc; // [k_index][atom_spec][l][r] for Calc
	double*** final_states_phase; // [k_index][atom_spec][l] for Calc
	complex<double>****** final_states_FP;  // [total_count_ext][FPIndex][ia][l][m+l][r] for FP_PAO and FP_AO
	int** final_states_EScale;    // [total_count_ext][FPIndex]=round(E/FPFS_energy_step)
	int** final_states_spin;      // [total_count_ext][FPIndex]=up(0) or down(1), always zero for spin_i==0 (off)
	double** final_states_EWidth; // [total_count_ext][FPIndex]=energy difference between final state bases
	int* final_states_FP_size;    // [total_count_ext]=FPIndex_size
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
				double* final_states_alloc=new double[5*wfn_length[j]];
				final_states_calc[i][j]=new double*[5];
				for(k=0; k<5; k++){
					final_states_calc[i][j][k]=&final_states_alloc[k*wfn_length[j]];
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
					// interpolate final states so that they correspond to wfn_r[j]
					int x1, x2;
					for(x1=0; x1<wfn_length[j]; x1++){
						bool r_found=false;
						for(x2=0; x2<x_count-1; x2++){
							if(mu*x_coordinates[x2] < wfn_r[j][x1] && wfn_r[j][x1] < mu*x_coordinates[x2+1]){
								double dr1=wfn_r[j][x1]-mu*x_coordinates[x2];
								double dr2=mu*x_coordinates[x2+1]-wfn_r[j][x1];
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
			}
				
		// return;
	}else if(PA_FPFS){
		write_log((char*)"----Calculate First-principles final states----");
		sprintf(sprintf_buffer, "Fermi level = %.2f Eh = %.3f eV", EF_Eh, EF_Eh*Eh);
		write_log(sprintf_buffer);
		double E_min_tail=PA_E_min-PA_dE*PA_sigma_max;
		double E_max_tail=PA_E_max+PA_dE*PA_sigma_max;
		int E_min_scale=round(E_min_tail/PA_FPFS_energy_step);
		int E_max_scale=round(E_max_tail/PA_FPFS_energy_step);
		int scale_width=E_max_scale-E_min_scale+1;
		final_states_FP=new complex<double>*****[total_count_ext];
		final_states_EScale=new int*[total_count_ext];
		final_states_spin=new int*[total_count_ext];
		final_states_EWidth=new double*[total_count_ext];
		final_states_FP_size=new int[total_count_ext];
		
#pragma omp parallel firstprivate(scale_width, PA_ext_set, spin_i, num_bands, EF_Eh, Eh, PA_FPFS_energy_step, E_min_scale, E_max_scale, PA_excitation_energy, atom_above_surface, atom_length) private(j)
#pragma omp for
		
		for(i=0; i<total_count_ext; i++){
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
			final_states_FP[i]=new complex<double>****[FPIndex_size];
			for(j=0; j<FPIndex_size; j++){
				final_states_FP[i][j]=new complex<double>***[atom_length];
				for(int ia=0; ia<atom_length; ia++){
					int is=atom_spec_index[ia];
					final_states_FP[i][j][ia]=new complex<double>**[5];
					for(int l=0; l<5; l++){
						final_states_FP[i][j][ia][l]=new complex<double>*[2*l+1];
						for(int m=0; m<2*l+1; m++){
							final_states_FP[i][j][ia][l][m]=new complex<double>[wfn_length[is]];
							for(int ir=0; ir<wfn_length[is]; ir++){
								final_states_FP[i][j][ia][l][m][ir]=complex<double>(0, 0);
							}
						}
					}
				}
			}
			final_states_EScale[i]=new int[FPIndex_size];
			final_states_spin[i]=new int[FPIndex_size];
			final_states_EWidth[i]=new double[FPIndex_size];
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
			// for debug
			printf("k=%d, count=%d\n", i, FPIndex_size);
			for(j=0; j<FPIndex_size; j++){
				printf("Scale=%d, Energy=%f eV, spin=%d\n",
							 final_states_EScale[i][j], final_states_EScale[i][j]*PA_FPFS_energy_step, final_states_spin[i][j]);
							 }
			// obtain plane wave by linear combination
			for(j=0; j<FPIndex_size; j++){
				int eigen_scale=final_states_EScale[i][j];
				double kinetic_energy_Eh=(eigen_scale*PA_FPFS_energy_step+PA_excitation_energy)/Eh+EF_Eh;
				int eigen_spin=final_states_spin[i][j];
				// printf("Scale=%d, Ekin=%.2f, Spin=%d\n", eigen_scale, kinetic_energy_Eh, eigen_spin);
				double k_length_au=sqrt(2*kinetic_energy_Eh);
				double k_au[3]={0, k_point[1], k_point[2]};
				double kz_square=k_length_au*k_length_au-inner_product(k_au, k_au);
				if(kz_square<0){
					write_log((char*)"Error: negative kz^2");
					continue;
				}
				k_au[0]=sqrt(kz_square);
				// printf("k=(%.2f, %.2f, %.2f)\n", k_au[0], k_au[1], k_au[2]);
				double k_test_au[2*PA_FPFS_range+1][2*PA_FPFS_range+1][3];
				int k_count=0;
				for(int n1=-PA_FPFS_range; n1<=PA_FPFS_range; n1++){
					int n1p=n1+PA_FPFS_range;
					for(int n2=-PA_FPFS_range; n2<=PA_FPFS_range; n2++){
						int n2p=n2+PA_FPFS_range;
						k_test_au[n1p][n2p][0]=0.0;
						for(int p=1; p<=2; p++){
							k_test_au[n1p][n2p][p]=k_au[p]+rec_cell[1][p]*n1+rec_cell[2][p]*n2;
						}
						double kz_square=k_length_au*k_length_au-inner_product(k_test_au[n1p][n2p], k_test_au[n1p][n2p]);
						if(kz_square>0){
							k_test_au[n1p][n2p][0]=sqrt(kz_square);
							k_count++;
						}else{
							k_test_au[n1p][n2p][0]=-1;
						}
						/*
						if(k_test_au[n1p][n2p][0]>0){
							printf("k[%d][%d]=(%.2f, %.2f, %.2f)\n", n1, n2, k_test_au[n1p][n2p][0], k_test_au[n1p][n2p][1], k_test_au[n1p][n2p][2]);
							}
					}
				}
				// +- pm
				k_count*=2;
				// spin
				if(spin_i==2){
					k_count*=2;
				}
				// composite list
				// (0, 0, +, s), (0, 0, -, s), ..., (n1, n2, +, s), (n1, n2, -, s), {(n1, n2, +-, -s)}
				int n_list[k_count][2];
				double k_list[k_count][3];
				int spin_list[k_count];
				int list_index=0;
				int main_spin=final_states_spin[i][j];
				int sub_spin=main_spin==0 ? 1 : 0;
				n_list[list_index][0]=0; n_list[list_index][1]=0;
				k_list[list_index][0]=k_au[0], k_list[list_index][1]=k_au[1], k_list[list_index][2]=k_au[2];
				spin_list[list_index]=main_spin;
				list_index++;
				n_list[list_index][0]=0; n_list[list_index][1]=0;
				k_list[list_index][0]=-k_au[0], k_list[list_index][1]=k_au[1], k_list[list_index][2]=k_au[2];
				spin_list[list_index]=main_spin;
				list_index++;
				for(int n1=-PA_FPFS_range; n1<=PA_FPFS_range; n1++){
					int n1p=n1+PA_FPFS_range;
					for(int n2=-PA_FPFS_range; n2<=PA_FPFS_range; n2++){
						int n2p=n2+PA_FPFS_range;
						// (n1, n2)=(0, 0) is skipped because it is the first element
						if(k_test_au[n1p][n2p][0]>0 && (n1!=0 || n2!=0)){
							n_list[list_index][0]=n1; n_list[list_index][1]=n2;
							k_list[list_index][0]=k_test_au[n1p][n2p][0];
							k_list[list_index][1]=k_test_au[n1p][n2p][1];
							k_list[list_index][2]=k_test_au[n1p][n2p][2];
							spin_list[list_index]=main_spin;
							list_index++;
							n_list[list_index][0]=n1; n_list[list_index][1]=n2;
							k_list[list_index][0]=-k_test_au[n1p][n2p][0];
							k_list[list_index][1]=k_test_au[n1p][n2p][1];
							k_list[list_index][2]=k_test_au[n1p][n2p][2];
							spin_list[list_index]=main_spin;
							list_index++;
						}
					}
				}
				if(spin_i==2){
					for(int kp=0; kp<k_count/2; kp++){
						n_list[k_count/2+kp][0]=n_list[kp][0];
						n_list[k_count/2+kp][1]=n_list[kp][1];
						k_list[k_count/2+kp][0]=k_list[kp][0];
						k_list[k_count/2+kp][1]=k_list[kp][1];
						k_list[k_count/2+kp][2]=k_list[kp][2];
						spin_list[k_count/2+kp]=sub_spin;
					}
				}/*
				for(int list_index=0; list_index<k_count; list_index++){
					printf("n=(%2d, %2d), k=(%5.2f, %5.2f, %5.2f), s=%d\n",
								 n_list[list_index][0], n_list[list_index][1],
								 k_list[list_index][0], k_list[list_index][1], k_list[list_index][2],
								 spin_list[list_index]);
								 }
				// composite matrix
				complex<double> final_state_matrix[k_count][k_count];
				for(int im=0; im<k_count; im++){
					for(int ib=0; ib<k_count; ib++){
						final_state_matrix[im][ib]=complex<double>(0, 0);
					}
				}
				// calculate the integral part
				double** final_state_integral=new double*[atom_spec_length]; //[is][io]
				for(int is=0; is<atom_spec_length; is++){
					double* wfn_spBessel_rdr=new double[wfn_length[is]];
					final_state_integral[is]=new double[num_orbits[is]];
					for(int io=0; io<num_orbits[is]; io++){
						int final_state_l=l_list[is][io];
						for(int ir=0; ir<wfn_length[is]-1; ir++){
							wfn_spBessel_rdr[ir]=wfn_r[is][ir]*sp_bessel(final_state_l, wfn_r[is][ir]*k_length_au)*(wfn_r[is][ir+1]-wfn_r[is][ir]);
						}
						wfn_spBessel_rdr[wfn_length[is]-1]=0;
						final_state_integral[is][io]=ddot(&wfn_length[is], &wfn_spBessel_rdr[0], &wfn_phi[is][io][0]);
						// printf("is=%d, l=%d, mat=%e\n", is, final_state_l, final_state_integral[is][io]);
					}
				}
				// determine the band index
				int Ekin_band_index=-1;
				int final_state_band_index[k_count];
				for(int ib=0; ib<num_bands-1; ib++){
					double eigen_Eh1;
					double eigen_Eh2;
					if(spin_i==0 || spin_i==2){
						eigen_Eh1=band[k_index][ib];
						eigen_Eh2=band[k_index][ib+1];
					}else if(final_states_spin[i][j]==0){
						eigen_Eh1=band_up[k_index][ib];
						eigen_Eh2=band_up[k_index][ib+1];
					}else{
						eigen_Eh1=band_dn[k_index][ib];
						eigen_Eh2=band_dn[k_index][ib];
					}
					if(eigen_Eh1<kinetic_energy_Eh && kinetic_energy_Eh<eigen_Eh2){
						Ekin_band_index=ib;
						break;
					}
				}
				if(Ekin_band_index<0){
					write_log((char*)"Error: excited states not found");
					continue;
				}
				if(Ekin_band_index+k_count/2>=num_bands){
					write_log((char*)"Error: excited states are not enough");
					continue;
				}
				double Ekin_min=0.0;
				double Ekin_max=0.0;
				for(int ib=Ekin_band_index-k_count/2+1; ib<=Ekin_band_index+k_count/2; ib++){
					int ibp=ib-(Ekin_band_index-k_count/2+1);
					final_state_band_index[ibp]=ib;
					double eigen;
					if(spin_i==0 || spin_i==2){
						eigen=band[k_index][ib];
					}else if(final_states_spin[i][j]==0){
						eigen=band_up[k_index][ib];
					}else{
						eigen=band_dn[k_index][ib];
					}
					if(ibp==0){
						Ekin_min=eigen;
					}
					if(ibp==k_count-1){
						Ekin_max=eigen;
					}
					//printf("Band #%d: eigen=%.2f\n", ib, eigen);
				}
				final_states_EWidth[i][j]=(Ekin_max-Ekin_min)*Eh;
				sprintf(sprintf_buffer, "Energy width: %.2f Eh = %.1f eV", Ekin_max-Ekin_min, (Ekin_max-Ekin_min)*Eh);
				write_log(sprintf_buffer);
				// matrix element calculation
				complex<double> Ylm_k[5][9];
				complex<double> m1jlp[5]={complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1), complex<double>(1, 0)};
				for(int im=0; im<k_count; im++){
					double* k_im=k_list[im];
					int spin_im=spin_list[im];
					int digit=0;
					if(spin_i==2){
						digit=spin_im;
					}
					spherical_harmonics(k_im, &Ylm_k[0][0]);
					for(int ibi=0; ibi<k_count; ibi++){
						int ib=final_state_band_index[ibi];
						for(int ia=0; ia<atom_length; ia++){
							if(atom_above_surface[ia]==false){
								continue;
							}
							double kt=inner_product(k_im, atom_coordinates[ia]);
							complex<double> atom_phase(cos(kt), -sin(kt));
							int is=atom_spec_index[ia];
							for(int io=0; io<num_orbits[is]; io++){
								int io2=io;
								if(spin_i==1){
									if(spin_im==0){
										io2=io*2;
									}else{
										io2=io*2+1;
									}
								}
								int l=l_list[is][io];
								complex<double> coeff=4*M_PI*atom_phase*m1jlp[l];
								for(int mpl=0; mpl<2*l+1; mpl++){
									final_state_matrix[im][ibi]+=coeff*Ylm_k[l][mpl]*LCAO[ia][io2][k_index][ib][mpl][digit]*final_state_integral[is][io];
								}
							}
						}
					}
				}/*
				// matrix print for debug
				for(int im=0; im<k_count; im++){
					for(int ibi=0; ibi<k_count; ibi++){
						printf("(%6.3f, %6.3f) ", final_state_matrix[im][ibi].real(), final_state_matrix[im][ibi].imag());
					}
					printf("\n");
					}
				// index order is inverted due to the lapack implementation
				complex<double> fsm_inv[k_count][k_count];
				for(int im=0; im<k_count; im++){
					for(int ibi=0; ibi<k_count; ibi++){
						fsm_inv[ibi][im]=final_state_matrix[im][ibi];
					}
				}
				// solve
				complex<double> LCVector[k_count];
				if(zgesv_1_0(&k_count, &fsm_inv[0][0], LCVector)){
					/*
					write_log((char*)"zgesv succeeded");
					// check
					for(int im=0; im<k_count; im++){
						complex<double> MV=complex<double>(0, 0);
						for(int ibi=0; ibi<k_count; ibi++){
							MV+=final_state_matrix[im][ibi]*LCVector[ibi];
						}
						printf("MV[%2d]=(%6.3f, %6.3f)\n", im, MV.real(), MV.imag());
						}
				}else{
					write_log((char*)"zgesv failed");
					continue;
				}
				// composite final_states_FP
				// [total_count_ext][FPIndex][ia][l][m+l][r] for FP_PAO and FP_AO
				for(int ibi=0; ibi<k_count; ibi++){
					int digit=0;
					if(spin_i==2){
						digit=final_states_spin[i][j];
					}
					int ib=final_state_band_index[ibi];
					for(int ia=0; ia<atom_length; ia++){
						int is=atom_spec_index[ia];
						for(int io=0; io<num_orbits[is]; io++){
							int io2=io;						 
							if(spin_i==1){
								if(final_states_spin[i][j]==0){
									io2=io*2;
								}else{
									io2=io*2+1;
								}
							}
							int l=l_list[is][io];
							for(int mpl=0; mpl<2*l+1; mpl++){
								complex<double> coeff=LCVector[ibi]*LCAO[ia][io2][k_index][ib][mpl][digit];
								for(int ir=0; ir<wfn_length[is]; ir++){
									final_states_FP[i][j][ia][l][mpl][ir]+=coeff*wfn_phi[is][io][ir];
								}
							}
						}
					}
				}
			} // for(j=0; j<FPIndex_size; j++)
		} // omp for(i=0; i<total_count_ext; i++)
	} // if(FPFS)
  
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
	complex<double> m1jlp[5]={complex<double>(1, 0), complex<double>(0, -1), complex<double>(-1, 0), complex<double>(0, 1), complex<double>(1, 0)};
	complex<double> Ylm_k[5][9];
	complex<double> Ylm_k2[5][9];
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
#pragma omp parallel firstprivate(Y_coeff, sp_max, num_bands, EF_Eh, Eh, PA_E_min, PA_E_pixel, num_points_E, tail_index,  m1jlp, Gaunt_arr, atom_length, atom_spec_index, num_orbits, spin_i,  atom_coordinates, PA_ext_set, PA_final_state_step, k_index_min, atom_weighting_flag, atom_weighting, PA_weighting, axis_au, PA_reflection, PA_reflection_coef, PA_FPFS, PA_FPFS_energy_step) private(ib, sp, Ylm_k, Ylm_k2, ia, is, io, il, ir, j)
#pragma omp for
	for(ik=0; ik<total_count_ext; ik++){
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
			}
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
				// find appropriate FPIndex
				int FPIndex_1=-1; // for all types of spin_i
				int FPIndex_2=-1; // for spin_i==2
				int FPIndex_size=final_states_FP_size[ik];
				for(j=0; j<FPIndex_size; j++){
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
				complex<double> PAD_1(0, 0);
				complex<double> PAD_2(0, 0);
				if(strcmp(PA_output_data, "Band")==0){
					PAD_1=complex<double>(1, 0);
					PAD_2=complex<double>(1, 0);
				}else{
					for(ia=0; ia<atom_length; ia++){
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
							double final_states[5][wfn_length[is]];
							for(il=0; il<5; il++){
								for(ir=0; ir<wfn_length[is]; ir++){
									final_states[il][ir]=0.0;
								}
							}
							int k_index;
							if(strcmp(PA_final_state, "PW")==0 && k_length>1e-5){
								for(il=0; il<5; il++){
									for(ir=0; ir<wfn_length[is]; ir++){
										final_states[il][ir]=wfn_r[is][ir]*sp_bessel(il, wfn_r[is][ir]*k_length);
									}
								}
							}else if(strcmp(PA_final_state, "Calc")==0){
								k_index=round(k_length/PA_final_state_step);
								for(il=0; il<5; il++){
									for(ir=0; ir<wfn_length[is]; ir++){
										final_states[il][ir]=final_states_calc[k_index-k_index_min][is][il][ir];
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
								}
							}
						}else{
							// FPFS==true 
							is=atom_spec_index[ia];
							int num_orbits2=num_orbits[is];
							if(spin_i==1){
								num_orbits2*=2;
							}
							double final_states_re[wfn_length[is]];
							double final_states_im[wfn_length[is]];
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
									if(lp<0 || lp>=5){
										continue;
									}
									int mpl;
									for(mpl=0; mpl<=2*l; mpl++){
										// cout << mpl << endl;
										int m=mpl-l;
										int jp1St=max(-1, -(m+lp))+1; // include
										int jp1En=min(1, lp-m)+2; // not include
										int jp1;
										complex<double> coeff1(0, 0);
										complex<double> coeff2(0, 0);
										for(jp1=jp1St; jp1<jp1En; jp1++){
											int mpplp=jp1-1+m+lp;
											// load final state
											for(ir=0; ir<wfn_length[is]; ir++){
												final_states_re[ir]=final_states_FP[ik][FPIndex_1][ia][lp][mpplp][ir].real();
												final_states_im[ir]=final_states_FP[ik][FPIndex_1][ia][lp][mpplp][ir].imag();
											}
											// cout << jp1 << endl;
											double radial_part_re=ddot(&wfn_length[is], &final_states_re[0], &wfn_phi_rdr[is][io2][0]);
											double radial_part_im=ddot(&wfn_length[is], &final_states_im[0], &wfn_phi_rdr[is][io2][0]);
											complex<double> radial_part(radial_part_re, -radial_part_im);
											coeff1+=radial_part*Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1];
											if(spin_i==2){
												// load final state
												for(ir=0; ir<wfn_length[is]; ir++){
													final_states_re[ir]=final_states_FP[ik][FPIndex_2][ia][lp][mpplp][ir].real();
													final_states_im[ir]=final_states_FP[ik][FPIndex_2][ia][lp][mpplp][ir].imag();
												}
												radial_part_re=ddot(&wfn_length[is], &final_states_re[0], &wfn_phi_rdr[is][io2][0]);
												radial_part_im=ddot(&wfn_length[is], &final_states_im[0], &wfn_phi_rdr[is][io2][0]);
												complex<double> radial_part(radial_part_re, -radial_part_im);
												coeff2+=radial_part*Gaunt_arr[l][mpl][lp][mpplp]*Y_coeff[jp1];
											}
										}
										if(PA_weighting==true){
											coeff1*=atom_weighting[ia];
											coeff2*=atom_weighting[ia];
										}
										// cout << coeff1 << endl;
										PAD_1+=LCAO_use[mpl][0]*coeff1;
										if(spin_i==2){
											PAD_2+=LCAO_use[mpl][1]*coeff2;
										}
									}
								}
							}
						} // end of if(FPFS==false)
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

	w_data_2d(atomG_out, "Coordinates", atom_length, 3, (double**) atom_coordinates);
	if(PA_weighting==true){
		w_data_1d(atomG_out, "Weighting", atom_length, (double*) atom_weighting);
	}

	w_data_2d(atomG_out, "UnitCell", 3, 3, (double**) atom_cell);

	if(PA_FPFS){
		Group FPFSG(rootG.createGroup("FPFS"));
		w_att_double(FPFSG, "Excitation_energy", PA_excitation_energy);
		w_att_double(FPFSG, "FPFS_energy_step", PA_FPFS_energy_step);
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
			w_att_1i(FPFSG_k, "Energy_scale", final_states_FP_size[i], &final_states_EScale[i][0]);
			w_att_1d(FPFSG_k, "Energy_width", final_states_FP_size[i], &final_states_EWidth[i][0]);
			for(j=0; j<final_states_FP_size[i]; j++){
				sprintf(group_name, "%d", j);
				Group FPFSG_ke(FPFSG_k.createGroup(group_name));
				for(ia=0; ia<atom_length; ia++){
					int is=atom_spec_index[ia];
					sprintf(group_name, "%d_%s", (ia+1), atom_labels[ia]);
					Group FPFSG_kea(FPFSG_ke.createGroup(group_name));
					w_att_1d(FPFSG_kea, "r", wfn_length[is], &wfn_r[is][0]);
					double wfn_component[wfn_length[is]];
					for(int l=0; l<5; l++){
						for(int mpl=0; mpl<2*l+1; mpl++){
							int m=mpl-l;
							sprintf(group_name, "%d_%d_real", l, m);
							for(ir=0; ir<wfn_length[is]; ir++){
								wfn_component[ir]=final_states_FP[i][j][ia][l][mpl][ir].real();
							}
							w_data_1d(FPFSG_kea, group_name, wfn_length[is], &wfn_component[0]);
							sprintf(group_name, "%d_%d_imag", l, m);
							for(ir=0; ir<wfn_length[is]; ir++){
								wfn_component[ir]=final_states_FP[i][j][ia][l][mpl][ir].imag();
							}
							w_data_1d(FPFSG_kea, group_name, wfn_length[is], &wfn_component[0]);
						}
					}
				}
			}
		}
	}
	*/
	return 0;
}
