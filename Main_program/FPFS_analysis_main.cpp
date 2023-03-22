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
			double k_length=sqrt(inner_product(k, k));
			printf("    k vector = [%7.3f, %7.3f, %7.3f]\n", k[0], k[1], k[2]);
			printf("    spin = %s\n", spin[j]==0 ? "Up": "Dn");
			printf("    t kx ky kz  fourier components\n");
			double* k_t;
			double k_list[1024][3];
			int t_count=0;
			// prepare k_list
			for(double tau=1.1; tau>=-1.1; tau-=0.01){
				k_list[t_count][0]=tau*k[0];
				k_list[t_count][1]=k[1];
				k_list[t_count][2]=k[2];
				t_count++;
			}
			double k_test[3];
			for(int n1=-FPFS_range; n1<=FPFS_range; n1++){
					for(int n2=-FPFS_range; n2<=FPFS_range; n2++){
						// (n1, n2)=(0, 0) is skipped because it is the first element
						k_test[0]=0.0;
						for(int p=1; p<=2; p++){
							k_test[p]=k[p]+rec_cell[1][p]*n1+rec_cell[2][p]*n2;
						}
						double kz_square=k_length*k_length-inner_product(k_test, k_test);
						if(kz_square>0){
							// printf("%d %d\n", n1, n2);
							k_test[0]=sqrt(kz_square);
							for(int p=0; p<3; p++){
								k_list[t_count][p]=k_test[p];
							}
							t_count++;
							k_test[0]*=-1;
							for(int p=0; p<3; p++){
								k_list[t_count][p]=k_test[p];
							}
							t_count++;
						}
					}
				}

			
			for(int t=0; t<t_count; t++){
				k_t=k_list[t];
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
				printf("%3d %7.3f %7.3f %7.3f  ", t, k_t[0], k_t[1], k_t[2]);
				for(int id=0; id<digit; id++){
					printf("%7.3f %7.3f  ", final_state_fourier[id].real(), final_state_fourier[id].imag());
				}
				printf("\n");
			}
		}
	}
	return 0;
}
