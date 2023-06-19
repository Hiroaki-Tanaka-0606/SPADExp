#include <complex>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include "allocation_tools.hpp"
#include "log.hpp"
#include "variables_ext.hpp"

using namespace std;

int factorial(int n){
	int ret=1;
	int i;
	for(i=1; i<=n; i++){
		ret*=i;
	}
	return ret;
}

int** permutation_all(int n, int offset){
	int** ret=alloc_imatrix(factorial(n), n);
	if(n==1){
		ret[0][0]=0;
		return ret;
	}
	int** retm1=permutation_all(n-1, 0);
	int* listm1=new int[n-1];
	int index=0;
	for(int i=0; i<n; i++){
		int indexm1=0;
		for(int j=0; j<n; j++){
			if(j!=i){
				listm1[indexm1]=j;
				indexm1++;
			}
		}
		for(int j=0; j<factorial(n-1); j++){
			ret[index][0]=i+offset;
			for(int k=1; k<n; k++){
				ret[index][k]=listm1[retm1[j][k-1]]+offset;
			}
			index++;
		}
	}

	delete_imatrix(retm1);
	delete[] listm1;
	return ret;
}

void resolve_connection(int l_count, complex<double>* eigen_l, int r_count, complex<double>* eigen_r, int* connection){
	char* sprintf_buffer2=new char[Log_length+1];
	int offset=0;
	if(l_count>r_count){
		offset=r_count-l_count;
	}else{
	  sprintf(sprintf_buffer2, "Warning: resolve_connection encountered increasing eigenstates around E=%8.4f", eigen_l[0].real());
	  write_log(sprintf_buffer2);
		delete[] sprintf_buffer2;
		return;
	}

	if(l_count>=10){
	  sprintf(sprintf_buffer2, "Warning: permutation size = %d is quite large around E=%8.4f", l_count, eigen_l[0].real());
		write_log(sprintf_buffer2);
		// printf("%8.4f\n", eigen_l[0].real());
		delete[] sprintf_buffer2;
		return;
	}

	int** connection_candidates=permutation_all(l_count, offset);

	/*
	printf("Offset: %d\n", offset);
	
	for(int i=0; i<factorial(l_count); i++){
		for(int j=0; j<l_count; j++){
			printf("%3d ", connection_candidates[i][j]);
		}
		printf("\n");
	}
	printf("\n");*/

	double distance_min=-1;
	int distance_min_index=-1;
	for(int i=0; i<factorial(l_count); i++){
		double d=0;
		for(int j=0; j<l_count; j++){
			if(connection_candidates[i][j]>=0){
				d+=abs(eigen_l[j]-eigen_r[connection_candidates[i][j]]);
			}
		}
		if(distance_min_index<0 || distance_min>d){
			distance_min_index=i;
			distance_min=d;
		}
	}

	// printf("Averaged distance: %8.4f\n", distance_min/(l_count*1.0));
	for(int i=0; i<l_count; i++){
		connection[i]=connection_candidates[distance_min_index][i];
	}

	delete[] sprintf_buffer2;
	delete_imatrix(connection_candidates);
}
