#include <complex>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include "log.hpp"
#include "variables_ext.hpp"
#include "allocation_tools.hpp"
#include "resolve_connection.hpp"
using namespace std;


extern "C"{
	void zgesv_(
							int* N,
							int* NRHS,
							complex<double>* A,
							int* LDA,
							int* IPIV,
							complex<double>* B,
							int* LDB,
							int* INFO);
	void dgesv_(
							int* N,
							int* NRHS,
							double* A,
							int* LDA,
							int* IPIV,
							double* B,
							int* LDB,
							int* INFO);
	void dgetrf_(
							 int* M,
							 int* N,
							 double* A,
							 int* LDA,
							 int* IPIV,
							 int* INFO);
	void zgetrf_(
							 int* M,
							 int* N,
							 complex<double>* A,
							 int* LDA,
							 int* IPIV,
							 int* INFO);
	void zgetrs_(
							 char* TRANS,
							 int* N,
							 int* NRHS,
							 complex<double>* A,
							 int* LDA,
							 int* IPIV,
							 complex<double>* B,
							 int* LDB,
							 int* INFO);
	void dgetri_(
							 int* N,
							 double* A,
							 int* LDA,
							 int* IPIV,
							 double* WORK,
							 int* LWORK,
							 int* INFO);
	void zgetri_(
							 int* N,
							 complex<double>* A,
							 int* LDA,
							 int* IPIV,
							 complex<double>* WORK,
							 int* LWORK,
							 int* INFO);
	void dgemm_(
							char* TRANSA,
							char* TRANSB,
							int* M,
							int* N,
							int* K,
							double* ALPHA,
							double* A,
							int* LDA,
							double* B,
							int* LDB,
							double* BETA,
							double* C,
							int* LDC);
	void zgemm_(
							char* TRANSA,
							char* TRANSB,
							int* M,
							int* N,
							int* K,
							complex<double>* ALPHA,
							complex<double>* A,
							int* LDA,
							complex<double>* B,
							int* LDB,
							complex<double>* BETA,
							complex<double>* C,
							int* LDC);
	void dgemv_(
							char* TRANS,
							int* M,
							int* N,
							double* ALPHA,
							double* A,
							int* LDA,
							double* X,
							int* INCX,
							double* BETA,
							double* Y,
							int* INCY);
	void zgemv_(
							char* TRANS,
							int* M,
							int* N,
							complex<double>* ALPHA,
							complex<double>* A,
							int* LDA,
							complex<double>* X,
							int* INCX,
							complex<double>* BETA,
							complex<double>* Y,
							int* INCY);
	void dgeev_(
							char* JOBVL,
							char* JOBVR,
							int* N,
							double* A,
							int* LDA,
							double* WR,
							double* WI,
							double* VL,
							int* LDVL,
							double* VR,
							int* LDVR,
							double* WORK,
							int* LWORK,
							int* INFO);
	
	double ddot_(
							 int* N,
							 double* CX,
							 int* INCX,
							 double* CY,
							 int* INCY);

	void zgels_(
							char* TRANS,
							int* M,
							int* N,
							int* NRHS,
							complex<double>* A,
							int* LDA,
							complex<double>* B,
							int* LDB,
							complex<double>* WORK,
							int* LWORK,
							int* INFO);

	void zheev_(
							char* JOBZ,
							char* UPLO,
							int* N,
							complex<double>* A,
							int* LDA,
							double* W,
							complex<double>* WORK,
							int* LWORK,
							double* RWORK,
							int* INFO);
	void zgeev_(
							char* JOBVL,
							char* JOBVR,
							int* N,
							complex<double>* A,
							int* LDA,
							complex<double>* W,
							complex<double>* VL,
							int* LDVL,
							complex<double>* VR,
							int* LDVR,
							complex<double>* WORK,
							int* LWORK,
							double* RWORK,
							int* INFO);
  
}

// see Main_GUI/lib/physical_tools for the detail of the conversion
complex<double>**** convert_LCAO(int size_k, int size_b, int size_m, int size_d, double* LCAO_raw){
	int ik; // k point
	int ib; // band
	int im; // m+l (max 2l+1)
	int id; // digit (1 or 2)
	int sk=size_k;
	int sb=size_b;
	int sm=size_m;
	int sd=size_d;
	int sd2=sd/2;
	complex<double> I(0, 1);
	int idx1, idx2, idx3;
	if(size_d%2!=0){
		return NULL;
	}
  complex<double>* LCAO_alloc=new complex<double>[sk*sb*sm*sd2];
	complex<double>**** LCAO=new complex<double>***[sk];

#pragma omp parallel private(ib, im, id, idx1, idx2, idx3), firstprivate(I, sk, sb, sm, sd, sd2)
#pragma omp for	
	for(ik=0; ik<sk; ik++){
		idx1=ik*sb*sm*sd;
		LCAO[ik]=new complex<double>**[sb];
		for(ib=0; ib<sb; ib++){
			idx2=ib*sm*sd;
			LCAO[ik][ib]=new complex<double>*[sm];
			for(im=0; im<sm; im++){
				idx3=im*sd;
				LCAO[ik][ib][im]=&LCAO_alloc[(idx1+idx2+idx3)/2];
			}
			for(id=0; id<sd2; id++){
				complex<double> px, py, pz, pm1, pp0, pp1;
				complex<double> d3z2r2, dx2y2, dxy, dxz, dyz, dm2, dm1, dp0, dp1, dp2;
				complex<double> f5z23r2, f5xy2xr2, f5yz2yr2, fzx2zy2, fxyz, fx33xy2, f3yx2y3;
				complex<double> fm3, fm2, fm1, fp0, fp1, fp2, fp3;
				switch(sm){
				case 1:
					// s orbital: nothing to calculate
					LCAO[ik][ib][0][id]=complex<double>(LCAO_raw[idx1+idx2+2*id], LCAO_raw[idx1+idx2+2*id+1]);
					break;
				case 3:
					// p orbital
				  px=complex<double>(LCAO_raw[idx1+idx2+     2*id], LCAO_raw[idx1+idx2+     2*id+1]);
				  py=complex<double>(LCAO_raw[idx1+idx2+  sd+2*id], LCAO_raw[idx1+idx2+  sd+2*id+1]);
				  pz=complex<double>(LCAO_raw[idx1+idx2+2*sd+2*id], LCAO_raw[idx1+idx2+2*sd+2*id+1]);
					/*
						pm1= (px-I*py)/sqrt(2);
						pp0=  pz;
						pp1=-(px+I*py)/sqrt(2);*/
				  pm1= (px+I*py)/sqrt(2);
				  pp0=  pz;
					pp1=-(px-I*py)/sqrt(2);
					LCAO[ik][ib][0][id]=pm1;
					LCAO[ik][ib][1][id]=pp0;
					LCAO[ik][ib][2][id]=pp1;
					break;
				case 5:
					// d orbital
				  d3z2r2=complex<double>(LCAO_raw[idx1+idx2+     2*id], LCAO_raw[idx1+idx2+     2*id+1]);
				  dx2y2 =complex<double>(LCAO_raw[idx1+idx2+  sd+2*id], LCAO_raw[idx1+idx2+  sd+2*id+1]);
				  dxy   =complex<double>(LCAO_raw[idx1+idx2+2*sd+2*id], LCAO_raw[idx1+idx2+2*sd+2*id+1]);
				  dxz   =complex<double>(LCAO_raw[idx1+idx2+3*sd+2*id], LCAO_raw[idx1+idx2+3*sd+2*id+1]);
				  dyz   =complex<double>(LCAO_raw[idx1+idx2+4*sd+2*id], LCAO_raw[idx1+idx2+4*sd+2*id+1]);
					/*
						dm2= (dx2y2-I*dxy)/sqrt(2);
						dm1= (dxz-  I*dyz)/sqrt(2);
						dp0=  d3z2r2;
						dp1=-(dxz+  I*dyz)/sqrt(2);
						dp2= (dx2y2+I*dxy)/sqrt(2);*/
				  dm2= (dx2y2+I*dxy)/sqrt(2);
				  dm1= (dxz+  I*dyz)/sqrt(2);
				  dp0=  d3z2r2;
				  dp1=-(dxz-  I*dyz)/sqrt(2);
				  dp2= (dx2y2-I*dxy)/sqrt(2);
					LCAO[ik][ib][0][id]=dm2;
					LCAO[ik][ib][1][id]=dm1;
					LCAO[ik][ib][2][id]=dp0;
					LCAO[ik][ib][3][id]=dp1;
					LCAO[ik][ib][4][id]=dp2;
					break;
				case 7:
					// f orbital
					f5z23r2 =complex<double>(LCAO_raw[idx1+idx2+     2*id], LCAO_raw[idx1+idx2+     2*id+1]);
				  f5xy2xr2=complex<double>(LCAO_raw[idx1+idx2+  sd+2*id], LCAO_raw[idx1+idx2+  sd+2*id+1]);
				  f5yz2yr2=complex<double>(LCAO_raw[idx1+idx2+2*sd+2*id], LCAO_raw[idx1+idx2+2*sd+2*id+1]);
				  fzx2zy2 =complex<double>(LCAO_raw[idx1+idx2+3*sd+2*id], LCAO_raw[idx1+idx2+3*sd+2*id+1]);
				  fxyz    =complex<double>(LCAO_raw[idx1+idx2+4*sd+2*id], LCAO_raw[idx1+idx2+4*sd+2*id+1]);
					fx33xy2 =complex<double>(LCAO_raw[idx1+idx2+5*sd+2*id], LCAO_raw[idx1+idx2+5*sd+2*id+1]);
					f3yx2y3 =complex<double>(LCAO_raw[idx1+idx2+6*sd+2*id], LCAO_raw[idx1+idx2+6*sd+2*id+1]);
					/*
						fm3= (fx33xy2- I*f3yx2y3)/sqrt(2);
						fm2= (fzx2zy2- I*fxyz)/sqrt(2);
						fm1= (f5xy2xr2-I*f5yz2yr2)/sqrt(2);
						fp0=  f5z23r2;
						fp1=-(f5xy2xr2+I*f5yz2yr2)/sqrt(2);
						fp2= (fzx2zy2+ I*fxyz)/sqrt(2);
						fp3=-(fx33xy2+ I*f3yx2y3)/sqrt(2);*/
					fm3= (fx33xy2+ I*f3yx2y3)/sqrt(2);
					fm2= (fzx2zy2+ I*fxyz)/sqrt(2);
					fm1= (f5xy2xr2+I*f5yz2yr2)/sqrt(2);
					fp0=  f5z23r2;
					fp1=-(f5xy2xr2-I*f5yz2yr2)/sqrt(2);
					fp2= (fzx2zy2- I*fxyz)/sqrt(2);
					fp3=-(fx33xy2- I*f3yx2y3)/sqrt(2);
					LCAO[ik][ib][0][id]=fm3;
					LCAO[ik][ib][1][id]=fm2;
					LCAO[ik][ib][2][id]=fm1;
					LCAO[ik][ib][3][id]=fp0;
					LCAO[ik][ib][4][id]=fp1;
					LCAO[ik][ib][5][id]=fp2;
					LCAO[ik][ib][6][id]=fp3;
				}
			}
		}
	}
	return LCAO;
}

// outer product of three-dimensional vectors
void outer_product(double* a, double* b, double* c){
	// c=a\times b
	c[0]=a[1]*b[2]-a[2]*b[1];
	c[1]=a[2]*b[0]-a[0]*b[2];
	c[2]=a[0]*b[1]-a[1]*b[0];
}

// inner product of three-dimensional vectors
double inner_product(double* a, double* b){
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

// coefficient for Y_1m
/// theta, phi in degree
void operator_coefficient(char* polarization, double theta, double phi, complex<double>* Y_coeff){
	theta=M_PI*theta/180.0;
	phi=M_PI*phi/180.0;

	if(strcmp(polarization, "Linear")==0){
		Y_coeff[2]=-sqrt(2.0*M_PI/3.0)*sin(theta)*complex<double>(cos(phi), -sin(phi));
		Y_coeff[1]= sqrt(4.0*M_PI/3.0)*cos(theta);
		Y_coeff[0]= sqrt(2.0*M_PI/3.0)*sin(theta)*complex<double>(cos(phi), sin(phi));
	}else if(strcmp(polarization, "RCircular")==0){
		Y_coeff[2]=sqrt(2.0*M_PI/3.0)*(1.0+cos(theta))*complex<double>(cos(phi),-sin(phi));
		Y_coeff[1]=sqrt(4.0*M_PI/3.0)*sin(theta);
		Y_coeff[0]=sqrt(2.0*M_PI/3.0)*(1.0-cos(theta))*complex<double>(cos(phi), sin(phi));
	}else if(strcmp(polarization, "LCircular")==0){
		Y_coeff[2]= sqrt(2.0*M_PI/3.0)*(-1.0+cos(theta))*complex<double>(cos(phi), -sin(phi));
		Y_coeff[1]= sqrt(4.0*M_PI/3.0)*sin(theta);
		Y_coeff[0]=-sqrt(2.0*M_PI/3.0)*( 1.0+cos(theta))*complex<double>(cos(phi), sin(phi));
	}
}
// electric field vector
/// theta, phi in degree
void electric_field_vector(char* polarization, double theta, double phi, complex<double>* e_vec){
	theta=M_PI*theta/180.0;
	phi=M_PI*phi/180.0;

	if(strcmp(polarization, "Linear")==0){
		e_vec[0]=sin(theta)*cos(phi);
		e_vec[1]=sin(theta)*sin(phi);
		e_vec[2]=cos(theta);
	}else if(strcmp(polarization, "RCircular")==0){
		e_vec[0]=complex<double>(-cos(theta)*cos(phi), sin(phi));
		e_vec[1]=complex<double>(-cos(theta)*sin(phi), -cos(phi));
		e_vec[2]=sin(theta);
	}else if(strcmp(polarization, "LCircular")==0){
		e_vec[0]=complex<double>(-cos(theta)*cos(phi), -sin(phi));
		e_vec[1]=complex<double>(-cos(theta)*sin(phi), cos(phi));
		e_vec[2]=sin(theta);
	}
}

void spherical_harmonics(double* r, complex<double>* Ylm2){
	// r=[x, y, z]
  // x=r*sin(theta)*cos(phi)
  // y=r*sin(theta)*sin(phi)
  // sqrt(x^2+y^2)=r*sin(theta)
  // z=r*cos(theta)
	
	// Ylm2[5][9]

	complex<double>** Ylm=new complex<double>*[6];
	int i;
	for(i=0; i<6; i++){
		Ylm[i]=&Ylm2[11*i];
	}
	for(i=0; i<66; i++){
		Ylm2[i]=complex<double>(0, 0);
	}
    
	double r_length=sqrt(inner_product(r, r));
	double cosT=1;
	double sinT=0;
	if(r_length>PA_zero_threshold){
		cosT=r[2]/r_length;
		sinT=sqrt(r[0]*r[0]+r[1]*r[1])/r_length;	
	}
	double cosP=1;
	double sinP=0;
	if(abs(sinT)>PA_zero_threshold){
		cosP=r[0]/(r_length*sinT);
		sinP=r[1]/(r_length*sinT);
	}

	double cos2P=cosP*cosP-sinP*sinP;
	double sin2P=2*sinP*cosP;

	double cos3P=4*cosP*cosP*cosP-3*cosP;
	double sin3P=3*sinP-4*sinP*sinP*sinP;

	double cos4P=8*cosP*cosP*cosP*cosP-8*cosP*cosP+1;
	double sin4P=4*sinP*cosP*(2*cosP*cosP-1);

	double cos5P=16*cosP*cosP*cosP*cosP*cosP-20*cosP*cosP*cosP+5*cosP;
	double sin5P=16*sinP*sinP*sinP*sinP*sinP-20*sinP*sinP*sinP+5*sinP;

	// s(0): 1/2 1/sqrt(pi)
	Ylm[0][0]=1.0/(2.0*sqrt(M_PI));
    
	// p(-1): 1/2 sqrt(3/2pi) sin(T)(cos(P)-i*sin(P))
	Ylm[1][0]=(1.0/2.0)*sqrt(3.0/(2.0*M_PI))*sinT*complex<double>(cosP, -sinP);
	// p(0): 1/2 sqrt(3/pi) cos(T)
	Ylm[1][1]=(1.0/2.0)*sqrt(3.0/M_PI)*cosT;
	// p(1): -1/2 sqrt(3/2pi) sin(T)(cos(P)+i*sin(P))
	Ylm[1][2]=-(1.0/2.0)*sqrt(3.0/(2.0*M_PI))*sinT*complex<double>(cosP, sinP);

	// d(-2): 1/4 sqrt(15/2pi) sin^2(T)(cos(2P)-i*sin(2P))
	Ylm[2][0]=(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*sinT*sinT*complex<double>(cos2P, -sin2P);
	// d(-1): 1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)-i*sin(P))
	Ylm[2][1]=(1.0/2.0)*sqrt(15.0/(2.0*M_PI))*sinT*cosT*complex<double>(cosP, -sinP);
	// d(0): 1/4 sqrt(5/pi) (3cos^2(T)-1)
	Ylm[2][2]=(1.0/4.0)*sqrt(5.0/M_PI)*(3*cosT*cosT-1.0);
	// d(1): -1/2 sqrt(15/2pi) sin(T)cos(T)(cos(P)+i*sin(P))
	Ylm[2][3]=-(1.0/2.0)*sqrt(15.0/(2.0*M_PI))*sinT*cosT*complex<double>(cosP, sinP);
	// d(2): 1/4 sqrt(15/2pi) sin^2(T)(cos(2P)+i*sin(2P))
	Ylm[2][4]=(1.0/4.0)*sqrt(15.0/(2.0*M_PI))*sinT*sinT*complex<double>(cos2P, sin2P);
    
	// f(-3): 1/8 sqrt(35/pi) sin^3(T)(cos(3P)-i*sin(3P))
	Ylm[3][0]=(1.0/8.0)*sqrt(35.0/M_PI)*sinT*sinT*sinT*complex<double>(cos3P, -sin3P);
	// f(-2): 1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)-i*sin(2P))
	Ylm[3][1]=(1.0/4.0)*sqrt(105.0/(2.0*M_PI))*sinT*sinT*cosT*complex<double>(cos2P, -sin2P);
	// f(-1): 1/8 sqrt(21/pi) (5cos^2(T)-1)sin(T)(cos(P)-i*sin(P))
	Ylm[3][2]=(1.0/8.0)*sqrt(21.0/M_PI)*(5.0*cosT*cosT-1.0)*sinT*complex<double>(cosP, -sinP);
	// f(0): 1/4 sqrt(7/pi) (5cos^2(T)-3)cos(T)
	Ylm[3][3]=(1.0/4.0)*sqrt(7.0/M_PI)*(5.0*cosT*cosT-3.0)*cosT;
	// f(1): -1/8 sqrt(21/pi) (5cos^2(T)-1)sin(T)(cos(P)+i*sin(P))
	Ylm[3][4]=-(1.0/8.0)*sqrt(21.0/M_PI)*(5.0*cosT*cosT-1.0)*sinT*complex<double>(cosP, sinP);
	// f(2): 1/4 sqrt(105/2pi) sin^2(T)cos(T)(cos(2P)+i*sin(2P))
	Ylm[3][5]=(1.0/4.0)*sqrt(105.0/(2.0*M_PI))*sinT*sinT*cosT*complex<double>(cos2P, sin2P);
	// f(3): -1/8 sqrt(35/pi) sin^3(T)(cos(3P)+i*sin(3P))
	Ylm[3][6]=-(1.0/8.0)*sqrt(35.0/M_PI)*sinT*sinT*sinT*complex<double>(cos3P, sin3P);
    
	// g(-4): 3/16 sqrt(35/2pi) sin^4(T)(cos(4P)-i*sin(4P))
	Ylm[4][0]=(3.0/16.0)*sqrt(35.0/(2.0*M_PI))*sinT*sinT*sinT*sinT*complex<double>(cos4P, -sin4P);
	// g(-3): 3/8 sqrt(35/pi) sin^3(T)cos(T)(cos(3P)-i*sin(3P))
	Ylm[4][1]=(3.0/8.0)*sqrt(35.0/M_PI)*sinT*sinT*sinT*cosT*complex<double>(cos3P, -sin3P);
	// g(-2): 3/8 sqrt(5/2pi) (7cos^2(T)-1)sin^2(T)(cos(2P)-i*sin(2P))
	Ylm[4][2]=(3.0/8.0)*sqrt(5.0/(2.0*M_PI))*(7.0*cosT*cosT-1.0)*sinT*sinT*complex<double>(cos2P, -sin2P);
	// g(-1): 3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)-i*sin(P))
	Ylm[4][3]=(3.0/8.0)*sqrt(5.0/M_PI)*(7.0*cosT*cosT-3.0)*sinT*cosT*complex<double>(cosP, -sinP);
	// g(0): 3/16 sqrt(1/pi) (35cos^4(T)-30cos^2(T)+3)
	Ylm[4][4]=(3.0/16.0)*sqrt(1.0/M_PI)*(35.0*cosT*cosT*cosT*cosT-30.0*cosT*cosT+3.0);
	// g(1): -3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)+i*sin(P))
	Ylm[4][5]=-(3.0/8.0)*sqrt(5.0/M_PI)*(7.0*cosT*cosT-3.0)*sinT*cosT*complex<double>(cosP, sinP);
	// g(2): 3/8 sqrt(5/2pi) (7cos^2(T)-1)sin^2(T)(cos(2P)+i*sin(2P))
	Ylm[4][6]=(3.0/8.0)*sqrt(5.0/(2.0*M_PI))*(7.0*cosT*cosT-1.0)*sinT*sinT*complex<double>(cos2P, sin2P);
	// g(3): -3/8 sqrt(35/pi) sin^3(T)cos(T)(cos(3P)+i*sin(3P))
	Ylm[4][7]=-(3.0/8.0)*sqrt(35.0/M_PI)*sinT*sinT*sinT*cosT*complex<double>(cos3P, sin3P);
	// g(4): 3/16 sqrt(35/2pi) sin^4(T)(cos(4P)+i*sin(4P))
	Ylm[4][8]=(3.0/16.0)*sqrt(35.0/(2.0*M_PI))*sinT*sinT*sinT*sinT*complex<double>(cos4P, sin4P);

	// h(-5): 3/32 sqrt(77/pi) sin^5(T)(cos(5P)-i*sin(5P))
	Ylm[5][0]=(3.0/32.0)*sqrt(77.0/M_PI)*sinT*sinT*sinT*sinT*sinT*complex<double>(cos5P, -sin5P);
	// h(-4): 3/16 sqrt(385/2pi) sin^4(T)cos(T)(cos(4P)-i*sin(4P))
	Ylm[5][1]=(3.0/16.0)*sqrt(385.0/(2.0*M_PI))*sinT*sinT*sinT*sinT*cosT*complex<double>(cos4P, -sin4P);
	// h(-3): 1/32 sqrt(385/pi) sin^3(T)(9cos^2(T)-1)(cos(3P)-i*sin(3P))
	Ylm[5][2]=(1.0/32.0)*sqrt(385.0/M_PI)*sinT*sinT*sinT*(9.0*cosT*cosT-1.0)*complex<double>(cos3P, -sin3P);
	// h(-2): 1/8 sqrt(1155/2pi) sin^2(T)(3cos^3(T)-cos(T))(cos(2P)-i*sin(2P))
	Ylm[5][3]=(1.0/8.0)*sqrt(1155.0/(2.0*M_PI))*sinT*sinT*(3.0*cosT*cosT*cosT-cosT)*complex<double>(cos2P, -sin2P);
	// h(-1): 1/16 sqrt(165/2pi) sin(T)(21cos^4(T)-14cos^2(T)+1)(cos(P)-i*sin(P))
	Ylm[5][4]=(1.0/16.0)*sqrt(165.0/(2.0*M_PI))*sinT*(21.0*cosT*cosT*cosT*cosT-14.0*cosT*cosT+1.0)*complex<double>(cosP, -sinP);
	// h(0): 1/16 sqrt(11/pi) (63cos^5(T)-70cos^3(T)+15cos(T))
	Ylm[5][5]=(1.0/16.0)*sqrt(11.0/M_PI)*(63.0*cosT*cosT*cosT*cosT*cosT-70.0*cosT*cosT*cosT+15.0*cosT);
	// h(1): -1/16 sqrt(165/2pi) sin(T)(21cos^4(T)-14cos^2(T)+1)(cos(P)+i*sin(P))
	Ylm[5][6]=-(1.0/16.0)*sqrt(165.0/(2.0*M_PI))*sinT*(21.0*cosT*cosT*cosT*cosT-14.0*cosT*cosT+1.0)*complex<double>(cosP, sinP);
	// h(2): 1/8 sqrt(1155/2pi) sin^2(T)(3cos^3(T)-cos(T))(cos(2P)+i*sin(2P))
	Ylm[5][7]=(1.0/8.0)*sqrt(1155.0/(2.0*M_PI))*sinT*sinT*(3.0*cosT*cosT*cosT-cosT)*complex<double>(cos2P, sin2P);
	// h(3): -1/32 sqrt(385/pi) sin^3(T)(9cos^2(T)-1)(cos(3P)+i*sin(3P))
	Ylm[5][8]=-(1.0/32.0)*sqrt(385.0/M_PI)*sinT*sinT*sinT*(9.0*cosT*cosT-1.0)*complex<double>(cos3P, sin3P);
	// h(4): 3/16 sqrt(385/2pi) sin^4(T)cos(T)(cos(4P)+i*sin(4P))
	Ylm[5][9]=(3.0/16.0)*sqrt(385.0/(2.0*M_PI))*sinT*sinT*sinT*sinT*cosT*complex<double>(cos4P, sin4P);
	// h(5): -3/32 sqrt(77/pi) sin^5(T)(cos(5P)-i*sin(5P))
	Ylm[5][10]=-(3.0/32.0)*sqrt(77.0/M_PI)*sinT*sinT*sinT*sinT*sinT*complex<double>(cos5P, sin5P);
	
	delete[] Ylm;
}

double Gaunt(int lp, int mp, int l, int m){
	if(abs(mp)>lp || abs(m)>l){
		return 0;
	}

	if(mp==m){
		if(lp==l+1){
			return sqrt((3.0/(4.0*M_PI))*(((l+1.0)*(l+1.0)-m*m)/((2.0*l+3.0)*(2.0*l+1.0))));
		}else if(lp==l-1){
			return sqrt((3.0/(4.0*M_PI))*(1.0*(l*l-m*m)/((2.0*l-1.0)*(2.0*l+1.0))));
		}else{
			return 0;
		}
	}else if(mp==m+1){
		if(lp==l+1){
			return sqrt((3.0/(4.0*M_PI))*(((l+m+2.0)*(l+m+1.0))/(2.0*(2.0*l+3.0)*(2.0*l+1.0))));
		}else if(lp==l-1){
			return -sqrt((3.0/(4.0*M_PI))*(((l-m)*(l-m-1.0))/(2.0*(2.0*l-1.0)*(2.0*l+1.0))));
		}else{
			return 0;
		}
	}else if(mp==m-1){
		if(lp==l+1){
			return sqrt((3.0/(4.0*M_PI))*(((l-m+2.0)*(l-m+1.0))/(2.0*(2.0*l+3.0)*(2.0*l+1.0))));
		}else if(lp==l-1){
			return -sqrt((3.0/(4.0*M_PI))*(((l+m)*(l+m-1.0))/(2.0*(2.0*l-1.0)*(2.0*l+1.0))));
		}else{
			return 0;
		}	
	}else{
		return 0;
	}
}

double sp_bessel(int l, double x){
	// spherical Bessel function j_l(x)
	if(x<PA_zero_threshold){
		return l==0 ? 1.0 : 0.0;
	}
	// use the asymptotic form if x is small
	if(x<0.02){
		double coef=1.0;
		for(int il=0; il<l; il++){
			coef*=2.0;
		}
		for(int il=l+1; il<=2*l+1; il++){
			coef/=(il*1.0);
		}
		for(int il=0; il<l; il++){
			coef*=x;
		}
		return coef;
	}
	
	// j_l(x)=(-1)^l x^l (1/x d/dx)^l sin(x)/x
	if(l==0){
		// sin(x)/x
		return sin(x)/x;
	}else if(l==1){
		// {sin(x)-x*cos(x)}/x^2
		return (sin(x)-x*cos(x))/(x*x);
		//return -cos(x)/x+sin(x)/(x*x);
	}else if(l==2){
		// {(3-x^2)sin(x)-3x*cos(x)}/x^3
		return ((3.0-x*x)*sin(x)-3.0*x*cos(x))/(x*x*x);
		// return -sin(x)/x-3.0*cos(x)/(x*x)+3.0*sin(x)/(x*x*x);
	}else if(l==3){
		// {(15-6x^2)sin(x)+(x^3-15x)cos(x)}/x^4
		return ((15.0-6.0*x*x)*sin(x)+(x*x*x-15.0*x)*cos(x))/(x*x*x*x);
		// return cos(x)/x-6.0*sin(x)/(x*x)-15.0*cos(x)/(x*x*x)+15.0*sin(x)/(x*x*x*x);
	}else if(l==4){
		// {(x^4-45x^2+105)sin(x)+(10x^3-105x)cos(x)}/x^5
		return ((x*x*x*x-45.0*x*x+105.0)*sin(x)+(10.0*x*x*x-105.0*x)*cos(x))/(x*x*x*x*x);
		// return sin(x)/x+10.0*cos(x)/(x*x)-45.0*sin(x)/(x*x)-105.0*cos(x)/(x*x*x*x)+105.0*sin(x)/(x*x*x*x*x);
	}else if(l==5){
		// {(15x^4-420x^2+945)*sin(x)-(x^5-105x^3+945x)*cos(x)}/x^6
		return ((15.0*x*x*x*x-420.0*x*x+945.0)*sin(x)-(x*x*x*x*x-105.0*x*x*x+945.0*x)*cos(x))/(x*x*x*x*x*x);
	}else{
		write_log((char*)"Invalid l");
		return 0.0;
	}
}

double cubeValue(int n1, int n2, int n3, int* count, double* cube){
	n1=n1%count[0];
	n2=n2%count[1];
	n3=n3%count[2];
	return cube[n1*count[1]*count[2]+n2*count[2]+n3];
}
// since the unit cell vectors a_i and reciprocal vectors b_j satisfy (a_i, b_j)=2pi d_{ij},
// r=p_1 a_1 + p_2 a_2 + p_3 a_3 --> p_i=(b_i, r)/2pi
double interpolate_potential(double* r, int* count, double* cube, double* atom_cell_buffer){
	int i, j;
	double atom_cell[3][3];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			atom_cell[i][j]=atom_cell_buffer[i*3+j];
		}
	}

	double rec_cell[3][3];
	double op[3]; // for outer product
	double det; // for determinant
	/// reciprocal unit cell
	outer_product(atom_cell[1], atom_cell[2], op);
	det=inner_product(atom_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[0][i]=2.0*M_PI*op[i]/det;
	}

	outer_product(atom_cell[2], atom_cell[0], op);
	for(i=0; i<3; i++){
		rec_cell[1][i]=2.0*M_PI*op[i]/det;
	}

	outer_product(atom_cell[0], atom_cell[1], op);
	for(i=0; i<3; i++){
		rec_cell[2][i]=2.0*M_PI*op[i]/det;
	}
	
	double p[3]; // fractional
	double q[3]; // non-integer index
	int qf[3]; // integer index (floored)
	double** rc=new double*[3];
	for(int i=0; i<3; i++){
		rc[i]=&rec_cell[i][0];
	}
	for(int i=0; i<3; i++){
		p[i]=inner_product(r, rc[i])/(2*M_PI);
		q[i]=p[i]*count[i];
		qf[i]=floor(q[i]);
	}
	double v000=cubeValue(qf[0],   qf[1],   qf[2],   count, cube);
	double v100=cubeValue(qf[0]+1, qf[1],   qf[2],   count, cube);
	double v010=cubeValue(qf[0],   qf[1]+1, qf[2],   count, cube);
	double v001=cubeValue(qf[0],   qf[1],   qf[2]+1, count, cube);
	double v011=cubeValue(qf[0],   qf[1]+1, qf[2]+1, count, cube);
	double v101=cubeValue(qf[0]+1, qf[1],   qf[2]+1, count, cube);
	double v110=cubeValue(qf[0]+1, qf[1]+1, qf[2],   count, cube);
	double v111=cubeValue(qf[0]+1, qf[1]+1, qf[2]+1, count, cube);

	double w000=(qf[0]+1-q[0])*(qf[1]+1-q[1])*(qf[2]+1-q[2]);
	double w100=(q[0]-qf[0])  *(qf[1]+1-q[1])*(qf[2]+1-q[2]);
	double w010=(qf[0]+1-q[0])*(q[1]-qf[1])  *(qf[2]+1-q[2]);
	double w001=(qf[0]+1-q[0])*(qf[1]+1-q[1])*(q[2]-qf[2])  ;
	double w011=(qf[0]+1-q[0])*(q[1]-qf[1])  *(q[2]-qf[2])  ;
	double w101=(q[0]-qf[0])  *(qf[1]+1-q[1])*(q[2]-qf[2])  ;
	double w110=(q[0]-qf[0])  *(q[1]-qf[1])  *(qf[2]+1-q[2]);
	double w111=(q[0]-qf[0])  *(q[1]-qf[1])  *(q[2]-qf[2]);

	delete[] rc;
	
	return
		v000*w000
		+v100*w100
		+v010*w010
		+v001*w001
		+v011*w011
		+v101*w101
		+v110*w110
		+v111*w111;
}


void Fourier_expansion_z(double* V_buffer, int* V_count, double* g, double* atom_cell_buffer, complex<double>* Vg, double* Vg_abs){
	int iz, ix, iy;
	int i,j;
	// buffer preparation
	double*** V=new double**[V_count[0]];
	for(iz=0; iz<V_count[0]; iz++){
		V[iz]=new double*[V_count[1]];
		for(ix=0; ix<V_count[1]; ix++){
			V[iz][ix]=&V_buffer[iz*V_count[1]*V_count[2]+ix*V_count[2]];
		}
	}
	double atom_cell[3][3];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			atom_cell[i][j]=atom_cell_buffer[i*3+j];
		}
	}

	int n1n2=V_count[1]*V_count[2];
	double r12[3];
	double norm_sum=0.0;
	for(iz=0; iz<V_count[0]; iz++){
		Vg[iz]=0;
		for(ix=0; ix<V_count[1]; ix++){
			double ix_frac=ix*1.0/V_count[1];
			for(iy=0; iy<V_count[2]; iy++){
				double iy_frac=iy*1.0/V_count[2];
				for(i=0; i<3; i++){
					r12[i]=atom_cell[1][i]*ix_frac+atom_cell[2][i]*iy_frac;
				}
				double gr=inner_product(g, r12);
				Vg[iz]+=complex<double>(cos(gr), -sin(gr))*V[iz][ix][iy];
			}
		}
		Vg[iz]/=n1n2;
		Vg_abs[iz]=abs(Vg[iz]);
	}
	for(iz=0; iz<V_count[0]; iz++){
		delete[] V[iz];
	}
	delete[] V;
}

// In the Numerov method, the solving direction is downward
// f_{i-1}=(1+h^2/12 a_{i-1})^-1 * [ 2(1-5h^2/12 a_i)*f_i - (1+h^2/12 a_{i+1})*f_{i+1} ]
void solve_final_state_Numerov(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int z_start, int V00_index, complex<double>*** Vgg, double** g_vec, complex<double>** FP_loc){
	// composite matrix
	/*
	complex<double>*** Vgg=alloc_zpmatrix(g_count);
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			Vgg[ig1][ig2]=Vgg_buffer[ig1*g_count+ig2];
		}
		}*/
	//double g_vec[g_count][3];
	//for(int ig=0; ig<g_count; ig++){
	//	for(int ix=0; ix<3; ix++){
	//		g_vec[ig][ix]=g_vec_buffer[ig*3+ix];
	//	}
	//}

	
	double kzz0=kz*dz*(z_count-1);
	double kzz0h=kz*dz*z_count;

	// Numerov and Euler
	//complex<double>* FP_loc[g_count];
	//for(int ig=0; ig<g_count; ig++){
	//	FP_loc[ig]=&FP_loc_buffer[ig*z_count];
	//}

	// i=z_count-1 (boundary condition)
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			FP_loc[ig][z_count-1]=complex<double>(cos(kzz0), sin(kzz0));
		}else{
			FP_loc[ig][z_count-1]=complex<double>(0, 0);
		}
	}

	// i=z_count-2 (Euler)
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			FP_loc[ig][z_count-2]=FP_loc[ig][z_count-1]-dz*kz*complex<double>(-sin(kzz0), cos(kzz0));
		}else{
			FP_loc[ig][z_count-2]=complex<double>(0, 0);
		}
	}

	// i=z_count-3 ~ 0 (Numerov)
	double* diagonal_part=new double[g_count];
	for(int ig=0; ig<g_count; ig++){
		double kpg[3];
		for(int p=0; p<3; p++){
			kpg[p]=k_para[p]+g_vec[ig][p];
		}
		//printf("g= (%f, %f, %f)\n", g_vec[ig][0], g_vec[ig][1], g_vec[ig][2]);
		//printf("k+g= (%f, %f, %f)\n", kpg[0], kpg[1], kpg[2]);
		diagonal_part[ig]=2*Ekin-inner_product(kpg, kpg);
		//printf("diag[%3d] = %f\n", ig, diagonal_part[ig]);
	}
	complex<double>*** a_matrix=alloc_zcube(z_count, g_count, g_count);
	for(int iz=0; iz<z_count; iz++){
		for(int ig1=0; ig1<g_count; ig1++){
			for(int ig2=0; ig2<g_count; ig2++){
				if(ig1==ig2){
					a_matrix[iz][ig1][ig2]=complex<double>(diagonal_part[ig1], 0);
					// printf("%d %f %f\n", iz, Vgg[ig1][ig2][iz].real(), Vgg[ig1][ig2][iz].imag());
				}else{
					a_matrix[iz][ig1][ig2]=complex<double>(0, 0);
				}
				a_matrix[iz][ig1][ig2]-=complex<double>(2, 0)*Vgg[ig1][ig2][iz];
				
				//printf("%d %f %f\n", iz, Vgg[ig1][ig2][iz].real(), Vgg[ig1][ig2][iz].imag());
			}
		}
	}
	complex<double>** left_matrix=alloc_zmatrix(g_count); // transposed!!
	complex<double>* right_vector=new complex<double>[g_count];
	int INFO=0;
	int NRHS=1;
	int* IPIV=new int[g_count];
	for(int im1=z_count-3; im1>=z_start; im1--){
		int i=im1+1;
		int ip1=i+1;
		for(int ig1=0; ig1<g_count; ig1++){
			right_vector[ig1]=complex<double>(2, 0)*FP_loc[ig1][i]-FP_loc[ig1][ip1];
			
			for(int ig2=0; ig2<g_count; ig2++){
				if(ig1==ig2){
					left_matrix[ig2][ig1]=complex<double>(1, 0);
				}else{
					left_matrix[ig2][ig1]=complex<double>(0, 0);
				}
				left_matrix[ig2][ig1]+=complex<double>(dz*dz/12.0, 0)*a_matrix[im1][ig1][ig2];

				right_vector[ig1]-=complex<double>(5.0/6.0*dz*dz, 0)*a_matrix[i][ig1][ig2]*FP_loc[ig2][i];
				right_vector[ig1]-=complex<double>(1.0/12.0*dz*dz, 0)*a_matrix[ip1][ig1][ig2]*FP_loc[ig2][ip1];
			}
		}
		zgesv_(&g_count, &NRHS, &left_matrix[0][0], &g_count, &IPIV[0], &right_vector[0], &g_count, &INFO);
		if(INFO==0){
			// ok
			for(int ig=0; ig<g_count; ig++){
				FP_loc[ig][im1]=right_vector[ig];
			}
		}else{
			write_log((char*)"zgesv failed");
			return;
		}
	}
	
	//delete_zpmatrix(Vgg);
	delete_zcube(a_matrix);
	delete_zmatrix(left_matrix);
	delete[] diagonal_part;
	delete[] right_vector;
	delete[] IPIV;
}

double solve_final_state_Matrix(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int z_start, int V00_index, complex<double>*** Vgg, double** g_vec, complex<double>** FP_loc, complex<double>** left_matrix, complex<double>** right_matrix, complex<double>* loc_edge){
  //cout << "Matrix " << endl;
	char* sprintf_buffer2=new char[Log_length+1];
	// composite matrix
	/*
	complex<double>*** Vgg=alloc_zpmatrix(g_count);
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			Vgg[ig1][ig2]=Vgg_buffer[ig1*g_count+ig2];
		}
		}*/
	//double g_vec[g_count][3];
	//for(int ig=0; ig<g_count; ig++){
	//	for(int ix=0; ix<3; ix++){
	//		g_vec[ig][ix]=g_vec_buffer[ig*3+ix];
	//	}
	//}

	
	double kzz0=kz*dz*(z_count-1);
	double kzz0h=kz*dz*z_count;

	// Linear equations
	int z_count_new=z_count-z_start;
	int eq_dim=z_count_new*g_count;
	//complex<double>** left_matrix=alloc_zmatrix(eq_dim); // transposed
	for(int i=0; i<eq_dim; i++){
		for(int j=0; j<eq_dim; j++){
			left_matrix[j][i]=complex<double>(0.0, 0.0);
		}
	}
	
	complex<double>* right_vector=new complex<double>[g_count];
	complex<double>* right_vector2=new complex<double>[g_count];

	double kpg[3];
	
	// left matrix
	/// M
	for(int ig1=0; ig1<g_count; ig1++){
		for(int p=0; p<3; p++){
			kpg[p]=k_para[p]+g_vec[ig1][p];
		}
		for(int ig2=0; ig2<g_count; ig2++){
			for(int iz=z_start; iz<z_count; iz++){
				int index1=ig1*z_count_new+(iz-z_start);
				int index2=ig2*z_count_new+(iz-z_start);
				complex<double> value=-2.0*Vgg[ig1][ig2][iz];
				if(ig1==ig2){
					value+=2.0*Ekin-inner_product(kpg, kpg);
				}
				left_matrix[index2][index1]+=value;
			}
		}
	}
	/// A
	for(int ig=0; ig<g_count; ig++){
		for(int iz1=z_start; iz1<z_count; iz1++){
			for(int iz2=z_start; iz2<z_count; iz2++){
				int index1=ig*z_count_new+(iz1-z_start);
				int index2=ig*z_count_new+(iz2-z_start);
				if(iz1==iz2){
					left_matrix[index2][index1]-=2.0/dz/dz;
				}
				if(iz1-1==iz2 || iz1+1==iz2){
					left_matrix[index2][index1]+=1.0/dz/dz;
				}
			}
		}
	}

	int info;
	int* ipiv=new int[eq_dim];
	int nrhs=2*g_count;
	char notrans='N';
	complex<double> alpha(1.0, 0.0);
	complex<double> beta(1.0, 0.0);
	int inc=1;
	// 0~g: for RL, g~2g: for RR
	//complex<double>** right_matrix=alloc_zmatrix(2*g_count, eq_dim); //transposed
	for(int igz1=0; igz1<eq_dim; igz1++){
		for(int ig2=0; ig2<2*g_count; ig2++){
			right_matrix[ig2][igz1]=complex<double>(0.0, 0.0);
		}
	}
	for(int ig1=0; ig1<g_count; ig1++){
		int index1=ig1*z_count_new;
		int index2=ig1;
		right_matrix[index2][index1]=complex<double>(1.0, 0.0);
		index1=(ig1+1)*z_count_new-1;
		index2=ig1+g_count;
		right_matrix[index2][index1]=complex<double>(1.0, 0.0);
	}
	/*
		for(int igz1=0; igz1<eq_dim; igz1++){
		for(int ig2=0; ig2<2*g_count; ig2++){
		printf("%4.1f ", abs(right_matrix[ig2][igz1]));
		}
		printf("\n");
		}*/
	zgesv_(&eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &right_matrix[0][0], &eq_dim, &info);
	if(info!=0){
		write_log((char*)"zgesv failed");
		return -1;
	}

	// extract the RL and RR matrices
	complex<double>** H_RL=alloc_zmatrix(g_count); // inverted
	complex<double>** H_RR=alloc_zmatrix(g_count); // inverted
	int index1, index2;
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			index1=(ig1+1)*z_count_new-1;
			index2=ig2;
			H_RL[ig2][ig1]=right_matrix[index2][index1];
				
			index1=(ig1+1)*z_count_new-1;
			index2=ig2+g_count;
			H_RR[ig2][ig1]=right_matrix[index2][index1];
		}
	}

	/*
		printf("H_RL\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
		printf("(%10.2e, %10.2e) ", H_RL[ig2][ig1].real(), H_RL[ig2][ig1].imag());
		}
		printf("\n");
		}
		
		printf("H_RR\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
		printf("(%10.2e, %10.2e) ", H_RR[ig2][ig1].real(), H_RR[ig2][ig1].imag());
		}
		printf("\n");
		}*/

	// compute the right vector
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			right_vector[ig]=-dz*dz*complex<double>(cos(kzz0), sin(kzz0));
		}else{
			right_vector[ig]=complex<double>(0.0, 0.0);
		}
	}
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			right_vector2[ig]=-complex<double>(cos(kzz0h), sin(kzz0h));
		}else{
			right_vector2[ig]=complex<double>(0.0, 0.0);
		}
	}

	// alpha=1
	// beta=1
	zgemv_(&notrans, &g_count, &g_count, &alpha, &H_RR[0][0], &g_count, &right_vector2[0], &inc, &beta, &right_vector[0], &inc);
		
	/*
		printf("Right\n");
		for(int ig=0; ig<g_count; ig++){
		printf("(%10.2e, %10.2e)\n", right_vector[ig].real(), right_vector[ig].imag());
		}
		printf("\n");*/

	// obtain f_G(-h)
	// by zgesv
	/*
		int ipiv2[g_count];
		zgesv_(&g_count, &nrhs, &H_RL[0][0], &g_count, &ipiv2[0], &right_vector[0], &g_count, &info);
		if(info==0){
		// ok
		for(int ig=0; ig<g_count; ig++){
		printf("(%10.2e, %10.2e)\n", right_vector[ig].real(), right_vector[ig].imag());
		}
		}else{
		write_log((char*)"zgesv failed");
		return;
		}*/

	
	bool* nonzero_edge=new bool[g_count];
	int nonzero_edge_count=0;
	int* nonzero_edge_index=new int[g_count];
		
	for(int ig=0; ig<g_count; ig++){
		for(int p=0; p<3; p++){
			kpg[p]=k_para[p]+g_vec[ig][p];
		}
		double Ekpg=inner_product(kpg, kpg)/2.0;
		nonzero_edge[ig]=Ekpg<Ekin;
		// printf("nonzero_edge[%2d]=%s\n", ig, nonzero_edge[ig]?"true":"false");
		if(nonzero_edge[ig]){
			nonzero_edge_index[nonzero_edge_count]=ig;
			nonzero_edge_count++;
		}
	}
		
	// by zgels
	complex<double>** A_complex=alloc_zmatrix(nonzero_edge_count, g_count);
	for(int ig=0; ig<g_count; ig++){
		for(int iv=0; iv<nonzero_edge_count; iv++){
			int ig2=nonzero_edge_index[iv];
			A_complex[iv][ig]=H_RL[ig2][ig];
		}
	}
	complex<double>* B_complex=new complex<double>[g_count];
	for(int ig=0; ig<g_count; ig++){
		B_complex[ig]=right_vector[ig];
	}

	int lwork=-1;
	complex<double> work_dummy;
	nrhs=1;
	zgels_(&notrans, &g_count, &nonzero_edge_count, &nrhs, &A_complex[0][0], &g_count, &B_complex[0], &g_count, &work_dummy, &lwork, &info);
	if(info!=0){
		write_log((char*)"Work space calculation failed");
		return -1;
	}
	lwork=round(abs(work_dummy));
	complex<double>* work=new complex<double>[lwork];
	zgels_(&notrans, &g_count, &nonzero_edge_count, &nrhs, &A_complex[0][0], &g_count, &B_complex[0], &g_count, &work[0], &lwork, &info);
	if(info!=0){
		write_log((char*)"zgels failed");
		return -1;
	}
	delete[] work;
	// by zgels finished
		
	for(int ig=0; ig<g_count; ig++){
		right_vector[ig]=complex<double>(0.0, 0.0);
	}
	for(int iv=0; iv<nonzero_edge_count; iv++){
		right_vector[nonzero_edge_index[iv]]=B_complex[iv];
	}
	for(int ig=0; ig<g_count; ig++){
		loc_edge[ig]=right_vector[ig];
	}
		
	//printf("Right(solution) by zgels\n");
	//for(int iv=0; iv<nonzero_edge_count; iv++){
	//	printf("%10.6f %10.6f  ", B_complex[iv].real(), B_complex[iv].imag());
	//}
	//printf("\n");

	double zgels_norm=0.0;
	for(int iv=nonzero_edge_count; iv<g_count; iv++){
		zgels_norm+=norm(B_complex[iv]);
	}
	sprintf(sprintf_buffer2, "zgels norm: %10.2e", zgels_norm);
	write_log(sprintf_buffer2);

	// obtain f_G(z)
	complex<double>* right_vector_all=new complex<double>[eq_dim];
	for(int ig=0; ig<g_count; ig++){
		int index=ig*z_count_new;
		right_vector_all[index]=-right_vector[ig]/dz/dz;
		if(ig==V00_index){
			index=(ig+1)*z_count_new-1;
			right_vector_all[index]=-1.0/dz/dz*complex<double>(cos(kzz0h), sin(kzz0h));
		}
	}
		
	// beta=complex<double>(0.0, 0.0);
	//for(int igz=0; igz<eq_dim; igz++){
	//printf("(%8.4f, %8.4f)\n", right_vector_all[igz].real(), right_vector_all[igz].imag());
	//}
	// zgemv_(&notrans, &eq_dim, &eq_dim, &alpha, &left_matrix[0][0], &eq_dim, &right_vector_all[0], &inc, &beta, &FP_loc_buffer[0], &inc);
	nrhs=1;
	zgetrs_(&notrans, &eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &right_vector_all[0], &eq_dim, &info);
	if(info!=0){
		write_log((char*)"zgetrs failed");
	}

	for(int ig=0; ig<g_count; ig++){
		for(int iz=z_start; iz<z_count; iz++){
			int index_sol=ig*z_count_new+iz-z_start;
			FP_loc[ig][iz]=right_vector_all[index_sol];
		}
	}
		
	//for(int ig=0; ig<g_count; ig++){
	//	for(int iz=0; iz<z_count; iz++){
	//		int index=ig*z_count+iz;
	//printf("(%8.4f, %8.4f)\n", FP_loc_buffer[index].real(), FP_loc_buffer[index].imag());
	//	}
	//printf("\n");
	//}
	delete[] right_vector_all;
		
	//delete_zmatrix(left_matrix);
	//delete_zmatrix(right_matrix);
	delete[] ipiv;
	delete[] right_vector;
	delete[] right_vector2;
		
	delete[] sprintf_buffer2;
	//delete_zpmatrix(Vgg);
	//delete_dmatrix(A_real);
	delete_zmatrix(A_complex);
	delete_zmatrix(H_RR);
	delete_zmatrix(H_RL);
	delete[] nonzero_edge;
	delete[] nonzero_edge_index;
	//delete[] B_real;
	//delete[] x_real;
	//delete[] nabla;
	//delete[] Anabla;
	delete[] B_complex;

	return zgels_norm;
}

double solve_final_state_from_bulk(double Ekin, double* k_para, double kz, int g_count, int z_count, int bulk_count, double dz, int z_start, int V00_index, complex<double>*** Vgg, double** g_vec, complex<double>*** bulk_z, complex<double>** FP_loc, complex<double>** left_matrix, complex<double>** right_matrix, complex<double>* bulk_coefs){
	char* sprintf_buffer2=new char[Log_length+1];
	
	// composite matrix
	//complex<double>*** bulk_z=alloc_zpmatrix(bulk_count, g_count);
	//for(int in=0; in<bulk_count; in++){
	//	for(int ig=0; ig<g_count; ig++){
	//		bulk_z[in][ig]=&bulk_z_buffer[in*g_count*z_count+ig*z_count];
	//	}
	//}
	/*
	complex<double>*** Vgg=alloc_zpmatrix(g_count);
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			Vgg[ig1][ig2]=Vgg_buffer[ig1*g_count+ig2];
		}
		}*/
	//double g_vec[g_count][3];
	//for(int ig=0; ig<g_count; ig++){
	//	for(int ix=0; ix<3; ix++){
	//g_vec[ig][ix]=g_vec_buffer[ig*3+ix];
	//	}
	//}

	double kzz0=kz*dz*z_start;
	double kzz0mh=kz*dz*(z_start-1);

	double kzz1=kz*dz*(z_count-1);
	double kzz1ph=kz*dz*z_count;

	int z_count_new=z_count-z_start;
	int eq_dim=z_count_new*g_count;

	//complex<double>** left_matrix=alloc_zmatrix(eq_dim); // transposed
	// buffer reset
	for(int i=0; i<eq_dim; i++){
		for(int j=0; j<eq_dim; j++){
			left_matrix[j][i]=complex<double>(0.0, 0.0);
		}
	}
	
	double kpg[3];

	// left matrix
	/// M
	for(int ig1=0; ig1<g_count; ig1++){
		for(int p=0; p<3; p++){
			kpg[p]=k_para[p]+g_vec[ig1][p];
		}
		for(int ig2=0; ig2<g_count; ig2++){
			for(int iz=z_start; iz<z_count; iz++){
				int index1=ig1*z_count_new+(iz-z_start);
				int index2=ig2*z_count_new+(iz-z_start);
				complex<double> value=-2.0*Vgg[ig1][ig2][iz];
				if(ig1==ig2){
					value+=2.0*Ekin-inner_product(kpg, kpg);
				}
				left_matrix[index2][index1]+=value;
			}
		}
	}
	/// A
	for(int ig=0; ig<g_count; ig++){
		for(int iz1=z_start; iz1<z_count; iz1++){
			for(int iz2=z_start; iz2<z_count; iz2++){
				int index1=ig*z_count_new+(iz1-z_start);
				int index2=ig*z_count_new+(iz2-z_start);
				if(iz1==iz2){
					left_matrix[index2][index1]-=2.0/dz/dz;
				}
				if(iz1-1==iz2 || iz1+1==iz2){
					left_matrix[index2][index1]+=1.0/dz/dz;
				}
			}
		}
	}
	
	int info;
	int* ipiv=new int[eq_dim];
	int nrhs=2*g_count;
	
	//complex<double>** right_matrix=alloc_zmatrix(2*g_count, eq_dim); //transposed
	// buffer reset
	for(int igz1=0; igz1<eq_dim; igz1++){
		for(int ig2=0; ig2<2*g_count; ig2++){
			right_matrix[ig2][igz1]=complex<double>(0.0, 0.0);
		}
	}
	// 0~g: for LL and RL, g~2g: for LR and RR
	for(int ig1=0; ig1<g_count; ig1++){
		int index1=ig1*z_count_new;
		int index2=ig1;
		right_matrix[index2][index1]=complex<double>(1.0, 0.0);
		index1=(ig1+1)*z_count_new-1;
		index2=ig1+g_count;
		right_matrix[index2][index1]=complex<double>(1.0, 0.0);
	}
	
	zgesv_(&eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &right_matrix[0][0], &eq_dim, &info);
	if(info!=0){
		write_log((char*)"zgesv failed");
		return -1;
	}

	// extract the RL and RR matrices
	complex<double>** H_RL=alloc_zmatrix(g_count); // inverted
	complex<double>** H_RR=alloc_zmatrix(g_count); // inverted
	complex<double>** H_LL=alloc_zmatrix(g_count); // inverted
	complex<double>** H_LR=alloc_zmatrix(g_count); // inverted
	int index1, index2;
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			index1=(ig1+1)*z_count_new-1;
			index2=ig2;
			H_RL[ig2][ig1]=right_matrix[index2][index1];

			index1=ig1*z_count_new;
			H_LL[ig2][ig1]=right_matrix[index2][index1];
				
			index1=(ig1+1)*z_count_new-1;
			index2=ig2+g_count;
			H_RR[ig2][ig1]=right_matrix[index2][index1];

			index1=ig1*z_count_new;
			H_LR[ig2][ig1]=right_matrix[index2][index1];
		}
	}

	/*
		printf("H_RL\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			printf("%10.2e ", abs(H_RL[ig2][ig1]));
		}
		printf("\n");
		}
		
		printf("H_RR\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			printf("%10.2e ", abs(H_RR[ig2][ig1]));
		}
		printf("\n");
		}

				printf("H_LL\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			printf("%10.2e ", abs(H_LL[ig2][ig1]));
		}
		printf("\n");
		}
		
		printf("H_LR\n");
		for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			printf("%10.2e ", abs(H_LR[ig2][ig1]));
		}
		printf("\n");
		}*/

	complex<double>* right_vector=new complex<double>[g_count*2];
	complex<double>* right_vector2=new complex<double>[g_count*2];

	char notrans='N';
	complex<double> alpha(1.0, 0.0);
	complex<double> beta(1.0, 0.0);
	int inc=1;
  
	// compute the right vector
	// upper part (left)
	for(int ig=0; ig<g_count; ig++){
		right_vector[ig]=complex<double>(0.0, 0.0);
	}
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			right_vector2[ig]=-complex<double>(cos(kzz1ph), sin(kzz1ph));
		}else{
			right_vector2[ig]=complex<double>(0.0, 0.0);
		}
	}
	zgemv_(&notrans, &g_count, &g_count, &alpha, &H_LR[0][0], &g_count, &right_vector2[0], &inc, &beta, &right_vector[0], &inc);
	
	// lower part (right)
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			right_vector[ig+g_count]=-dz*dz*complex<double>(cos(kzz1), sin(kzz1));
		}else{
			right_vector[ig+g_count]=complex<double>(0.0, 0.0);
		}
	}
	for(int ig=0; ig<g_count; ig++){
		if(ig==V00_index){
			right_vector2[ig+g_count]=-complex<double>(cos(kzz1ph), sin(kzz1ph));
		}else{
			right_vector2[ig+g_count]=complex<double>(0.0, 0.0);
		}
	}
	zgemv_(&notrans, &g_count, &g_count, &alpha, &H_RR[0][0], &g_count, &right_vector2[g_count], &inc, &beta, &right_vector[g_count], &inc);

	/*
	printf("Right vector\n");
	for(int ig=0; ig<2*g_count; ig++){
		printf("%10.2e %10.2e\n", right_vector[ig].real(), right_vector[ig].imag());
		}*/

	// left matrix
	// upper (left)
	complex<double>** left_matrix_upper=alloc_zmatrix(bulk_count, g_count); // transposed
	complex<double>** zgemm_buffer=alloc_zmatrix(bulk_count, g_count); // transposed
	// fill h^2 B(z0)
	for(int in=0; in<bulk_count; in++){
		for(int ig=0; ig<g_count; ig++){
			left_matrix_upper[in][ig]=dz*dz*bulk_z[in][ig][z_start];
			//printf("%10.2e ", abs(left_matrix_upper[in][ig]));
		}
		//printf("\n");
	}
	// add H_LL*B(z0mh)
	for(int in=0; in<bulk_count; in++){
		for(int ig=0; ig<g_count; ig++){
		  zgemm_buffer[in][ig]=bulk_z[in][ig][z_start-1];
		}
	}
	zgemm_(&notrans, &notrans, &g_count, &bulk_count, &g_count, &alpha, &H_LL[0][0], &g_count, &zgemm_buffer[0][0], &g_count, &beta, &left_matrix_upper[0][0], &g_count);
	
	// lower (right)
	complex<double>** left_matrix_lower=alloc_zmatrix(bulk_count, g_count); // transposed
	beta=complex<double>(0.0, 0.0);
	zgemm_(&notrans, &notrans, &g_count, &bulk_count, &g_count, &alpha, &H_RL[0][0], &g_count, &zgemm_buffer[0][0], &g_count, &beta, &left_matrix_lower[0][0], &g_count);

	// combine left_matrix_lower and left_matrix_upper
	complex<double>** left_matrix_all=alloc_zmatrix(bulk_count, g_count*2); // transposed
	//complex<double>** left_matrix_all=alloc_zmatrix(bulk_count, g_count); // transposed
	for(int in=0; in<bulk_count; in++){
		for(int ig=0; ig<g_count; ig++){
			left_matrix_all[in][ig]=left_matrix_upper[in][ig];
			//left_matrix_all[in][ig]=left_matrix_lower[in][ig];
			left_matrix_all[in][ig+g_count]=left_matrix_lower[in][ig];
		}
	}
	/*
	for(int ig=0; ig<g_count; ig++){
		for(int in=0; in<bulk_count; in++){
			printf("(%10.2e %10.2e) ", left_matrix_all[in][ig].real(), left_matrix_all[in][ig].imag());
		}
		printf("\n");
		}*/

	// put high weight
	/*
	for(int ig=0; ig<g_count; ig++){
		for(int in=0; in<bulk_count; in++){
			left_matrix_all[in][ig+g_count]*=100;
		}
		right_vector[ig+g_count]*=100;
		}*/
	
	int lwork=-1;
	complex<double> work_dummy;
	nrhs=1;
	int numRows=g_count*2;
	
	zgels_(&notrans, &numRows, &bulk_count, &nrhs, &left_matrix_all[0][0], &numRows, &right_vector[0], &numRows, &work_dummy, &lwork, &info);
	//zgels_(&notrans, &g_count, &bulk_count, &nrhs, &left_matrix_all[0][0], &g_count, &right_vector[g_count], &numRows, &work_dummy, &lwork, &info);
	if(info!=0){
		write_log((char*)"Work space calculation failed");
		return -1;
	}
	lwork=round(abs(work_dummy));
	complex<double>* work=new complex<double>[lwork];
	zgels_(&notrans, &numRows, &bulk_count, &nrhs, &left_matrix_all[0][0], &numRows, &right_vector[0], &numRows, &work[0], &lwork, &info);
	//zgels_(&notrans, &g_count, &bulk_count, &nrhs, &left_matrix_all[0][0], &g_count, &right_vector[g_count], &numRows, &work[0], &lwork, &info);
	if(info!=0){
		write_log((char*)"zgels failed");
		return -1;
	}
	delete[] work;
	
	double zgels_norm=0.0;
	for(int iv=bulk_count; iv<numRows; iv++){
	//for(int iv=bulk_count+g_count; iv<numRows; iv++){
		zgels_norm+=norm(right_vector[iv]);
	}
	sprintf(sprintf_buffer2, "zgels norm: %10.2e", zgels_norm);
	write_log(sprintf_buffer2);

	/*
	printf("Solution by zgels\n");
	for(int in=0; in<bulk_count; in++){
		printf("%10.2e %10.2e\n", right_vector[in].real(), right_vector[in].imag());
		}*/

	// export the solution
	for(int in=0; in<bulk_count; in++){
		bulk_coefs[in]=right_vector[in];
		//bulk_coefs[in]=right_vector[in+g_count];
	}

	complex<double>* left_edge=new complex<double>[g_count];
	beta=complex<double>(0.0, 0.0);

	zgemv_(&notrans, &g_count, &bulk_count, &alpha, &zgemm_buffer[0][0], &g_count, &bulk_coefs[0], &inc, &beta, &left_edge[0], &inc);
	/*
	for(int ig=0; ig<g_count; ig++){
	printf("%8.4f %8.4f = %8.4f %8.4f\n", g_vec[ig][0], g_vec[ig][1], left_edge[ig].real(), left_edge[ig].imag());
	}*/

	complex<double>* right_vector_sol=new complex<double>[eq_dim];
	for(int ig=0; ig<g_count; ig++){
		for(int iz=z_start; iz<z_count; iz++){
			int index=ig*z_count_new+(iz-z_start);
			if(iz==z_start){
				right_vector_sol[index]=-1.0/dz/dz*left_edge[ig];
			}else	if(iz==z_count-1 && ig==V00_index){
				right_vector_sol[index]=-1.0/dz/dz*complex<double>(cos(kzz1ph), sin(kzz1ph));
			}else{
				right_vector_sol[index]=complex<double>(0.0, 0.0);
			}
		}
	}
	zgetrs_(&notrans, &eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &right_vector_sol[0], &eq_dim, &info);
	if(info!=0){
		write_log((char*)"zgetrs failed");
		return -1;
	}
	
	for(int ig=0; ig<g_count; ig++){
		for(int iz=0; iz<=z_start; iz++){
			for(int in=0; in<bulk_count; in++){
				FP_loc[ig][iz]+=bulk_coefs[in]*bulk_z[in][ig][iz];
			}
		}
		//printf("%8.4f %8.4f\n", FP_loc_buffer[ig*z_count+z_start].real(), FP_loc_buffer[ig*z_count+z_start].imag());
		for(int iz=z_start; iz<z_count; iz++){
			int index_sol=ig*z_count_new+iz-z_start;
			FP_loc[ig][iz]=right_vector_sol[index_sol];
		}
		//printf("%8.4f %8.4f\n", FP_loc_buffer[ig*z_count+z_start].real(), FP_loc_buffer[ig*z_count+z_start].imag());
		//printf("\n");
	}

	/*
	for(int ig=0; ig<g_count; ig++){
		for(int iz=z_start-84; iz<z_count; iz++){
			int index=ig*z_count+iz;
			printf("%d %8.4f %8.4f\n", iz, FP_loc_buffer[index].real(), FP_loc_buffer[index].imag());
		}
		printf("\n");
		}*/
	
	//delete_zpmatrix(bulk_z);
	//delete_zpmatrix(Vgg);
	delete_zmatrix(H_RL);
	delete_zmatrix(H_RR);
	delete_zmatrix(H_LL);
	delete_zmatrix(H_LR);
	delete_zmatrix(left_matrix_upper);
	delete_zmatrix(left_matrix_lower);
	delete_zmatrix(zgemm_buffer);
	delete[] left_edge;
	delete[] sprintf_buffer2;
	delete[] right_vector;
	delete[] right_vector2;
	delete[] right_vector_sol;
	return zgels_norm;
}

complex<double> interpolate_fgz(double z, complex<double>* fgz, double dz, int z_count){
	double index_d=z/dz;
	int index_floor=floor(z/dz);
	if(index_floor<0){
		write_log((char*)"Out of range");
		return fgz[0];
	}
	if(index_floor>=z_count-1){
		write_log((char*)"Out of range");
		return fgz[z_count-1];
	}
	return fgz[index_floor]*(index_floor+1.0-index_d)+fgz[index_floor+1]*(index_d-index_floor);
}

double spherical_harmonic_theta(int l, int m, double theta){
	double sinT=sin(theta);
	double cosT=cos(theta);
	
	// s(0): 1/sqrt(2)
	if(l==0){
		if(m==0){
			return 1.0/(sqrt(2.0));
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
	if(l==1){
		if(m==-1){
			// p(-1): 1/2 sqrt(3) sin(T)
			return (1.0/2.0)*sqrt(3.0)*sinT;
		}else if(m==0){
			// p(0): sqrt(3/2) cos(T)
			return sqrt(3.0/2.0)*cosT;
		}else if(m==1){
			// p(1): -1/2 sqrt(3/2pi) sin(T)
			return -(1.0/2.0)*sqrt(3.0)*sinT;
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
	if(l==2){
		if(m==-2){
			// d(-2): 1/4 sqrt(15) sin^2(T)
		  return (1.0/4.0)*sqrt(15.0)*sinT*sinT;
		}else if(m==-1){
			// d(-1): 1/2 sqrt(15) sin(T)cos(T)
			return (1.0/2.0)*sqrt(15.0)*sinT*cosT;
		}else if(m==0){
			// d(0): 1/2 sqrt(5/2) (3cos^2(T)-1)
			return (1.0/2.0)*sqrt(5.0/2.0)*(3*cosT*cosT-1.0);
		}else if(m==1){
			// d(1): -1/2 sqrt(15) sin(T)cos(T)
			return -(1.0/2.0)*sqrt(15.0)*sinT*cosT;
		}else if(m==2){
			// d(2): 1/4 sqrt(15) sin^2(T)
			return (1.0/4.0)*sqrt(15.0)*sinT*sinT;
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
	if(l==3){
		if(m==-3){
			// f(-3): 1/8 sqrt(70) sin^3(T)
			return (1.0/8.0)*sqrt(70.0)*sinT*sinT*sinT;
		}else if(m==-2){
			// f(-2): 1/4 sqrt(105) sin^2(T)cos(T)
			return (1.0/4.0)*sqrt(105.0)*sinT*sinT*cosT;
		}else if(m==-1){
			// f(-1): 1/8 sqrt(42) (5cos^2(T)-1)sin(T)
			return (1.0/8.0)*sqrt(42.0)*(5.0*cosT*cosT-1.0)*sinT;
		}else if(m==0){
			// f(0): 1/2 sqrt(7/2) (5cos^2(T)-3)cos(T)
			return (1.0/2.0)*sqrt(7.0/2.0)*(5.0*cosT*cosT-3.0)*cosT;
		}else if(m==1){
			// f(1): -1/8 sqrt(42) (5cos^2(T)-1)sin(T)
			return -(1.0/8.0)*sqrt(42.0)*(5.0*cosT*cosT-1.0)*sinT;
		}else if(m==2){
			// f(2): 1/4 sqrt(105) sin^2(T)cos(T)
			return (1.0/4.0)*sqrt(105.0)*sinT*sinT*cosT;
		}else if(m==3){
			// f(3): -1/8 sqrt(70) sin^3(T)
			return -(1.0/8.0)*sqrt(70.0)*sinT*sinT*sinT;
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
	if(l==4){
		if(m==-4){
			// g(-4): 3/16 sqrt(35) sin^4(T)
			return (3.0/16.0)*sqrt(35.0)*sinT*sinT*sinT*sinT;
		}else if(m==-3){
			// g(-3): 3/8 sqrt(35/pi) sin^3(T)cos(T)
			return (3.0/8.0)*sqrt(70.0)*sinT*sinT*sinT*cosT;
		}else if(m==-2){
			// g(-2): 3/8 sqrt(5) (7cos^2(T)-1)sin^2(T)
			return (3.0/8.0)*sqrt(5.0)*(7.0*cosT*cosT-1.0)*sinT*sinT;
		}else if(m==-1){			
			// g(-1): 3/8 sqrt(10) (7cos^2(T)-3)sin(T)cos(T)
			return (3.0/8.0)*sqrt(10.0)*(7.0*cosT*cosT-3.0)*sinT*cosT;
		}else if(m==0){
			// g(0): 3/16 sqrt(2) (35cos^4(T)-30cos^2(T)+3)
			return (3.0/16.0)*sqrt(2.0)*(35.0*cosT*cosT*cosT*cosT-30.0*cosT*cosT+3.0);
		}else if(m==1){
			// g(1): -3/8 sqrt(10) (7cos^2(T)-3)sin(T)cos(T)
			return -(3.0/8.0)*sqrt(10.0)*(7.0*cosT*cosT-3.0)*sinT*cosT;
		}else if(m==2){
			// g(2): 3/8 sqrt(5) (7cos^2(T)-1)sin^2(T)
			return (3.0/8.0)*sqrt(5.0)*(7.0*cosT*cosT-1.0)*sinT*sinT;
		}else if(m==3){
			// g(3): -3/8 sqrt(70) sin^3(T)cos(T)
			return -(3.0/8.0)*sqrt(70.0)*sinT*sinT*sinT*cosT;
		}else if(m==4){
			// g(4): 3/16 sqrt(35) sin^4(T)
			return (3.0/16.0)*sqrt(35.0)*sinT*sinT*sinT*sinT;
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
	if(l==5){
		if(m==-5){
			// h(-5): 3/32 sqrt(154) sin^5(T)
			return (3.0/32.0)*sqrt(154.0)*sinT*sinT*sinT*sinT*sinT;
		}else if(m==-4){
			// h(-4): 3/16 sqrt(385) sin^4(T)cos(T)
			return (3.0/16.0)*sqrt(385.0)*sinT*sinT*sinT*sinT*cosT;
		}else if(m==-3){
			// h(-3): 1/32 sqrt(770) sin^3(T)(9cos^2(T)-1)
			return (1.0/32.0)*sqrt(770)*sinT*sinT*sinT*(9*cosT*cosT-1.0);
		}else if(m==-2){
			// h(-2): 1/8 sqrt(1155) sin^2(T)(3cos^3(T)-cos(T))
			return (1.0/8.0)*sqrt(1155.0)*sinT*sinT*(3*cosT*cosT*cosT-cosT);
		}else if(m==-1){
			// h(-1): 1/16 sqrt(165) sin(T)(21cos^4(T)-14cos^2(T)+1)
			return (1.0/16.0)*sqrt(165.0)*sinT*(21*cosT*cosT*cosT*cosT-14*cosT*cosT+1.0);
		}else if(m==0){
			// h(0): 1/16 sqrt(22) (63cos^5(T)-70cos^3(T)+15cos(T))
			return (1.0/16.0)*sqrt(22.0)*(63.0*cosT*cosT*cosT*cosT*cosT-70*cosT*cosT*cosT+15.0*cosT);
		}else if(m==1){
			// h(1): -1/16 sqrt(165) sin(T)(21cos^4(T)-14cos^2(T)+1)
			return -(1.0/16.0)*sqrt(165.0)*sinT*(21*cosT*cosT*cosT*cosT-14*cosT*cosT+1.0);
		}else if(m==2){
			// h(2): 1/8 sqrt(1155) sin^2(T)(3cos^3(T)-cos(T))
			return (1.0/8.0)*sqrt(1155.0)*sinT*sinT*(3*cosT*cosT*cosT-cosT);
		}else if(m==3){
			// h(3): -1/32 sqrt(770) sin^3(T)(9cos^2(T)-1)
			return -(1.0/32.0)*sqrt(770)*sinT*sinT*sinT*(9*cosT*cosT-1.0);
		}else if(m==4){
			// h(4): 3/16 sqrt(385) sin^4(T)cos(T)(cos(4P)+i*sin(4P))
			return (3.0/16.0)*sqrt(385.0)*sinT*sinT*sinT*sinT*cosT;
		}else if(m==5){
			// h(5): -3/32 sqrt(154) sin^5(T)(cos(5P)-i*sin(5P))
			return -(3.0/32.0)*sqrt(154)*sinT*sinT*sinT*sinT*sinT;
		}else{
			write_log((char*)"Invalid m");
			return 0.0;
		}
	}
  
	write_log((char*)"Invalid l");
	return 0.0;
}

void psi_normalize(int count, double* dr, double* psi);
double expectation_value(int count, double* matrix, double* psi, double* dr);
double matrix_element(int count, double* left, double* center, double* right);

void solve_nonlocal_wfn(double Ekin, int l, int r_count, double* r_arr, double* V_loc, int VPS_l_length, int* VPS_l, double* VPS_E_buffer, double** VPS_nonloc_buffer, double** psi){
	char* sprintf_buffer2=new char[Log_length+1];
	int N=2;
	// printf("Ekin %f, l %d\n", Ekin, l);
	double** VPS_nonloc=new double*[VPS_l_length];
	double* VPS_E=new double[VPS_l_length];
	int l_count=0;
	for(int il=0; il<VPS_l_length; il++){
		if(VPS_l[il]==l){
			VPS_nonloc[l_count]=&VPS_nonloc_buffer[il][0];
			VPS_E[l_count]=VPS_E_buffer[il];
			//printf("VPS_E[%d] %f\n", count, VPS_E[count]);
			l_count++;
		}
	}
	/*
		for(int ir=0; ir<r_count; ir++){
		printf("%8.4f ", r_arr[ir]);
		printf("%8.4f ", V_loc[ir]);
		for(int il=0; il<l_count; il++){
		printf("%13.4e ", VPS_nonloc[il][ir]);
		}
		printf("\n");
		}*/

	// delta r
	double dr[r_count];
	for(int ix=0; ix<r_count; ix++){
		if(ix==0){
			dr[ix]=r_arr[ix];
		}else{
			dr[ix]=r_arr[ix]-r_arr[ix-1];
		}
	}
	// delta x
	double dx=(log(r_arr[r_count-1])-log(r_arr[0]))/(r_count-1);
	double sum_dx2=0.0;
	for(int ir=0; ir<r_count-1; ir++){
		double dxi=log(r_arr[ir+1])-log(r_arr[ir]);
		sum_dx2+=dxi*dxi;
	}
	double ave_dx2=sum_dx2/(r_count-1);
	double ddx=ave_dx2-dx*dx;
	// printf("dx = %13.4e, V[dx] = %13.4e\n", dx, ddx);

	// prepare the left matrix
	double** matrix=alloc_dmatrix(r_count);
	for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
			matrix[iy][ix]=0.0;
			// (-1) * second derivative
			/*if(ix==0){
					if(iy==0 || iy==2){
					matrix[iy][ix]-=1.0/dx/dx;
					}else if(iy==1){
					matrix[iy][ix]-=-2.0/dx/dx;
					}
				if(iy==-1 || iy==1){
					matrix[iy][ix]-=1.0/dx/dx;
				}else if(iy==0){
					matrix[iy][ix]-=-2.0/dx/dx;
				}
			}else if(ix==r_count-1){
				if(iy==r_count-3 || iy==r_count-1){
					matrix[iy][ix]-=1.0/dx/dx;
				}else if(iy==r_count-2){
					matrix[iy][ix]-=-2.0/dx/dx;
				}
				}else{*/
			if(iy==ix-1 || iy==ix+1){
				matrix[iy][ix]-=1.0/dx/dx;
			}else if(iy==ix){
				matrix[iy][ix]-=-2.0/dx/dx;
			}
				//}
			// first derivative
			/*if(ix==0){
									if(iy==0){
									matrix[iy][ix]+=-1.0/dx;
									}else if(iy==1){
									matrix[iy][ix]+=1.0/dx;
									}
				if(iy==-1){
					matrix[iy][ix]+=-0.5/dx;
				}else if(iy==1){
					matrix[iy][ix]+=0.5/dx;
				}
			}else if(ix==r_count-1){
				if(iy==r_count-2){
					matrix[iy][ix]+=-1.0/dx;
				}else if(iy==r_count-1){
					matrix[iy][ix]+=1.0/dx;
				}
			}else{*/
			if(iy==ix-1){
				matrix[iy][ix]+=-0.5/dx;
			}else if(iy==ix+1){
				matrix[iy][ix]+=0.5/dx;
			}
				//}
			// l(l+1)+2r^2(V-E)
			if(iy==ix){
				matrix[iy][ix]+=l*(l+1)*1.0+2.0*r_arr[ix]*r_arr[ix]*(V_loc[ix]-Ekin);
			}
			// 2r^2 * projector
			for(int il=0; il<l_count; il++){
				matrix[iy][ix]+=2.0*r_arr[ix]*r_arr[ix]*VPS_E[il]*VPS_nonloc[il][ix]
					*VPS_nonloc[il][iy]*dr[iy];
			}
		}
	}

	// debug
	/*
		for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
		printf("%8.1e ", matrix[iy][ix]);
		}
		printf("\n");
		}*/

	// test
	/*
		for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
		matrix[iy][ix]=(ix==iy)?((ix+0.1)/10000.0-0.1):0.0;
		}
		}*/

	// copy
	//double** inverse=alloc_dmatrix(r_count);
	//for(int ix=0; ix<r_count; ix++){
	//	for(int iy=0; iy<r_count; iy++){
	//		inverse[iy][ix]=matrix[iy][ix];
	//	}
	//}
	
	//double** right_vector=alloc_dmatrix(2, r_count);
	for(int ir=0; ir<r_count; ir++){
		if(ir==r_count-1){
			psi[0][ir]=1.0/dx/dx-0.5/dx;
			psi[1][ir]=1.0/dx/dx-0.5/dx;
		}else if(ir==0){
			psi[0][ir]=0.0;
			psi[1][ir]=1.0/dx/dx+0.5/dx;
		}else{
			psi[0][ir]=0.0;
			psi[1][ir]=0.0;
		}
	}

	int nrhs=2;
	int ipiv[r_count];
	int info;
	dgesv_(&r_count, &nrhs, &matrix[0][0], &r_count, &ipiv[0], &psi[0][0], &r_count, &info);

	//for(int ir=0; ir<r_count; ir++){
	//	printf("%10.6f %10.6f %10.6f\n", r_arr[ir], psi[0][ir], psi[1][ir]);
	//}
	//printf("\n");

	/*
	// LU decomposition
	//int info;
	//int ipiv[r_count];
	dgetrf_(&r_count, &r_count, &inverse[0][0], &r_count, &ipiv[0], &info);
	if(info!=0){
	write_log((char*)"LU decomposition failed");
	return;
	}
	// work space calculation for inverse matrix
	int lwork=-1;
	double work_dummy;
	dgetri_(&r_count, &inverse[0][0], &r_count, &ipiv[0], &work_dummy, &lwork, &info);
	if(info!=0){
	write_log((char*)"Work space calculation failed");
	return;
	}
	lwork=round(work_dummy);
	// printf("Work space size: %d\n", lwork);
	double work[lwork];
	
	dgetri_(&r_count, &inverse[0][0], &r_count, &ipiv[0], &work[0], &lwork, &info);
	if(info!=0){
	write_log((char*)"Inverse matrix calculation failed");
	return;
	}
	// check
	//double mm[r_count][r_count];
	//double alpha=1.0;
	//double beta=0.0;
	//char trans='N';
	//dgemm_(&trans, &trans, &r_count, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &inverse[0][0], &r_count, &beta, &mm[0][0], &r_count);
	//for(int ix=0; ix<r_count; ix++){
	//for(int iy=0; iy<r_count; iy++){
	//printf("%4.1f ", mm[iy][ix]);
	//}
	//printf("\n");
	//}
	// initial vector
	double psi_vector[r_count];
	double hpsi_vector[r_count];
	for(int ir=0; ir<r_count; ir++){
	  //psi_vector[ir]=r_arr[ir];
	  psi_vector[ir]=1.0;
	}
	psi_normalize(r_count, dr, psi_vector);
	//for(int ir=0; ir<r_count; ir++){
	//printf("%10.6f %10.6f\n", r_arr[ir], psi_vector[ir]);
	//	}
	double expect=expectation_value(r_count, &matrix[0][0], psi_vector, dr);
	// printf("Expect: %f\n", expect);
	char no_trans='N';
	double alpha=1.0;
	int inc=1;
	double beta=0.0;
	int trial;
	for(trial=0; trial<100; trial++){
		// multiply H^-1
		dgemv_(&no_trans, &r_count, &r_count, &alpha, &inverse[0][0], &r_count, &psi_vector[0], &inc, &beta, &hpsi_vector[0], &inc);
		psi_normalize(r_count, dr, hpsi_vector);
		expect=expectation_value(r_count, &matrix[0][0], hpsi_vector, dr);
		// copy
		for(int ir=0; ir<r_count; ir++){
			psi_vector[ir]=hpsi_vector[ir];
		}
		if(abs(expect)<0.1 && trial>10){
			break;
		}
	}
	// printf("Power iteration trials: %3d, Expectation value: %8.4f\n", trial, expect);
	//double hpsi_norm=0.0;
	//dgemv_(&trans, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &hpsi_vector[0], &inc);
	//	for(int ir=0; ir<r_count; ir++){
	//hpsi_norm+=hpsi_vector[ir]*hpsi_vector[ir]*dr[ir];
	//}
	//printf("Norm: %f\n", hpsi_norm);
	// prepare H^t * diag(dr) * H
	double** dr_matrix=alloc_dmatrix(r_count);
	double** drH_matrix=alloc_dmatrix(r_count);
	double** A_matrix=alloc_dmatrix(r_count);
	for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
			if(ix==iy){
				dr_matrix[iy][ix]=dr[ix];
			}else{
				dr_matrix[iy][ix]=0.0;
			}
		}
	}
	// approx solution search by conjugate gradient method
	char trans='T';
	dgemm_(&no_trans, &no_trans, &r_count, &r_count, &r_count, &alpha, &dr_matrix[0][0], &r_count, &matrix[0][0], &r_count, &beta, &drH_matrix[0][0], &r_count);
	dgemm_(&trans, &no_trans, &r_count, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &drH_matrix[0][0], &r_count, &beta, &A_matrix[0][0], &r_count);
	// printf("Norm: %f\n", matrix_element(r_count, psi_vector, &A_matrix[0][0], psi_vector));
	double mv[r_count];
	double cg_p[r_count];
	double cg_r[r_count];
	double cg_alpha;
	double cg_beta;
	double cg_norm;
	double cg_norm_prev;
	dgemv_(&no_trans, &r_count, &r_count, &alpha, &A_matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &mv[0], &inc);
	for(int ir=0; ir<r_count; ir++){
		cg_p[ir]=-mv[ir];
		cg_r[ir]=-mv[ir];
	}
	cg_norm_prev=matrix_element(r_count, psi_vector, &A_matrix[0][0], psi_vector);
	for(trial=0; trial<100000; trial++){
		cg_alpha=-matrix_element(r_count, psi_vector, &A_matrix[0][0], &cg_p[0])/matrix_element(r_count, &cg_p[0], &A_matrix[0][0], &cg_p[0]);
		//printf("Alpha: %e\n", sd_alpha);
		for(int ir=0; ir<r_count; ir++){
			psi_vector[ir]+=cg_alpha*cg_p[ir];
		}
		dgemv_(&no_trans, &r_count, &r_count, &alpha, &A_matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &mv[0], &inc);
		for(int ir=0; ir<r_count; ir++){
			cg_r[ir]=-mv[ir];
		}
		cg_beta=-matrix_element(r_count, &cg_p[0], &A_matrix[0][0], &cg_r[0])/matrix_element(r_count, &cg_p[0], &A_matrix[0][0], &cg_p[0]);
		for(int ir=0; ir<r_count; ir++){
			cg_p[ir]=cg_r[ir]+cg_beta*cg_p[ir];
		}
		//printf("Norm: %f\n", matrix_element(r_count, psi_vector, &A_matrix[0][0], psi_vector));
		psi_normalize(r_count, dr, psi_vector);
		cg_norm=matrix_element(r_count, psi_vector, &A_matrix[0][0], psi_vector);
		if(cg_norm<0.01 && trial>10){
			break;
		}
		if(cg_norm>cg_norm_prev){
			// printf("!! norm increasing !!\n");
			// printf("trial: %d, cg_norm: %e\n", trial, cg_norm);
			dgemv_(&no_trans, &r_count, &r_count, &alpha, &A_matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &mv[0], &inc);
			for(int ir=0; ir<r_count; ir++){
				cg_p[ir]=-mv[ir];
				cg_r[ir]=-mv[ir];
			}
		}
		cg_norm_prev=cg_norm;
	}
	sprintf(sprintf_buffer2, "CG trials: %6d, CG norm: %10.4e", trial, cg_norm);
	write_log(sprintf_buffer2);
	dgemv_(&no_trans, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &hpsi_vector[0], &inc);*/
	// copy & normalization so that the edge is 1.0
	//double psi_edge=right_vector[0][r_count-1];
	//for(int ir=0; ir<r_count; ir++){
	//	psi[ir]=right_vector[0][ir]/psi_edge;
	//printf("%10.6f %10.6f %10.6f\n", r_arr[ir], psi[ir], hpsi_vector[ir]/psi_edge);
	//}
	delete[] VPS_E;
	delete[] VPS_nonloc;
	delete[] sprintf_buffer2;
	delete_dmatrix(matrix);
	//delete_dmatrix(right_vector);
	//delete_dmatrix(inverse);
	//delete_dmatrix(dr_matrix);
	//delete_dmatrix(drH_matrix);
	//delete_dmatrix(A_matrix);
}

void psi_normalize(int count, double* dr, double* psi){
	double norm=0.0;
	for(int i=0; i<count; i++){
		norm+=psi[i]*psi[i]*dr[i];
	}
	double coef=sqrt(norm);
	//printf("Coef: %e\n", coef);
	for(int i=0; i<count; i++){
		psi[i]/=coef;
	}
}

double expectation_value(int count, double* matrix, double* psi, double* dr){
	char trans='N';
	double alpha=1.0;
	int inc=1;
	double beta=0.0;
	double y[count];

	dgemv_(&trans, &count, &count, &alpha, matrix, &count, &psi[0], &inc, &beta, &y[0], &inc);
	double ret=0.0;
	for(int i=0; i<count; i++){
		ret+=psi[i]*y[i]*dr[i];
	}
	return ret;
}

double matrix_element(int count, double* left, double* center, double* right){
	double cr[count];
	char trans='N';
	double alpha=1.0;
	double beta=0.0;
	int inc=1;

	/*
		for(int i=0; i<count; i++){
		printf("%f %f\n", left[i], right[i]);
		}*/
	 
	dgemv_(&trans, &count, &count, &alpha, center, &count, right, &inc, &beta, &cr[0], &inc);
	return ddot_(&count, left, &inc, &cr[0], &inc);
}


void calc_bulk_potential(complex<double>* Vgg, double dz_slab, int slab_count, double z_offset, double bulk_height, double dz_bulk, int z_count_bulk, int bulk_count, complex<double>* Vgg_average, double* Vgg_stdev){
	for(int iz=0; iz<z_count_bulk; iz++){
		for(int ib=0; ib<bulk_count; ib++){
			double z=z_offset+ib*bulk_height+iz*dz_bulk;
			complex<double> V=interpolate_fgz(z, Vgg, dz_slab, slab_count);
			Vgg_average[iz]+=V;
			Vgg_stdev[iz]+=norm(V);
		}
		Vgg_average[iz]/=(1.0*bulk_count);
		Vgg_stdev[iz]=sqrt(Vgg_stdev[iz]/(1.0*bulk_count)-norm(Vgg_average[iz]));			
	}
	
}


complex<double> Fourier_expansion_1D(complex<double>* Vgg_bulk, double g, double dz, int count){
	complex<double> ret(0.0, 0.0);
	for(int iz=0; iz<count; iz++){
		double gz=g*dz*iz;
		ret+=Vgg_bulk[iz]*complex<double>(cos(gz), -sin(gz));
	}
	ret/=(count*1.0);
	return ret;
}

void prepare_matrix_bulk(int g_count, complex<double>** mat, double Ekin, double* k, double** g_vec, complex<double>** Vgg){
	double kpg[3];
	for(int i=0; i<g_count; i++){
		for(int j=0; j<g_count; j++){
			mat[j][i]=2.0*Vgg[i][j];
			if(i==j){
				mat[j][i]-=2.0*Ekin;
				for(int ix=0; ix<3; ix++){
					kpg[ix]=g_vec[i][ix]+k[ix];
				}
				//printf("kpg %.3f %.3f %.3f\n", kpg[0], kpg[1], kpg[2]);
				mat[j][i]+=inner_product(kpg, kpg);
			}
			//printf("(%4.1f %4.1f) ", mat[j][i].real(), mat[j][i].imag());
			
		}
		//printf("\n");
	}
}
void prepare_matrix_bulk_complex(int g_count, complex<double>** mat, double Ekin, double* k, complex<double> kz, double** g_vec, complex<double>** Vgg){
	double kpg[3];
	kpg[2]=0;
	for(int i=0; i<g_count; i++){
		for(int j=0; j<g_count; j++){
			mat[j][i]=2.0*Vgg[i][j];
			if(i==j){
				mat[j][i]-=2.0*Ekin;
				for(int ix=0; ix<2; ix++){
					kpg[ix]=g_vec[i][ix]+k[ix];
				}
				//printf("kpg %.3f %.3f %.3f\n", kpg[0], kpg[1], kpg[2]);
				mat[j][i]+=inner_product(kpg, kpg)+(kz+g_vec[i][2])*(kz+g_vec[i][2]);
			}
			//printf("(%4.1f %4.1f) ", mat[j][i].real(), mat[j][i].imag());
			
		}
		//printf("\n");
	}
}
// export_flag:
// 0: do nothing
// 1: calculate -> export to export_mat
// 2: import from export_mat, do not calculate
void add_nonlocal_term(int g_count, double** g_vec, complex<double>** mat, double* k_para, double kz_real,
											 int atom_length, int* VPS_l_length, int* vps_cutoff_index,
											 int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
											 double** atom_coordinates, int* atom_spec_index,
											 double FPFS_bulk_min, double FPFS_bulk_height, int export_flag, complex<double>** export_mat){
	if(export_flag==2){
		for(int ig1=0; ig1<g_count; ig1++){
			for(int ig2=0; ig2<g_count; ig2++){
				mat[ig2][ig1]+=export_mat[ig2][ig1];
			}
		}
		return;
	}
	complex<double> p1jlp[12]={
		complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1),
		complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1),
		complex<double>(1, 0), complex<double>(0, 1), complex<double>(-1, 0), complex<double>(0, -1)
	};
	double FPFS_bulk_max=FPFS_bulk_min+FPFS_bulk_height;
	
	complex<double>* Q=new complex<double>[g_count];
	double kpg[3];
	double atom_pos_bulk[3];
	double kpg_length;
	int inc=1;
	complex<double> Ylm_k[6][11];
	double* rj;

	if(export_flag==1){
		for(int ig1=0; ig1<g_count; ig1++){
			for(int ig2=0; ig2<g_count; ig2++){
				export_mat[ig2][ig1]=complex<double>(0, 0);
			}
		}
	}
	
	for(int ia=0; ia<atom_length; ia++){
		if(!(FPFS_bulk_min<=atom_coordinates[ia][2] && atom_coordinates[ia][2]<=FPFS_bulk_max)){
			continue;
		}
		atom_pos_bulk[0]=atom_coordinates[ia][0];
		atom_pos_bulk[1]=atom_coordinates[ia][1];
		atom_pos_bulk[2]=atom_coordinates[ia][2]-FPFS_bulk_min;
		
		int is=atom_spec_index[ia];
		rj=new double[vps_cutoff_index[is]];
		for(int il=0; il<VPS_l_length[is]; il++){
			// calculate Q = sum_m <Ylm bpl|phig>
			int l=VPS_l[is][il];
			for(int ig=0; ig<g_count; ig++){
				Q[ig]=complex<double>(0, 0);
				kpg[0]=k_para[0]+g_vec[ig][0];
				kpg[1]=k_para[1]+g_vec[ig][1];
				kpg[2]=kz_real  +g_vec[ig][2];
				kpg_length=sqrt(inner_product(kpg, kpg));
				spherical_harmonics(kpg, &Ylm_k[0][0]);
				for(int m=-l; m<=l; m++){
					Q[ig]+=conj(Ylm_k[l][l+m]);
				}
				double kpgt=inner_product(kpg, atom_pos_bulk);
				//printf("kpg %f %f %f\n", kpg[0], kpg[1], kpg[2]);
				for(int ir=0; ir<vps_cutoff_index[is]-1; ir++){
					rj[ir]=VPS_r[is][ir]*sp_bessel(l, kpg_length*VPS_r[is][ir])*(VPS_r[is][ir+1]-VPS_r[is][ir]);
					//printf("%f ", rj[ir]);
				}
				//printf("\n");
				rj[vps_cutoff_index[is]-1]=0.0;
				double integral=ddot_(&vps_cutoff_index[is], &VPS_nonloc_ave[is][il][0], &inc, &rj[0], &inc);
				Q[ig]*=4*M_PI*p1jlp[l]*complex<double>(cos(kpgt), sin(kpgt))*integral;
				//Q[ig]=integral;
				//printf("Q[%4d] = (%8.3f %8.3f)\n", ig, Q[ig].real(), Q[ig].imag());
			}
			// mat+=2*E*Q*Q/V
			for(int ig1=0; ig1<g_count; ig1++){
				for(int ig2=0; ig2<g_count; ig2++){
					complex<double> nonlocal_term=2.0*conj(Q[ig1])*Q[ig2]*VPS_E_ave[is][il];
					mat[ig2][ig1]+=nonlocal_term/PA_FPFS_bulk_volume;
					if(export_flag==1){
						export_mat[ig2][ig1]+=nonlocal_term/PA_FPFS_bulk_volume;
					}
					// printf("Add (%f %f)\n", nonlocal_term.real(), nonlocal_term.imag());
				}					
			}
		}
		delete[] rj;
	}
	delete[] Q;
}
double determinant(int g_count, complex<double>** mat, double Ekin, double* k, double** g_vec, complex<double>** Vgg,
									 int atom_length, int* VPS_l_length, int* vps_cutoff_index,
									 int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
									 double** atom_coordinates, int* atom_spec_index,
									 double FPFS_bulk_min, double FPFS_bulk_height, int export_flag, complex<double>** export_mat){
	prepare_matrix_bulk(g_count, mat, Ekin, k, g_vec, Vgg);
	if(PA_FPFS_bulk_include_nonlocal){
		add_nonlocal_term(g_count, g_vec, mat, k, k[2],
											atom_length, VPS_l_length, vps_cutoff_index,
											VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
											atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
											export_flag, export_mat);
	}
	int* ipiv=new int[g_count];
	int info;
	zgetrf_(&g_count, &g_count, &mat[0][0], &g_count, &ipiv[0], &info);
	if(info!=0){
		write_log((char*)"zgetrf failed");
		return -1;
	}
	complex<double> zdet=1.0;
	int i;
	for(i=0; i<g_count; i++){
		zdet*=mat[i][i];
		if(ipiv[i]!=i+1){
			zdet*=-1;
		}
	}
	//printf("%8.4f\n", zdet.imag());
	delete[] ipiv;
	return zdet.real();
}

double determinant_sign(int g_count, complex<double>** mat, double Ekin, double* k, double** g_vec, complex<double>** Vgg,
												int atom_length, int* VPS_l_length, int* vps_cutoff_index,
												int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
												double** atom_coordinates, int* atom_spec_index,
												double FPFS_bulk_min, double FPFS_bulk_height, int export_flag, complex<double>** export_mat){
	prepare_matrix_bulk(g_count, mat, Ekin, k, g_vec, Vgg);
	if(PA_FPFS_bulk_include_nonlocal){
		add_nonlocal_term(g_count, g_vec, mat, k, k[2],
											atom_length, VPS_l_length, vps_cutoff_index,
											VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
											atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
											export_flag, export_mat);
	}
	int* ipiv=new int[g_count];
	int info;
	zgetrf_(&g_count, &g_count, &mat[0][0], &g_count, &ipiv[0], &info);
	if(info!=0){
		write_log((char*)"zgetrf failed");
		return -1;
	}
	complex<double> zdet=1.0;
	int i;
	for(i=0; i<g_count; i++){
		zdet*=mat[i][i]/abs(mat[i][i]);
		if(ipiv[i]!=i+1){
			zdet*=-1;
		}
	}
	delete[] ipiv;
	return zdet.real();
}

complex<double> determinant_complex(int g_count, complex<double>** mat, double Ekin, double* k, complex<double> kz, double** g_vec, complex<double>** Vgg,
																		int atom_length, int* VPS_l_length, int* vps_cutoff_index,
																		int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
																		double** atom_coordinates, int* atom_spec_index,
																		double FPFS_bulk_min, double FPFS_bulk_height, int export_flag, complex<double>** export_mat){
	prepare_matrix_bulk_complex(g_count, mat, Ekin, k, kz, g_vec, Vgg);
	if(PA_FPFS_bulk_include_nonlocal){
		add_nonlocal_term(g_count, g_vec, mat, k, kz.real(),
											atom_length, VPS_l_length, vps_cutoff_index,
											VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
											atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
											export_flag, export_mat);
	}
	int* ipiv=new int[g_count];
	int info;
	zgetrf_(&g_count, &g_count, &mat[0][0], &g_count, &ipiv[0], &info);
	if(info!=0){
		write_log((char*)"zgetrf failed");
		return -1;
	}
	complex<double> zdet=1.0;
	int i;
	for(i=0; i<g_count; i++){
		zdet*=mat[i][i]/abs(mat[i][i]);
		if(ipiv[i]!=i+1){
			zdet*=-1;
		}
	}
	delete[] ipiv;
	return zdet;
}

bool encloseOrigin(complex<double> lb, complex<double> lt, complex<double> rt, complex<double> rb){
	/*
	complex<double> lblt=lt/lb;
	complex<double> ltrt=rt/lt;
	complex<double> rtrb=rb/rt;
	complex<double> rblb=lb/rb;
	return (lblt.imag()>0 && ltrt.imag()>0 && rtrb.imag()>0 && rblb.imag()>0) ||
	(lblt.imag()<0 && ltrt.imag()<0 && rtrb.imag()<0 && rblb.imag()<0);*/

	return !((lb.real()>0 && lt.real()>0 && rt.real()>0 && rb.real()>0) ||
					 (lb.real()<0 && lt.real()<0 && rt.real()<0 && rb.real()<0) ||
					 (lb.imag()>0 && lt.imag()>0 && rt.imag()>0 && rb.imag()>0) ||
					 (lb.imag()<0 && lt.imag()<0 && rt.imag()<0 && rb.imag()<0));
}

void sort_ascend(int count, complex<double>* array){
	complex<double>* array2=new complex<double>[count];
	bool* picked=new bool[count];
	int i, j;
	for(i=0; i<count; i++){
		picked[i]=false;
	}

	complex<double> v_min;
	int v_min_index;
	for(i=0; i<count; i++){
		v_min=-1;
		v_min_index=-1;
		for(j=0; j<count; j++){
			if(picked[j]){
				continue;
			}
			if(v_min_index<0 || v_min.real()>array[j].real()){
				v_min=array[j];
				v_min_index=j;
			}
		}
		array2[i]=v_min;
		picked[v_min_index]=true;
	}
	for(i=0; i<count; i++){
		array[i]=array2[i];
	}

	delete[] array2;
	delete[] picked;
}

int find_eigenstate(int g_count, int band_index, int band_index2, double Ekin, complex<double>* w, complex<double>** vr, int** g_index){
	// band_order
	// -3: initial value -> not used (complex eigenvalue)
	// -1: will be used
	// 0~eigen_count-1: order index
	int* band_order=new int[g_count];
	int i, j, t;
	for(i=0; i<g_count; i++){
		band_order[i]=-3;
	}
	int eigen_real_count=0;
	for(i=0; i<g_count; i++){
		if(abs(w[i].imag())*0.5<PA_FPFS_real_eigenvalue_criterion){
			band_order[i]=-1;
			eigen_real_count++;
		}
	}
	double eigen_min=-1;
	int eigen_min_index=-1;
	for(t=0; t<eigen_real_count; t++){
		eigen_min_index=-1;
		for(i=0; i<g_count; i++){
			if(band_order[i]==-1){
				double eigen=(0.5*w[i]).real();
				if(eigen_min_index<0 || eigen_min>eigen){
					eigen_min_index=i;
					eigen_min=eigen;
				}
			}
		}
		if(eigen_min_index==-1){
			break;
		}
		band_order[eigen_min_index]=t;
	}
	if(band_index2<0){
		for(i=0; i<g_count; i++){
			if(band_order[i]==band_index){
				return i;
			}
		}
	}
	int nearest_index=-1;
	double nearest_distance=-1;
	for(i=0; i<g_count; i++){
	  if((band_index>=band_order[i] && band_order[i]>=band_index2) ||
	     (band_index<=band_order[i] && band_order[i]<=band_index2)){
			double distance=abs(0.5*w[i]-Ekin);
			if(nearest_index<0 || nearest_distance>distance){
				nearest_index=i;
				nearest_distance=distance;
			}
		}
	}
	if(nearest_index<0){
	  write_log((char*)"Error eigenstate not found");
	}
	return nearest_index;
}

int solve_final_states_bulk(double Ekin, double* k_para, double gz, int g_count, int** g_index, double** g_vec, complex<double>** Vgg, int kz_count, double* dispersion_kz, int kappaz_count, int kappaz_border_index, int g_para_count, double* dispersion_kappaz, double* dispersion_mkappaz, double** dispersion_r,
														complex<double>** dispersion_c, int* dispersion_c_count, int** connection_c,
														complex<double>** dispersion_c_BZ, int* dispersion_c_BZ_count, int** connection_c_BZ,
														complex<double>** dispersion_mc, int* dispersion_mc_count, int** connection_mc,
														complex<double>** dispersion_mc_BZ, int* dispersion_mc_BZ_count, int** connection_mc_BZ,
														complex<double>*** final_states_pointer, double** kz_pointer, double** kappaz_pointer, complex<double>** mat, complex<double>** vr,
														int atom_length, int* VPS_l_length, int* vps_cutoff_index,
														int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
														double** atom_coordinates, int* atom_spec_index,
														double FPFS_bulk_min, double FPFS_bulk_height){
	// printf("Ekin %10.6f\n", Ekin);
	char* sprintf_buffer2=new char[Log_length+1];
	int i, j;
	// prepare matrix
	//double** g_vec=new double*[g_count];
	//for(i=0; i<g_count; i++){
	//	g_vec[i]=&g_vec_buffer[3*i];
		//}
	//complex<double>** Vgg=new complex<double>*[g_count];
	//for(i=0; i<g_count; i++){
	//	Vgg[i]=&Vgg_buffer[i*g_count];
	//}
	//double** dispersion=new double*[kz_count];
	//for(i=0; i<kz_count; i++){
	//	dispersion[i]=&dispersion_buffer[i*g_count];
	//}
	
	//complex<double>** mat=alloc_zmatrix(g_count); // transposed

	int buffer_size=PA_FPFS_bulk_buffer_size;
	double* solution_kz=new double[buffer_size];
	double* solution_kappaz=new double[buffer_size];
	
	double* solution_kz_alt=new double[buffer_size];
	double* solution_kappaz_alt=new double[buffer_size];
	bool* solution_alt_use=new bool[buffer_size];
	for(int is=0; is<buffer_size; is++){
		solution_alt_use[is]=false;
	}
	int* solution_band_indices=new int[buffer_size];
	int* solution_band_indices2=new int[buffer_size];

	double k_bloch[3];
	k_bloch[0]=k_para[0];
	k_bloch[1]=k_para[1];
	
	int solution_count_real=0;
	int solution_count_complex_sum=0;
	int solution_count_cspace=0;
	int solution_index=0;

	int cspace_count=0;
	bool* cspace_search_flag=new bool[kz_count];
	for(int ikz=0; ikz<kz_count; ikz++){
		cspace_search_flag[ikz]=false;
	}
	int* cspace_search_indices=new int[buffer_size];
	double* cspace_search_scale=new double[buffer_size];

	int lmax_count=0;
	bool* lmax_exist_flag=new bool[kz_count];
	for(int ikz=0; ikz<kz_count; ikz++){
		lmax_exist_flag[ikz]=false;
	}
	int* lmax_indices=new int[buffer_size];
	double* lmax_eigen=new double[buffer_size];
	int* lmax_band=new int[buffer_size];
	bool* lmax_listed=new bool[buffer_size];
	
	// real region: use dispersion_r
	for(int ib=0; ib<g_count; ib++){
		for(int ikz=0; ikz<kz_count; ikz++){
			// find local maxima just below EKin
			if(ib<g_count-1 && dispersion_r[ikz][ib]<Ekin /*&& Ekin<dispersion_r[ikz][ib+1]*/){
				int index_next=(ikz+1)%kz_count;
				int index_prev=(ikz+kz_count-1)%kz_count;
				double eigen_prev=dispersion_r[index_prev][ib];
				double eigen_curr=dispersion_r[ikz][ib];
				double eigen_next=dispersion_r[index_next][ib];
				if(eigen_prev < eigen_curr && eigen_curr > eigen_next){
					bool below_flag=false;
					bool local_minimum_flag=false;
					for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
						int index2_curr=(ikz2+kz_count)%kz_count;
						int index2_next=(ikz2+kz_count+1)%kz_count;
						int index2_prev=(ikz2+kz_count-1)%kz_count;
						if(Ekin < dispersion_r[index2_curr][ib+1]){
							below_flag=true;
						}
						if(dispersion_r[index2_prev][ib+1] > dispersion_r[index2_curr][ib+1] && dispersion_r[index2_curr][ib+1] < dispersion_r[index2_next][ib+1] && dispersion_r[index2_curr][ib+1]<Ekin){
							local_minimum_flag=true;
							break;
						}
					}
					if(below_flag && !local_minimum_flag){
						sprintf(sprintf_buffer2, "Local maximum below Ekin at kz=%8.4f", dispersion_kz[ikz]);
						write_log(sprintf_buffer2);
						for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
							lmax_exist_flag[(ikz2+kz_count)%kz_count]=true;
						}
						lmax_indices[lmax_count]=ikz;
						lmax_eigen[lmax_count]=eigen_curr;
						lmax_band[lmax_count]=ib;
						lmax_listed[lmax_count]=false;
						lmax_count++;
					}
				}
			}
			// local minima just above Ekin
			if(ib>0 && /*dispersion_r[ikz][ib-1]<Ekin &&*/ Ekin<dispersion_r[ikz][ib]){
				int index_next=(ikz+1)%kz_count;
				int index_prev=(ikz+kz_count-1)%kz_count;
				double eigen_prev=dispersion_r[index_prev][ib];
				double eigen_curr=dispersion_r[ikz][ib];
				double eigen_next=dispersion_r[index_next][ib];
				if(eigen_prev > eigen_curr && eigen_curr < eigen_next){
					bool below_flag=false;
					bool local_maximum_flag=false;
					for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
						int index2_curr=(ikz2+kz_count)%kz_count;
						int index2_next=(ikz2+kz_count+1)%kz_count;
						int index2_prev=(ikz2+kz_count-1)%kz_count;
						if(dispersion_r[index2_curr][ib-1] < Ekin){
							below_flag=true;
						}
						if(dispersion_r[index2_prev][ib-1] < dispersion_r[index2_curr][ib-1] && dispersion_r[index2_curr][ib-1] > dispersion_r[index2_next][ib-1] && dispersion_r[index2_curr][ib-1]>Ekin){
							local_maximum_flag=true;
							break;
						}
					}
					if(below_flag && !local_maximum_flag){
						if(lmax_exist_flag[ikz]){
							int lmax_index=-1;
							int lmax_distance_min=-1;
							for(int il=0; il<lmax_count; il++){
								int lmax_distance=min(min(abs(lmax_indices[il]-ikz), abs(lmax_indices[il]-ikz-kz_count)), abs(lmax_indices[il]-ikz+kz_count));
								if(lmax_distance<=PA_FPFS_kz_margin_index_size && (lmax_index<0 || lmax_distance<lmax_distance_min)){
									lmax_index=il;
									lmax_distance_min=lmax_distance;
								}
							}
							if(lmax_index>=0){
								double eigen_dn=dispersion_r[lmax_indices[lmax_index]][lmax_band[lmax_index]];
								double eigen_gap=eigen_curr-eigen_dn;
								sprintf(sprintf_buffer2, "Local minimum above Ekin at kz=%8.4f, gap=%10.6f", dispersion_kz[ikz], eigen_gap);
								write_log(sprintf_buffer2);
								if(eigen_gap<PA_FPFS_negligible_gap_size){
									solution_kz[solution_index]=dispersion_kz[ikz];
									solution_kappaz[solution_index]=0.0;
									solution_band_indices[solution_index]=ib;
									solution_index++;
									solution_count_real++;
									if(lmax_listed[lmax_index]==false){
										solution_kz[solution_index]=dispersion_kz[lmax_indices[lmax_index]];
										solution_kappaz[solution_index]=0.0;
										solution_band_indices[solution_index]=lmax_band[lmax_index];
										solution_index++;
										solution_count_real++;
										sprintf(sprintf_buffer2, "Local maximum above Ekin at kz=%8.4f is also listed", dispersion_kz[lmax_indices[lmax_index]]);
										write_log(sprintf_buffer2);
										lmax_listed[lmax_index]=true;
									}else{
										sprintf(sprintf_buffer2, "Local maximum above Ekin at kz=%8.4f was already listed", dispersion_kz[lmax_indices[lmax_index]]);
										write_log(sprintf_buffer2);
									}
								}else{
									if(!(ikz==0 || (kz_count%2==0 && ikz==kz_count/2))){
										for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
											if(ikz2>=0 && ikz2<kz_count){
												cspace_search_flag[ikz2]=true;
											}
										}
										cspace_search_indices[cspace_count]=ikz;
										cspace_search_scale[cspace_count]=eigen_gap;
										cspace_count++;
									}
								}
							}
						}else{
							double eigen_dn=dispersion_r[ikz][ib-1];
							double eigen_gap=eigen_curr-eigen_dn;
							double d_Gamma=abs(dispersion_kz[ikz]);
							double d_pBZ=abs(dispersion_kz[ikz]-gz*0.5);
							double d_mBZ=abs(dispersion_kz[ikz]+gz*0.5);
							if(d_Gamma<PA_FPFS_kz_exclude_criterion){
								sprintf(sprintf_buffer2, "Unpaired local minimum above Ekin at kz=%8.4f is not used because it is too close to kz=0", dispersion_kz[ikz]);
								write_log(sprintf_buffer2);
							}else if(d_pBZ<PA_FPFS_kz_exclude_criterion){
								sprintf(sprintf_buffer2, "Unpaired local minimum above Ekin at kz=%8.4f is not used because it is too close to kz=pi/c", dispersion_kz[ikz]);
								write_log(sprintf_buffer2);
							}else if(d_mBZ<PA_FPFS_kz_exclude_criterion){
								sprintf(sprintf_buffer2, "Unpaired local minimum above Ekin at kz=%8.4f is not used because it is too close to kz=-pi/c", dispersion_kz[ikz]);
								write_log(sprintf_buffer2);
							}else{
								sprintf(sprintf_buffer2, "Unpaired local minimum above Ekin at kz=%8.4f, gap=%10.6f", dispersion_kz[ikz], eigen_gap);
								write_log(sprintf_buffer2);
								
								for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
									cspace_search_flag[(ikz2+kz_count)%kz_count]=true;
								}
								cspace_search_indices[cspace_count]=ikz;
								cspace_search_scale[cspace_count]=eigen_gap;
								cspace_count++;
							}
						}
					}
				}
			}
			// find bands crossing Ekin
			int index_next=(ikz+1)%kz_count;
			double eigen_curr=dispersion_r[ikz][ib];
			double eigen_next=dispersion_r[index_next][ib];
			if((eigen_curr-Ekin)*(eigen_next-Ekin)<0){
				// determine kz by interpolation
				double kz_curr=dispersion_kz[ikz];
				double kz_next=dispersion_kz[index_next];
				if(kz_curr==kz_count-1){
					kz_next+=gz;
				}
				double kz_interp=(kz_curr*abs(eigen_next-Ekin)+kz_next*abs(eigen_curr-Ekin))/abs(eigen_next-eigen_curr);

				// determine kz by the bisection
				double kz_left=kz_curr;
				double kz_right=kz_next;
				k_bloch[2]=kz_left;
				double det_left=determinant_sign(g_count, mat, Ekin, k_bloch, g_vec, Vgg,
																				 atom_length, VPS_l_length, vps_cutoff_index,
																				 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																				 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
				k_bloch[2]=kz_right;
				double det_right=determinant_sign(g_count, mat, Ekin, k_bloch, g_vec, Vgg,
																					atom_length, VPS_l_length, vps_cutoff_index,
																					VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																					atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);

				//printf("kz          %10.5f -- %10.5f\n", kz_left, kz_right);
				//printf("Determinant %10.3e -- %10.3e\n", det_left, det_right);

				while(kz_right-kz_left>1e-8){
					double kz_center=(kz_left+kz_right)/2.0;
					k_bloch[2]=kz_center;
					double det_center=determinant_sign(g_count, mat, Ekin, k_bloch, g_vec, Vgg,
																						 atom_length, VPS_l_length, vps_cutoff_index,
																						 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																						 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
					if(det_left*det_center<0){
						kz_right=kz_center;
						det_right=det_center;
					}else if(det_center*det_right<0){
						kz_left=kz_center;
						det_left=det_center;
					}else{
						break;
					}
					//printf("kz          %10.5f -- %10.5f\n", kz_left, kz_right);
					//printf("Determinant %10.3e -- %10.3e\n", det_left, det_right);
				}

				solution_kz[solution_index]=(kz_left+kz_right)/2.0;
				solution_kappaz[solution_index]=0.0;
				solution_kz_alt[solution_index]=kz_interp;
				solution_alt_use[solution_index]=true;
			  
				solution_band_indices[solution_index]=ib;
				solution_index++;
				solution_count_real++;
			}
		}	
	}

	for(int il=0; il<lmax_count; il++){
		if(lmax_listed[il]==false){
			int ikz=lmax_indices[il];
			int ib=lmax_band[il];
			double eigen_curr=lmax_eigen[il];
			double eigen_up=dispersion_r[ikz][ib+1];
			double eigen_gap=eigen_up-eigen_curr;
			double d_Gamma=abs(dispersion_kz[ikz]);
			double d_pBZ=abs(dispersion_kz[ikz]-gz*0.5);
			double d_mBZ=abs(dispersion_kz[ikz]+gz*0.5);
			if(d_Gamma<PA_FPFS_kz_exclude_criterion){
				sprintf(sprintf_buffer2, "Unpaired local maximum above Ekin at kz=%8.4f is not used because it is too close to kz=0", dispersion_kz[ikz]);
				write_log(sprintf_buffer2);
			}else if(d_pBZ<PA_FPFS_kz_exclude_criterion){
				sprintf(sprintf_buffer2, "Unpaired local maximum above Ekin at kz=%8.4f is not used because it is too close to kz=pi/c", dispersion_kz[ikz]);
				write_log(sprintf_buffer2);
			}else if(d_mBZ<PA_FPFS_kz_exclude_criterion){
				sprintf(sprintf_buffer2, "Unpaired local maximum above Ekin at kz=%8.4f is not used because it is too close to kz=-pi/c", dispersion_kz[ikz]);
				write_log(sprintf_buffer2);
			}else{
				sprintf(sprintf_buffer2, "Unpaired local maximum below Ekin at kz=%8.4f, gap=%10.6f", dispersion_kz[ikz], eigen_gap);
				write_log(sprintf_buffer2);
				for(int ikz2=ikz-PA_FPFS_kz_margin_index_size; ikz2<=ikz+PA_FPFS_kz_margin_index_size; ikz2++){
					cspace_search_flag[(ikz2+kz_count)%kz_count]=true;
				}
				cspace_search_indices[cspace_count]=ikz;
				cspace_search_scale[cspace_count]=eigen_gap;
				cspace_count++;
			}
		}
	}

	// complex axis: use dispersion_c
	complex<double> kz;
	double kz_r[4];
	kz_r[0]=0.0;
	kz_r[1]=gz*0.5;
	kz_r[2]=0.0;
	kz_r[3]=gz*0.5;
	int solution_count_complex[4];
	solution_count_complex[0]=0;
	solution_count_complex[1]=0;
	solution_count_complex[2]=0;
	solution_count_complex[3]=0;

	complex<double>** nonlocal_mat_left;
	complex<double>** nonlocal_mat_right;
	complex<double>** nonlocal_mat_center;
	if(PA_FPFS_bulk_include_nonlocal){
		nonlocal_mat_left=alloc_zmatrix(g_count);
		nonlocal_mat_right=alloc_zmatrix(g_count);
		nonlocal_mat_center=alloc_zmatrix(g_count);
	}
	bool nonlocal_mat_set;

	complex<double>** dispersion_use;
	int** connection_use;
	int* dispersion_count_use;
	double* dispersion_kappaz_use;
	int* solution_ikappaz=new int[buffer_size];
	for(int ikzr=0; ikzr<4; ikzr++){
		int kappaz_count_temp;
		if(ikzr==0){
			dispersion_use=dispersion_c;
			connection_use=connection_c;
			dispersion_count_use=dispersion_c_count;
			dispersion_kappaz_use=dispersion_kappaz;
			kappaz_count_temp=kappaz_count;
		}else if(ikzr==1){
			dispersion_use=dispersion_c_BZ;
			connection_use=connection_c_BZ;
			dispersion_count_use=dispersion_c_BZ_count;
			dispersion_kappaz_use=dispersion_kappaz;
			kappaz_count_temp=kappaz_border_index;
		}else if(ikzr==2){
			dispersion_use=dispersion_mc;
			connection_use=connection_mc;
			dispersion_count_use=dispersion_mc_count;
			dispersion_kappaz_use=dispersion_mkappaz;
			kappaz_count_temp=kappaz_border_index;
		}else if(ikzr==3){
			dispersion_use=dispersion_mc_BZ;
			connection_use=connection_mc_BZ;
			dispersion_count_use=dispersion_mc_BZ_count;
			dispersion_kappaz_use=dispersion_mkappaz;
			kappaz_count_temp=kappaz_border_index;
		}else{
			break;
		}
		nonlocal_mat_set=false;
		for(int ikappaz=0; ikappaz<kappaz_count_temp-1; ikappaz++){
			for(int ib=0; ib<dispersion_count_use[ikappaz]; ib++){
				if(connection_use[ikappaz][ib]<0){
					continue;
				}
				double eigen_curr=dispersion_use[ikappaz][ib].real();
				double eigen_next=dispersion_use[ikappaz+1][connection_use[ikappaz][ib]].real();
				double kappaz_curr=dispersion_kappaz_use[ikappaz];
				double kappaz_next=dispersion_kappaz_use[ikappaz+1];
				if((eigen_curr-Ekin)*(eigen_next-Ekin)<0){
					// interpolation
					double kappaz_interp=(kappaz_curr*abs(eigen_next-Ekin)+kappaz_next*abs(eigen_curr-Ekin))/abs(eigen_next-eigen_curr);

					// bisection
					k_bloch[2]=0.0;
					double kappaz_left=dispersion_kappaz_use[ikappaz];
					double kappaz_right=dispersion_kappaz_use[ikappaz+1];
					double kappaz_center;
					complex<double> kz_left(kz_r[ikzr], -kappaz_left);
					complex<double> kz_right(kz_r[ikzr], -kappaz_right);
					complex<double> kz_center;
				
					complex<double> det_left=determinant_complex(g_count, mat, Ekin, k_bloch, kz_left, g_vec, Vgg,
																											 atom_length, VPS_l_length, vps_cutoff_index,
																											 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																											 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																											 nonlocal_mat_set?2:1, nonlocal_mat_center);
					nonlocal_mat_set=true;
					complex<double> det_right=determinant_complex(g_count, mat, Ekin, k_bloch, kz_right, g_vec, Vgg,
																												atom_length, VPS_l_length, vps_cutoff_index,
																												VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																												atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																												2, nonlocal_mat_center);
					complex<double> det_center;

					while(kappaz_right-kappaz_left>1e-8){
						kappaz_center=(kappaz_left+kappaz_right)*0.5;
						kz_center=complex<double>(kz_r[ikzr], -kappaz_center);
						det_center=determinant_complex(g_count, mat, Ekin, k_bloch, kz_center, g_vec, Vgg,
																					 atom_length, VPS_l_length, vps_cutoff_index,
																					 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																					 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																					 2, nonlocal_mat_center);
						if(det_left.real()*det_center.real()<0){
							kappaz_right=kappaz_center;
							kz_right=complex<double>(kz_r[ikzr], -kappaz_right);
							det_right=determinant_complex(g_count, mat, Ekin, k_bloch, kz_right, g_vec, Vgg,
																						atom_length, VPS_l_length, vps_cutoff_index,
																						VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																						atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																						2, nonlocal_mat_center);
						}else if(det_center.real()*det_right.real()<0){
							kappaz_left=kappaz_center;
							kz_left=complex<double>(kz_r[ikzr], -kappaz_center);
							det_center=determinant_complex(g_count, mat, Ekin, k_bloch, kz_left, g_vec, Vgg,
																						 atom_length, VPS_l_length, vps_cutoff_index,
																						 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																						 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																						 2, nonlocal_mat_center);
						}else{
							break;
						}
					}

					solution_kz[solution_index]=kz_r[ikzr];
					solution_kappaz[solution_index]=kappaz_interp;
					solution_ikappaz[solution_index]=ikappaz;

					solution_kappaz_alt[solution_index]=(kappaz_left+kappaz_right)*0.5;
					solution_alt_use[solution_index]=true;

					solution_band_indices[solution_index]=ib;
					if(connection_use[ikappaz][ib]!=ib){
						solution_band_indices2[solution_index]=connection_use[ikappaz][ib];
					}else{
						solution_band_indices2[solution_index]=-1;
					}
					
					solution_index++;
					solution_count_complex[ikzr]++;
				}
			} // end of for(ib)
		} // end of for(ikappaz)
		// check
		int solution_index_offset=solution_index-solution_count_complex[ikzr];
		bool* passed=new bool[solution_count_complex[ikzr]];
		for(int is=0; is<solution_count_complex[ikzr]; is++){
			passed[is]=true;
		}
		int removed_count=0;
		// connected to another
		for(int is=0; is<solution_count_complex[ikzr]; is++){
			int is_offset=is+solution_index_offset;
			for(int is_prev=0; is_prev<solution_count_complex[ikzr]; is_prev++){
				int is_prev_offset=is_prev+solution_index_offset;
				if(is_prev!=is && solution_ikappaz[is_prev_offset] < solution_ikappaz[is_offset]){
					int band_index=solution_band_indices[is_prev_offset];
					bool terminated=false;
					for(int ikappaz_path=solution_ikappaz[is_prev_offset]; ikappaz_path < solution_ikappaz[is_offset]; ikappaz_path++){
						band_index=connection_use[ikappaz_path][band_index];
						if(band_index<0){
							terminated=true;
							break;
						}
					}
					if(!terminated && band_index==solution_band_indices[is_offset]){
						passed[is]=false;
						removed_count++;
						sprintf(sprintf_buffer2, "Complex axis: kz=%8.4f, kappaz=%8.4f is removed because it is connected to another point",
										solution_kz[is_offset], solution_kappaz[is_offset]);
						write_log(sprintf_buffer2);
						break;
					}
				}
			}
		}
		// close to another removed point
		for(int is=0; is<solution_count_complex[ikzr]; is++){
		  if(passed[is]==false){
		    continue;
		  }
		      
			int is_offset=is+solution_index_offset;
			for(int is_prev=0; is_prev<solution_count_complex[ikzr]; is_prev++){
			  if(is==is_prev){
			    continue;
			  }
				int is_prev_offset=is_prev+solution_index_offset;
				if(passed[is_prev]==false){
					if(abs(solution_band_indices[is_prev_offset]-solution_band_indices[is_offset])<=1 &&
						 abs(solution_ikappaz[is_prev_offset]-solution_ikappaz[is_offset])<=2){
						passed[is]=false;
						sprintf(sprintf_buffer2, "Complex axis: kz=%8.4f, kappaz=%8.4f is removed because it is close to another removed point",
										solution_kz[is_offset], solution_kappaz[is_offset]);
						write_log(sprintf_buffer2);
						removed_count++;
						break;
					}
				}
			}
		}
		// re-fill
		if(removed_count>0){
			solution_index=solution_index_offset;
			for(int is=0; is<solution_count_complex[ikzr]; is++){
				int is_offset=is+solution_index_offset;
				if(passed[is]){
					solution_kz[solution_index]=solution_kz[is_offset];
					solution_kappaz[solution_index]=solution_kappaz[is_offset];

					solution_kappaz_alt[solution_index]=solution_kappaz_alt[is_offset];
					solution_alt_use[solution_index]=solution_alt_use[is_offset];

					solution_band_indices[solution_index]=solution_band_indices[is_offset];
					solution_band_indices2[solution_index]=solution_band_indices2[is_offset];

					solution_index++;
				}
			}
			sprintf(sprintf_buffer2, "Complex axis: %2d points are removed from %2d candidates", removed_count, solution_count_complex[ikzr]);
			write_log(sprintf_buffer2);
			solution_count_complex[ikzr]-=removed_count;
		}
		delete[] passed;
	}
	delete[] solution_ikappaz;
	solution_count_complex_sum=solution_count_complex[0]+solution_count_complex[1]
		+solution_count_complex[2]+solution_count_complex[3];

	// complex plane
	// cout << "C plane" << endl;
	complex<double>** kz_corner1=alloc_zmatrix(buffer_size, 4);
	complex<double>** kz_corner2;
	int kz_corner_index=0;
	double* dispersion_kappaz_temp=new double[kappaz_border_index];
	double dkz=dispersion_kappaz[1]-dispersion_kappaz[0];
	for(int ikz=1; ikz<kz_count-1; ikz++){
		// printf("ikz %d / %d\n", ikz, kz_count);
		if(kz_count%2==0){
			if(ikz==kz_count/2 || ikz==kz_count/2-1){
				continue;
			}
		}else{
			if(ikz==(kz_count-1)/2){
				continue;
			}
		}
		if(cspace_search_flag[ikz]==false){
			continue;
		}
		int area_min=ikz;
		int area_max=ikz;
		for(; area_min>=0; area_min--){
			if(cspace_search_flag[area_min]==false){
				break;
			}
			if(kz_count%2==0){
				if(area_min==kz_count/2){
					break;
				}
			}else{
				if(area_min==(kz_count+1)/2){
					break;
				}
			}
		}
		area_min++;
		for(; area_max<kz_count; area_max++){
			if(cspace_search_flag[area_max]==false){
				break;
			}
			if(kz_count%2==0){
				if(area_max==kz_count/2){
					break;
				}
			}else{
				if(area_max==(kz_count-1)/2){
					break;
				}
			}
		}
		area_max--;
		int count=0;
		double gap_sum=0.0;
		for(int ic=0; ic<cspace_count; ic++){
			if(area_min<=cspace_search_indices[ic] && cspace_search_indices[ic]<=area_max){
				gap_sum+=cspace_search_scale[ic];
				count++;
			}
		}
		if(count==0){
			continue;
		}
		double gap_ave=gap_sum/(1.0*count);
		// printf("Area: %8.4f to %8.4f, gap=%10.6f\n", dispersion_kz[area_min], dispersion_kz[area_max], gap_ave);
		nonlocal_mat_set=false;
		double dkappaz_temp=gap_ave*PA_FPFS_cspace_size/(kappaz_border_index*sqrt(2.0*PA_excitation_energy/Eh));
		for(int inp=0; inp<2; inp++){
			for(int ikappaz=0; ikappaz<kappaz_border_index; ikappaz++){
				if(inp==0){
					dispersion_kappaz_temp[ikappaz]=dkappaz_temp*(1.0+ikappaz);
				}else{
					dispersion_kappaz_temp[ikappaz]=-dkappaz_temp*(1.0+ikappaz);
				}
			}
			for(int ikappaz=0; ikappaz<kappaz_border_index-1; ikappaz++){
				// printf("%d %d\n", ikz, ikappaz);
				int ikz_next=(ikz+1)%kz_count;
				complex<double> kz_lb(dispersion_kz[ikz], -dispersion_kappaz_temp[ikappaz]);
				complex<double> kz_lt(dispersion_kz[ikz], -dispersion_kappaz_temp[ikappaz+1]);
				complex<double> kz_rt(dispersion_kz[ikz_next], -dispersion_kappaz_temp[ikappaz+1]);				
				complex<double> kz_rb(dispersion_kz[ikz_next], -dispersion_kappaz_temp[ikappaz]);
				if(ikz==kz_count-1){
					kz_rt+=gz;
					kz_rb+=gz;
				}
				
				complex<double> det_lb=determinant_complex(g_count, mat, Ekin, k_bloch, kz_lb, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 nonlocal_mat_set?2:1, nonlocal_mat_left);
				complex<double> det_lt=determinant_complex(g_count, mat, Ekin, k_bloch, kz_lt, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_left);
				complex<double> det_rt=determinant_complex(g_count, mat, Ekin, k_bloch, kz_rt, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 nonlocal_mat_set?2:1, nonlocal_mat_right);
				complex<double> det_rb=determinant_complex(g_count, mat, Ekin, k_bloch, kz_rb, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_right);
				nonlocal_mat_set=true;

				if(encloseOrigin(det_lb, det_lt, det_rt, det_rb)){
					// printf("Det: (%f %f) (%f %f) (%f %f) (%f %f)\n", det_lb.real(), det_lb.imag(), det_lt.real(), det_lt.imag(), det_rt.real(), det_rt.imag(), det_rb.real(), det_rt.imag());
					if(kz_corner_index>=buffer_size){
						write_log((char*)"Candidates are more than buffer size");
						break;
					}
					kz_corner1[kz_corner_index][0]=kz_lb;
					kz_corner1[kz_corner_index][1]=kz_lt;
					kz_corner1[kz_corner_index][2]=kz_rt;
					kz_corner1[kz_corner_index][3]=kz_rb;
					kz_corner_index++;
				}
			}
		} // end of for(inp)
	} // end of for(kz)
	int kz_corner_count=kz_corner_index;
	if(kz_corner_count>0){
		kz_corner2=alloc_zmatrix(kz_corner_count*4, 4);
		while(dkz>1e-6){
			// printf("Corner count: %d\n", kz_corner_count);
			int is1, is2;
			is2=0;
			for(is1=0; is1<kz_corner_count; is1++){
				// printf("is = %d / %d kz = %f %f\n", is1, kz_corner_count, kz_corner1[is1][0].real(), kz_corner1[is1][0].imag());
				// find kz by bisection
				complex<double> kz_lb=kz_corner1[is1][0];
				complex<double> kz_lt=kz_corner1[is1][1];
				complex<double> kz_rt=kz_corner1[is1][2];
				complex<double> kz_rb=kz_corner1[is1][3];
				complex<double> det_lb=determinant_complex(g_count, mat, Ekin, k_bloch, kz_lb, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 1, nonlocal_mat_left);
				complex<double> det_lt=determinant_complex(g_count, mat, Ekin, k_bloch, kz_lt, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_left);
				complex<double> det_rt=determinant_complex(g_count, mat, Ekin, k_bloch, kz_rt, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 1, nonlocal_mat_right);
				complex<double> det_rb=determinant_complex(g_count, mat, Ekin, k_bloch, kz_rb, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_right);
				// printf("Det: (%f %f) (%f %f) (%f %f) (%f %f)\n", det_lb.real(), det_lb.imag(), det_lt.real(), det_lt.imag(), det_rt.real(), det_rt.imag(), det_rb.real(), det_rt.imag());
			
				complex<double> kz_lc=(kz_lb+kz_lt)*0.5;
				complex<double> kz_cb=(kz_lb+kz_rb)*0.5;
				complex<double> kz_cc=(kz_lb+kz_rt)*0.5;
				complex<double> kz_ct=(kz_lt+kz_rt)*0.5;
				complex<double> kz_rc=(kz_rb+kz_rt)*0.5;
				complex<double> det_lc=determinant_complex(g_count, mat, Ekin, k_bloch, kz_lc, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_left);
				complex<double> det_cb=determinant_complex(g_count, mat, Ekin, k_bloch, kz_cb, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 1, nonlocal_mat_center);
				complex<double> det_cc=determinant_complex(g_count, mat, Ekin, k_bloch, kz_cc, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_center);
				complex<double> det_ct=determinant_complex(g_count, mat, Ekin, k_bloch, kz_ct, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_center);
				complex<double> det_rc=determinant_complex(g_count, mat, Ekin, k_bloch, kz_rc, g_vec, Vgg,
																									 atom_length, VPS_l_length, vps_cutoff_index,
																									 VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
																									 atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
																									 2, nonlocal_mat_right);

				if(encloseOrigin(det_lb, det_lc, det_cc, det_cb)){
					kz_corner2[is2][0]=kz_lb;			
					kz_corner2[is2][1]=kz_lc;
					kz_corner2[is2][2]=kz_cc;
					kz_corner2[is2][3]=kz_cb;
					is2++;
				}
				if(encloseOrigin(det_lc, det_lt, det_ct, det_cc)){
					kz_corner2[is2][0]=kz_lc;
					kz_corner2[is2][1]=kz_lt;
					kz_corner2[is2][2]=kz_ct;
					kz_corner2[is2][3]=kz_cc;
					is2++;
				}
				if(encloseOrigin(det_cb, det_cc, det_rc, det_rb)){
					kz_corner2[is2][0]=kz_cb;
					kz_corner2[is2][1]=kz_cc;
					kz_corner2[is2][2]=kz_rc;
					kz_corner2[is2][3]=kz_rb;
					is2++;
				}
				if(encloseOrigin(det_cc, det_ct, det_rt, det_rc)){
					kz_corner2[is2][0]=kz_cc;
					kz_corner2[is2][1]=kz_ct;
					kz_corner2[is2][2]=kz_rt;
					kz_corner2[is2][3]=kz_rc;
					is2++;
				}
			}
			kz_corner_count=is2;
			delete_zmatrix(kz_corner1);
			kz_corner1=kz_corner2;
			if(kz_corner_count>0){
				kz_corner2=alloc_zmatrix(kz_corner_count*4, 4);
			}else{
				break;
			}
			dkz*=0.5;
		} // end of while

		int index_offset=solution_count_real+solution_count_complex_sum;
		for(int is=0; is<kz_corner_count; is++){
			complex<double> kz_sol=(kz_corner1[is][0]+kz_corner1[is][2])*0.5;
			//printf("kz_sol %10.6f %10.6f\n", kz_sol.real(), kz_sol.imag());

			double nearest_distance=-1;
			for(int is2=index_offset; is2<solution_index; is2++){
				double distance=abs(complex<double>(solution_kz[is2], -solution_kappaz[is2])-kz_sol);
				if(kz_sol.imag()*(-solution_kappaz[is2])>0 && (nearest_distance<0 || nearest_distance>distance)){
					nearest_distance=distance;
				}
			}
			if(nearest_distance<0 || nearest_distance>0.01){
				solution_kz[solution_index]=kz_sol.real();
				solution_kappaz[solution_index]=-kz_sol.imag();
				solution_alt_use[solution_index]=false;
				solution_count_cspace++;
				solution_index++;
			}
		}
		if(kz_corner_count>0){
			delete_zmatrix(kz_corner2);
		}
	}
	delete_zmatrix(kz_corner1);

	//int solution_count=solution_count_cross+solution_count_local;
	int solution_count=solution_count_real+solution_count_complex_sum+solution_count_cspace;
	sprintf(sprintf_buffer2, "%3d solutions are found (real: %3d, complex: %3d %3d %3d %3d, complex space: %3d)", solution_count, solution_count_real, solution_count_complex[0], solution_count_complex[1], solution_count_complex[2], solution_count_complex[3], solution_count_cspace);
	write_log(sprintf_buffer2);

	// for debug
	//for(int is=0; is<solution_count; is++){
	//	printf("%3d | %8.4f %8.4f | %8.4f %8.4f\n", is, solution_kz[is], solution_kappaz[is], solution_kz[is], solution_kappaz_alt[is]);
	//}
	//printf("\n");

	if(solution_count==0){
		write_log((char*)"Warning: No bulk solution found");
		return 0;
	}
	complex<double>** solution=alloc_zmatrix(solution_count, g_count);
	*final_states_pointer=solution;

	*kz_pointer=solution_kz;
	*kappaz_pointer=solution_kappaz;
			
	// for zheev and zgeev
	char compute='V';
	char nocompute='N';
	char uplo='U';
	int ld=1;
	
	complex<double> work_dummy;
	int lwork=-1;
	double* rwork=new double[3*g_count-2];
	int info;
	double* w=new double[g_count];
	double* rwork2=new double[2*g_count];
	complex<double>* wc=new complex<double>[g_count];
	complex<double>* work;
		
	double eigen_diff_max=0.0;
	// solution (real)
	for(int is=0; is<solution_count_real; is++){
		k_bloch[2]=solution_kz[is];
		prepare_matrix_bulk(g_count, mat, 0.0, k_bloch, g_vec, Vgg);
		if(PA_FPFS_bulk_include_nonlocal){
			add_nonlocal_term(g_count, g_vec, mat, k_bloch, k_bloch[2],
												atom_length, VPS_l_length, vps_cutoff_index,
												VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
												atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
		}

		if(lwork==-1){
			zheev_(&compute, &uplo, &g_count, &mat[0][0], &g_count, &w[0], &work_dummy, &lwork, &rwork[0], &info);
			if(info!=0){
				write_log((char*)"zheev preparation failed");
				return 0;
			}
			lwork=round(abs(work_dummy));
			work=new complex<double>[lwork];
		}
		
		zheev_(&compute, &uplo, &g_count, &mat[0][0], &g_count, &w[0], &work[0], &lwork, &rwork[0], &info);
		if(info!=0){
			write_log((char*)"zheev failed");
			return 0;
		}
		for(i=0; i<g_count; i++){
			solution[is][i]=mat[solution_band_indices[is]][i];
		}
		double eigen_diff=abs(w[solution_band_indices[is]]*0.5-Ekin);

		if(solution_alt_use[is]){
			k_bloch[2]=solution_kz_alt[is];
			
			prepare_matrix_bulk(g_count, mat, 0.0, k_bloch, g_vec, Vgg);
			if(PA_FPFS_bulk_include_nonlocal){
				add_nonlocal_term(g_count, g_vec, mat, k_bloch, k_bloch[2],
													atom_length, VPS_l_length, vps_cutoff_index,
													VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
													atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
			}
			zheev_(&compute, &uplo, &g_count, &mat[0][0], &g_count, &w[0], &work[0], &lwork, &rwork[0], &info);
			if(info!=0){
				write_log((char*)"zheev failed");
				return 0;
			}
			double eigen_diff_2=abs(w[solution_band_indices[is]]*0.5-Ekin);
			if(eigen_diff_2<eigen_diff){
				eigen_diff=eigen_diff_2;
				solution_kz[is]=solution_kz_alt[is];
				for(i=0; i<g_count; i++){
					solution[is][i]=mat[solution_band_indices[is]][i];
				}
				//cout << "alt" << endl;
			}
		}
		//printf("Diff: %8.4f\n", eigen_diff);
		if(eigen_diff>eigen_diff_max){
			eigen_diff_max=eigen_diff;
		}
	}
	if(lwork>0){
	  lwork=-1;
	  delete[] work;
	}
	
	// solution (complex)
	for(int is=solution_count_real; is<solution_count; is++){
		k_bloch[2]=0.0;
		complex<double> kz(solution_kz[is], -solution_kappaz[is]);
		prepare_matrix_bulk_complex(g_count, mat, 0.0, k_bloch, kz, g_vec, Vgg);
		if(PA_FPFS_bulk_include_nonlocal){
			add_nonlocal_term(g_count, g_vec, mat, k_bloch, kz.real(),
												atom_length, VPS_l_length, vps_cutoff_index,
												VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
												atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
		}
		if(lwork==-1){
			zgeev_(&nocompute, &compute, &g_count, &mat[0][0], &g_count, &wc[0], NULL, &ld, &vr[0][0], &g_count, &work_dummy, &lwork, &rwork2[0], &info);
			if(info!=0){
				write_log((char*)"zgeev preparation failed");
				return 0;
			}
			lwork=round(abs(work_dummy));
			work=new complex<double>[lwork];
		}
		zgeev_(&nocompute, &compute, &g_count, &mat[0][0], &g_count, &wc[0], NULL, &ld, &vr[0][0], &g_count, &work[0], &lwork, &rwork2[0], &info);
		if(info!=0){
			write_log((char*)"zgeev failed");
			return 0;
		}
		double eigen_diff;
		if(is<solution_count_real+solution_count_complex_sum){
			// solution (complex)
			int eigen_index=find_eigenstate(g_count, solution_band_indices[is], solution_band_indices2[is], Ekin, wc, vr, g_index);
			eigen_diff=abs(0.5*wc[eigen_index]-Ekin);
			for(int ig=0; ig<g_count; ig++){
				solution[is][ig]=vr[eigen_index][ig];
			}
			if(solution_alt_use[is]){
				complex<double> kz(0.0, -solution_kappaz_alt[is]);
				prepare_matrix_bulk_complex(g_count, mat, 0.0, k_bloch, kz, g_vec, Vgg);
				if(PA_FPFS_bulk_include_nonlocal){
					add_nonlocal_term(g_count, g_vec, mat, k_bloch, kz.real(),
														atom_length, VPS_l_length, vps_cutoff_index,
														VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
														atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
				}
				zgeev_(&nocompute, &compute, &g_count, &mat[0][0], &g_count, &wc[0], NULL, &ld, &vr[0][0], &g_count, &work[0], &lwork, &rwork2[0], &info);
				if(info!=0){
					write_log((char*)"zgeev failed");
					return 0;
				}
				int eigen_index2=find_eigenstate(g_count, solution_band_indices[is], solution_band_indices2[is], Ekin, wc, vr, g_index);
				double eigen_diff2=abs((0.5*wc[eigen_index2]).real()-Ekin);
				if(eigen_diff2<eigen_diff){
					eigen_index=eigen_index2;
					eigen_diff=eigen_diff2;
					for(int ig=0; ig<g_count; ig++){
						solution[is][ig]=vr[eigen_index][ig];
					}
				
				}
			}
		}else{
			// solution (cspace)
			int nearest_index=-1;
			double nearest_distance=-1;
			for(int ig=0; ig<g_count; ig++){
				double distance=abs(0.5*wc[ig]-Ekin);
				if(nearest_index<0 || nearest_distance>distance){
					nearest_index=ig;
					nearest_distance=distance;
				}
			}
			eigen_diff=nearest_distance;
			for(int ig=0; ig<g_count; ig++){
				solution[is][ig]=vr[nearest_index][ig];
			}
		}
		
		//printf("eigen diff: %12.4e\n", eigen_diff);
		if(eigen_diff_max<eigen_diff){
			eigen_diff_max=eigen_diff;
		}
	}
	
	sprintf(sprintf_buffer2, "Diff max: %12.8f", eigen_diff_max);
	write_log(sprintf_buffer2);

	//delete[] g_vec;
	//delete[] Vgg;
	delete[] sprintf_buffer2;
	//delete[] dispersion;
	//delete_zmatrix(mat);
	//delete_bmatrix(isSolution);
	delete[] solution_kz_alt;
	delete[] solution_kappaz_alt;
	delete[] solution_alt_use;
	delete[] solution_band_indices;
	delete[] solution_band_indices2;
	delete[] rwork;
	delete[] w;
	delete[] rwork2;
	delete[] wc;
	if(lwork>0){
	  delete[] work;
	}
	delete[] cspace_search_flag;
	delete[] cspace_search_indices;
	delete[] cspace_search_scale;
	delete[] lmax_exist_flag;
	delete[] lmax_indices;
	delete[] lmax_eigen;
	delete[] lmax_listed;
	delete[] lmax_band;
	delete[] dispersion_kappaz_temp;
	if(PA_FPFS_bulk_include_nonlocal){
		delete_zmatrix(nonlocal_mat_left);
		delete_zmatrix(nonlocal_mat_right);
		delete_zmatrix(nonlocal_mat_center);
	}
	return solution_count;
	
}

double* interpolate_wfn(int wfn_length, double* wfn, double* r, int wfn_length_reduced, double dr){
	double* wfn_reduced=new double[wfn_length_reduced];
	int* wfn_count=new int[wfn_length_reduced];

	for(int ir=0; ir<wfn_length_reduced; ir++){
		wfn_reduced[ir]=0.0;
		wfn_count[ir]=0;
	}
	for(int ir=0; ir<wfn_length; ir++){
		int index=floor(r[ir]/dr);
		if(index>=0 && index<wfn_length_reduced){
			wfn_reduced[index]+=wfn[ir];
			wfn_count[index]++;
		}
	}
	for(int ir=0; ir<wfn_length_reduced; ir++){
		if(wfn_count[ir]>0){
			wfn_reduced[ir]/=wfn_count[ir];
		}
	}
	for(int ir=0; ir<wfn_length_reduced; ir++){
		if(wfn_count[ir]==0){
			int ir_left, ir_right;
			bool ir_left_found=false;
			bool ir_right_found=false;
			for(ir_left=ir-1; ir_left>=0; ir_left--){
				if(wfn_count[ir_left]>0){
					ir_left_found=true;
					break;
				}
			}
			for(ir_right=ir+1; ir_right<wfn_length_reduced; ir_right++){
				if(wfn_count[ir_right]>0){
					ir_right_found=true;
					break;
				}
			}
			if(ir_left_found && ir_right_found){
				wfn_reduced[ir]=(wfn_reduced[ir_right]*(ir-ir_left)+wfn_reduced[ir_left]*(ir_right-ir))/(ir_right-ir_left);
			}
		}
	}
  
	delete[] wfn_count;
	/*
	for(int ir=0; ir<wfn_length; ir++){
		printf("%8.4f %8.4f ", r[ir], wfn[ir]);
		if(ir<wfn_length_reduced){
			printf("%8.4f %8.4f", (ir+0.5)*dr, wfn_reduced[ir]);
		}
		printf("\n");
		}*/
	
	return wfn_reduced;
}


void calc_bulk_dispersion(double* k_para, int kz_count, double* kz, int g_count, double** g_vec, complex<double>** Vgg, double** eigen, complex<double>** mat,													
													int atom_length, int* VPS_l_length, int* vps_cutoff_index,
													int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
													double** atom_coordinates, int* atom_spec_index,
													double FPFS_bulk_min, double FPFS_bulk_height){
	int i, j;
	// prepare matrix
	//double** g_vec=new double*[g_count];
	//for(i=0; i<g_count; i++){
	//	g_vec[i]=&g_vec_buffer[3*i];
	//}
	//complex<double>** Vgg=new complex<double>*[g_count];
	//for(i=0; i<g_count; i++){
	//		Vgg[i]=&Vgg_buffer[i*g_count];
	//}
	//double** eigen=new double*[kz_count];
	//for(i=0; i<kz_count; i++){
	//	eigen[i]=&eigen_buffer[i*g_count];
	//}
	
	
	//complex<double>** mat=alloc_zmatrix(g_count); // transposed

	double k_bloch[3];
	k_bloch[0]=k_para[0];
	k_bloch[1]=k_para[1];
	
	// for zheev
	char nocompute='N';
	char uplo='U';
	
	complex<double> work_dummy;
	int lwork=-1;
	double rwork[3*g_count-2];
	int info;
	double w[g_count];
	complex<double>* work;

	for(int ikz=0; ikz<kz_count; ikz++){
		//printf("ikz=%4d/%4d start\n", ikz, kz_count);
		k_bloch[2]=kz[ikz];
		
		prepare_matrix_bulk(g_count, mat, 0.0, k_bloch, g_vec, Vgg);
		if(PA_FPFS_bulk_include_nonlocal){
			add_nonlocal_term(g_count, g_vec, mat, k_bloch, k_bloch[2],
												atom_length, VPS_l_length, vps_cutoff_index,
												VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
												atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height, 0, NULL);
		}
		if(lwork==-1){
			zheev_(&nocompute, &uplo, &g_count, &mat[0][0], &g_count, &w[0], &work_dummy, &lwork, &rwork[0], &info);
			if(info!=0){
				write_log((char*)"zheev preparation failed");
				return;
			}
			lwork=round(abs(work_dummy));
			work=new complex<double>[lwork];
		}
		
		zheev_(&nocompute, &uplo, &g_count, &mat[0][0], &g_count, &w[0], &work[0], &lwork, &rwork[0], &info);
		if(info!=0){
			write_log((char*)"zheev failed");
			return;
		}
		for(int ig=0; ig<g_count; ig++){
			eigen[ikz][ig]=0.5*w[ig];
		}
	}
	delete[] work;
	//delete[] g_vec;
	//delete[] Vgg;
	//delete[] eigen;
	//delete_zmatrix(mat);
	
}

void calc_bulk_dispersion_complex(double* k_para, double kz_r, int kappaz_count, int kappaz_border_index, double* kappaz_list, int g_count, int** g_index, double** g_vec, complex<double>** Vgg, complex<double>*** eigen_pointer, int* eigen_count, int*** connection_pointer, complex<double>** mat,
												int atom_length, int* VPS_l_length, int* vps_cutoff_index,
												int** VPS_l, double** VPS_r, double** VPS_E_ave, double*** VPS_nonloc_ave,
												double** atom_coordinates, int* atom_spec_index,
												double FPFS_bulk_min, double FPFS_bulk_height){
	char* sprintf_buffer2=new char[Log_length+1];
	int i,j;
	double k_bloch[3];
	k_bloch[0]=k_para[0];
	k_bloch[1]=k_para[1];
	k_bloch[2]=0;
	complex<double> kz;
	
	// for zgeev
	char nocompute='N';
	char compute='V';
	int ld=1;
	
	complex<double> work_dummy;
	int lwork=-1;
	double* rwork=new double[2*g_count];
	int info;
	complex<double>* w=new complex<double>[g_count];
	complex<double>* work;

	complex<double>* dispersion_real=new complex<double>[g_count];
	int dispersion_real_count;
	int row_size;

	complex<double>** nonlocal_matrix;
	if(PA_FPFS_bulk_include_nonlocal){
		nonlocal_matrix=alloc_zmatrix(g_count);
	}
	for(int ikappaz=0; ikappaz<kappaz_count; ikappaz++){
		//printf("ikappaz=%4d/%4d start\n", ikappaz, kappaz_count);
		
		kz=complex<double>(kz_r, -kappaz_list[ikappaz]);
		//printf("kz = (%f %f)\n", kz.real(), kz.imag());
					 
		prepare_matrix_bulk_complex(g_count, mat, 0.0, k_bloch, kz, g_vec, Vgg);
		if(PA_FPFS_bulk_include_nonlocal){
			add_nonlocal_term(g_count, g_vec, mat, k_bloch, kz_r,
												atom_length, VPS_l_length, vps_cutoff_index,
												VPS_l, VPS_r, VPS_E_ave, VPS_nonloc_ave,
												atom_coordinates, atom_spec_index, FPFS_bulk_min, FPFS_bulk_height,
												ikappaz==0?1:2, nonlocal_matrix);
		}
		
		//if(ikz==kz_count/2){
		/*
			for(int ig1=0; ig1<g_count; ig1++){
			for(int ig2=0; ig2<g_count; ig2++){
			printf("(%8.4f %8.4f) ", mat[ig2][ig1].real(), mat[ig2][ig1].imag());
			}
			printf("\n");
			}
			printf("\n");*/
			//}
		if(lwork==-1){
			zgeev_(&nocompute, &nocompute, &g_count, &mat[0][0], &g_count, &w[0], NULL, &ld, NULL, &ld, &work_dummy, &lwork, &rwork[0], &info);
			if(info!=0){
				write_log((char*)"zgeev preparation failed");
				return;
			}
			lwork=round(abs(work_dummy));
			work=new complex<double>[lwork];
		}
		
		zgeev_(&nocompute, &nocompute, &g_count, &mat[0][0], &g_count, &w[0], NULL, &ld, NULL, &ld, &work[0], &lwork, &rwork[0], &info);
		if(info!=0){
			write_log((char*)"zgeev failed");
			return;
		}

		dispersion_real_count=0;
		for(int ig=0; ig<g_count; ig++){
			if(abs(w[ig].imag())*0.5<PA_FPFS_real_eigenvalue_criterion){
				dispersion_real[dispersion_real_count]=0.5*w[ig];
				dispersion_real_count++;
			}
		}
		//printf("%d\n", dispersion_real_count);
		//printf("\n");

		// determine the matrix size using the count at the left edge
		if(ikappaz==0){
			row_size=dispersion_real_count;
			*eigen_pointer=alloc_zmatrix(kappaz_count, row_size);
			*connection_pointer=alloc_imatrix(kappaz_count, row_size);
		}
		eigen_count[ikappaz]=dispersion_real_count;

		// fill
		if(dispersion_real_count<=eigen_count[0]){
			// no problem
			for(int ib=0; ib<dispersion_real_count; ib++){
				(*eigen_pointer)[ikappaz][ib]=dispersion_real[ib];
			}
		}else{
			write_log((char*)"Warning: more eigenstates than expected");
			for(int ib=0; ib<eigen_count[0]; ib++){
				(*eigen_pointer)[ikappaz][ib]=dispersion_real[ib];
			}
		}

		// sort
		sort_ascend(dispersion_real_count, (*eigen_pointer)[ikappaz]);
	}
	if(PA_FPFS_bulk_include_nonlocal){
		delete_zmatrix(nonlocal_matrix);
	}

	// initialize connection_pointer
	// -1 means that that state has no connection with ikappaz+1
	for(int ikappaz=0; ikappaz<kappaz_count; ikappaz++){
		for(int ib=0; ib<row_size; ib++){
			(*connection_pointer)[ikappaz][ib]=-1;
		}
	}

	int* group_l=new int[row_size];
	int* group_r=new int[row_size];
	double crt;
	for(int ikappaz=0; ikappaz<kappaz_count-1; ikappaz++){
		if(eigen_count[ikappaz]==0){
			continue;
		}
		crt=ikappaz<kappaz_border_index ? PA_FPFS_dispersion_group_criterion_left: PA_FPFS_dispersion_group_criterion_right;
		int ik_l=ikappaz;
		int ik_r=ikappaz+1;
		// grouping
		for(int ib=0; ib<row_size; ib++){
			group_l[ib]=-1;
			group_r[ib]=-1;
		}
		int groupNo=0;
		int index_l=0;
		int index_r=0;
		while(index_l<eigen_count[ik_l] || index_r<eigen_count[ik_r]){
			double er_top=-1;
			bool er_top_set=false;
			if(index_l<eigen_count[ik_l]){
				er_top=(*eigen_pointer)[ik_l][index_l].real()+crt;
				er_top_set=true;
			}
			if(index_r<eigen_count[ik_r]){
				if(er_top_set){
					er_top=min(er_top, (*eigen_pointer)[ik_r][index_r].real()+crt);
				}else{
					er_top=(*eigen_pointer)[ik_r][index_r].real()+crt;
				}
			}
			while(true){
				bool added=false;
				if(index_l<eigen_count[ik_l] && (*eigen_pointer)[ik_l][index_l].real()<er_top){
					group_l[index_l]=groupNo;
					er_top=(*eigen_pointer)[ik_l][index_l].real()+crt;
					index_l++;
					added=true;
				}
				if(index_r<eigen_count[ik_r] && (*eigen_pointer)[ik_r][index_r].real()<er_top){
					group_r[index_r]=groupNo;
					er_top=(*eigen_pointer)[ik_r][index_r].real()+crt;
					index_r++;
					added=true;
				}
				if(!added){
					break;
				}
			}
			groupNo++;
		}

		int group_count=groupNo;
		for(int igr=0; igr<group_count; igr++){
			int count_l=0;
			int count_r=0;
			int offset_l=-1;
			int offset_r=-1;
			// count the eigenstates belonging to the group
			for(int ib=0; ib<row_size; ib++){
				if(group_l[ib]==igr){
					if(offset_l<0){
						offset_l=ib;
					}
					count_l++;
				}
				if(group_r[ib]==igr){
					if(offset_r<0){
						offset_r=ib;
					}
					count_r++;
				}
			}
			if(count_l==count_r){
				// simple one-by-one connection
				for(int ib=0; ib<count_l; ib++){
					(*connection_pointer)[ik_l][offset_l+ib]=offset_r+ib;
				}
			}else if(count_l!=0 && count_r!=0){
				complex<double>* eigen_l=new complex<double>[count_l];
				complex<double>* eigen_r=new complex<double>[count_r];
				int* connection=new int[count_l];

				for(int ib=0; ib<count_l; ib++){
					eigen_l[ib]=(*eigen_pointer)[ik_l][offset_l+ib];
					connection[ib]=-1;
				}
				for(int ib=0; ib<count_r; ib++){
					eigen_r[ib]=(*eigen_pointer)[ik_r][offset_r+ib];
				}

				if(eigen_l[0].real()<PA_excitation_energy/Eh+1){
					resolve_connection(count_l, eigen_l, count_r, eigen_r, connection);
					for(int ib=0; ib<count_l; ib++){
						if(connection[ib]>=0){
							(*connection_pointer)[ik_l][offset_l+ib]=offset_r+connection[ib];
						}
					}
				}
				
				delete[] eigen_l;
				delete[] eigen_r;
				delete[] connection;
			}
		}
		/*
		printf("ikappaz=%d\n", ikappaz);
		for(int ib=0; ib<row_size; ib++){
			printf("%d | %3d %8.4f -> [%d]| %3d %8.4f\n", ib, group_l[ib], (*eigen_pointer)[ik_l][ib].real(), (*connection_pointer)[ik_l][ib], group_r[ib], (*eigen_pointer)[ik_r][ib].real());
		}
		printf("\n");*/
	}
	
	delete[] work;
	delete[] rwork;
	delete[] w;
	delete[] dispersion_real;
	delete[] sprintf_buffer2;
	delete[] group_l;
	delete[] group_r;

}
