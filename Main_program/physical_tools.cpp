#include <complex>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>
#include "log.hpp"
#include "variables_ext.hpp"
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
double interpolate_potential(double* r, int* count, double* cube, double* rec_cell){
	double p[3]; // fractional
	double q[3]; // non-integer index
	int qf[3]; // integer index (floored)
	double** rc=new double*[3];
	for(int i=0; i<3; i++){
		rc[i]=&rec_cell[i*3];
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

void calc_residual_vector(int g_count, int v_count, double* A, double* B, double* x, double* res);
double calc_residual_error(int g_count, int v_count, double* A, double* B, double* x);
void calc_gradient_vector(int g_count, int v_count, double* A, double* B, double* x, double* gr);
double calc_norm(int count, double* vector);

complex<double>** alloc_zmatrix(int n);
complex<double>** alloc_zmatrix(int m, int n);
void delete_zmatrix(complex<double>** mat);

complex<double>*** alloc_zpmatrix(int n);
void delete_zpmatrix(complex<double>*** mat);

double** alloc_dmatrix(int n);
double** alloc_dmatrix(int m, int n);
void delete_dmatrix(double** mat);

// In the Numerov method, the solving direction is downward
// f_{i-1}=(1+h^2/12 a_{i-1})^-1 * [ 2(1-5h^2/12 a_i)*f_i - (1+h^2/12 a_{i+1})*f_{i+1} ]
void solve_final_state(double Ekin, double* k_para, double kz, int g_count, int z_count, double dz, int V00_index, complex<double>** Vgg_buffer, double* g_vec_buffer, complex<double>* FP_loc_buffer){
	char* sprintf_buffer2=new char[Log_length+1];
	// composite matrix
	complex<double>*** Vgg=alloc_zpmatrix(g_count);
	for(int ig1=0; ig1<g_count; ig1++){
		for(int ig2=0; ig2<g_count; ig2++){
			Vgg[ig1][ig2]=Vgg_buffer[ig1*g_count+ig2];
		}
	}
	double g_vec[g_count][3];
	for(int ig=0; ig<g_count; ig++){
		for(int ix=0; ix<3; ix++){
			g_vec[ig][ix]=g_vec_buffer[ig*3+ix];
		}
	}

	
	double kzz0=kz*dz*(z_count-1);
	double kzz0h=kz*dz*z_count;

	if(PA_FPFS_Numerov){
		// Numerov and Euler
		complex<double>* FP_loc[g_count];
		for(int ig=0; ig<g_count; ig++){
			FP_loc[ig]=&FP_loc_buffer[ig*z_count];
		}

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
		double diagonal_part[g_count];
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
		complex<double> a_matrix[z_count][g_count][g_count];
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
		complex<double> left_matrix[g_count][g_count]; // transposed!!
		complex<double> right_vector[g_count];
		int INFO=0;
		int NRHS=1;
		int IPIV[g_count];
		for(int im1=z_count-3; im1>=0; im1--){
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
	}else{
		// Linear equations
		int eq_dim=z_count*g_count;
		complex<double>** left_matrix=alloc_zmatrix(eq_dim); // transposed
		for(int i=0; i<eq_dim; i++){
			for(int j=0; j<eq_dim; j++){
				left_matrix[j][i]=complex<double>(0.0, 0.0);
			}
		}

		bool* vacuum_ok=new bool[g_count];
		int vacuum_ok_count=0;
		int* vacuum_ok_index=new int[g_count];
		
		double kpg[3];
		for(int ig=0; ig<g_count; ig++){
			for(int p=0; p<3; p++){
				kpg[p]=k_para[p]+g_vec[ig][p];
			}
			double Ekpg=inner_product(kpg, kpg)/2.0;
			vacuum_ok[ig]=Ekpg<Ekin;
			// printf("vacuum_ok[%2d]=%s\n", ig, vacuum_ok[ig]?"true":"false");
			if(vacuum_ok[ig]){
				vacuum_ok_index[vacuum_ok_count]=ig;
				vacuum_ok_count++;
			}
		}
		// left matrix
		/// M
		for(int ig1=0; ig1<g_count; ig1++){
			for(int p=0; p<3; p++){
				kpg[p]=k_para[p]+g_vec[ig1][p];
			}
			for(int ig2=0; ig2<g_count; ig2++){
				for(int iz=0; iz<z_count; iz++){
					int index1=ig1*z_count+iz;
					int index2=ig2*z_count+iz;
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
			for(int iz1=0; iz1<z_count; iz1++){
				for(int iz2=0; iz2<z_count; iz2++){
					int index1=ig*z_count+iz1;
					int index2=ig*z_count+iz2;
					if(iz1==iz2){
						left_matrix[index2][index1]-=2.0/dz/dz;
					}
					if(iz1-1==iz2 || iz1+1==iz2){
						left_matrix[index2][index1]+=1.0/dz/dz;
					}
				}
			}
		}

		/*
		// LU decomposition
		int info;
		int* ipiv=new int[eq_dim];
		int nrhs=1;

		cout << "LU start" << endl;
		zgetrf_(&eq_dim, &eq_dim, &left_matrix[0][0], &eq_dim, &ipiv[0], &info);
		if(info!=0){
			write_log((char*)"zgetrf failed");
			return;
		}
		cout << "LU end" << endl;

		// inverse matrix calculation
		cout << "Inverse 1 start" << endl;
		int lwork=-1;
		complex<double> work_dummy;
		zgetri_(&eq_dim, &left_matrix[0][0], &eq_dim, &ipiv[0], &work_dummy, &lwork, &info);
		if(info!=0){
			write_log((char*)"Work space calculation failed");
			return;
		}
		cout << "Inverse 1 end" << endl;
		lwork=round(abs(work_dummy));
		complex<double>* work=new complex<double>[lwork];
		cout << "Inverse 2 start" << endl;
		zgetri_(&eq_dim, &left_matrix[0][0], &eq_dim, &ipiv[0], &work[0], &lwork, &info);
		if(info!=0){
			write_log((char*)"zgetri failed");
			return;
		}
		cout << "Inverse 2 end" << endl;
		delete[] work;*/
		/*
		for(int i=0; i<eq_dim; i++){
			for(int j=0; j<eq_dim; j++){
				printf("(%10.2e, %10.2e) ", left_matrix[j][i].real(), left_matrix[j][i].imag());
			}
			printf("\n");
			}*/

		// 0~g: for RL, g~2g: for RR
		complex<double>** right_matrix=alloc_zmatrix(2*g_count, eq_dim); //transposed
		for(int igz1=0; igz1<eq_dim; igz1++){
			for(int ig2=0; ig2<2*g_count; ig2++){
				right_matrix[ig2][igz1]=complex<double>(0.0, 0.0);
			}
		}
		for(int ig1=0; ig1<g_count; ig1++){
			int index1=ig1*z_count;
			int index2=ig1;
			right_matrix[index2][index1]=complex<double>(1.0, 0.0);
			index1=(ig1+1)*z_count-1;
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
		int info;
		int* ipiv=new int[eq_dim];
		int nrhs=2*g_count;
		zgesv_(&eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &right_matrix[0][0], &eq_dim, &info);
		if(info!=0){
			write_log((char*)"zgesv failed");
			return;
		}

		// extract the RL and RR matrices
		complex<double>** H_RL=alloc_zmatrix(g_count); // inverted
		complex<double>** H_RR=alloc_zmatrix(g_count); // inverted
		int index1, index2;
		for(int ig1=0; ig1<g_count; ig1++){
			for(int ig2=0; ig2<g_count; ig2++){
				index1=(ig1+1)*z_count-1;
				index2=ig2;
				H_RL[ig2][ig1]=right_matrix[index2][index1];
				
				index1=(ig1+1)*z_count-1;
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
		complex<double>* right_vector=new complex<double>[g_count];
		for(int ig=0; ig<g_count; ig++){
			if(ig==V00_index){
				right_vector[ig]=-dz*dz*complex<double>(cos(kzz0), sin(kzz0));
			}else{
				right_vector[ig]=complex<double>(0.0, 0.0);
			}
		}
		complex<double>* right_vector2=new complex<double>[g_count];
		for(int ig=0; ig<g_count; ig++){
			if(ig==V00_index){
				right_vector2[ig]=-complex<double>(cos(kzz0h), sin(kzz0h));
			}else{
				right_vector2[ig]=complex<double>(0.0, 0.0);
			}
		}
		
		char notrans='N';
		complex<double> alpha(1.0, 0.0);
		complex<double> beta(1.0, 0.0);
		int inc=1;
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

		// by steepest-decent
		/*
		int g2_count=g_count*2;
		int v2_count=vacuum_ok_count*2;
		double** A_real=alloc_dmatrix(v2_count, g2_count); //transposed
		for(int ig=0; ig<g_count; ig++){
			for(int iv=0; iv<vacuum_ok_count; iv++){
				int ig2=vacuum_ok_index[iv];
				A_real[iv*2][ig*2]=H_RL[ig2][ig].real();
				A_real[iv*2+1][ig*2]=-H_RL[ig2][ig].imag();
				A_real[iv*2][ig*2+1]=H_RL[ig2][ig].imag();
				A_real[iv*2+1][ig*2+1]=H_RL[ig2][ig].real();
			}
		}
		
		printf("A_real\n");
		for(int igp=0; igp<g2_count; igp++){
			for(int ivp=0; ivp<v2_count; ivp++){
				printf("%10.2e ", A_real[ivp][igp]);
			}
			printf("\n");
			}
		double* B_real=new double[g2_count];
		for(int ig=0; ig<g_count; ig++){
			B_real[ig*2]=right_vector[ig].real();
			B_real[ig*2+1]=right_vector[ig].imag();
		}
		
		printf("B_real\n");
		for(int igp=0; igp<g2_count; igp++){
			printf("%10.2e\n", B_real[igp]);
			}

		double* x_real=new double[v2_count];
		for(int iv=0; iv<v2_count; iv++){
			x_real[iv]=0.0;
		}

		double* nabla=new double[v2_count];
		double* Anabla=new double[g2_count];

		double res_init=calc_residual_error(g2_count, v2_count, &A_real[0][0], &B_real[0], &x_real[0]);
		sprintf(sprintf_buffer2, "Residual error[%3d] = %10.2e", 0, res_init);
		write_log(sprintf_buffer2);

		int trial;
		double alpha2=1.0;
		double beta2=0.0;
		double sd_alpha;
		double res;
		double res_prev=res_init;
		for(trial=0; trial<100; trial++){
			calc_gradient_vector(g2_count, v2_count, &A_real[0][0], &B_real[0], &x_real[0], &nabla[0]);
			dgemv_(&notrans, &g2_count, &v2_count, &alpha2, &A_real[0][0], &g2_count, &nabla[0], &inc, &beta2, &Anabla[0], &inc);
			sd_alpha=-calc_norm(v2_count, nabla)/calc_norm(g2_count, Anabla);
			for(int iv=0; iv<v2_count; iv++){
				x_real[iv]+=sd_alpha*nabla[iv];
			}
			res=calc_residual_error(g2_count, v2_count, &A_real[0][0], &B_real[0], &x_real[0]);
			// printf("Residual error[%3d] = %10.2e\n", trial+1, res);
			if(res/res_init<1e-12){
				break;
			}
		}
		sprintf(sprintf_buffer2, "Last residual error[%3d] = %10.2e", trial+1, res);
		write_log(sprintf_buffer2);*/
		// by steepest decent finished
		
		// by zgels
		complex<double>** A_complex=alloc_zmatrix(vacuum_ok_count, g_count);
		for(int ig=0; ig<g_count; ig++){
			for(int iv=0; iv<vacuum_ok_count; iv++){
				int ig2=vacuum_ok_index[iv];
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
		zgels_(&notrans, &g_count, &vacuum_ok_count, &nrhs, &A_complex[0][0], &g_count, &B_complex[0], &g_count, &work_dummy, &lwork, &info);
		if(info!=0){
			write_log((char*)"Work space calculation failed");
			return;
		}
		lwork=round(abs(work_dummy));
		complex<double>* work=new complex<double>[lwork];
		zgels_(&notrans, &g_count, &vacuum_ok_count, &nrhs, &A_complex[0][0], &g_count, &B_complex[0], &g_count, &work[0], &lwork, &info);
		if(info!=0){
			write_log((char*)"zgels failed");
			return;
		}
		delete[] work;
		// by zgels finished
		
		for(int ig=0; ig<g_count; ig++){
			right_vector[ig]=complex<double>(0.0, 0.0);
		}
		for(int iv=0; iv<vacuum_ok_count; iv++){
			// right_vector[vacuum_ok_index[iv]]=complex<double>(x_real[2*iv], x_real[2*iv+1]);
			right_vector[vacuum_ok_index[iv]]=B_complex[iv];
		}

		/*
		printf("Right(solution)\n");
		for(int ig=0; ig<g_count; ig++){
			printf("(%10.2e, %10.2e)\n", right_vector[ig].real(), right_vector[ig].imag());
		}
		printf("\n");*/
		/*
		printf("Right(solution) by zgels\n");
		for(int iv=0; iv<vacuum_ok_count; iv++){
			printf("(%10.2e, %10.2e)\n", B_complex[iv].real(), B_complex[iv].imag());
			}*/
		
		double zgels_norm=0.0;
		for(int iv=vacuum_ok_count; iv<g_count; iv++){
			zgels_norm+=norm(B_complex[iv]);
		}
		sprintf(sprintf_buffer2, "zgels norm: %10.2e", zgels_norm);
		write_log(sprintf_buffer2);

		// printf("f_G(z)\n");

		// obtain f_G(z)
		//complex<double>* right_vector_all=new complex<double>[eq_dim];
		for(int igz=0; igz<eq_dim; igz++){
			//right_vector_all[igz]=complex<double>(0.0, 0.0);
			FP_loc_buffer[igz]=complex<double>(0.0, 0.0);
		}
		for(int ig=0; ig<g_count; ig++){
			int index=ig*z_count;
			FP_loc_buffer[index]=-right_vector[ig]/dz/dz;
			if(ig==V00_index){
				index=(ig+1)*z_count-1;
				FP_loc_buffer[index]=-1.0/dz/dz*complex<double>(cos(kzz0h), sin(kzz0h));
			}
		}
		
		// beta=complex<double>(0.0, 0.0);
		//for(int igz=0; igz<eq_dim; igz++){
			//printf("(%8.4f, %8.4f)\n", right_vector_all[igz].real(), right_vector_all[igz].imag());
		//}
		// zgemv_(&notrans, &eq_dim, &eq_dim, &alpha, &left_matrix[0][0], &eq_dim, &right_vector_all[0], &inc, &beta, &FP_loc_buffer[0], &inc);
		zgetrs_(&notrans, &eq_dim, &nrhs, &left_matrix[0][0], &eq_dim, &ipiv[0], &FP_loc_buffer[0], &eq_dim, &info);
		if(info!=0){
			write_log((char*)"zgetrs failed");
		}
		
		//for(int ig=0; ig<g_count; ig++){
		//for(int iz=0; iz<z_count; iz++){
		//int index=ig*z_count+iz;
		//printf("(%8.4f, %8.4f)\n", FP_loc_buffer[index].real(), FP_loc_buffer[index].imag());
		//}
		//printf("\n");
		//}
		
		delete_zmatrix(left_matrix);
		//delete_dmatrix(A_real);
		delete_zmatrix(A_complex);
		delete_zmatrix(right_matrix);
		delete[] vacuum_ok;
		delete[] vacuum_ok_index;
		delete[] ipiv;
		delete[] right_vector;
		delete[] right_vector2;
		//delete[] B_real;
		//delete[] x_real;
		//delete[] nabla;
		//delete[] Anabla;
		//delete[] right_vector_all;
		delete[] B_complex;
		
	}

	delete[] sprintf_buffer2;
	delete_zpmatrix(Vgg);
}

complex<double>** alloc_zmatrix(int n){
	complex<double>* buffer=new complex<double>[n*n];
	complex<double>** mat=new complex<double>*[n];
	for(int i=0; i<n; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

complex<double>** alloc_zmatrix(int m, int n){
	complex<double>* buffer=new complex<double>[m*n];
	complex<double>** mat=new complex<double>*[m];
	for(int i=0; i<m; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

complex<double>*** alloc_zpmatrix(int n){
	complex<double>** buffer=new complex<double>*[n*n];
	complex<double>*** mat=new complex<double>**[n];
	for(int i=0; i<n; i++){
		mat[i]=&buffer[i*n];
	}
	return mat;
}

void delete_zpmatrix(complex<double>*** mat){
	delete[] mat[0];
	delete[] mat;
}

void delete_zmatrix(complex<double>** mat){
	delete[] mat[0];
	delete[] mat;
}

// calculate res=Ax-B
void calc_residual_vector(int g_count, int v_count, double* A, double* B, double* x, double* res){
	char notrans='N';
	double alpha=1.0;
	double beta=-1.0;
	int inc=1;
	for(int i=0; i<g_count; i++){
		res[i]=B[i];
	}
	dgemv_(&notrans, &g_count, &v_count, &alpha, A, &g_count, x, &inc, &beta, res, &inc);
}

double calc_residual_error(int g_count, int v_count, double* A, double* B, double* x){
	double* res=new double[g_count];
	calc_residual_vector(g_count, v_count, A, B, x, &res[0]);
	double re=0.0;
	for(int i=0; i<g_count; i++){
		re+=res[i]*res[i]/2.0;
	}
	delete[] res;
	return re;
}

void calc_gradient_vector(int g_count, int v_count, double* A, double* B, double* x, double* gr){
	double* res=new double[g_count];
	calc_residual_vector(g_count, v_count, A, B, x, &res[0]);
	double alpha=1.0;
	double beta=0.0;
	char trans='T';
	int inc=1;
	dgemv_(&trans, &g_count, &v_count, &alpha, A, &g_count, &res[0], &inc, &beta, gr, &inc);
	delete[] res;
}

double calc_norm(int count, double* vector){
	double re=0.0;
	for(int i=0; i<count; i++){
		re+=vector[i]*vector[i];
	}
	return re;
}

/*
int conv_index_2D(int i, int j, int* size);

void solve_final_state_real_space(double Ekin, double* k_para, double kz, double* atom_cell_buffer, int* VKS_count, double* VKS_buffer, complex<double>* FP_loc_buffer){
	int z_size=VKS_count[0];
	int x_size=VKS_count[1];
	int y_size=VKS_count[2];
	int xy_size=x_size*y_size;
	// prepare the multi-dimension array
	double** VKS=new double*[z_size];
	for(int i=0; i<z_size; i++){
		VKS[i]=&VKS_buffer[i*xy_size];
	}
	complex<double>** FP_loc=new complex<double>*[z_size];
	for(int i=0; i<z_size; i++){
		FP_loc[i]=&FP_loc_buffer[i*xy_size];
		if(i==z_size-1 || i==z_size-2){
			for(int j=0; j<xy_size; j++){
				FP_loc[i][j]=complex<double>(1.0, 0.0);
			}
		}
	}
	double atom_cell[3][3];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			atom_cell[i][j]=atom_cell_buffer[i*3+j];
		}
	}

	// coefficients for the differential parts
	double x1=atom_cell[1][0]/x_size;
	double y1=atom_cell[1][1]/x_size;
	double x2=atom_cell[2][0]/y_size;
	double y2=atom_cell[2][1]/y_size;

	double det2=x1*y2-x2*y1;
	double t1=(y2*k_para[0]-x2*k_para[1])/det2;
	double t2=(-y1*k_para[0]+x1*k_para[1])/det2;
	// debug
	// printf("t1=%10.2e, t2=%10.2e\n", t1, t2);

	double s_matrix[3][3]; // inverted
	s_matrix[0][0]=x1*x1;
	s_matrix[1][0]=x2*x2;
	s_matrix[2][0]=(x1+x2)*(x1+x2);
	s_matrix[0][1]=2.0*x1*y1;
	s_matrix[1][1]=2.0*x2*y2;
	s_matrix[2][1]=2.0*(x1+x2)*(y1+y2);
	s_matrix[0][2]=y1*y1;
	s_matrix[1][2]=y2*y2;
	s_matrix[2][2]=(y1+y2)*(y1+y2);

	// LU decomposition
	int info;
	int ipiv[3];
	int s_size=3;
	int nrhs=1;
	double s_vec[3];
	double s1, s2, s3;
	double s_matrix2[3][3];
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			s_matrix2[i][j]=s_matrix[i][j];
		}
	}
		
	s_vec[0]=1.0;
	s_vec[1]=0.0;
	s_vec[2]=1.0;
	dgesv_(&s_size, &nrhs, &s_matrix2[0][0], &s_size, &ipiv[0], &s_vec[0], &s_size, &info);
	if(info==0){
		s1=s_vec[0];
		s2=s_vec[1];
		s3=s_vec[2];
	}else{
		write_log((char*)"dgesv failed");
		return;
	}

	// debug
	
	// double test1=s_matrix[0][0]*s1+s_matrix[1][0]*s2+s_matrix[2][0]*s3;
	// double test2=s_matrix[0][1]*s1+s_matrix[1][1]*s2+s_matrix[2][1]*s3;
	// double test3=s_matrix[0][2]*s1+s_matrix[1][2]*s2+s_matrix[2][2]*s3;
	// printf("s1=%10.2e, s2=%10.2e, s3=%10.2e\n", s1, s2, s3);
	// printf("te1=%10.2e, te2=%10.2e, te3=%10.2e\n", test1, test2, test3);

	// prepare the A matrix
	double k_para_norm=inner_product(k_para, k_para);
	complex<double> a_matrix_wo_V[xy_size][xy_size]; // transposed

	/// diagonal
	for(int i=0; i<xy_size; i++){
		for(int j=0; j<xy_size; j++){
			a_matrix_wo_V[j][i]=complex<double>(0.0, 0.0);
			if(i==j){
				a_matrix_wo_V[j][i]=-2.0*(-Ekin+k_para_norm/2.0);
			}
		}
	}

	/// derivative
	
	complex<double> im(0.0, 1.0);
	for(int i=0; i<x_size; i++){
		for(int j=0; j<y_size; j++){
			int index_ij  =conv_index_2D(i  , j  , VKS_count);
			int index_ipj =conv_index_2D(i+1, j  , VKS_count);
			int index_imj =conv_index_2D(i-1, j  , VKS_count);
			int index_ijp =conv_index_2D(i  , j+1, VKS_count);
			int index_ijm =conv_index_2D(i  , j-1, VKS_count);
			int index_ipjp=conv_index_2D(i+1, j+1, VKS_count);
			int index_imjm=conv_index_2D(i-1, j-1, VKS_count);

			/// second derivative
			a_matrix_wo_V[index_ij][index_ij]-=2.0*(s1+s2+s3);
			a_matrix_wo_V[index_ipj][index_ij]+=s1;
			a_matrix_wo_V[index_imj][index_ij]+=s1;
			a_matrix_wo_V[index_ijp][index_ij]+=s2;
			a_matrix_wo_V[index_ijm][index_ij]+=s2;
			a_matrix_wo_V[index_ipjp][index_ij]+=s3;
			a_matrix_wo_V[index_imjm][index_ij]+=s3;

			/// first derivative
			a_matrix_wo_V[index_ipj][index_ij]+=t1*im;
			a_matrix_wo_V[index_imj][index_ij]-=t1*im;
			a_matrix_wo_V[index_ijp][index_ij]+=t2*im;
			a_matrix_wo_V[index_ijm][index_ij]-=t2*im;
		}
		}

	complex<double> left_matrix[xy_size][xy_size]; // transposed
	complex<double> right_vector[xy_size];
	double dz=atom_cell[0][2]/z_size;
	// dz=0.0;
	complex<double> alpha;
	complex<double> beta(1.0, 0.0);
	int inc=1;
	char notrans='N';
	int ipiv2[xy_size];
	
	for(int im1=z_size-3; im1>=0; im1--){
		int i=im1+1;
		int ip1=i+1;
		// left matrix
		for(int ixy1=0; ixy1<xy_size; ixy1++){
			for(int ixy2=0; ixy2<xy_size; ixy2++){
				left_matrix[ixy1][ixy2]=dz*dz/12.0*a_matrix_wo_V[ixy1][ixy2];
				if(ixy1==ixy2){
					left_matrix[ixy1][ixy2]+=1.0+dz*dz/12.0*(-2.0*VKS[im1][ixy1]);
				}
			}
		}
		// right vector
		/// diagonal
		for(int ixy=0; ixy<xy_size; ixy++){
			right_vector[ixy]=2.0*(1.0-5.0/12.0*dz*dz*(-2.0*VKS[i][ixy]))*FP_loc[i][ixy];
			right_vector[ixy]-=(1.0+1.0/12.0*dz*dz*(-2.0*VKS[ip1][ixy]))*FP_loc[ip1][ixy];
		}
		/// off-diagonal
		alpha=complex<double>(-5.0/6.0*dz*dz, 0.0);
		zgemv_(&notrans, &xy_size, &xy_size, &alpha, &a_matrix_wo_V[0][0], &xy_size, &FP_loc[i][0], &inc, &beta, &right_vector[0], &inc);
		
		alpha=complex<double>(-1.0/12.0*dz*dz, 0.0);
		zgemv_(&notrans, &xy_size, &xy_size, &alpha, &a_matrix_wo_V[0][0], &xy_size, &FP_loc[ip1][0], &inc, &beta, &right_vector[0], &inc);

		// solve
		zgesv_(&xy_size, &nrhs, &left_matrix[0][0], &xy_size, &ipiv2[0], &right_vector[0], &xy_size, &info);
		if(info==0){
			// ok
			for(int ixy=0; ixy<xy_size; ixy++){
				FP_loc[im1][ixy]=right_vector[ixy];
				printf("(%8.4f %8.4f) ", FP_loc[im1][ixy].real(), FP_loc[im1][ixy].imag());
			}
			printf("\n");
		}else{
			write_log((char*)"zgesv failed");
			return;
		}
	}
}


int conv_index_2D(int i, int j, int* size){
	// i is for size[1]
	// j is for size[2]
	while(i<0){
		i+=size[1];
	}
	while(j<0){
		j+=size[2];
	}
	i=i%size[1];
	j=j%size[2];
	return i*size[2]+j;
}*/


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

void solve_nonlocal_wfn(double Ekin, int l, int r_count, double* r_arr, double* V_loc, int VPS_l_length, int* VPS_l, double* VPS_E_buffer, double** VPS_nonloc_buffer, double* psi){
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
			if(ix==0){
				/*
				if(iy==0 || iy==2){
					matrix[iy][ix]-=1.0/dx/dx;
				}else if(iy==1){
					matrix[iy][ix]-=-2.0/dx/dx;
					}*/
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
			}else{
				if(iy==ix-1 || iy==ix+1){
					matrix[iy][ix]-=1.0/dx/dx;
				}else if(iy==ix){
					matrix[iy][ix]-=-2.0/dx/dx;
				}
			}
			// first derivative
			if(ix==0){/*
				if(iy==0){
					matrix[iy][ix]+=-1.0/dx;
				}else if(iy==1){
				matrix[iy][ix]+=1.0/dx;
				}*/
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
			}else{
				if(iy==ix-1){
					matrix[iy][ix]+=-0.5/dx;
				}else if(iy==ix+1){
					matrix[iy][ix]+=0.5/dx;
				}
			}
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
	double** inverse=alloc_dmatrix(r_count);
	for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
			inverse[iy][ix]=matrix[iy][ix];
		}
	}

	// LU decomposition
	int info;
	int ipiv[r_count];
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
	/*
	double mm[r_count][r_count];
	double alpha=1.0;
	double beta=0.0;
	char trans='N';
	dgemm_(&trans, &trans, &r_count, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &inverse[0][0], &r_count, &beta, &mm[0][0], &r_count);
	for(int ix=0; ix<r_count; ix++){
		for(int iy=0; iy<r_count; iy++){
			printf("%4.1f ", mm[iy][ix]);
		}
		printf("\n");
		}*/

	
	double dr_dummy[r_count];
	for(int ir=0; ir<r_count; ir++){
		dr_dummy[ir]=1.0;
	}
	
	// initial vector
	double psi_vector[r_count];
	double hpsi_vector[r_count];
	for(int ir=0; ir<r_count; ir++){
		psi_vector[ir]=1.0;
	}
	psi_normalize(r_count, dr, psi_vector);
	/*
	for(int ir=0; ir<r_count; ir++){
		printf("%10.6f %10.6f\n", r_arr[ir], psi_vector[ir]);
		}*/

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

	/*
	double hpsi_norm=0.0;
	dgemv_(&trans, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &hpsi_vector[0], &inc);
	for(int ir=0; ir<r_count; ir++){
		hpsi_norm+=hpsi_vector[ir]*hpsi_vector[ir]*dr[ir];
	}
	printf("Norm: %f\n", hpsi_norm);*/

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
	sprintf(sprintf_buffer2, "CG trials: %6d, CG norm: %10.4e\n", trial, cg_norm);
	write_log(sprintf_buffer2);

	dgemv_(&no_trans, &r_count, &r_count, &alpha, &matrix[0][0], &r_count, &psi_vector[0], &inc, &beta, &hpsi_vector[0], &inc);
	// copy & normalization so that the edge is 1.0
	double psi_edge=psi_vector[r_count-1];
	for(int ir=0; ir<r_count; ir++){
		psi[ir]=psi_vector[ir]/psi_edge;
		//printf("%10.6f %10.6f %10.6f\n", r_arr[ir], psi[ir], hpsi_vector[ir]/psi_edge);
	}
	delete[] VPS_E;
	delete[] VPS_nonloc;
	delete[] sprintf_buffer2;
	delete_dmatrix(matrix);
	delete_dmatrix(inverse);
	delete_dmatrix(dr_matrix);
	delete_dmatrix(drH_matrix);
	delete_dmatrix(A_matrix);
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

double** alloc_dmatrix(int n){
	double* buffer=new double[n*n];
	double** ret=new double*[n];
	for(int i=0; i<n; i++){
		ret[i]=&buffer[i*n];
	}
	return ret;
}

double** alloc_dmatrix(int m, int n){
	double* buffer=new double[m*n];
	double** ret=new double*[m];
	for(int i=0; i<m; i++){
		ret[i]=&buffer[i*n];
	}
	return ret;
}

void delete_dmatrix(double** mat){
	delete[] mat[0];
	delete[] mat;
}
