#include <complex>
#include <omp.h>
#define _USE_MATH_DEFINES
#include <cmath>
#include <cstring>
#include <iostream>

// see Main_GUI/lib/physical_tools for the detail of the conversion
using namespace std;
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
				  pm1= (px-I*py)/sqrt(2);
				  pp0=  pz;
					pp1=-(px+I*py)/sqrt(2);
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
				  dm2= (dx2y2-I*dxy)/sqrt(2);
				  dm1= (dxz-  I*dyz)/sqrt(2);
				  dp0=  d3z2r2;
				  dp1=-(dxz+  I*dyz)/sqrt(2);
				  dp2= (dx2y2+I*dxy)/sqrt(2);
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
					fm3= (fx33xy2- I*f3yx2y3)/sqrt(2);
					fm2= (fzx2zy2- I*fxyz)/sqrt(2);
					fm1= (f5xy2xr2-I*f5yz2yr2)/sqrt(2);
					fp0=  f5z23r2;
					fp1=-(f5xy2xr2+I*f5yz2yr2)/sqrt(2);
					fp2= (fzx2zy2+ I*fxyz)/sqrt(2);
					fp3=-(fx33xy2+ I*f3yx2y3)/sqrt(2);
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

void spherical_harmonics(double* r, complex<double>* Ylm2){
	// r=[x, y, z]
  // x=r*sin(theta)*cos(phi)
  // y=r*sin(theta)*sin(phi)
  // sqrt(x^2+y^2)=r*sin(theta)
  // z=r*cos(theta)
	
	// Ylm2[5][9]

	complex<double>** Ylm=new complex<double>*[5];
	int i;
	for(i=0; i<5; i++){
		Ylm[i]=&Ylm2[9*i];
	}
	for(i=0; i<45; i++){
		Ylm2[i]=complex<double>(0, 0);
	}
    
	double r_length=sqrt(inner_product(r, r));
	double cosT=1;
	double sinT=0;
	if(r_length>1e-5){
		cosT=r[2]/r_length;
		sinT=sqrt(r[0]*r[0]+r[1]*r[1])/r_length;	
	}
	double cosP=1;
	double sinP=0;
	if(abs(sinT)>1e-5){
		cosP=r[0]/(r_length*sinT);
		sinP=r[1]/(r_length*sinT);
	}

	double cos2P=cosP*cosP-sinP*sinP;
	double sin2P=2*sinP*cosP;

	double cos3P=4*cosP*cosP*cosP-3*cosP;
	double sin3P=3*sinP-4*sinP*sinP*sinP;

	double cos4P=8*cosP*cosP*cosP*cosP-8*cosP*cosP+1;
	double sin4P=4*sinP*cosP*(2*cosP*cosP-1);

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
	Ylm[4][2]=(3.0/8.0)*sqrt(5.0/(2.0*M_PI))*(7.0*cosT*cosT-1.0)*complex<double>(cos2P, -sin2P);
	// g(-1): 3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)-i*sin(P))
	Ylm[4][3]=(3.0/8.0)*sqrt(5.0/M_PI)*(7.0*cosT*cosT-3.0)*sinT*cosT*complex<double>(cosP, -sinP);
	// g(0): 3/16 sqrt(1/pi) (35cos^4(T)-30cos^2(T)+3)
	Ylm[4][4]=(3.0/16.0)*sqrt(1.0/M_PI)*(35.0*cosT*cosT*cosT*cosT-30.0*cosT*cosT+3.0);
	// g(1): -3/8 sqrt(5/pi) (7cos^2(T)-3)sin(T)cos(T)(cos(P)+i*sin(P))
	Ylm[4][5]=-(3.0/8.0)*sqrt(5.0/M_PI)*(7.0*cosT*cosT-3.0)*sinT*cosT*complex<double>(cosP, sinP);
	// g(2): 3/8 sqrt(5/2pi) (7cos^2(T)-1)sin^2(T)(cos(2P)+i*sin(2P))
	Ylm[4][6]=(3.0/8.0)*sqrt(5.0/(2.0*M_PI))*(7.0*cosT*cosT-1.0)*complex<double>(cos2P, sin2P);
	// g(3): -3/8 sqrt(35/pi) sin^3(T)cos(T)(cos(3P)+i*sin(3P))
	Ylm[4][7]=-(3.0/8.0)*sqrt(35.0/M_PI)*sinT*sinT*sinT*cosT*complex<double>(cos3P, sin3P);
	// g(4): 3/16 sqrt(35/2pi) sin^4(T)(cos(4P)+i*sin(4P))
	Ylm[4][8]=(3.0/16.0)*sqrt(35.0/(2.0*M_PI))*sinT*sinT*sinT*sinT*complex<double>(cos4P, sin4P);

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
	// j_l(x)=(-1)^l x^l (1/x d/dx)^l sin(x)/x
	if(l==0){
		// sin(x)/x
		return sin(x)/x;
	}else if(l==1){
		// {sin(x)-x*cos(x)}/x^2
		return (sin(x)-x*cos(x))/(x*x);
	}else if(l==2){
		// {(3-x^2)sin(x)-3x*cos(x)}/x^3
		return ((3-x*x)*sin(x)-3*x*cos(x))/(x*x*x);
	}else if(l==3){
		// {(15-6x^2)sin(x)+(x^3-15x)cos(x)}/x^4
		return ((15-6*x*x)*sin(x)+(x*x*x-15*x)*cos(x))/(x*x*x*x);
	}else if(l==4){
		// {(x^4-45x^2+105)sin(x)+(10x^3-105x)cos(x)}/x^5
		return ((x*x*x*x-45*x*x+105)*sin(x)+(10*x*x*x-105*x)*cos(x))/(x*x*x*x*x);
	}else{
		return 0;
	}
}
