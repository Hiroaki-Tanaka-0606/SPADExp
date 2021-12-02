#include <complex>
#include <omp.h>
#include <cmath>

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
				  pm1=(px-I*py)/sqrt(2);
				  pp0=pz;
					pp1=-(px+I*py)/sqrt(2);
					LCAO[ik][ib][0][id]=pm1;
					LCAO[ik][ib][1][id]=pp0;
					LCAO[ik][ib][2][id]=pp1;
					break;
				case 5:
					// d orbital
				  d3z2r2=complex<double>(LCAO_raw[idx1+idx2+2*id], LCAO_raw[idx1+idx2+2*id+1]);
				  dx2y2=complex<double>(LCAO_raw[idx1+idx2+sd+2*id], LCAO_raw[idx1+idx2+sd+2*id+1]);
				  dxy=complex<double>(LCAO_raw[idx1+idx2+2*sd+2*id], LCAO_raw[idx1+idx2+2*sd+2*id+1]);
				  dxz=complex<double>(LCAO_raw[idx1+idx2+3*sd+2*id], LCAO_raw[idx1+idx2+3*sd+2*id+1]);
				  dyz=complex<double>(LCAO_raw[idx1+idx2+4*sd+2*id], LCAO_raw[idx1+idx2+4*sd+2*id+1]);
				  dm2=(dx2y2-I*dxy)/sqrt(2);
				  dm1=(dxz-I*dyz)/sqrt(2);
				  dp0=d3z2r2;
				  dp1=-(dxz+I*dyz)/sqrt(2);
				  dp2=(dx2y2+I*dxy)/sqrt(2);
					LCAO[ik][ib][0][id]=dm2;
					LCAO[ik][ib][1][id]=dm1;
					LCAO[ik][ib][2][id]=dp0;
					LCAO[ik][ib][3][id]=dp1;
					LCAO[ik][ib][4][id]=dp2;
					break;
				case 7:
					// f orbital
					f5z23r2=complex<double>(LCAO_raw[idx1+idx2+2*id], LCAO_raw[idx1+idx2+2*id+1]);
				  f5xy2xr2=complex<double>(LCAO_raw[idx1+idx2+sd+2*id], LCAO_raw[idx1+idx2+sd+2*id+1]);
				  f5yz2yr2=complex<double>(LCAO_raw[idx1+idx2+2*sd+2*id], LCAO_raw[idx1+idx2+2*sd+2*id+1]);
				  fzx2zy2=complex<double>(LCAO_raw[idx1+idx2+3*sd+2*id], LCAO_raw[idx1+idx2+3*sd+2*id+1]);
				  fxyz=complex<double>(LCAO_raw[idx1+idx2+4*sd+2*id], LCAO_raw[idx1+idx2+4*sd+2*id+1]);
					fx33xy2=complex<double>(LCAO_raw[idx1+idx2+5*sd+2*id], LCAO_raw[idx1+idx2+5*sd+2*id+1]);
					f3yx2y3=complex<double>(LCAO_raw[idx1+idx2+6*sd+2*id], LCAO_raw[idx1+idx2+6*sd+2*id+1]);
					fm3=(fx33xy2-I*f3yx2y3)/sqrt(2);
					fm2=(fzx2zy2-I*fxyz)/sqrt(2);
					fm1=(f5xy2xr2-I*f5yz2yr2)/sqrt(2);
					fp0=f5z23r2;
					fp1=-(f5xy2xr2+I*f5yz2yr2)/sqrt(2);
					fp2=(fzx2zy2+I+fxyz)/sqrt(2);
					fp3=-(fx33xy2+I*f3yx2y3)/sqrt(2);
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
