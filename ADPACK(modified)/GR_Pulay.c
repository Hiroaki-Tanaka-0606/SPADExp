/**********************************************************************
  GR_Pulay.c:
  
     GR_Pulay.c is a subroutine to achieve self-consistent field
     using the guaranteed-reduction Pulay method proposed by D.Bowler.

   Ref:
     D.R.Bowler and M.J.Gillan, Chem.Phys.Lett. 325, 475(2000)

  Log of GR_Pulay.c:

     31/May/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

static void Inverse(int n, double a[ASIZE9][ASIZE9],
                    double ia[ASIZE9][ASIZE9]);
static void GR(int state_num,int SCF_iter);

void GR_Pulay(int state_num, int SCF_iter)
{
  static int i,n,l;
  static double r,p,p2,sum;

  if (SCF_iter<=(Pulay_SCF-1))
    Simple_Mixing(state_num,rho[0],rho[1],rho[2],Rrho[2]);
  else if (SCF_iter==Pulay_SCF)
    GR(state_num,2);
  else 
    GR(state_num,SCF_iter-(Pulay_SCF-2));
}


void GR(int state_num,int SCF_iter)
{
  static int i,j,k,l,n,NumMix,NumSlide;
  static int SCFi,SCFj; 
  static int pSCF_iter;
  static double r,p,p2,s0,s1,s2;
  static double alden[ASIZE9];
  static double A[ASIZE9][ASIZE9];
  static double IA[ASIZE9][ASIZE9];
  static double B[ASIZE9][ASIZE9];
  static double x[ASIZE9];
  static double Av_dia,IAv_dia;
  static double dum1,dum2,bunsi,bunbo,sum;
  static double OptRrho[ASIZE1],rp,OptNorm_Rrho,coef_OptRDM;

  /****************************************************
                       GR_Pulay
  ****************************************************/

  if (SCF_iter==1){
    Simple_Mixing(state_num,rho[0],rho[1],rho[2],Rrho[2]);
  }
  else if (SCF_iter==2){
    for (i=0; i<Grid_Num; i++){
      Rrho[2][i] = rho[0][i] - rho[1][i];
      rho[2][i] = rho[1][i];
    }
    Simple_Mixing(state_num,rho[0],rho[1],rho[2],Rrho[2]);
  }
  else {

    if (SCF_iter%2==1){
      /****************************************************
                             Prho
      ****************************************************/
      for (i=0; i<Grid_Num; i++){
        Prho[i] = rho[0][i];
      }
      /****************************************************
                          Calc of Rrho1
      ****************************************************/
      for (i=0; i<Grid_Num; i++){
        Rrho[1][i] = rho[0][i] - rho[1][i];
      }

      /****************************************************
                              NormRD
      ****************************************************/

      s0 = Rrho[1][0]*Rrho[1][0]*MRV[0]*MRV[0]*MRV[0]
 	 + Rrho[1][Grid_Num-1]*Rrho[1][Grid_Num-1]*
	   MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

      s1 = 0.0;
      for (i=1; i<=(Grid_Num-2); i=i+2){
	s1 = s1 + Rrho[1][i]*Rrho[1][i]*MRV[i]*MRV[i]*MRV[i];
      }

      s2 = 0.0; 
      for (i=2; i<=(Grid_Num-3); i=i+2){
	s2 = s2 + Rrho[1][i]*Rrho[1][i]*MRV[i]*MRV[i]*MRV[i];
      }
      NormRD[0] = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

    }

    else {

      /****************************************************
                          Calc of RDM0
      ****************************************************/
      for (i=0; i<Grid_Num; i++){
        Rrho[0][i] = rho[0][i] - Prho[i];
      }

      if (((SCF_iter-2)/2 + 2 - 1)<Num_Mixing_pDM){
	NumMix   = (SCF_iter-2)/2 + 2 - 1;
	NumSlide = NumMix + 1;
      }
      else{
	NumMix   = Num_Mixing_pDM;
	NumSlide = NumMix;
      }

      /****************************************************
                         alpha from RDM
      ****************************************************/

      for (SCFi=0; SCFi<=NumMix; SCFi++){
	for (SCFj=SCFi; SCFj<=NumMix; SCFj++){
	  A[SCFi][SCFj] = 0.0;
          for (i=0; i<Grid_Num; i++){
            dum1 = Rrho[SCFi][i];
            dum2 = Rrho[SCFj][i];
            A[SCFi][SCFj] = A[SCFi][SCFj] + dum1*dum2;
          }
          A[SCFj][SCFi] = A[SCFi][SCFj];  
        }
      }

      Av_dia = A[0][0];
      IAv_dia = 1.0/Av_dia;
      for (SCFi=0; SCFi<=NumMix; SCFi++){
	for (SCFj=0; SCFj<=NumMix; SCFj++){
	  A[SCFi][SCFj] = A[SCFi][SCFj]*IAv_dia;
	}
      }

      Inverse(NumMix,A,IA);

      bunbo = 0.0;
      for (SCFi=0; SCFi<=NumMix; SCFi++){
	for (SCFj=0; SCFj<=NumMix; SCFj++){
	  bunbo = bunbo + IA[SCFi][SCFj];
	}
      } 

      sum = 0.0;
      for (SCFi=0; SCFi<=NumMix; SCFi++){
	bunsi = 0.0;
	for (SCFj=0; SCFj<=NumMix; SCFj++){
	  bunsi = bunsi + IA[SCFj][SCFi];
	}
	alden[SCFi] = bunsi/bunbo;
	sum = sum + alden[SCFi];
      }

      /****************************************************
                Calculate an optimized residual rho
      ****************************************************/

      OptNorm_Rrho = 0.0;

      for (i=0; i<Grid_Num; i++){
        OptRrho[i] = 0.0;
	for (pSCF_iter=0; pSCF_iter<=NumMix; pSCF_iter++){
          if (pSCF_iter==0){
            OptRrho[i] = OptRrho[i] + alden[pSCF_iter]*Rrho[pSCF_iter][i];
          }
          else{
            OptRrho[i] = OptRrho[i] + alden[pSCF_iter]*Rrho[pSCF_iter][i];
          } 
        }
        rp = MRV[i];
        OptNorm_Rrho = OptNorm_Rrho + OptRrho[i]*OptRrho[i]*rp*rp*rp*dx;
      }

      if (1.0e-2<=OptNorm_Rrho)
	coef_OptRDM = 0.2;
      else if (1.0e-4<=OptNorm_Rrho && OptNorm_Rrho<1.0e-2)
	coef_OptRDM = 0.3;
      else if (1.0e-6<=OptNorm_Rrho && OptNorm_Rrho<1.0e-4)
	coef_OptRDM = 0.4;
      else if (1.0e-8<=OptNorm_Rrho && OptNorm_Rrho<1.0e-6)
	coef_OptRDM = 0.5;
      else
	coef_OptRDM = 1.0;

      /*
      printf("OptNorm_Rrho=%15.12f\n",OptNorm_Rrho);
      */

      /****************************************************
                          Mixing of DM 
      ****************************************************/

      for (i=0; i<Grid_Num; i++){
        rho[0][i] = 0.0;
	for (pSCF_iter=0; pSCF_iter<=NumMix; pSCF_iter++){
          if (pSCF_iter==0){
            rho[0][i] = rho[0][i] + alden[pSCF_iter]*Prho[i];
          }
          else{
            rho[0][i] = rho[0][i] + alden[pSCF_iter]*rho[pSCF_iter][i];
          } 
        }

        /* Correction by the optimized rho */ 
        rho[0][i] = rho[0][i] + coef_OptRDM*OptRrho[i];
        if (rho[0][i]<0.0) rho[0][i] = 1.0e-15;
      }

      /****************************************************
                         Shift of rho
      ****************************************************/

      for (pSCF_iter=NumSlide; 0<pSCF_iter; pSCF_iter--){
        for (i=0; i<Grid_Num; i++){
          rho[pSCF_iter][i] = rho[pSCF_iter-1][i];
        }
      }

      /****************************************************
                     Shift of residual rho
      ****************************************************/

      for (pSCF_iter=NumSlide; 1<pSCF_iter; pSCF_iter--){
        for (i=0; i<Grid_Num; i++){
          Rrho[pSCF_iter][i] = Rrho[pSCF_iter-1][i];
        }
      }

    }
  }

}



void Inverse(int n, double a[ASIZE9][ASIZE9],
             double ia[ASIZE9][ASIZE9])
{
  /****************************************************
                  LU decomposition
                      0 to n
   NOTE:
   This routine does not take into account
   the reduction of rank.
  ****************************************************/

  static int i,j,k,l,m;
  static double w,sum;
  static double x[ASIZE9],y[ASIZE9];
  static double da[ASIZE9][ASIZE9];

  if (n==-1){
    for (i=0; i<=(ASIZE9-1); i++){
      for (j=0; j<=(ASIZE9-1); j++){
	a[i][j] = 0.0;
      }
    }
  }
  else{
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
      }
    }

    /****************************************************
                     LU factorization
    ****************************************************/

    for (k=0; k<=n-1; k++){
      w = 1.0/a[k][k];
      for (i=k+1; i<=n; i++){
	a[i][k] = w*a[i][k];
	for (j=k+1; j<=n; j++){
	  a[i][j] = a[i][j] - a[i][k]*a[k][j];
	}
      }
    }
    for (k=0; k<=n; k++){

      /****************************************************
                             Ly = b
      ****************************************************/

      for (i=0; i<=n; i++){
	if (i==k)
	  y[i] = 1.0;
	else
	  y[i] = 0.0;
	for (j=0; j<=i-1; j++){
	  y[i] = y[i] - a[i][j]*y[j];
	}
      }

      /****************************************************
                             Ux = y 
      ****************************************************/

      for (i=n; 0<=i; i--){
	x[i] = y[i];
	for (j=n; (i+1)<=j; j--){
	  x[i] = x[i] - a[i][j]*x[j];
	}
	x[i] = x[i]/a[i][i];
	ia[i][k] = x[i];
      }
    }

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	a[i][j] = da[i][j];
      }
    }
  }
}

