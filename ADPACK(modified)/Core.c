/**********************************************************************
  Core.c:

     Core.c is a subroutine to calculate a nuclear Coulomb potential
     or a potential modified by an error function.

  Log of Core.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

#define LSIZE1  10

static void RIH(int n, double a[LSIZE1][LSIZE1], double b[LSIZE1]);

void Core(int SCF, double Z, double Vcore[ASIZE1])
{
  static int i,j,po;
  static double a,b,c,LA,r0,rc,r,dr,Gnew,Gold,r1,SCF2;
  static double cd,d,AC[LSIZE1][LSIZE1],Coes[LSIZE1];

  if (Calc_Type==0 || Calc_Type==1 || Calc_Type==3){

    Vinf = 0.0;
    for (i=0; i<Grid_Num; i++){
      Vcore[i] = -Z/MRV[i];
    }

    /* consider a large wall height in case of GVPS */
    if (Calc_Type==1 && VPP_switch==4){
      Vinf = 1.0e+30;
    }
  }
  else {

    SCF2 = 5.0; 
    a = 5.0;

    c = height_wall;
    d = rising_edge;

    if (SCF<=SCF2){
      cd = c*SCF/SCF2;
    }
    else{
      cd = c; 
    }

    Vinf = cd;

    if (SCF<=SCF2){
      rc = PAO_Rcut + 4.0/(double)SCF;
    }
    else {
      rc = PAO_Rcut;
    }

    rc = PAO_Rcut;
    Vinf = c;

    r0 = rc - d;
     
    /*************************************************************
               f(r) = c0 + c1*r + c2*r^2 + c3*r^3 
    *************************************************************/

    AC[0][0] = 1.0;
    AC[0][1] = r0;
    AC[0][2] = r0*r0;
    AC[0][3] = r0*r0*r0;

    AC[1][0] = 0.0;
    AC[1][1] = 1.0;
    AC[1][2] = 2.0*r0;
    AC[1][3] = 3.0*r0*r0;

    AC[2][0] = 1.0;
    AC[2][1] = rc;
    AC[2][2] = rc*rc;
    AC[2][3] = rc*rc*rc;

    AC[3][0] = 0.0;
    AC[3][1] = 1.0;
    AC[3][2] = 2.0*rc;
    AC[3][3] = 3.0*rc*rc;

    Coes[0] = -Z/r0;
    Coes[1] = Z/r0/r0;
    Coes[2] = Vinf;
    Coes[3] = 0.0;

    RIH(3,AC,Coes);

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      if (MRV[i]<=r0)
        Vcore[i] = -Z/r;
      else if (MRV[i]<=rc)
        Vcore[i] = Coes[0] + Coes[1]*r + Coes[2]*r*r + Coes[3]*r*r*r;
      else
        Vcore[i] = Vinf;
    }
  }

}



void RIH(int n, double a[LSIZE1][LSIZE1], double b[LSIZE1])
{
  /****************************************************
                  LU decomposition
                      0 to n
   NOTE:
   This routine does not consider the reduction of rank
  ****************************************************/

  static int i,j,k,l,m;
  static double w,sum;
  static double x[LSIZE1],y[LSIZE1];
  static double da[LSIZE1][LSIZE1];
  static double ia[LSIZE1][LSIZE1];

  if (n==-1){
    for (i=0; i<=(LSIZE1-1); i++){
      for (j=0; j<=(LSIZE1-1); j++){
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

    for (i=0; i<=n; i++){
      sum = 0.0;
      for (j=0; j<=n; j++){
        sum = sum + ia[i][j]*b[j];
      }
      x[i] = sum;
    }

    for (i=0; i<=n; i++){
      b[i] = x[i];
    }        

    /*
    printf("A\n");
    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
        sum = 0.0;
        for (k=0; k<=n; k++){
          sum = sum + ia[i][k]*da[k][j];
        }
        printf("%10.7f ",sum);
      }
      printf("\n");
    }
    */

  }
}
