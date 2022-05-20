/**********************************************************************
  Density_V.c:

     Density_V.c is a subroutine to calculate densities of valence 
     electrons on grids.

  Log of Density_V.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Density_V(int state_num, double OcpN0[ASIZE15][2][ASIZE3+1][ASIZE3])
{
  static int i,j,k,l,n,m,so;
  static double sum,r,p,p2;

  /* calculate rho_V */

  for (i=0; i<Grid_Num; i++){

    r = MRV[i];
    sum = 0.0;
  
    for (so=0; so<SOI_switch; so++){
      for (m=0; m<Number_VPS; m++){
        n = NVPS[m];
        l = LVPS[m];
        p = W2[so][m][i];
        p2 = p*p;
        sum = sum + OcpN0[state_num][so][n][l]*p2/4.0/PI/r/r;
      }
    }

    rho_V[state_num][i] = sum;
  }

  /* normalization of rho_V */

  sum = 0.0;
  dx = MXV[1] - MXV[0];
  for (i=0; i<Grid_Num; i++){
    r = MRV[i];
    sum += rho_V[state_num][i]*r*r*r;
  }   
  sum *= (4.0*PI*dx); 

  for (i=0; i<Grid_Num; i++){
    rho_V[state_num][i] *= (valence_electron/sum);
  }

}


