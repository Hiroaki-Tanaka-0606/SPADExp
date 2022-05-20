/**********************************************************************
  Hartree.c:

     Hartree.c is a subroutine to calculate Hartree potentials on 
     grids.

  Log of Hartree.c:

     10/Dec/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Hartree(double rh[ASIZE1], double h[ASIZE1])
{
  static int i,j,k;
  static double r,rp,Inside,Outside;
  static double dum;

  for (i=0; i<Grid_Num; i++){
    r = MRV[i];

    /****************************************************
                          Inside
    ****************************************************/

    Inside = 0.0;
    for (j=0; j<=i; j++){
      rp = MRV[j];
      Inside = Inside + rh[j]*rp*rp*rp*dx;
    }     
    Inside = 4.0*PI*Inside/r;

    /****************************************************
                          Outside
    ****************************************************/

    Outside = 0.0;
    for (j=(i+1); j<Grid_Num; j++){
      rp = MRV[j];
      Outside = Outside + rh[j]*rp*rp*dx;
    }     
    Outside = 4.0*PI*Outside;
    h[i] = Inside + Outside;
  }

  if (Check_switch==1){ 
    dum = 0.0;
    for (j=0; j<Grid_Num; j++){
      rp = MRV[j];
      dum = dum + 4.0*PI*rh[j]*rp*rp*rp*dx; 
    }     
    printf("<ALL>  total_electron=%15.12f  Calcd total_electron=%15.12f\n",
            total_electron,dum);
  }

}
