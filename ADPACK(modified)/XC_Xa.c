/**********************************************************************
  XC_Xa.c:

     XC_Xa.c is a subroutine to calculate the exchange-correlation
     potential within LDA by X alpha potential.

  Log of XC_Xa.c:

     10/Nov/2002  Released by T.Ozaki,

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void XC_Xa(double rh[ASIZE1], double xc[ASIZE1])
{

  /****************************************************
     Xalpha exchange-correlation potential by Slater
  ****************************************************/

  static int i,j;
  static double dum,alpha;

  alpha = 0.70;

  for (i=0; i<Grid_Num; i++){
    dum = 3.0/PI*rh[i];
    xc[i] = -3.0/2.0*alpha*pow(dum,0.33333333333333333);
  }
} 
