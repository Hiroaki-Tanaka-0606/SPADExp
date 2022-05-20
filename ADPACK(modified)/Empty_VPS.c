/**********************************************************************
  Empty_VPS.c:

     Empty_VPS.c is a subroutine to set zero VPS potentials of an empty atom

  Log of Empty_VPS.c:

     18/Feb/2005  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Empty_VPS()
{
  static int so,m,i;

  for (so=0; so<SOI_switch; so++){
    for (m=0; m<Number_VPS; m++){
      for (i=0; i<Grid_Num; i++){
        VPS[0][so][m][i] = V[i];
        Vh_V[i] = 0.0;
        Vxc_V[i] = 0.0;
        rho_V[0][i] = 0.0;
      }
    }
  }

} 

