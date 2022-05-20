/**********************************************************************
  VP.c:

     VP.c is a subroutine to calculate the sum of potetials
     on radial grids. 

  Log of VP.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void VP()
{
  static int i;

  for (i=0; i<Grid_Num; i++){
    V[i] = Vcore[i] + Vh[i] + Vxc[i];
    V_all_ele[i] = V[i];

    /*
    printf("VP %4d %18.15f %18.15f %18.15f %18.15f\n",
           i,V[i],Vcore[i],Vh[i],Vxc[i]);
    */
  }

} 

