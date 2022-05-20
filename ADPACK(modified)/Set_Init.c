/**********************************************************************
  Set_Init.c:

     Set_Init.c is a subroutine to initialize several quantities.

  Log of Set_Init.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Set_Init()
{

  int i,j,n,l;

  dx = (Grid_Xmax - Grid_Xmin)/(double)(Grid_Num-1);

  for (i=0; i<Grid_Num; i++){
    MXV[i] = Grid_Xmin + dx*(double)i;
    MRV[i] = exp(MXV[i]);
  }

  for (i=0; i<Grid_Num; i++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
	PF[0][n][l][i] = 0.0;
	PF[1][n][l][i] = 0.0;
      }
    }
  }

  for (i=0; i<Grid_Num; i++){
    for (j=0; j<ASIZE9; j++){
      rho[j][i] = 0.0;
    }
  } 

}


