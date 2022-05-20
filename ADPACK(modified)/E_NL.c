/**********************************************************************
  E_NL.c:

     E_NL.c is a subroutine to calculate projection energies of
     the non-local part in the pseudo potential.

  Log of E_NL.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void E_NL(int NumVPS)
{
  static int i,n,j,l,nf,fg;
  static double r,rmin,rmax,Sr,Dr,sum,dum;

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  rmax = MRV[Grid_Num-1];
  Sr = rmax + rmin;
  Dr = rmax - rmin;

  for (l=0; l<=(NumVPS-1); l++){
    sum = 0.0;
    for (i=0; i<=(n-1); i++){
      r = 0.50*(Dr*x[i] + Sr);
      sum += 0.5*Dr*VNLF(l,0,r)*PAO_RadialF(l,r,W2)*PAO_RadialF(l,r,W2)*w[i];
    }

    if (l!=Local_Part_VPS) ProjectEnl[l] = 1.0/sum;
  }
}

