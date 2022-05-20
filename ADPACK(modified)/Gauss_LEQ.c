/**********************************************************************
  Gauss_LEQ.c:

     Gauss_LEQ.c is a subroutine to solve a linear equation
     using Gauss's method.

  Log of Gauss_LEQ.c:

     22/Nov/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Gauss_LEQ(int n, double a[ASIZE6][ASIZE6], double x[ASIZE6])
{
  /****************************************************
                   From 0 to n,  Ax = b
              The (n+1) column of a[][] is b.
  ****************************************************/

  static int i,j,k,max_i,po;
  static double max,dum1,dum2,w;

  for (i=0; i<=n; i++){

    /****************************************************
         choose the maximum element of the subspace.
    ****************************************************/

    po = 0;
    max = fabs(a[i][i]);
    for (j=i+1; j<=n; j++){
      if (max<fabs(a[j][i])){
	po = 1;
	max = fabs(a[j][i]);
	max_i = j;
      }         
    }  
  
    if (po==1){
      for (j=i; j<=(n+1); j++){
	dum1 = a[i][j];
	dum2 = a[max_i][j];
	a[i][j]     = dum2;
	a[max_i][j] = dum1;
      }          
    }

    /****************************************************
                       Gauss's method 
    ****************************************************/

    w = 1.0/a[i][i];
    for (j=i; j<=(n+1); j++){
      a[i][j] = a[i][j]*w;
    }

    for (j=(i+1); j<=n; j++){
      w = a[j][i];
      for (k=i; k<=(n+1); k++){
	a[j][k] = a[j][k] - a[i][k]*w;
      } 
    }
  }

  /****************************************************
                       Inverting
  ****************************************************/
 
  x[n] = a[n][n+1];
  for (i=(n-1); 0<=i; i--){
    dum1 = a[i][n+1];
    for (j=n; (i+1)<=j; j--){
      dum1 = dum1 - a[i][j]*x[j];
    }
    x[i] = dum1;
  }

}

