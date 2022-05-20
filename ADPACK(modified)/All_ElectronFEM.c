/**********************************************************************
  All_Electron.c:

   All_Electron.c is a subroutine to perform the self-consistent
   calculation of an atomic Kohn-Sham equation including all electrons.

  Log of All_Electron.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

#ifndef ___INTEGER_definition___
typedef int INTEGER; /* for fortran integer */
#define ___INTEGER_definition___ 
#endif




double find_Emin(long int n, long double *H0, long double *H1, long double *S0, long double *S1);
void diagonalize(INTEGER N0, long double **H, long double **S);
long double TH(long int n, long double x);
long double TH1(long int n, long double x);


void All_ElectronFEM()
{
  static long int i,j,l,N;
  static long double l2,d,q,fac,fac0,fac1;
  static long double tmp0,tmp1;
  static long double **H,**S;
 
  /* allocation of arrays */

  N = (long int)Grid_Num;
  
  H = (long double**)malloc(sizeof(long double*)*2*N);
  for (i=0; i<2*N; i++){
    H[i] = (long double*)malloc(sizeof(long double)*6);
    for (j=0; j<6; j++){
      H[i][j] = 0.0;
    }
  }

  S = (long double**)malloc(sizeof(long double*)*2*N);
  for (i=0; i<2*N; i++){
    S[i] = (long double*)malloc(sizeof(long double)*6);
    for (j=0; j<6; j++){
      S[i][j] = 0.0;
    }
  }

  /* step for the x-coordinate */

  d = (long double)Grid_Xmax/(long double)N;



  /*

  {
    long double *v0,*v1,*v2;
    long double d0,d5,xi,xi0,xj,c0,c1,c2;
    long int num,k;

    num = 40000;

    d0 = (long double)Grid_Xmax/(long double)num;
    d5 = (long double)Grid_Xmax/(long double)(num*5);

    v0 = (long double*)malloc(sizeof(long double)*5*num);
    v1 = (long double*)malloc(sizeof(long double)*5*num);
    v2 = (long double*)malloc(sizeof(long double)*5*num);

    for (i=0; i<5*num; i++){
      v0[i] = 0.0;
      v1[i] = 0.0;
      v2[i] = 0.0;
    }

    for (i=0; i<num; i++){
      xi0 = (long double)i*d0; 
      xi = (long double)(i-1)*d0; 

      for (j=0; j<5*2; j++){
        xj = xi + (long double)j*d5;  
        if (0.0<=xj) {

          c0 = 2.0*exp(-xi0*xi0);
          c1 = -4.0*xi0*exp(-xi0*xi0);
          c2 = -4.0*exp(-xi0*xi0) + 8.0*xi0*xi0*exp(-xi0*xi0);

          v0[5*(i-1)+j] += TH(0,(xj-(long double)i*d0)/d0)*c0 
                         + TH(1,(xj-(long double)i*d0)/d0)*c1 
                         + TH(2,(xj-(long double)i*d0)/d0)*c2;

          v1[5*(i-1)+j] += TH1(0,(xj-(long double)i*d0)/d0)*c0 
                         + TH1(1,(xj-(long double)i*d0)/d0)*c1 
                         + TH1(2,(xj-(long double)i*d0)/d0)*c2;

        }
      }
    }


    for (j=0; j<5*num; j++){
      xj = (long double)j*d5;
      c0 = 2.0*exp(-xj*xj);
      c1 = -4.0*xj*exp(-xj*xj);
      c2 = -4.0*exp(-xj*xj) + 8.0*xj*xj*exp(-xj*xj);
      printf("%18.15Lf %18.15Lf %18.15Lf %18.15Lf %18.15Lf %18.15Lf\n",xj,c0,c1,c2,v1[j],c0-v0[j]); 
    }


  }

  exit(0);
  */

  
  /**************************************************************
                          kinetic terms
  **************************************************************/
  
  /* diagonal element, i=0, for the kinetic operator */

  H[0][0] = d*d*15.0/280.0;
  H[0][1] = d*d*7.0/560.0;
  H[1][0] = d*d*7.0/560.0;
  H[1][1] = d*d*11.0/3360.0;

  /* diagonal element for the kinetic operator */

  for (i=1; i<N; i++){
    q = (long double)i; 
    H[i*2+0][2] = 3.0*d*d*q*(6.0 + 7.0*q*q)/35.0;
    H[i*2+0][3] = d*d*(1.0 + 6.0*q*q)/40.0;
    H[i*2+1][2] = d*d*(1.0 + 6.0*q*q)/40.0;
    H[i*2+1][3] = d*d*q*(3.0 + 7.0*q*q)/105.0;
  }
  
  /* off-diagonal element for the kinetic operator */

  for (i=0; i<(N-1); i++){

    q = (long double)i; 

    if (i==0) j = 2;
    else      j = 4;

    H[i*2+0][j+0] = -3.0*d*d*(5.0 + 24.0*q + 42.0*q*q + 28.0*q*q*q)/280.0;
    H[i*2+0][j+1] = d*d*(-5.0 - 12.0*q + 14.0*q*q*q)/560.0;
    H[i*2+1][j+0] = -d*d*(7.0 + 30.0*q + 42.0*q*q + 14.0*q*q*q)/560.0; 
    H[i*2+1][j+1] = -d*d*(11.0 + 36.0*q + 42.0*q*q + 28.0*q*q*q)/3360.0;

    H[(i+1)*2+0][0] = H[i*2+0][j+0];
    H[(i+1)*2+1][0] = H[i*2+0][j+1];
    H[(i+1)*2+0][1] = H[i*2+1][j+0];
    H[(i+1)*2+1][1] = H[i*2+1][j+1];
  }


  /*
   for (i=0; i<3*N; i++){
     for (j=0; j<9; j++){
       printf("%9.5Lf ",H[i][j]);
     }
     printf("\n");
   }
  */


  /*

  {
    long int k,l;
    long double sum,c0,c1,c2,c3,c4,c5,c6,c7,c8,x;
    long double *v0;

    v0 = (long double*)malloc(sizeof(long double)*3*N);

    sum = 0.0;
 
    j = 0;
    x = d*(long double)j;
    c0 = 2.0*exp(-x*x);
    c1 = -4.0*x*exp(-x*x);
    c2 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);
    j = 1;
    x = d*(long double)j;
    c3 = 2.0*exp(-x*x);
    c4 = -4.0*x*exp(-x*x);
    c5 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

    for (k=0; k<3; k++){  
     v0[k] = H[k][0]*c0 + H[k][1]*c1 + H[k][2]*c2 
           + H[k][3]*c3 + H[k][4]*c4 + H[k][5]*c5; 
    }

    for (j=1; j<(N-1); j++){

      x = d*(long double)(j-1);
      c0 = 2.0*exp(-x*x);
      c1 = -4.0*x*exp(-x*x);
      c2 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

      x = d*(long double)j;
      c3 = 2.0*exp(-x*x);
      c4 = -4.0*x*exp(-x*x);
      c5 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

      x = d*(long double)(j+1);
      c6 = 2.0*exp(-x*x);
      c7 = -4.0*x*exp(-x*x);
      c8 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

      for (k=0; k<3; k++){  
       v0[j*3+k] = H[j*3+k][0]*c0 + H[j*3+k][1]*c1 + H[j*3+k][2]*c2 
                 + H[j*3+k][3]*c3 + H[j*3+k][4]*c4 + H[j*3+k][5]*c5
                 + H[j*3+k][6]*c6 + H[j*3+k][7]*c7 + H[j*3+k][8]*c8;
      }
    }

    x = d*(long double)(N-2);
    c0 = 2.0*exp(-x*x);
    c1 = -4.0*x*exp(-x*x);
    c2 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

    x = d*(long double)(N-1);
    c3 = 2.0*exp(-x*x);
    c4 = -4.0*x*exp(-x*x);
    c5 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);

    j = N-1;
    for (k=0; k<3; k++){  
     v0[j*3+k] = H[j*3+k][0]*c0 + H[j*3+k][1]*c1 + H[j*3+k][2]*c2 
               + H[j*3+k][3]*c3 + H[j*3+k][4]*c4 + H[j*3+k][5]*c5;
    }

    sum = 0.0;
    for (j=0; j<N; j++){

      x = d*(long double)j;
      c0 = 2.0*exp(-x*x);
      c1 = -4.0*x*exp(-x*x);
      c2 = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);
      sum += c0*v0[j*3+0] + c1*v0[j*3+1] + c2*v0[j*3+2];
    }

    printf("kin  %18.15Lf\n",sum); 

    free(v0);
  }
  */



  /**************************************************************
                           l*(l+1)/(2*x^4) 
  **************************************************************/

  l = 0;

  l2 = (long double)l*((long double)l+1.0);

  /* diagonal element, i=0, for l*(l+1)/(2*x^4) */

  H[0][0] += d*d*3.0*l2/35.0;
  H[0][1] += d*d*7.0*l2/420.0;
  H[1][0] += d*d*7.0*l2/420.0;
  H[1][1] += d*d*3.0*l2/840.0;

  /* diagonal element for for l*(l+1)/(2*x^4) */

  for (i=1; i<N; i++){
    q = (long double)i;  
    H[i*2+0][2] += 26.0*d*d*q*l2/35.0;
    H[i*2+0][3] += d*d*l2/30.0;
    H[i*2+1][2] += d*d*l2/30.0;
    H[i*2+1][3] += 2.0*d*d*q*l2/105.0;
  }

  /* off-diagonal element for l*(l+1)/(2*x^4) */

  for (i=0; i<(N-1); i++){

    q = (long double)i; 

    if (i==0) j = 2;
    else      j = 4;

    H[i*2+0][j+0] += 9.0*d*d*(1.0 + 2.0*q)*l2/140.0;
    H[i*2+0][j+1] += -d*d*(6.0 + 13.0*q)*l2/420.0;
    H[i*2+1][j+0] += d*d*(7.0 + 13.0*q)*l2/420.0;
    H[i*2+1][j+1] += -d*d*(1.0 + 2.0*q)*l2/280.0;

    H[(i+1)*2+0][0] = H[i*2+0][j+0];
    H[(i+1)*2+1][0] = H[i*2+0][j+1];
    H[(i+1)*2+0][1] = H[i*2+1][j+0];
    H[(i+1)*2+1][1] = H[i*2+1][j+1];

  }

  /**************************************************************
                            -Z/(x^2)
  **************************************************************/

  AtomNum = 1;

  /* diagonal element, i=0, for -Z/(x^2) */

  H[0][0] +=-d*d*d*d*11.0/420.0*(long double)AtomNum;
  H[0][1] +=-d*d*d*d*8.0/1260.0*(long double)AtomNum;
  H[1][0] +=-d*d*d*d*8.0/1260.0*(long double)AtomNum; 
  H[1][1] +=-d*d*d*d*2.0/1260.0*(long double)AtomNum; 

  /* diagonal element for -Z/(x^2) */

  for (i=1; i<N; i++){
    q = (long double)i; 
    H[i*2+0][2] += -2.0*d*d*d*d*q*(19.0 + 78.0*q*q)*(long double)AtomNum/105.0;
    H[i*2+0][3] += -d*d*d*d*(4.0 + 63.0*q*q)*(long double)AtomNum/315.0;
    H[i*2+1][2] += -d*d*d*d*(4.0 + 63.0*q*q)*(long double)AtomNum/315.0;
    H[i*2+1][3] += -2.0*d*d*d*d*q*(1.0 + 2.0*q*q)*(long double)AtomNum/105.0;
  }

  /* off-diagonal element for -Z/(x^2) */

  for (i=0; i<(N-1); i++){

    q = (long double)i; 

    if (i==0) j = 2;
    else      j = 4;

    H[i*2+0][j+0] += -d*d*d*d*(19.0 + 92.0*q + 162.0*q*q + 108.0*q*q*q)*(long double)AtomNum/420.0;
    H[i*2+0][j+1] += d*d*d*d*(11.0 + 57.0*q + 108.0*q*q + 78.0*q*q*q)*(long double)AtomNum/1260.0;
    H[i*2+1][j+0] += -d*d*d*d*(16.0 + 75.0*q + 126.0*q*q + 78.0*q*q*q)*(long double)AtomNum/1260.0;
    H[i*2+1][j+1] += d*d*d*d*(1.0 + 5.0*q + 9.0*q*q + 6.0*q*q*q)*(long double)AtomNum/420.0;

    H[(i+1)*2+0][0] = H[i*2+0][j+0];
    H[(i+1)*2+1][0] = H[i*2+0][j+1];
    H[(i+1)*2+0][1] = H[i*2+1][j+0];
    H[(i+1)*2+1][1] = H[i*2+1][j+1];

  }


  /*
  {
    long int k,l;
    long double sum,c0,c1,c2,c3,c4,c5,c6,c7,c8,x;
    long double *v0,*v1;

    v0 = (long double*)malloc(sizeof(long double)*3*N);
    v1 = (long double*)malloc(sizeof(long double)*3*N);


    for (j=0; j<N; j++){
      x = d*(long double)j;
      v0[j*3+0] = 2.0*exp(-x*x);
      v0[j*3+1] = -4.0*x*exp(-x*x);
      v0[j*3+2] = -4.0*exp(-x*x) + 8.0*x*x*exp(-x*x);
    }

    for (k=0; k<3; k++){  
     v1[k] = H[k][0]*v0[0] + H[k][1]*v0[1] + H[k][2]*v0[2]
           + H[k][3]*v0[3] + H[k][4]*v0[4] + H[k][5]*v0[5]; 
    }

    for (j=1; j<(N-1); j++){

      for (k=0; k<3; k++){  
       v1[j*3+k] = H[j*3+k][0]*v0[j*3-3] + H[j*3+k][1]*v0[j*3-2] + H[j*3+k][2]*v0[j*3-1] 
                 + H[j*3+k][3]*v0[j*3+0] + H[j*3+k][4]*v0[j*3+1] + H[j*3+k][5]*v0[j*3+2]
                 + H[j*3+k][6]*v0[j*3+3] + H[j*3+k][7]*v0[j*3+4] + H[j*3+k][8]*v0[j*3+5];
      }
    }

    j = N-1;
    for (k=0; k<3; k++){  
     v1[j*3+k] = H[j*3+k][0]*v0[j*3-3] + H[j*3+k][1]*v0[j*3-2] + H[j*3+k][2]*v0[j*3-1] 
               + H[j*3+k][3]*v0[j*3+0] + H[j*3+k][4]*v0[j*3+1] + H[j*3+k][5]*v0[j*3+2];
    }

    sum = 0.0;
    for (j=0; j<N; j++){
      sum += v1[j*3+0]*v0[j*3+0] + v1[j*3+1]*v0[j*3+1] + v1[j*3+2]*v0[j*3+2];
    }

    printf("-Z/x^2  %18.15Lf\n",sum); 

    free(v0);
    free(v1);

  }

  exit(0);
  */



  /*
  {
    long int j;
    long double sum,ci,cj,x;
    sum = 0.0;
 
    for (i=0; i<N; i++){
      x = d*(long double)i;
      ci = 2.0*exp(-x*x);
      for (j=0; j<N; j++){
        x = d*(long double)j;
        cj = 2.0*exp(-x*x);

        if (i==j)      sum += H0[i  ]*ci*cj;
        if (j==(i+1))  sum += H1[i  ]*ci*cj;
        if (j==(i-1))  sum += H1[i-1]*ci*cj;
      }
    } 

    printf("-Z/x^2  %18.15Lf\n",sum); 
  }
  */

  /**************************************************************
                         overlap integral
  **************************************************************/

  /* diagonal element, i=0, for overlap integral */

  S[0][0] = d*d*d*d*d*d*49.0/6930.0;
  S[0][1] = d*d*d*d*d*d*13.0/6930.0;
  S[1][0] = d*d*d*d*d*d*13.0/6930.0;
  S[1][1] = d*d*d*d*d*d*7.0/13860.0;

  /* diagonal element for overlap integral */

  for (i=1; i<N; i++){
    q = (long double)i; 
    S[i*2+0][2] = (2.0*d*d*d*d*d*d*q*(225.0 + 2090.0*q*q + 2574.0*q*q*q*q))/3465.0;
    S[i*2+0][3] = (d*d*d*d*d*d*(13.0 + 440.0*q*q + 1155.0*q*q*q*q))/3465.0;
    S[i*2+1][2] = (d*d*d*d*d*d*(13.0 + 440.0*q*q + 1155.0*q*q*q*q))/3465.0;
    S[i*2+1][3] = (2.0*d*d*d*d*d*d*q*(15.0 + 110.0*q*q + 66.0*q*q*q*q))/3465.0;
  }

  /* off-diagonal element for overlap integral */

  for (i=0; i<(N-1); i++){

    q = (long double)i; 

    if (i==0) j = 2;
    else      j = 4;

    S[i*2+0][j+0] = (d*d*d*d*d*d*(1.0 + 2.0*q)*(287.0 + 22.0*q*(1.0 + q)*(68.0 + 81.0*q*(1.0 + q))))/13860.0;
    S[i*2+0][j+1] = -(d*d*d*d*d*d*(49.0 + 375.0*q + 1210.0*q*q + 2090.0*q*q*q + 1980.0*q*q*q*q + 858.0*q*q*q*q*q))/13860.0;
    S[i*2+1][j+0] = (d*d*d*d*d*d*(84.0 + 595.0*q + 1760.0*q*q + 2750.0*q*q*q + 2310.0*q*q*q*q + 858.0*q*q*q*q*q))/13860.0;
    S[i*2+1][j+1] = -(d*d*d*d*d*d*(1.0 + 2.0*q)*(14.0 + 11.0*q*(1.0 + q)*(7.0 + 9.0*q*(1.0 + q))))/13860.0;

    S[(i+1)*2+0][0] = S[i*2+0][j+0];
    S[(i+1)*2+1][0] = S[i*2+0][j+1];
    S[(i+1)*2+0][1] = S[i*2+1][j+0];
    S[(i+1)*2+1][1] = S[i*2+1][j+1];

  }

  /*
  {
    long int j;
    long double sum,ci,cj,x;
    sum = 0.0;
 
    for (i=0; i<N; i++){
      x = d*(long double)i;
      ci = 2.0*exp(-x*x);
      for (j=0; j<N; j++){
        x = d*(long double)j;
        cj = 2.0*exp(-x*x);

        if (i==j)      sum += S0[i  ]*ci*cj;
        if (j==(i+1))  sum += S1[i  ]*ci*cj;
        if (j==(i-1))  sum += S1[i-1]*ci*cj;
      }
    } 

    printf("overlap  %18.15Lf\n",sum); 
  }
  */


  /*
  for (i=0; i<N; i++){
    printf("%6d %18.14Lf %18.14Lf\n",i,H0[i],H1[i]);
  }


  for (i=0; i<N; i++){
    printf("%6d %25.20Lf %25.20Lf\n",i,S0[i],S1[i]);
  }
  exit(0);
  */

  /*
  H0[0] = 1.0;
  H0[1] =-1.0;
  H0[2] =-0.5;
  H1[0] = 1.0;
  H1[1] = 0.5;

  S0[0] = 1.0;
  S0[1] = 1.0;
  S0[2] = 1.0;
  S1[0] = 0.1;
  S1[1] = 0.2;

  diagonalize(3,H0,H1,S0,S1);  
  */

  /*
  printf("%15.12f\n",find_Emin(3,H0,H1,S0,S1));
  */



  /*
  for (i=0; i<N; i++){
    H0[i] = H0[i] - 300.0*S0[i];  
  }

  for (i=0; i<(N-1); i++){
    H1[i] = H1[i] - 300.0*S1[i];  
  }
  */


  diagonalize(N,H,S);  


  /*
  printf("%15.12f\n",find_Emin(N,H0,H1,S0,S1));
  */

  exit(0);


  /* freeing of arrays */

  for (i=0; i<3*N; i++){
    free(H[i]);
  }
  free(H);

  for (i=0; i<3*N; i++){
    free(S[i]);
  }
  free(S);

}



double find_Emin(long int n, long double *H0, long double *H1, long double *S0, long double *S1)
{
  static long int i,j,k,k1,k2,po,num;
  static long double *v1,*v2;
  static long double tmp,e0,e1;
  /* !!! change below for accuracy and efficiency !!! */
  static long double diff=1.0e-10;
  static long int nummax=20000; 

  /* allocation of arrays */

  v1 = (long double*)malloc(sizeof(long double)*(n+3)); 
  v2 = (long double*)malloc(sizeof(long double)*(n+3)); 

  /* initial vector */

  for (i=0; i<n; i++) v1[i] = (long double)i+1.0;

  i = 0;
  v2[i] = S0[i]*v1[i] + S1[i]*v1[i+1];        
  for (i=1; i<n; i++){
    v2[i] = S1[i-1]*v1[i-1]*S0[i]*v1[i] + S1[i]*v1[i+1];        
  }

  tmp = 0.0;
  for (i=0; i<n; i++)  tmp += v1[i]*v2[i];

  tmp = 1.0/sqrt(fabs(tmp));
  for (i=0; i<n; i++)  v1[i] = v1[i]*tmp;

  /* steepest decent method */

  e0 = 1.0e+100;
  e1 = 0.0;
  po = 0;
  num = 0;

  do {

    num++;

    i = 0;
    v2[i] = H0[i]*v1[i] + H1[i]*v1[i+1] - e1*(S0[i]*v1[i] + S1[i]*v1[i+1]);
    for (i=1; i<n; i++){
      v2[i] = H1[i-1]*v1[i-1] + H0[i]*v1[i] + H1[i]*v1[i+1] 
         -e1*(S1[i-1]*v1[i-1] + S0[i]*v1[i] + S1[i]*v1[i+1]);
    }

    for (i=0; i<n; i++) v1[i] -= 0.01*v2[i]; 

    i = 0;
    v2[i] = S0[i]*v1[i] + S1[i]*v1[i+1];
    for (i=1; i<n; i++){
      v2[i] = S1[i-1]*v1[i-1] + S0[i]*v1[i] + S1[i]*v1[i+1];
    }

    tmp = 0.0;
    for (i=0; i<n; i++) tmp += v1[i]*v2[i];   

    tmp = 1.0/sqrt(fabs(tmp));

    for (i=0; i<n; i++) v1[i] = v1[i]*tmp;

    i = 0;     
    v2[i] = H0[i]*v1[i] + H1[i]*v1[i+1];
    for (i=1; i<n; i++){
      v2[i] = H1[i-1]*v1[i-1] + H0[i]*v1[i] + H1[i]*v1[i+1];
    }

    tmp = 0.0;
    for (i=0; i<n; i++)  tmp += v1[i]*v2[i];   
    e1 = tmp;

    /* converge? */

    if (fabs(e1-e0)<diff) po = 1;
    else                  e0 = e1;


    printf("num=%5d  e1=%18.15Lf\n",num,e1);

  } while(po==0 && num<nummax);

  /* freeing of arrays */

  free(v1);
  free(v2);

  /* return */
  return e1;
}





void diagonalize(INTEGER N0, long double **H, long double **S)
{
  int i,j,i1,j1,ii,jj;

  char  *JOBZ="V";
  char  *RANGE="A";
  char  *UPLO="L";

  INTEGER ITYPE;
  double VL,VU; /* dummy */
  INTEGER IL,IU;  /* dummy */
  double ABSTOL=1.0e-14;
  double *Z;
  double *W;
  double *A;
  double *B;
  double *WORK;
  INTEGER N;
  INTEGER LDZ;
  INTEGER LDA;
  INTEGER LDB;
  INTEGER M;
  INTEGER LWORK;
  INTEGER *IWORK;
  INTEGER *IFAIL;
  INTEGER INFO;

  ITYPE = 1;

  N = 2*N0; 

  LDA = N;
  LDB = N;
  LDZ = N;
  LWORK = 8*N;

  A = (double*)malloc(sizeof(double)*(N+4)*(N+4));
  B = (double*)malloc(sizeof(double)*(N+4)*(N+4));
  Z = (double*)malloc(sizeof(double)*LDZ*N);
  W = (double*)malloc(sizeof(double)*N);
  WORK = (double*)malloc(sizeof(double)*LWORK);
  IWORK = (INTEGER*)malloc(sizeof(INTEGER)*5*N);
  IFAIL = (INTEGER*)malloc(sizeof(INTEGER)*N);

  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      A[i*N+j] = 0.0;
      B[i*N+j] = 0.0;
    }
  }

  i = 0;
  for (i1=0; i1<2; i1++) {
    ii = i*2 + i1;
    for (j=0; j<4; j++) {
      jj = j;
      A[ii*N+jj] = (double)H[ii][j];
      B[ii*N+jj] = (double)S[ii][j];
    }
  }

  for (i=1; i<(N0-1); i++) {
    for (i1=0; i1<2; i1++) {
      ii = i*2 + i1;
      for (j=0; j<6; j++) {
        jj = (i-1)*2 + j;
        A[ii*N+jj] = (double)H[ii][j];
        B[ii*N+jj] = (double)S[ii][j];
      }
    }
  }

  i = N0-1;
  for (i1=0; i1<2; i1++) {
    ii = i*2 + i1;
    for (j=0; j<4; j++) {
      jj = (i-1)*2 + j;
      A[ii*N+jj] = (double)H[ii][j];
      B[ii*N+jj] = (double)S[ii][j];
    }
  }

  dsygvx_( &ITYPE, JOBZ, RANGE, UPLO, &N, A, &LDA, B, &LDB, &VL, &VU, &IL, &IU, &ABSTOL,
           &M, W, Z, &LDZ, WORK, &LWORK, IWORK, IFAIL, &INFO );

  /* store eigenvectors */

  /*
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      ev[i+1][j+1]= Z[i*N+j];
    }
  }
  */

  /* shift ko by 1 */
  for (i=0; i<N; i++) {
    printf("%4d %20.16f\n",i,W[i]);
  }

  if (INFO>0) {
    /*
    printf("\n error in dstevx_, info=%d\n\n",INFO);fflush(stdout);
    */
  }
  if (INFO<0) {
    printf("info=%d in dstevx_\n",INFO);fflush(stdout);
    exit(0);
  }

  free(A);
  free(B);
  free(Z);
  free(W);
  free(WORK);
  free(IWORK);
  free(IFAIL);
}



long double TH(long int n, long double x)
{

  long double t;

  if (n==0){
  
    if (0.0<=x && x<=1.0)
      t = 1.0 - 10.0*x*x*x + 15.0*x*x*x*x - 6.0*x*x*x*x*x;  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = 1.0 - 10.0*x*x*x + 15.0*x*x*x*x - 6.0*x*x*x*x*x;  
    }
    else 
      t = 0.0;
  }

  else if (n==1){
    if (0.0<=x && x<=1.0)
      t = x - 6.0*x*x*x + 8.0*x*x*x*x - 3.0*x*x*x*x*x;  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = -(x - 6.0*x*x*x + 8.0*x*x*x*x - 3.0*x*x*x*x*x);
    }
    else 
      t = 0.0;
  }

  else if (n==2){
    if (0.0<=x && x<=1.0)
      t = 0.5*(x*x - 3.0*x*x*x + 3.0*x*x*x*x - x*x*x*x*x);  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = 0.5*(x*x - 3.0*x*x*x + 3.0*x*x*x*x - x*x*x*x*x);
    }
    else 
      t = 0.0;
  }

  return t;
}



long double TH1(long int n, long double x)
{

  long double t;

  if (n==0){
  
    if (0.0<=x && x<=1.0)
      t = -30.0*x*x + 60.0*x*x*x - 30.0*x*x*x*x;  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = -30.0*x*x + 60.0*x*x*x - 30.0*x*x*x*x;  
    }
    else 
      t = 0.0;
  }

  else if (n==1){
    if (0.0<=x && x<=1.0)
      t = 1 - 18.0*x*x + 32.0*x*x*x - 15.0*x*x*x*x;  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = -(1 - 18.0*x*x + 32.0*x*x*x - 15.0*x*x*x*x);
    }
    else 
      t = 0.0;
  }

  else if (n==2){
    if (0.0<=x && x<=1.0)
      t = 0.5*(2*x - 9.0*x*x + 12.0*x*x*x - 5.0*x*x*x*x);  
    else if (-1.0<=x && x<0.0){
      x = -x;
      t = 0.5*(2*x - 9.0*x*x + 12.0*x*x*x - 5.0*x*x*x*x);
    }
    else 
      t = 0.0;
  }

  return t;
}
