/*----------------------------------------------------------------------
  FEMHF_ERI.c
----------------------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "FEMHF_ERI.h"
#include "FEMHF_JKLM.h"


#define Sqrt(x) sqrtl(x)

#define PI  3.141592653589793238463L /* 22 digits */
#define LMAX_MAX 512
#define MAX_ITR 50
#define ERR_MAX 1e-12
#define CG_JMAX 30


static long double CG(int j1, int m1, int j2, int m2, int j, int m);
static void swap_sk(int *k1, int *s1, int *k2, int *s2);


long double FEMHF_ERI(int k1, int k2, int k3, int k4, 
                      int s1, int s2, int s3, int s4, int l)
{
  int s, ktype, kinf, ksup;
  long double j1, j2, int_l, int_m, int_k;

  assert( s1==0 || s1==1 );
  assert( s2==0 || s2==1 );
  assert( s3==0 || s3==1 );
  assert( s4==0 || s4==1 );

#if 0
  fprintf(stderr, "IN:\n");
  fprintf(stderr, "  K1= %3d K2= %3d K3= %3d K4= %3d ", k1, k2, k3, k4);
  fprintf(stderr, "S1= %1d S2= %1d S3= %1d S4= %1d ", s1, s2, s3, s4);
  fprintf(stderr, "L= %1d\n", l);
#endif

  /* normalize s's and k's */
  if (k1>k3) { swap_sk(&k1, &s1, &k3, &s3); }
  if (k2>k4) { swap_sk(&k2, &s2, &k4, &s4); }
  if (k1>k2) { 
    swap_sk(&k1, &s1, &k2, &s2); 
    swap_sk(&k3, &s3, &k4, &s4); 
  } else if (k1==k2 && k3>k4) {
    swap_sk(&k1, &s1, &k2, &s2); 
    swap_sk(&k3, &s3, &k4, &s4); 
  }

#if 0
  fprintf(stderr, "NORMALIZED:\n");
  fprintf(stderr, "  K1= %3d K2= %3d K3= %3d K4= %3d ", k1, k2, k3, k4);
  fprintf(stderr, "S1= %1d S2= %1d S3= %1d S4= %1d\n", s1, s2, s3, s4);
#endif

  /* single index for s's */
  s = s1*8 + s2*4 + s3*2 + s4;

  j1 = FEMHF_J(k1, s1, k3, s3, 2*l+5);
  j2 = 0.0L;
  if ((l<=1 && k4>0) || k4>1) { j2 = FEMHF_J(k2, s2, k4, s4, 3-2*l); } 
    
  if (isnan(j2) || isinf(j2)) {
    fprintf(stderr, "*** NAN FOUND!\n");
    fprintf(stderr, "  J2(K2= %3d, S2= %1d, K4= %3d, S4= %1d, L= %1d)\n",
      k2, s2, k4, s4, l);
    abort();
  }

  if (k1+2<=k4) { return j1*j2; } 

  /* overlap in integration range */
  int_l = 0.0L;
  int_m = 0.0L;

  if (k1==k3 && k4==k1+1 && k1>0) {
    int_l = FEMHF_L(k1, l, s1, s3);
    if (isnan(int_l)) {
      fprintf(stderr, "*** NAN FOUND!\n");
      fprintf(stderr, "  L(K1= %3d, S1= %1d, S3= %1d, L= %1d)\n",
        k1, s1,s3, l);
      abort();
    }
  }

  if (k1+1==k2 && k2==k4 && k4>0) {
    int_m = FEMHF_M(k4,l,s2,s4);
    if (isnan(int_m)) {
      fprintf(stderr, "*** NAN FOUND!\n");
      fprintf(stderr, "  M(K4= %3d, S2= %1d, S4= %1d, L= %1d)\n",
        k4, s2,s4, l);
      abort();
    }
  }

  int_k = 0.0L;
  if (l<=1 || k4>1) { int_k = FEMHF_K(k1, k2, k3, k4, s1, s2, s3, s4, l); } 

  if (isnan(int_k)) {
    fprintf(stderr, "*** NAN FOUND!\n");
    fprintf(stderr, "  K(K1= %3d, S1= %1d, K2= %3d, S2= %1d",
      k1, s1, k2, s2);
    fprintf(stderr, " , K3= %3d, S3= %1d, K4= %3d, S4= %1d",
      k3, s3, k4, s4);
    fprintf(stderr, ", L= %1d)\n", l);
    abort();  
  }

  return int_l*j2 + j1*int_m - int_l*int_m + int_k;
}


static void swap_sk(int *k1, int *s1, int *k2, int *s2)
{
  int tmp;
  tmp = *k1; *k1 = *k2; *k2 = tmp;
  tmp = *s1; *s1 = *s2; *s2 = tmp;
}




/*----------------------------------------------------------------------
  FEMHF_Gaunt

  Routine to calculate the Gaunt coefficients.
  
  See Eq. (3.7.73) in Modern Quantum Mechanics by J. J. Sakurai.
----------------------------------------------------------------------*/
long double FEMHF_Gaunt(int l,  int m, int l1, int m1, int l2, int m2)
{
  long double tmp0, tmp1, tmp2, tmp3;
  long double result, cg1, cg2;

  cg1 = CG(l1, 0,  l2,  0, l, 0);
  cg2 = CG(l1, m1, l2, m2, l, m);
  
  return cg1*cg2*Sqrt( (2*l1+1)*(2*l2+1)/(8.0L*l+4.0L)/PI );
}



/*----------------------------------------------------------------------
  Clebsch-Gordan coefficients

  <j1 j2; m1 m2 | j1 j2; j m>
  =
  delta(m,m1+m2)*sqrt(2j+1)
    *sqrt{ (j1+j2-j)! (j+j1-j2)! (j+j2-j1)! / (j+j1+j2+1)! } 
    *sqrt{ (j1+m1)! (j1-m1)! (j2+m2)! (j2-m2)! (j+m)! (j-m!) }
    *sum_k (-1)^k / ( k! (j1+j2-j-k)! (j1-m1-k)! (j2+m2-k)! 
                       (j-j2+m1+k)!  (j-j1-m2+k)! )

  See M. Avbramowitz and I.A. Stegun (eds.), "Handbook of Mathematical 
  Functions with Formulas, Graphs, and Tables," 
  Dover Publications, Inc. (Mineola, NY), 1972; an electronic copy of 
  this book is available at http://www.math.sfu.ca/~cbm/aands/.
----------------------------------------------------------------------*/
static long double CG(int j1, int m1, int j2, int m2, int j, int m)
{
  int kmin, kmax, k;
  long double sgn, cg;
  const int fmax = 3*CG_JMAX;
  long double f[3*CG_JMAX];

  /* this routine safely calculates when j1, j2, j < CG_JMAX */
  if (j1>CG_JMAX || j2>CG_JMAX || j>CG_JMAX) {
    fprintf(stderr, "*** out of bound (%s, %d)\n", __FILE__, __LINE__);
    abort();
  }

  /* factorials */
  f[0] = f[1] = 1.0L;
  for (k=2; k<fmax; k++) { f[k] = f[k-1]*k; }

  /* delta(m,m1+m2) */
  if (m != m1+m2) return 0.0L;

  /* conditions for j1, j2, and j */
  if (j1+j2 < j || j2+j < j1 || j+j1 < j2) return 0.0L;

  /* conditions for m1, m2 and m */
  if (abs(m1)>j1 || abs(m2)>j2 || abs(m)>j) return 0.0L;

  /* determin the range of k 
     max(0, -j+j2-m1, -j+j1+m2) <= k <= min(j1+j2-j, j1-m1, j2+m2) */
  kmin = 0;
  if (kmin < -j+j2-m1) kmin = -j+j2-m1;
  if (kmin < -j+j1+m2) kmin = -j+j1+m2;

  kmax = j1+j2-j;
  if (kmax > j1-m1) kmax = j1-m1;
  if (kmax > j2+m2) kmax = j2+m2;

  if (kmin>kmax) return 0.0L;

  cg = 0.0L;
  sgn = (kmin%2==0) ? 1.0L : -1.0L;
  for (k=kmin; k<=kmax; k++) {
    cg += sgn/(f[k]*f[j1+j2-j-k]*f[j1-m1-k]*f[j2+m2-k]
          *f[j-j2+m1+k]*f[j-j1-m2+k]);
    sgn = -sgn;
  }
  
  cg *= Sqrt( (2*j+1)*f[j1+j2-j]*f[j+j1-j2]*f[j+j2-j1]*f[j1+m1]*f[j1-m1]
                *f[j2+m2]*f[j2-m2]*f[j+m]*f[j-m]/f[j+j1+j2+1] );

  return cg;
}

