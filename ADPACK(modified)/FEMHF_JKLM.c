#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "FEMHF_JKLM.h"

#define Power(a,b) powl(a,b)

#define JKLM_MEMORY 1

static long double K1(int l, int s);
static long double K2(int ik, int l, int s);
static long double K3(int ik, int l, int s);
static long double K4(int ik, int l, int s);
static long double K5(int ik, int l, int s);
static long double K6(int ik, int l, int s);

static long double J1(int l, int s);
static long double J2(int ik, int l, int s);
static long double J3(int ik, int l, int s);

static long double Log10_p(int ik);
static long double Log20_p(int ik);



#if JKLM_MEMORY
#define LOG10( i ) Log10_saved[N_saved+(i)]
#define LOG20( i ) Log20_saved[N_saved+(i)]

static int N_saved = 0;
static long double *Log10_saved = NULL;
static long double *Log20_saved = NULL;
static long double **J1_saved = NULL;
static long double ***J2_saved = NULL;
static long double ***J3_saved = NULL;

#else
#define LOG10( i ) Log10_p(i)
#define LOG20( i ) Log20_p(i)
#endif


void FEMHF_JKLM_Init(int N)
{
#if JKLM_MEMORY
  int i, j, k, kmin;
  /* Logarithms */
  Log10_saved = (long double*)malloc(sizeof(long double)*2*N);
  Log10_saved[N] = 0.0L;
  Log10_saved[N-1] = 0.0L;
  for (i=1; i<N; i++) { Log10_saved[N+i] = Log10_p(i); }
  for (i=2; i<N; i++) { Log10_saved[N-i] = Log10_p(-i); }
  Log20_saved = (long double*)malloc(sizeof(long double)*2*N);
  Log20_saved[N] = 0.0L;
  Log20_saved[N-1] = 0.0L;
  for (i=1; i<N; i++) { Log20_saved[N+i] = Log20_p(i); }
  for (i=2; i<N; i++) { Log20_saved[N-i] = Log20_p(-i);  }
  N_saved = N;

  /* J integrals */
  J1_saved = (long double**)malloc(sizeof(long double*)*3);
  for (i=0; i<3; i++) {
    J1_saved[i] = (long double*)malloc(sizeof(long double)*16);
    for (j=0; j<16; j++) { J1_saved[i][j] = 0.0L; }
  }
  J2_saved = (long double***)malloc(sizeof(long double**)*3);
  for (i=0; i<3; i++) {
    J2_saved[i] = (long double**)malloc(sizeof(long double*)*16);
    for (j=0; j<16; j++) {
      J2_saved[i][j] = (long double*)malloc(sizeof(long double)*N);
      for (k=0; k<N; k++) { J2_saved[i][j][k] = 0.0L; }
    }
  }
  J3_saved = (long double***)malloc(sizeof(long double**)*4);
  for (i=0; i<4; i++) {
    J3_saved[i] = (long double**)malloc(sizeof(long double*)*16);
    for (j=0; j<16; j++) {
      J3_saved[i][j] = (long double*)malloc(sizeof(long double)*N);
      for (k=0; k<N; k++) { J3_saved[i][j][k] = 0.0L; }
    }
  }

  for (j=0; j<16; j++) {
    J1_saved[0][j] = J1(2*j-11, 0);
    J1_saved[1][j] = J1(2*j-11, 1);
    J1_saved[2][j] = J1(2*j-11, 3);
    kmin = (j==7 || j==8) ? 1 : 2;
    for (k=kmin; k<N; k++) {
      J2_saved[0][j][k] = J2(k, 2*j-11, 0);
      J2_saved[1][j][k] = J2(k, 2*j-11, 1);
      J2_saved[2][j][k] = J2(k, 2*j-11, 3);
    }
    kmin = (j==7 || j==8) ? 0 : 1;
    for (k=kmin; k<N; k++) {
      J3_saved[0][j][k] = J3(k, 2*j-11, 0);
      J3_saved[1][j][k] = J3(k, 2*j-11, 1);
      J3_saved[2][j][k] = J3(k, 2*j-11, 2);
      J3_saved[3][j][k] = J3(k, 2*j-11, 3);
    }
  } 
#endif 
}


void FEMHF_JKLM_Free(void)
{
#if JKLM_MEMORY
  int i, j;

  free(Log10_saved); 
  free(Log20_saved); 

  for (i=0; i<3; i++) { free(J1_saved[i]); }
  free(J1_saved);
  for (i=0; i<3; i++) {
    for (j=0; j<16; j++) { free(J2_saved[i][j]); }
    free(J2_saved[i]);
  }
  free(J2_saved);
  for (i=0; i<4; i++) {
    for (j=0; j<16; j++) { free(J3_saved[i][j]); }
    free(J3_saved[i]);
  }
  free(J3_saved); 
#endif
}


long double FEMHF_J(int k1, int s1, int k2, int s2, int l)
{
  int  ktype;
  const int s = s1*2 + s2;
#if JKLM_MEMORY
  const int il = (l+11)/2;
#endif

  /* k-type */
  ktype = 0;
  if (k1==k2) {
    if (k1==0) { ktype = 1; } else { ktype = 2; }
  } else if (k1+1==k2) {
    ktype = 3; 
  }

#if JKLM_MEMORY
  switch (ktype) {
  case 1:
    switch (s) {
    case 0:         return J1_saved[0][il];
    case 1: case 2: return J1_saved[1][il];
    case 3:         return J1_saved[2][il];
    }
    break;
  case 2: 
    switch (s) {
    case 0:         return J2_saved[0][il][k1];
    case 1: case 2: return J2_saved[1][il][k1];
    case 3:         return J2_saved[2][il][k1];
    }
    break;
  case 3:
    return J3_saved[s][il][k1];
  }
#else
  switch (ktype) {
  case 1: return J1(l, s);
  case 2: return J2(k1, l, s);
  case 3: return J3(k1, l, s);
  }
#endif
 
  fprintf(stderr, "***ERROR in J\n");
  fprintf(stderr, "  k1= %4d  s1= %1d k2= %4d s2= %1d l= %1d\n",
    k1, s1, k2, s2, l);
  abort();

  return 0.0L;
}



long double FEMHF_K(int k1, int k2, int k3, int k4, 
                    int s1, int s2, int s3, int s4, int l)
{
  int ktype;
  const int s = s1*8 + s2*4 + s3*2 + s4;

  /* normalize s's and k's */
  if (k1>k3 || k2>k4 || k1>k2) {
    fprintf(stderr, "***ERROR in K\n");
    abort();
  }

  /* k-type */
  ktype = 0;
  if (k1==k4) {
    if (k1==0) { ktype = 1; } else { ktype = 2; }
  } else {
    if (k2==k1 && k3==k4) { ktype=3; }
    if (k2==k1 && k3==k1) { ktype=4; }
    if (k2==k4 && k3==k4) { ktype=5; }
    if (k2==k4 && k3==k1) { ktype=6; }
  }

  switch (ktype) {
  case 1: return K1(l, s);
  case 2: return K2(k1, l, s);
  case 3: return K3(k1, l, s);
  case 4: return K4(k1, l, s);
  case 5: return K5(k1, l, s);
  case 6: return K6(k1, l, s);
  }

  fprintf(stderr, "***ERROR in K\n");
  abort();

  return 0.0L;
}


long double FEMHF_L(int ik, int l, int s1, int s2)
{
  const int s = (s1==0) ? (s2==0 ? 0 : 1) : (s2==0 ? 1 : 2);
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;

  switch (l) {
  case 0:
    switch (s) {
    case 0: return -(49-(450-(1815-(4180-(5940-5148*k)*k)*k)*k)*k)
      /13860.0L;
    case 1: return (13-(115-(440-(935-(1155-726*k)*k)*k)*k)*k)/13860.0L;
    case 2: return -(7-(60-(220-(440-(495-264*k)*k)*k)*k)*k)/27720.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0: return -(153-(1736-(8918-(27300-(55055-(76076-(72072-44616*k
      )*k)*k)*k)*k)*k)*k)/120120.0L;
    case 1: return (128-(1421-(7098-(20930-(40040-(51051-(42042-18876*k)
      *k)*k)*k)*k)*k)*k)/360360.0L;
    case 2: return -(36-(392-(1911-(5460-(10010-(12012-(9009-3432*k)*k)
      *k)*k)*k)*k)*k)/360360.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0: return -(11-(148-(918-(3472-(8918-(16380-(22022-(21736
      -(15444-7436*k)*k)*k)*k)*k)*k)*k)*k)*k)/20020.0L;
    case 1: return (19-(252-(1536-(5684-(14196-(25116-(32032-(29172
      -(18018-6292*k)*k)*k)*k)*k)*k)*k)*k)*k)/120120.0L;
    case 2: return -(11-(144-(864-(3136-(7644-(13104-(16016-(13728
      -(7722-2288*k)*k)*k)*k)*k)*k)*k)*k)*k)/240240.0L;
    }
    break;

  case 3:
    switch (s) {
    case 0: return -(299-(4644-(33660-(150960-(468180-(1062432-(1819272
      -(2386800-(2406690-(1847560-(1050192-413712*k)*k)*k)*k)*k)*k)*k)
      *k)*k)*k)*k)/1.11384e6L;
    case 1: return (88-(1353-(9690-(42840-(130560-(289884-(482664-(
      609960-(583440-(413270-(204204-58344*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      *k)/1.11384e6L;
    case 2: return -(26-(396-(2805-(12240-(36720-(79968-(129948-(159120
      -(145860-(97240-(43758-10608*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      /1.11384e6L;
    }
    break;

  case 4:
    switch (s) {
    case 0: return -(39-(686-(5681-(29412-(106590-(286824-(593028-(
      961248-(1234506-(1259700-(1016158-(638248-(302328-100776*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)/271320.0L;
    case 1: return (35-(611-(5016-(25707-(92055-(244188-(496128-(
      786828-(982566-(965770-(739024-(428298-(176358-42636*k)*k)*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)/813960.0L;
    case 2: return -(21-(364-(2964-(15048-(53295-(139536-(279072-(
      434112-(529074-(503880-(369512-(201552-(75582-15504*k)*k)*k)*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)/1.62792e6L;
    }
    break;

  case 5:
    switch (s) {
    case 0: return -(493-(9680-(90090-(528220-(2187185-(6794172-(
      16414860-(31550640-(48924810-(61680080-(63371308-(52907400-(
      35565530-(18901960-(7674480-2217072*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      *k)*k)*k)*k)*k)/5.96904e6L;
    case 1: return (448-(8745-(80850-(470470-(1931160-(5938317-(
      14176470-(26860680-(40930560-(50488130-(50438388-(40562340-(
      25865840-(12684210-(4476780-937992*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      *k)*k)*k)*k)*k)/1.790712e7L;
    case 2: return -(136-(2640-(24255-(140140-(570570-(1738044-(4103715
      -(7674480-(11511720-(13927760-(13579566-(10581480-(6466460-(
      2984520-(959310-170544*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      *k)/1.790712e7L;
    }
    break;

  case 6:
    switch (s) {
    case 0: return -(152-(3294-(34017-(222640-(1036035-(3644718-(
      10061051-(22323708-(40450905-(60472060-(75018042-(77380464-(
      66251822-(46802700-(26967270-(12421288-(4412826-1124838*k)*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)/3.02841e6L;
    case 1: return (93-(2006-(20608-(134090-(619850-(2164162-(5922224-(
      13007742-(23289915-(34321980-(41840128-(42226436-(35154028-(
      23921380-(13075040-(5556892-(1716099-317262*k)*k)*k)*k)*k)*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)*k)/6.05682e6L;
    case 2: return -(57-(1224-(12512-(80960-(371910-(1289288-(3499496-(
      7614288-(13483635-(19612560-(23535072-(23297344-(18929092-(
      12480720-(6537520-(2615008-(735471-115368*k)*k)*k)*k)*k)*k)*k)*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)/1.211364e7L;
    }
    break;

  case 7:
    switch (s) {
    case 0: return -(147-(3484-(39520-(285480-(1474070-(5788640-(
      17957940-(45125080-(93424045-(161226780-(233716340-(285867920-(
      295525620-(257934880-(189290920-(115892400-(58429085-(23746580-(
      7498920-1710280*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      *k)*k)/4.6046e6L;
    case 1: return (136-(3211-(36270-(260780-(1339520-(5229510-(16116100
      -(40191580-(82488120-(140917205-(201845930-(243374040-(247237120
      -(211132180-(150660120-(88850840-(42493880-(15935205-(4374370-
      723580*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      /1.38138e7L;
    case 2: return -(42-(988-(11115-(79560-(406640-(1578720-(4834830-(
      11971960-(24371490-(41244060-(58429085-(69535440-(69535440-(
      58243360-(40562340-(23178480-(10623470-(3749460-(937365-131560*k)
      *k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
      /13813800.0L;
    }
    break;
  }

  fprintf(stderr, "***ERROR in L\n");
  abort();

  return 0.0L;
}


long double FEMHF_M(int ik, int l, int s1, int s2)
{
  const int s = (s1==0) ? ((s2==0) ? 0 : 1) : ((s2==0) ? 1 : 2);
  const long double k   = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
 
  switch (l) {
  case 0:
    switch (s) {
    case 0: return 11/840.0L+19*k/210.0L+9*k2/35.0L+13*k3/35.0L;
    case 1: return 1/315.0L+17*k/840.0L+k2/20.0L+11*k3/210.0L;
    case 2: return 1/1260.0L+k/210.0L+3*k2/280.0L+k3/105.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0: return 3/35.0L + (13*k)/35.0L;
    case 1: return 1/60.0L+ (11*k)/210.0L;
    case 2: return 1/280.0L+ k/105.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0: return -1/10.0L/k10+1/9.0L/k9+19/40.0L/k8-13/105.0L/k7-137
      /180.0L/k6-5/14.0L/k5-11/840.0L/k4+19/630.0L/k3-3/35.0L/k2+13
      /35.0L/k+Power(1-2*k,2)*Power(1+k,4)*LOG10(ik);
    case 1: return 1/10.0L/k9+4/45.0L/k8-107/360.0L/k7-593/1260.0L/k6
      -227/1260.0L/k5-1/315.0L/k4+17/2520.0L/k3-1/60.0L/k2+11/210.0L/k+
      (-k-2*k2+2*k3+8*k4+7*k5+2*k6)*LOG10(ik);
    case 2: return -1/10.0L/k8-13/45.0L/k7-101/360.0L/k6-19/210.0L/k5-1
      /1260.0L/k4+1/630.0L/k3-1/280.0L/k2+1/105.0L/k+k2*Power(1+k,4)
      *LOG10(ik);
    }
    break;

  case 3:
    switch (s) {
    case 0: return 3/5.0L/k10+8/15.0L/k9-359/60.0L/k8-75/14.0L/k7-11
      /84.0L/k6+19/105.0L/k5-9/35.0L/k4+13/35.0L/k3-(6+12*k-54*k2-120*k3
      -60*k4)*LOG10(ik);
    case 1: return 1/5.0L/k10-37/45.0L/k9-233/60.0L/k8-227/84.0L/k7-2
       /63.0L/k6+17/420.0L/k5-1/20.0L/k4+11/210.0L/k3-(2-6*k-48*k2-70*k3
       -30*k4)*LOG10(ik);
    case 2: return -1/10.0L/k10-49/45.0L/k9-287/120.0L/k8-19/14.0L/k7-1
       /126.0L/k6+1/105.0L/k5-3/280.0L/k4+1/105.0L/k3+Power(1+k,2)*(1+10
       *k+15*k2)*LOG10(ik);
    }
    break;

  case 4:
    switch (s) {
    case 0: return -9/10.0L/k10-5/k9-11/24.0L/k8+19/42.0L/k7-3/7.0L/k6
      +13/35.0L/k5+(9+60*k*(1+k))*LOG10(ik);
    case 1: return -4/5.0L/k10-47/18.0L/k9-1/9.0L/k8+17/168.0L/k7-1
      /12.0L/k6+11/210.0L/k5+(8+35*k+30*k2)*LOG10(ik);
    case 2: return -3/5.0L/k10-4/3.0L/k9-1/36.0L/k8+1/42.0L/k7-1/56.0L
      /k6+1/105.0L/k5+(6+20*k+15*k2)*LOG10(ik);
    }
    break;

  case 5:
    switch (s) {
    case 0: return (-2/5.0L/k10-16/45.0L/k9-1/90.0L/k8+1/63.0L/k7+1/7.0L
      /k6+13/35.0L/k5)/Power(1+k,2)+4*LOG10(ik);
    case 1: return (-1/5.0L/k10-8/45.0L/k9-1/180.0L/k8+1/126.0L/k7-1
      /84.0L/k6+11/210.0L/k5)/Power(1+k,2)+2*LOG10(ik);
    case 2: return (-1/10.0L/k10-4/45.0L/k9-1/360.0L/k8+1/252.0L/k7-1
      /168.0L/k6+1/105.0L/k5)/Power(1+k,2)+LOG10(ik);
    }
    break;

  case 6:
    switch (s) {
    case 0: return (35.0L+140*k+200*k2+104*k3)/(280.0L*k8*Power(1+k,4));
    case 1: return (15+50*k+44*k2)/(840.0L*k7*Power(1+k,4));
    case 2: return (5+8*k)/(840.0L*k6*Power(1+k,4));
    }
    break;

  case 7:
    switch (s) {
    case 0: return (42.0L+252*k+623*k2+800*k3+540*k4+156*k5)/(420.0L*k10
      *Power(1+k,6));
    case 1: return (14.0L+77*k+165*k2+165*k3+66*k4)/(1260.0L*k9
      *Power(1+k,6));
    case 2: return (7+30*k+45*k2+24*k3)/(2520.0L*k8*Power(1+k,6));
    }
    break;
  }

  fprintf(stderr, "***ERROR in M\n");
  abort();

  return 0.0L;
}


/*----------------------------------------------------------------------
  J1

  l : order 5+2*l or 3-2*l
  s = s1*2 + s2 
----------------------------------------------------------------------*/
static long double J1(int l, int s)
{
  switch (s) {
  case 0: 
    return 72*(13+3*l)/(1.0L+l)/(3.0L+l)/(4.0L+l)/(5.0L+l)/(6.0L+l)
      /(7.0L+l);
  case 1: case 2: 
    return 24*(11+3*l)/(7.0L+l)/(6.0L+l)/(5.0L+l)/(4.0L+l)/(3.0L+l)
      /(2.0L+l);
  case 3:
    return 24/(7.0L+l)/(6.0L+l)/(5.0L+l)/(4.0L+l)/(3.0L+l);
  }

  fprintf(stderr, "***ERROR in J1\n");
  abort();

  return 0.0L; 
}


/*----------------------------------------------------------------------
  J2
----------------------------------------------------------------------*/
long double J2(int ik, int l, int s)
{
  const long double k   = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;


  switch (l) {

  case 5: /* l=0 high */
    switch (s) {
    case 0: 
      return 5*k/77.0L+38*k3/63.0L+26*k5/35.0L;
    case 1: case 2:
      return 13/6930.0L+4*k2/63.0L+k4/6.0L;
    case 3: 
      return k/231.0L+2*k3/63.0L+2*k5/105.0L;
    }
    break;

  case 7: /* l=1 high */
    switch (s) {
    case 0: 
      return 26*k7/35.0L+19*k5/15.0L+5*k3/11.0L+62*k/2145.0L;
    case 1: case 2:
      return 7*k6/30.0L+2*k4/9.0L+13*k2/330.0L+32/45045.0L;
    case 3: 
      return 2*k7/105.0L+k5/15.0L+k3/33.0L+14*k/6435.0L;
    }
    break;
  
  case 9: /* l=2 high */
    switch (s) {
    case 0:  
      return 26*k9/35.0L+76*k7/35.0L+18*k5/11.0L+248*k3/715.0L+74*k
        /5005.0L;
    case 1: case 2:
      return 3*k8/10.0L+8*k6/15.0L+13*k4/55.0L+128*k2/5005.0L+19
        /60060.0L;
    case 3: 
      return 2*k9/105.0L+4*k7/35.0L+6*k5/55.0L+56*k3/2145.0L+6*k
        /5005.0L;
    }
    break;

  case 11: /* l=3 high */
    switch (s) {
    case 0: 
      return 26*k11/35.0L+209*k9/63.0L+30*k7/7.0L+124*k5/65.0L+74*k3
        /273.0L+129*k/15470.0L;
    case 1: case 2:
      return 11*k10/30.0L+22*k8/21.0L+13*k6/15.0L+64*k4/273.0L+19*k2
        /1092.0L+11/69615.0L;
    case 3:
      return 2*k11/105.0L+11*k9/63.0L+2*k7/7.0L+28*k5/195.0L+2*k3/91.0L
        +11*k/15470.0;
    }
    break;

  case 13: /* l=4 high */
    switch (s) {
    case 0: 
      return 26*k13/35.0L+494*k11/105.0L+65*k9/7.0L+248*k7/35.0L+74*k5
        /35.0L+129*k3/595.0L+49*k/9690.0L;
    case 1: case 2:
      return 13*k12/30.0L+572*k10/315.0L+169*k8/70.0L+128*k6/105.0L+19
        *k4/84.0L+22*k2/1785.0L+1/11628.0L;
    case 3:
      return 2*k13/105.0L+26*k11/105.0L+13*k9/21.0L+8*k7/15.0L+6*k5
        /35.0L+11*k3/595.0L+13*k/29070.0L;
    }
    break;

  case 15: /* l=5 high */
    switch (s) {
    case 0: 
      return 26*k15/35.0L+19*k13/3.0L+195*k11/11.0L+62*k9/3.0L+74*k7
        /7.0L+387*k5/170.0L+343*k3/1938.0L+22*k/6783.0L;
    case 1: case 2:
      return k14/2.0L+26*k12/9.0L+169*k10/30.0L+32*k8/7.0L+19*k6/12.0L
        +11*k4/51.0L+35*k2/3876.0L+8/159885.0L;
    case 3:
      return 2*k15/105.0L+k13/3.0L+13*k11/11.0L+14*k9/9.0L+6*k7/7.0L
        +33*k5/170.0L+91*k3/5814.0L+2*k/6783.0L;
    }
    break;

  case 17: /* l=6 high */
    switch (s) {
    case 0: 
      return 26*k17/35.0L+2584*k15/315.0L+340*k13/11.0L+8432*k11/165.0L
        +2516*k9/63.0L+516*k7/35.0L+686*k5/285.0L+176*k3/1197.0L+366*k
        /168245.0L;
    case 1: case 2:
      return 17*k16/30.0L+272*k14/63.0L+5746*k12/495.0L+4352*k10/315.0L
        +323*k8/42.0L+88*k6/45.0L+35*k4/171.0L+64*k2/9405.0L+31
        /1009470.0L;
    case 3:
      return 2*k17/105.0L+136*k15/315.0L+68*k13/33.0L+1904*k11/495.0L
        +68*k9/21.0L+44*k7/35.0L+182*k5/855.0L+16*k3/1197.0L+34*k
        /168245.0L;
    }
    break;
  
  case 19: /* l=7 high */
    switch (s) {
    case 0: 
      return 26*k19/35.0L+361*k17/35.0L+3876*k15/77.0L+80104*k13/715.0L
        +47804*k11/385.0L+2451*k9/35.0L+98*k7/5.0L+88*k5/35.0L+1098*k3
        /8855.0L+67*k/44275.0L;
    case 1: case 2:
      return 19*k18/30.0L+646*k16/105.0L+8398*k14/385.0L+41344*k12
        /1155.0L+6137*k10/210.0L+418*k8/35.0L+7*k6/3.0L+32*k4/165.0L+93
        *k2 /17710.0L+34/1726725.0L;
    case 3: 
      return 2*k19/105.0L+19*k17/35.0L+1292*k15/385.0L+18088*k13/2145.0L
        +3876*k11/385.0L+209*k9/35.0L+26*k7/15.0L+8*k5/35.0L+102*k3
        /8855.0L+19*k/132825.0L;
    }
    break;

  case 3: /* l = 0 low */
    switch (s) {
    case 0: 
      return 19*k/105.0L+26*k3/35.0L;
    case 1: case 2:
      return 2/315.0L+k2/10.0L;
    case 3: 
      return k/105.0L+2*k3/105.0L;
    }
    break;

  case 1: /* l = -1 low */
    switch (s) {
    case 0: 
      return 26*k/35.0L;
    case 1: case 2:
      return 1/30.0L;
    case 3: 
      return 2*k/105.0L;
    }
    break;

  case -1: /* l = -2 low */
    switch (s) {
    case 0: 
      return 2/9.0L/k9-26/105.0L/k7-5/7.0L/k5+19/315.0L/k3+26/35.0L/k
        -Power(k-1,4)*Power(1+2*k,2)*LOG10(-ik)+Power(1+k,4)
        *Power(1-2*k,2)*LOG10(ik);
    case 1: case 2:
      return 8/45.0L/k8-593/630.0L/k6-2/315.0L/k4-1/30.0L/k2
        +Power(k-1,4)*k*(1+2*k)*LOG10(-ik)+k*Power(1+k,4)*(-1+2*k)
        *LOG10(ik);
    case 3:
      return -26/45.0L/k7-19/105.0L/k5+1/315.0L/k3+2/105.0L/k
        -Power(k-1,4)*k2*LOG10(-ik)+k2*Power(1+k,4)*LOG10(ik);
    }
    break;

  case -3: /* l = -3 low */
    switch (s) {
    case 0: 
      return 16/15.0L/k9-75/7.0L/k7+38/105.0L/k5+26/35.0L/k3-6
        *Power(k-1,2)*(-1+10*k2)*LOG10(-ik)+6*Power(1+k,2)*(-1+10*k2)
        *LOG10(ik);
    case 1: case 2:
      return 2/5.0L/k10-233/30.0L/k8-4/63.0L/k6-1/10.0L/k4+2
        *Power(k-1,2)*(-1+5*k*(-1+3*k))*LOG10(-ik)+2*Power(1+k,2)
        *(-1+5*k*(1+3*k))*LOG10(ik);
    case 3:
      return -98/45.0L/k9-19/7.0L/k7+2/105.0L/k5+2/105.0L/k3
        -Power(k-1,2)*(1+5*k*(-2+3*k))*LOG10(-ik)+Power(1+k,2)*(1+5*k*(2
        +3*k))*LOG10(ik);
    }
    break;
 
  case -5: /* l = -4 low */
    switch (s) {
    case 0:
      return -10/k9+19/21.0L/k7+26/35.0L/k5-3*(3+20*(-1+k)*k)*LOG10(-ik)
        +3*(3+20*k*(1+k))*LOG10(ik);
    case 1: case 2:
      return -8/5.0L/k10-2/9.0L/k8-1/6.0L/k6-(-8+5*(7-6*k)*k)*LOG10(-ik)
        +(8+5*k*(7+6*k))*LOG10(ik);
    case 3:
      return -8/3.0L/k9+1/21.0L/k7+2/105.0L/k5-(6+5*k*(-4+3*k))
        *LOG10(-ik)+(6+5*k*(4+3*k))*LOG10(ik);
    }
    break;
 
  case -7: /* l = -5 low */
    switch (s) {
    case 0:
      return (2*(140-100*k2+32*k4+117*k6))/(315.0L*k9*Power(k2-1,2))
        -4*LOG10(-ik)+4*LOG10(ik);
    case 1: case 2:
      return (-12+9*k2-2*k4-7*k6)/(30.0L*k10*Power(k2-1,2))+2
        *LOG10(-ik)+2*LOG10(ik);
    case 3:
      return (2*(35-25*k2+8*k4+3*k6))/(315.0L*k9*Power(k2-1,2))
        -LOG10(-ik)+LOG10(ik);
    }
    break;

  case -9: /* l=6 low */
    switch (s) {
    case 0:
      return (1-9*k2+26*k4)/(35.0L*k5*Power(k2-1,4));
    case 1: case 2:
      return (-5+32*k2-63*k4)/(210.0L*k6*Power(k2-1,4));
    case 3:
      return (-3+7*k2+2*k4)/(105.0L*k5*Power(k2-1,4));
    }
    break;

  case -11: /* l=7 low */
    switch (s) {
    case 0: 
      return (1-8*k2+27*k4-50*k6+78*k8)/(105.0L*k7*Power(k2-1,6));
    case 1: case 2:
      return (-7+50*k2-150*k4+242*k6-231*k8)/(630.0L*k8*Power(k2-1,6));
    case 3:
      return (-3+16*k2-33*k4+30*k6+6*k8)/(315.0L*k7*Power(k2-1,6));
    }
    break;
  }

  fprintf(stderr, "***ERROR in J2\n");
  abort();

  return 0.0L;
}


/*----------------------------------------------------------------------
  J3
----------------------------------------------------------------------*/
long double J3(int ik, int l, int s)
{
  const long double k   = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;

  switch (l) {
  case 5: /* l=0 high */
    switch (s) {
    case 0: 
      return 41/3960.0L+23*k/308.0L+19*k2/84.0L+23*k3/63.0L+9*k4/28.0L+9
        *k5/70.0L;
    case 1: 
      return -7/3960.0L-25*k/1848.0L-11*k2/252.0L-19*k3/252.0L-k4/14.0L
        -13*k5/420.0L;
    case 2:
      return 1/330.0L+17*k/792.0L+4*k2/63.0L+25*k3/252.0L+k4/12.0L+13*k5
        /420.0L;
    case 3: 
      return -1/1980.0L-k/264.0L-k2/84.0L-5*k3/252.0L-k4/56.0L-k5
        /140.0L;
    }
    break;

  case 7: /* l=1 high */
    switch (s) {
    case 0: return 111/20020.0L+112*k/2145.0L+287*k2/1320.0L+23*k3/44.0L
      +19*k4/24.0L+23*k5/30.0L+9*k6/20.0L+9*k7/70.0L;
    case 1: return -17/20020.0L-217*k/25740.0L-49*k2/1320.0L-25*k3
      /264.0L-11*k4/72.0L-19*k5/120.0L -k6/10.0L-13*k7/420.0L;
    case 2: return 5/3003.0L+133*k/8580.0L+7*k2/110.0L+119*k3/792.0L+2
      *k4/9.0L+5*k5/24.0L+7*k6/60.0L+13*k7/420.0L;
    case 3: return -1/4004.0L-7*k/2860.0L-7*k2/660.0L-7*k3/264.0L-k4
      /24.0L-k5/24.0L-k6/40.0L-k7/140.0L;
    }
    break;

  case 9: /* l=2 high */
    switch (s) {
    case 0: return 3/910.0L+381*k/10010.0L+999*k2/5005.0L+448*k3/715.0L
      +287*k4/220.0L+207*k5/110.0L+19*k6/10.0L+46*k7/35.0L+81*k8/140.0L
      +9*k9/70.0L;
    case 1: return -1/2184.0L-111*k/20020.0L-153*k2/5005.0L-217*k3
      /2145.0L-49*k4/220.0L-15*k5/44.0L-11*k6/30.0L-19*k7/70.0L-9*k8
      /70.0L-13*k9/420.0L;
    case 2: return 11/10920.0L+3*k/260.0L+60*k2/1001.0L+133*k3/715.0L
      +21*k4/55.0L+119*k5/220.0L+8*k6/15.0L+5*k7/14.0L+3*k8/20.0L+13*k9
      /420.0L;
    case 3: return -1/7280.0L-3*k/1820.0L-9*k2/1001.0L-21*k3/715.0L-7*k4
      /110.0L-21*k5/220.0L-k6/10.0L-k7/14.0L-9*k8/280.0L-k9/140.0;
    }
    break;
 
  case 11: /* l=3 high */
    switch (s) {
    case 0: return 181/85680.0L+891*k/30940.0L+33*k2/182.0L+127*k3
      /182.0L+333*k4/182.0L+224*k5/65.0L+287*k6/60.0L+69*k7/14.0L+209
      *k8/56.0L+253*k9/126.0L+99*k10/140.0L+9*k11/70.0L;
    case 1: return -23/85680.0L-473*k/123760.0L-55*k2/2184.0L-37*k3
      /364.0L-51*k4/182.0L-217*k5/390.0L-49*k6/60.0L-25*k7/28.0L-121*k8
      /168.0L-209*k9/504.0L-11*k10/70.0L-13*k11/420.0L;
    case 2: return 1/1530.0L+253*k/28560.0L+121*k2/2184.0L+11*k3/52.0L
      +50*k4/91.0L+133*k5/130.0L+7*k6/5.0L+17*k7/12.0L+22*k8/21.0L+275
      *k9/504.0L+11*k10/60.0L+13*k11/420.0L;
    case 3: return -1/12240.0L-11*k/9520.0L-11*k2/1456.0L-11*k3/364.0L
      -15*k4/182.0L-21*k5/130.0L-7*k6/30.0L-k7/4.0L-11*k8/56.0L-55*k9
      /504.0L-11*k10/280.0L-k11/140.0L;
    }
    break;

  case 13: /* l=4 high */
    switch (s) {
    case 0: return 37/25840.0L+871*k/38760.0L+2353*k2/14280.0L+891*k3
      /1190.0L+33*k4/14.0L+381*k5/70.0L+333*k6/35.0L+64*k7/5.0L+533*k8
      /40.0L+299*k9/28.0L+2717*k10/420.0L+299*k11/105.0L+117*k12/140.0L
      +9*k13/70.0L;
    case 1: return -13/77520.0L-637*k/232560.0L-299*k2/14280.0L-473*k3
      /4760.0L-55*k4/168.0L-111*k5/140.0L-51*k6/35.0L-31*k7/15.0L-91*k8
      /40.0L-325*k9/168.0L-1573*k10/1260.0L-247*k11/420.0L-13*k12/70.0L
      -13*k13/420.0L;
    case 2: return 13/29070.0L+325*k/46512.0L+13*k2/255.0L+3289*k3
      /14280.0L+121*k4/168.0L+33*k5/20.0L+20*k6/7.0L+19*k7/5.0L+39*k8
      /10.0L+221*k9/72.0L+572*k10/315.0L+65*k11/84.0L+13*k12/60.0L+13
      *k13/420.0L;
    case 3: return -1/19380.0L-13*k/15504.0L-13*k2/2040.0L-143*k3
      /4760.0L-11*k4/112.0L-33*k5/140.0L-3*k6/7.0L-3*k7/5.0L-13*k8/20.0L
      -13*k9/24.0L-143*k10/420.0L-13*k11/84.0L-13*k12/280.0L-k13/140.0L;
    }
    break;
 
  case 15: /* l=5 high */
    switch (s) {
    case 0: return 89/87780.0L+122*k/6783.0L+777*k2/5168.0L+6097*k3
      /7752.0L+2353*k4/816.0L+2673*k5/340.0L+33*k6/2.0L+381*k7/14.0L
      +999*k8/28.0L+112*k9/3.0L+3731*k10/120.0L+897*k11/44.0L+247*k12 
      /24.0L+23*k13/6.0L+27*k14/28.0L+9*k15/70.0L;
    case 1: return -29/263340.0L-55*k/27132.0L-91*k2/5168.0L-4459*k3
      /46512.0L-299*k4/816.0L-1419*k5/1360.0L-55*k6/24.0L-111*k7/28.0L
      -153*k8/28.0L-217*k9/36.0L-637*k10/120.0L-325*k11/88.0L-143*k12
      /72.0L-19*k13/24.0L-3*k14/14.0L-13*k15/420.0L;
    case 2: return 1/3135.0L+3*k/532.0L+91*k2/1938.0L+11375*k3/46512.0L
      +91*k4/102.0L+3289*k5/1360.0L+121*k6/24.0L+33*k7/4.0L+75*k8/7.0L
      +133*k9/12.0L+91*k10/10.0L+1547*k11/264.0L+26*k12/9.0L+25*k13
      /24.0L+k14/4.0L+13*k15/420.0L;
    case 3: return -1/29260.0L-k/1596.0L-7*k2/1292.0L-455*k3/15504.0L
      -91*k4/816.0L-429*k5/1360.0L-11*k6/16.0L-33*k7/28.0L-45*k8/28.0L
      -7*k9/4.0L-91*k10/60.0L-91*k11/88.0L-13*k12/24.0L-5*k13/24.0L-3
      *k14/56.0L-k15/140.0L;
    }
    break;

  case 17: /* l=6 high */
    switch (s) {
    case 0: return 79/106260.0L+4947*k/336490.0L+3026*k2/21945.0L+976*k3
      /1197.0L+259*k4/76.0L+6097*k5/570.0L+2353*k6/90.0L+1782*k7/35.0L
      +561*k8/7.0L+2159*k9/21.0L+3774*k10/35.0L+15232*k11/165.0L+63427
      *k12/990.0L+391*k13/11.0L+323*k14/21.0L+1564*k15/315.0L+153*k16
      /140.0L+9*k17/70.0L;
    case 1: return -2/26565.0L-1037*k/672980.0L-986*k2/65835.0L-110*k3
      /1197.0L-91*k4/228.0L-4459*k5/3420.0L-299*k6/90.0L-473*k7/70.0L
      -935*k8/84.0L-629*k9/42.0L-578*k10/35.0L-7378*k11/495.0L-10829*k12
      /990.0L-425*k13/66.0L-187*k14/63.0L-323*k15/315.0L-17*k16/70.0L-13
      *k17/420.0L;
    case 2: return 5/21252.0L+493*k/106260.0L+136*k2/3135.0L+34*k3
      /133.0L+182*k4/171.0L+2275*k5/684.0L+364*k6/45.0L+3289*k7/210.0L
      +2057*k8/84.0L+187*k9/6.0L+680*k10/21.0L+4522*k11/165.0L+3094*k12
      /165.0L+2023*k13/198.0L+272*k14/63.0L+85*k15/63.0L+17*k16/60.0L+13
      *k17/420.0L;
    case 3: return -1/42504.0L-17*k/35420.0L-34*k2/7315.0L-34*k3/1197.0L
      -7*k4/57.0L-91*k5/228.0L-91*k6/90.0L-143*k7/70.0L-187*k8/56.0L-187
      *k9/42.0L-34*k10/7.0L-238*k11/55.0L-1547*k12/495.0L-119*k13/66.0L
      -17*k14/21.0L-17*k15/63.0L-17*k16/280.0L-k17/140.0L;
    }
    break;

  case 19: /* l=7 high */
    switch (s) {
    case 0: return 369/657800.0L+1083*k/88550.0L+4503*Power(k,2)
      /35420.0L+14841*k3/17710.0L+1513*k4/385.0L+488*k5/35.0L+777*k6
      /20.0L+871*k7/10.0L+44707*k8/280.0L+16929*k9/70.0L+10659*k10/35.0L
      +123063*k11/385.0L+107559*k12/385.0L+144704*k13/715.0L+13243*k14
      /110.0L+22287*k15/385.0L+6137*k16/280.0L+437*k17/70.0L+171*k18
      /140.0L+9*k19/70.0L;
    case 1: return -7/131560.0L-1273*k/1.0626e6L-114*k2/8855.0L-3111*k3
      /35420.0L-493*k4/1155.0L-11*k5/7.0L-91*k6/20.0L-637*k7/60.0L-5681
      *k8/280.0L-8987*k9/280.0L-3553*k10/84.0L-35853*k11/770.0L-16473
      *k12/385.0L-70091*k13/2145.0L-2261*k14/110.0L-1615*k15/154.0L-3553
      *k16/840.0L-361*k17/280.0L-19*k18/70.0L-13*k19/420.0L;
    case 2: return 4/22425.0L+589*k/151800.0L+285*k2/7084.0L+9367*k3
      /35420.0L+68*k4/55.0L+153*k5/35.0L+182*k6/15.0L+325*k7/12.0L+247
      *k8/5.0L+62491*k9/840.0L+39083*k10/420.0L+969*k11/10.0L+6460*k12
      /77.0L+42959*k13/715.0L+1938*k14/55.0L+5491*k15/330.0L+646*k16
      /105.0L+95*k17/56.0L+19*k18/60.0L+13*k19/420.0L;
    case 3: return -1/59800.0L-19*k/50600.0L-57*k2/14168.0L-969*k3
      /35420.0L-51*k4/385.0L-17*k5/35.0L-7*k6/5.0L-13*k7/4.0L-247*k8
      /40.0L-2717*k9/280.0L-3553*k10/280.0L-969*k11/70.0L-969*k12/77.0L
      -6783*k13/715.0L-323*k14/55.0L-323*k15/110.0L-323*k16/280.0L-19
      *k17/56.0L-19*k18/280.0L-k19/140.0L;
    }
    break;

  case 3: /* l=0 low */
    switch (s) {
    case 0: return 19/840.0L+23*k/210.0L+27*k2/140.0L+9*k3/70.0L;
    case 1: return -11/2520.0L-19*k/840.0L-3*k2/70.0L-13*k3/420.0L;
    case 2: return 2/315.0L+5*k/168.0L+k2/20.0L+13*k3/420.0L;
    case 3: return -1/840.0L-k/168.0L-3*k2/280.0L-k3/140.0l;
    }
    break;

  case 1: /* l=1 low */
    switch (s) {
    case 0: return 9/140.0L+9*k/70.0L;
    case 1: return -1/70.0L-13*k/420.0L;
    case 2: return 1/60.0L+13*k/420.0L;
    case 3: return -1/280.0L-k/140.0L;
    }
    break;

  case -1: /* l = -1 low */
    switch (s) {
    case 0: return -3/10.0L/k8+2/15.0L/k7+269/360.0L/k6+53/140.0L/k5-19
      /840.0L/k4+23/630.0L/k3-9/140.0L/k2+9/70.0L/k-k2*Power(1+k,2)*(-3
      +4*k*(1+k))*LOG10(ik);
    case 1: return 1/10.0L/k8-1/90.0L/k7-103/360.0L/k6-31/168.0L/k5+11
      /2520.0L/k4-19/2520.0L/k3+1/70.0L/k2-13/420.0L/k+k2*Power(1+k,3)
      *(-1+2*k)*LOG10(ik);
    case 2: return 3/10.0L/k7+7/15.0L/k6+67/360.0L/k5-2/315.0L/k4+5
      /504.0L/k3-1/60.0L/k2+13/420.0L/k-k3*Power(1+k,2)*(3+2*k)
      *LOG10(ik);
    case 3: return -1/10.0L/k7-17/90.0L/k6-11/120.0L/k5+1/840.0L/k4-1
      /504.0L/k3+1/280.0L/k2-1/140.0L/k+k3*Power(1+k,3)*LOG10(ik);
    }
    break;
  
  case -3: /* l = -3 low */
    switch (s) {
    case 0: return -3/10.0L/k10-4/15.0L/k9+683/120.0L/k8+159/28.0L/k7-19
      /84.0L/k6+23/105.0L/k5-27/140.0L/k4+9/70.0L/k3-3*(-1+2*k*(1+k)*(-1
      +10*k*(1+k)))*LOG10(ik);
    case 1: return 1/10.0L/k10+17/90.0L/k9-241/120.0L/k8-155/56.0L/k7+11
      /252.0L/k6-19/420.0L/k5+3/70.0L/k4-13/420.0L/k3+(1+k)*(-1+2*k*(-1+
      5*k*(2+3*k)))*LOG10(ik);
    case 2: return 9/10.0L/k9+19/5.0L/k8+67/24.0L/k7-4/63.0L/k6+5/84.0L
      /k5-1/20.0L/k4+13/420.0L/k3-(9*k+48*k2+70*k3+30*k4)*LOG10(ik);
    case 3: return -3/10.0L/k9-22/15.0L/k8-11/8.0L/k7+1/84.0L/k6-1/84.0L
      /k5+3/280.0L/k4-1/140.0L/k3+(3*k+18*k2+30*k3+15*k4)*LOG10(ik);
    }
    break;

  case -5: /* l = -4 low */
    switch (s) {
    case 0: return (9/10.0L/k10+34/5.0L/k9+1363/120.0L/k8+153/28.0L/k7
      -1/56.0L/k6+1/30.0L/k5-9/140.0L/k4+9/70.0L/k3)/Power(1+k,2)-3*(3
      +20*k*(1+k))*LOG10(ik);
    case 1: return (-3/10.0L/k10-37/15.0L/k9-199/72.0L/k8+5/126.0L/k7-1
      /24.0L/k6+17/420.0L/k5-13/420.0L/k4)/(1+k)+(3+5*k*(5+6*k))
      *LOG10(ik);
    case 2: return (4/5.0L/k10+379/90.0L/k9+92/15.0L/k8+153/56.0L/k7-1
      /126.0L/k6+11/840.0L/k5-3/140.0L/k4+13/420.0L/k3)/Power(1+k,2)-(8
      +35*k+30*k2)*LOG10(ik);
    case 3: return (-3/10.0L/k10-22/15.0L/k9-11/8.0L/k8+1/84.0L/k7-1
      /84.0L/k6+3/280.0L/k5-1/140.0L/k4)/(1+k)+(3+15*k+15*k2)*LOG10(ik);
    }
    break;
 
  case -7: /* l = -5 low */
    switch (s) {
    case 0: return (2/5.0L/k10+52/45.0L/k9+101/90.0L/k8+38/105.0L/k7+1
      /315.0L/k6-2/315.0L/k5+9/140.0L/k4+9/70.0L/k3)/Power(1+k,4)-4
      *LOG10(ik);
    case 1: return (-1/5.0L/k10-17/45.0L/k9-11/60.0L/k8+1/420.0L/k7-1
      /252.0L/k6+1/140.0L/k5-13/420.0L/k4)/Power(1+k,3)+2*LOG10(ik);
    case 2: return (1/5.0L/k10+26/45.0L/k9+101/180.0L/k8+19/105.0L/k7+1
      /630.0L/k6-1/315.0L/k5+1/140.0L/k4+13/420.0L/k3)/Power(1+k,4)-2
      *LOG10(ik);
    case 3: return (-1/10.0L/k10-17/90.0L/k9-11/120.0L/k8+1/840.0L/k7-1
      /504.0L/k6+1/280.0L/k5-1/140.0L/k4)/Power(1+k,3)+LOG10(ik);
    }
    break;

  case -9: /* l=6 low */
    switch (s) {
    case 0: return (5+28*k+54*k2+36*k3)/(280.0L*k6*Power(1+k,6));
    case 1: return -(5+22*k+26*k2)/(840.0L*k6*Power(1+k,5));
    case 2: return  (9+30*k+26*k2)/(840.0L*k5*Power(1+k,6));
    case 3: return -(1+2*k)/(280.0L*k5*Power(1+k,5));
    }
    break;

  case -11: /* l=7 low */
    switch (s) {
    case 0: return (7+54*k+174*k2+296*k3+270*k4+108*k5)/(840.0L*k8
      *Power(1+k,8));
    case 1: return -(7+2*k*(23+3*k*(20+k*(25+13*k))))/(2520.0L*k8
      *Power(1+k,7));
    case 2: return (9+2*k*(28+3*k*(23+k*(27+13*k))))/(2520.0L*k7
      *Power(1+k,8));
    case 3: return -(1+5*k+9*k2+6*k3)/(840.0L*k7*Power(1+k,7));
    }
    break;
  }

  fprintf(stderr, "***ERROR in J3\n");
  abort();

  return 0.0L;
}




/*----------------------------------------------------------------------
  K1
  
  k1 = k2 = k3 = k4 = 0
----------------------------------------------------------------------*/
static long double K1(int il, int s)
{
  const long double l = (long double)il;

  switch (s) {
  case 0: 
    return (14419218+11762731*l+3970469*Power(l,2)+692106*Power(l,3)
      +61820*Power(l,4)+2248*Power(l,5))/(1.293292e7L*(3+l)*(4+l)*(5+l)
      *(6+l)*(9+2*l)*(11+2*l));
  case 1: case 2: case 4: case 8: 
    return (78116973+87228356*l+40546448*Power(l,2)+10131560*Power(l,3)
      +1441420*Power(l,4)+110608*Power(l,5)+3568*Power(l,6))
      /(3.879876e7L*(3+l)*(4+l)*(5+l)*(6+l)*(7+2*l)*(9+2*l)*(11+2*l));
  case 3: case 6: case 9: case 12:
    return (10103280+8108438*l+2677539*Power(l,2)+456697*Power(l,3)
      +39984*Power(l,4)+1428*Power(l,5))/(5.819814e7L*(4+l)*(5+l)*(6+l)
      *(7+2*l)*(9+2*l)*(11+2*l));
  case 5: case 10:
    return (5853837+4916959*l+1674971*Power(l,2)+292854*Power(l,3)+26180
      *Power(l,4)+952*Power(l,5))/(7.759752e7L*(3+l)*(4+l)*(5+l)*(6+l)
      *(9+2*l)*(11+2*l));
  case 7: case 11: case 13: case 14:
    return (1770714+1443507*l+479012*Power(l,2)+81836*Power(l,3)+7168
      *Power(l,4)+256*Power(l,5))/(3.879876e7L*(4+l)*(5+l)*(6+l)*(7+2*l)
      *(9+2*l)*(11+2*l));
  case 15:
    return (15378+8387*l+1870*Power(l,2)+196*Power(l,3)+8*Power(l,4))
      /(8.95356e6L*(11880+12126*l+4925*Power(l,2)+995*Power(l,3)+100
      *Power(l,4)+4*Power(l,5)));
  }

  fprintf(stderr, "***ERROR in K1\n");
  abort();
  
  return 0.0L;
}


/*----------------------------------------------------------------------
  K2
  
  k1 = k2 = k3 = k4 > 0
----------------------------------------------------------------------*/
static long double K2(int ik, int l, int s)
{
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;
  const long double k20 = k10*k10;
  const long double k21 = k10*k11;
  const long double k22 = k11*k11;

  switch (l) {
  case 0:
    switch (s) {
    case 0:
      return -3256159*k/1.22216094e10L+19*k2/1617.0L-2421031*k3
        /1.38881925e8L+1636*k4/10395.0L-40944*k5/425425.0L+6422*k6
        /11025.0L-93944*k7/525525.0L+676*k8/1225.0L;
    case 1: case 2: case 4: case 8:
      return -184967/29331862560.0L+547*k/1.4553e6L-10077437*k2
        /1.22216094e10L+50657*k3/4.3659e6L-1583039*k4/1.7153136e8L
        +3139*k5/44100.0L-197*k6/8580.0L+52*k7/525.0L;
    case 3: case 6: case 9: case 12:
      return 13/1091475.0L -3004297*k/7.33296564e10L+2579*k2/4.3659e6L
        -1869317*k3/3.3331662e9L+k4/135.0L-261649*k5/8.576568e7L+k6
        /60.0L+554*k7/225225.0L;
    case 5: case 10: 
      return -287*k/2.078505e7L+17*k2/24255.0L-9601621*k3/1.04756652e10L
        +116*k4/14553.0L-25141*k5/6.8068e6L+754*k6/33075.0L-881*k7
        /150150.0L+52*k8/3675.0L;
    case 7: case 11: case 13: case 14:
      return -7417/14665931280.0L+k/44100.0L-243739*k2/6.1108047e9L+253
        *k3/396900.0L-33151*k4/8.576568e7L+403*k5/132300.0L-191*k6
        /540540.0L+4*k7/1575.0L;
    case 15:
      return -2143*k/4.0738698e9L+k2/24255.0L-229267*k3/5.2378326e9L+4
        *k4/10395.0L-5507*k5/7.14714e7L+26*k6/33075.0L-89*k7/525525.0L
        +4*k8/11025.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0:
      return -56377*k/4.526522e7L+124*k2/5775.0L-1415543*k3/1.616615e7L
        +26*k4/77.0L-656168*k5/1.276275e6L+494*k6/525.0L-93944*k7
        /175175.0L+676*k8/1225.0L;
    case 1: case 2: case 4: case 8:
      return -174269/5431826400.0L+2351*k/3.15315e6L-95659*k2/2.28228e7L
        +171*k3/7700.0L-204823*k4/4.900896e6L+653*k5/6300.0L-197*k6
        /2860.0L+52*k7/525.0L;
    case 3: case 6: case 9: case 12:
      return 16/675675.0L-599099*k/2.7159132e9L+13*k2/9900.0L-4459207*k3
        /1.7459442e9L+k4/135.0L-7253*k5/4.7124e6L+7*k6/900.0L+554*k7
        /75075.0L;
    case 5: case 10:
      return -505*k/9.053044e6L+244*k2/225225.0L-119771*k3/2.98452e7L+6
        *k4/385.0L-1295951*k5/6.12612e7L+58*k6/1575.0L-881*k7/50050.0L
        +52*k8/3675.0L;
    case 7: case 11: case 13: case 14:
      return -10627/4073869800.0L+37*k/859950.0L-18491*k2/1.058148e8L+61
        *k3/69300.0L-3133*k4/2.82744e6L+61*k5/18900.0L-191*k6/180180.0L
        +4*k7/1575.0L;
    case 15:
      return -29*k/3.133746e7L+4*k2/96525.0L-70183*k3/5.819814e8L+2*k4
        /3465.0L-4507*k5/7.0686e6L+2*k6/1575.0L-89*k7/175175.0L+4*k8
        /11025.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0:
      return (-33693660+k*(453332880+k*(-2775420648+k*(10144278144+k*(
        -24151042916+k*(37022457472+k*(-27718294989+k*(-26190987284+k*(
        124496397157+k*(-230717277012+k*(288930949056+k*(-265198213020
        +k*(177727185372+k*(-81592834868+k*(20742015240+k*(1819784340+13
        *k*(-427882309+2*k*(182427854+k*(-317760090+k*(709637460+11*k*(
        -159512152+k*(394789980+k*(-813680753+76*k*(18174377+k*(
        -21551460+17*k*(1306877+6*k*(-117430+72501*k))))))))))))))))))))
        )))))))/(2.912816907e11L*Power(k,19))
        +(Power(k-1,18)*Power(1+2*k,2)*(11+6*k+k2)*LOG20(-ik))/10010.0L
        -(Power(1-2*k,2)*Power(1+k,4)*(11-148*k+918*k2-3472*k3+8918*k4
        -16380*k5+22022*k6-21736*k7+15444*k8-7436*k9+2002*k10-182*k12
        -28*k13+18*k14+8*k15+k16)*LOG20(ik))/10010.0L;
    case 1: case 2: case 4: case 8:
      return (116396280+k*(-1907673768+k*(14179733568+k*(-61655818224
        +k*(166158624632+k*(-250170951030+k*(22049152818+k*(895288088948
        +k*(-2471076246894+k*(3929473346202+k*(-4263922743120+k
        *(3241870034520+k*(-1658791287744+k*(489085041792+k*(-628456176
        +k*(-76858798752+13*k*(-75428176+k*(716917656+k2*(-536147957+11
        *k*(363413760+k*(-1314762186+19*k*(182583128+3*k*(-100068175+136
        *k*(1063205+21*k*(-34475+29744*k)))))))))))))))))))))))))
        /(6.9907605768e12L*k19)
        -(Power(k-1,18)*(1+2*k)*(19+118*k+67*k2+12*k3)*LOG20(-ik))
        /120120.0L-((-19+186*k-666*k2+648*k3+2409*k4-8274*k5+7112*k6
        +10560*k7-29142*k8+16900*k9+22308*k10-40560*k11+13650*k12+19740
        *k13-18456*k14+9216*k15-2601*k16-17406*k17-2126*k18+2544*k19
        +1347*k20+286*k21+24*k22)*LOG20(ik))/120120.0L;
    case 3: case 6: case 9: case 12:
      return (104756652+k*(-1389404016+k*(7932564640+k*(-24227291088+k*(
        34923449561+k*(21927565660+k*(-224013019780+k*(547969052708+k*(
        -804692328845+k*(788731063904+k*(-508791350864+k*(186192650504
        +k*(-711687608+k*(-37156456000+13*k*(-116697392+k*(167641952+k*(
        -293101828+k*(657862396+11*k*(-87184160+k*(166979055+k*(
        -289109548+3*k*(148143198+19*k*(-8440432+k*(9782335+476*k*(-9009
        +11080*k)))))))))))))))))))))))))/(3.4953802884e12L*k18)
        +(Power(k-1,18)*k*(1+2*k)*(19+k*(14+3*k))*LOG20(-ik))/60060.0L
        -(k*Power(1+k,4)*(-1+2*k)*(-19+k*(252+k*(-1536+k*(5684+k*(-14196
        +k*(25116+k*(-32032+k*(29172+k*(-18018+k*(6292+k2*(-1092+k*(-364
        +k*(84+k*(96+k*(28+3*k)))))))))))))))*LOG20(ik))/60060.0L;
    case 5: case 10:
      return (-33693660+k*(441080640+k*(-1477908432+k*(-6103633536+k*(
        74417483140+k*(-334108879104+k*(934209115455+k*(-1837386194336
        +k*(2677106860131+2*k*(-1482633126192+k*(1262250579204+k*(
        -823590879408+k*(402028726944+k*(-135585019792+k*(20410174908+k*
        (1574871120+13*k*(2146252+k*(303437904+k*(-385910166+k*(
        1069931040+11*k*(-221014332+k*(506577696+k*(-972089243+19*k*(
        81617952+k*(-88558365+136*k*(613756+9*k*(-30835+14872*k)))))))))
        ))))))))))))))))))/(6.9907605768e12L*k19)+(Power(k-1,18)*(11+k*(
        54+k*(219+4*k*(31+6*k))))*LOG20(-ik))/240240.0L-(Power(1+k,4)*(
        11+k*(-188+k*(1616+k*(-8944+k*(34660+k*(-97888+k*(206024+k*(
        -326768+k*(390962+k*(-348920+k*(225368+k*(-98384+k*(24388+k*(
        -1120+k*(-1928+k*(464+k*(667+4*k*(53+6*k))))))))))))))))))
        *LOG20(ik))/240240.0L;
    case 7: case 11: case 13: case 14:
      return (-60648588+k*(142125984+k*(4191931744+k*(-39678358720+k*(
        180026829983+k*(-520285593828+k*(1052046728788+k*(-1550352479764
        +k*(1688566040733+4*k*(-336990910952+k*(189276585184+k*(
        -66063981432+k*(6008402498+k*(16600248+k*(-1332684952+13*k*(
        1854248+k*(257811509+k2*(20994286+11*k*(20672+k*(-4630857+19*k*(
        1323348+k*(-2154925+68*k*(67925+7*k*(-4775+6864*k)))))))))))))))
        )))))))))/(1.39815211536e13L*k18)-(Power(k-1,18)*k*(11+k*(70+3*k
        *(17+4*k)))*LOG20(-ik))/240240.0L-(k*Power(1+k,4)*(-11+k*(128+k
        *(-648+k*(1792+k*(-2548+k2*(8008+k*(-18304+k*(23166+k*(-18304+k
        *(8008+k2*(-2548+k*(-336+k*(696+k*(448+3*k*(39+4*k))))))))))))))
        )*LOG20(ik))/240240.0L;
    case 15:
      return (94342248+k*(-1235025792+k*(7440633200+k*(-27295107840+k*(
        67953218333+k*(-121076331376+k*(158541178125+2*k*(-77021929688+k
        *(55297779460+k*(-28820358312+k*(10681005884+k*(-3143834680+k*(
        2012012196+k*(55446840+13*k*(192468012+k*(10320648+k*(9285706+k*
        (32026096+11*k*(-5744694+k*(11591824+k*(-19923321+19*k*(1580592
        +k*(-1555055+136*k*(10582+3*k*(-1335+572*k))))))))))))))))))))))
        )))/(3.4953802884e12L*k17)+(Power(k-1,18)*Power(k,2)*(11+k*(10+3
        *k))*LOG20(-ik))/120120.0L-(k2*Power(1+k,4)*(11+k*(-144+k*(864+k
        *(-3136+k*(7644+k*(-13104+k*(16016+k*(-13728+k*(7722+k*(-2288+k3
        *(364+k*(336+k*(144+k*(32+3*k))))))))))))))*LOG20(ik))/120120.0L;
    }
    break;

  case 3:
    switch (s) {
    case 0:
      return (-1648538892+k*(25604731152+k*(-169087318400+k*(
        576085910400+k*(-724222295797+2*k*(-1234957149426+k*(
        7894807462971+k*(-22703043707124+17*k*(2556822073140+k*(
        -3560707952136+k*(3709495731776+k*(-2888523460856+k*(
        1627925595516+k*(-591001636120+k*(51090766408+k*(124798836456+13
        *k*(-12898584906+k*(15323939736+k*(-19871425800+k*(26461209104+k
        *(-35172077095+8*k*(5652966735+2*k*(-3333299879+19*k*(182077480+
        3*k*(-51315936+17*k*(2154581+6*k*(-164402+72501*k)))))))))))))))
        ))))))))))))/(5.4019513548e12L*k19)+(Power(k-1,18)*(-1+10*k2)*(
        299+20*k*(7+k))*LOG20(-ik))/92820.0L-(Power(1+k,2)*(-1+10*k2)*(
        299+k*(-4644+k*(33660+k*(-150960+k*(468180+k*(-1062432+k*(
        1819272+k*(-2386800+k*(2406690+k*(-1847560+k*(1050192+k*(-413712
        +k*(92820+k2*(-6120+k*(-816+k*(459+20*k*(9+k)))))))))))))))))
        *LOG20(ik))/92820.0L;
    case 1: case 2: case 4: case 8:
      return (-3480248772+k*(56965564656+k*(-367879768256+k*(
        981993444432+k*(1311839423895+k*(-21194455622340+k*(
        90954902631450+k*(-243535587603660+k*(465770121728175+34*k*(
        -19632653434350+k*(21451617712242+k*(-17842652297216+k*(
        11052705633996+k*(-4805130597072+k*(1227292648840+k*(
        -21110722416+k*(-112321395990+13*k*(-1583991696+k*(1740594240+k2
        *(-2709822649+k*(8925040392+k*(-18933730032+19*k*(1643775016+9*k
        *(-224784525+136*k*(1677533+21*k*(-48265+29744*k))))))))))))))))
        ))))))))))/(6.48234162576e13L*k20)-(Power(k-1,18)*(-563-1800*k
        +6375*k2+3650*k3+600*k4)*LOG20(-ik))/556920.0L-((-35+2010*k
        -21840*k2+111820*k3-319725*k4+481950*k5-85680*k6-1187280*k7
        +2467890*k8-2011100*k9-583440*k10+3155880*k11-3264170*k12
        +1135260*k13+930240*k14-1405968*k15+755361*k16-560934*k17-135188
        *k18+47700*k19+32475*k20+7150*k21+600*k22)*LOG20(ik))/556920.0L;
    case 3: case 6: case 9: case 12:
      return (1024287264+k*(-15748416684+k*(92577789304+k*(-187909746024
        +k*(-705795690600+k*(6465264859053+k*(-24371722611522+k*(
        59498951040672+17*k*(-6125525472264+k*(7968142232985+2*k*(
        -3881589750885+k*(2788439731528+k*(-1402865716608+k*(
        426399310500+k*(-36218270600+k*(-7686523392+13*k*(-4052722800+k
        *(4453957494+k*(-6155138388+k*(7877085540+k*(-9468091248+k*(
        11249556141+k*(-12292001436+7*k*(1689523994+57*k*(-23532080+3*k
        *(5221179+476*k*(-4719+2216*k)))))))))))))))))))))))))))
        /(3.24117081288e13L*k20)+(Power(k-1,18)*(88+5*k*(11+2*k))*(-1+5
        *k*(-1+3*k))*LOG20(-ik))/278460.0L-(Power(1+k,2)*(-1+5*k*(1+3*k)
        )*(-88+k*(1353+k*(-9690+k*(42840+k*(-130560+k*(289884+k*(-482664
        +k*(609960+k*(-583440+k*(413270+k*(-204204+k*(58344+k2*(-7140+k
        *(-2040+k*(408+k*(408+5*k*(21+2*k)))))))))))))))))*LOG20(ik))
        /278460.0L;
    case 5: case 10:
      return (9094837752+k*(-141388399152+k*(1042583461920+k*(
        -4864780000080+k*(16168075400477+k*(-40812194070648+k*(
        81307215062073+2*k*(-65351625981576+17*k*(5034640752405+k*(
        -5365407505824+k*(4598012005964+k*(-3090911565024+k*(
        1553559746364+k*(-520246261360+k*(61877211972+k*(54394505424+13
        *k*(-5112759714+k*(6372195984+k*(-7793986590+k*(10120504736+k*(
        -12654917280+k*(15171981840+k*(-16586446961+57*k*(278919680+k*(
        -214701879+136*k*(1011868+9*k*(-43169+14872*k)))))))))))))))))))
        ))))))))/(3.24117081288e13L*k19)+(Power(k-1,18)*(143+5*k*(-594+k
        *(927+20*k*(31+6*k))))*LOG20(-ik))/1.11384e6L-(Power(1+k,2)*(143
        +k*(722+k*(-23565+k*(165660+k*(-588540+k*(1100376+k*(-358836+k*(
        -3973920+k*(12749490+k*(-22546420+k*(27163994+k*(-23395944+k*(
        14334060+k*(-5913960+k*(1380060+k*(-66912+k*(-40137+5*k*(8226+k
        *(5407+20*k*(65 + 6*k))))))))))))))))))))*LOG20(ik))/1.11384e6L;
    case 7: case 11: case 13: case 14:
      return (-302630328+k*(-997332336+k*(59523900436+k*(-573147006432+k
        *(3058516791330+k*(-11002941289340+k*(28946700120195+k*(
        -58106577171060+k*(91082257883430+17*k*(-6627724610100+k*(
        6469485401061+4*k*(-1230535180472+k*(707086551732+k*(
        -284918570544+k*(63733545730+k*(573318648+k*(-9194595270+13*k*(
        38939208+k*(321598785+k2*(64349357+k*(-76826196+k*(69481146+19*k
        *(1576172+3*k*(-2552935+68*k*(78221+7*k*(-6685+6864*k)))))))))))
        )))))))))))))))/(6.48234162576e13L*k20)-(Power(k-1,18)*(12+25*k
        *(-15+k*(2+k)*(9+4*k)))*LOG20(-ik))/371280.0L-(Power(1+k,2)*(
        -140+k*(1525+k*(-6270+k*(6885+k*(40800+k*(-213180+k*(499800+k*(
        -649740+k*(291720+k*(607750+k*(-1507220+3*k*(585650+k*(-424320+k
        *(183260+k*(-25160+k*(-39508+k*(4012+25*k*(503+k*(226+k*(47+4*k)
        )))))))))))))))))))*LOG20(ik))/1.11384e6L;
    case 15:
      return (828251424+k*(-12614906304+k*(90434624280+k*(-406349863920
        +k*(1286155084214+k*(-3055388331996+k* (5663275490301+2*k*(
        -4193261676252+17*k*(295171311375+k*(-285073519368+k*(
        218542898108+k*(-127849154808+k*(52337558388+k*(-11601131080+k*(
        -596033196+k*(1937868408+13*k*(70484082+k*(216733608+k*(
        -222806910+k*(311573552+k*(-358394685+k*(389034120+k*(-381708047
        +57*k*(5799040+k*(-4007213+136*k*(17446+3*k*(-1869+572*k))))))))
        )))))))))))))))))))/(1.62058540644e13L*k19)+(Power(k-1,18)*(26+5
        *k*(4+k))*(1+5*k*(-2+3*k))*LOG20(-ik))/556920.0L-(Power(1+k,2)
        *(1+5*k*(2+3*k))*(26+k*(-396+k*(2805+k*(-12240+k*(36720+k*(
        -79968+k*(129948+k*(-159120+k*(145860+k*(-97240+k*(43758+k*(
        -10608+k3*(1020+k*(816+k*(306+5*k*(12+k))))))))))))))))
        *LOG20(ik))/556920.0L;
    }
    break;

  case 4:
    switch (s) {
    case 0:
      return (1911349440+k*(-33620146560+k*(278152094220+k*(
        -1436740865560+19*k*(272871966367+k*(-728943613398+k*(
        1488927638913+k*(-2364120573276+17*k*(172136028090+k*(
        -163427070696+k*(111671515392+k*(-40367682072+k*(-21275522456+k
        *(58747655880+k*(-75296729320+k*(81124012296+k*(-83624971014+k*(
        85376235672+k*(-86450172120+k*(85986315792+k*(-82837532199+8*k*(
        9464743353+2*k*(-3990278088+19*k*(158532140+3*k*(-34074768+17*k
        *(1094951+6*k*(-70458+24167*k)))))))))))))))))))))))))))
        /(1.3158599454e12L*k19)+(Power(k-1,18)*(3+20*(-1+k)*k)*(39+2*k*(
        8+k))*LOG20(-ik))/45220.0L-((3+20*k*(1+k))*(39+k*(-686+k*(5681+k
        *(-29412+k*(106590+k*(-286824+k*(593028+k*(-961248+k*(1234506+k
        *(-1259700+k*(1016158+k*(-638248+k*(302328+k*(-100776+k*(19380
        +k2*(-969+k*(-114+k*(57+2*k*(10+k)))))))))))))))))))*LOG20(ik))
        /45220.0L;
    case 1: case 2: case 4: case 8:
      return (1815781968+k*(-35369766432+k*(324852419892+k*(
        -1868727724864+k*(7542188213560+19*k*(-1192984134342+k*(
        2764461463335+k*(-5044205625690+k*(7333458648945+17*k*(
        -502035642690+k*(466834560198+k*(-343890388944+k*(197694110464+k
        *(-87604003728+k*(31692975320+k*(-13569248304+k*(10113455970+k*(
        -8825096592+k*(5085465840+k2*(-6572900469+k*(14320168824+k*(
        -21602455704+19*k*(1365778232+3*k*(-434189955+136*k*(2496637+21
        *k*(-62055+29744*k))))))))))))))))))))))))))/(1.57903193448e13L
        *k20)-(Power(k-1,18)*(417-1880*k+955*k2+730*k3+120*k4)
        *LOG20(-ik))/271320.0L-((207-2990*k+19080*k2-68020*k3+131385*k4
        -40698*k5-542640*k6+1782960*k7-3110490*k8+3275220*k9-1478048*k10
        -1511640*k11+3821090*k12-4205460*k13+3023280*k14-1503888*k15
        +493221*k16-131214*k17-48108*k18+5300*k19+6175*k20+1430*k21+120
        *k22)*LOG20(ik))/271320.0L; 
    case 3: case 6: case 9: case 12:
      return (-1629547920+k*(28447250832+k*(-233954412692+k*(
        1204157838884+19*k*(-228743855340+k*(614890460184+k*(
        -1276787140629+k*(2098579746549+17*k*(-163711051284+k*(
        178790556945+k*(-166007644641+k*(138994922240+k*(-115183324080+k
        *(102973948000+k*(-100347363000+k*(101701042920+k*(-103441809744
        +k*(103382696898+k*(-102878741628+k*(99766919460+k*(-93163582512
        +k*(82224855099+k*(-66607712444+3*k*(15963759966+19*k*(
        -513232720+3*k*(86167155+476*k*(-65065+19944*k))))))))))))))))))
        )))))))))/(2.36854790172e13L*k20)+(Power(k-1,18)*(35+k*(19+3*k))
        *(8+5*k*(-7+6*k))*LOG20(-ik))/406980.0L-((8+5*k*(7+6*k))*(-35+k*
        (611+k*(-5016+k*(25707+k*(-92055+k*(244188+k*(-496128+k*(786828
        +k*(-982566+k*(965770+k*(-739024+k*(428298+k*(-176358+k*(42636
        +k2*(-3876+k*(-969+k*(171+k*(152+k*(35+3*k)))))))))))))))))))
        *LOG20(ik))/406980.0L;
    case 5: case 10:
      return (4135131000+k*(-72471999600+k*(597408591780+k*(
        -3074846494720+19*k*(582071521031+k*(-1550804719464+k*(
        3163560416214+k*(-5031498273648+17*k*(369438348375+2*k*(
        -179858813664+k*(132726184908+k*(-67178770128+k*(11798937476+3*k
        *(6704729520+k*(-10847000960+k*(11832017328+k*(-11914076382+k*(
        11834078256+k*(-11555037390+k*(10996894176+k*(-10027547472+k*(
        8590451152+k*(-6727005759+19*k*(245964160+k*(-143855187+136*k*(
        514228+3*k*(-55503+14872*k)))))))))))))))))))))))))))
        /(1.57903193448e13L*k19)+(Power(k-1,18)*(531+k*(-1746+k*(703+20
        *k*(31+6*k))))*LOG20(-ik))/542640.0L+((-531-k*(-7344+k*(43934+k*(
        -139428+k*(191235+k*(285684+k*(-2046528+k*(5069808+k*(-7358586+k
        *(5814000+k*(453492+k*(-8095672+k*(12168702+k*(-10829544+k*(
        6472920+k*(-2558160+k*(572679+k*(-18696+k*(798+k*(17460+k*(7903
        +20*k*(77+6*k))))))))))))))))))))))*LOG20(ik))/542640.0L;
    case 7: case 11: case 13: case 14:
      return (977728752+k*(-22522067568+k*(235569121788+k*(
        -1503637431296+k*(6607042181740+19*k*(-1121370037788+k*(
        2756190411975+k*(-5283230729700+k*(8000682677460+17*k*(
        -565759559460+k*(538091909457+4*k*(-99807985224+k*(55870817864+3
        *k*(-7234199376+k*(1393277690+k*(374115672+k*(-475690830+k*(
        216947016+k*(-71823765+k2*(105192657+k*(-157987052+k*(139199382
        +19*k*(-2562716+k*(-2486985+68*k*(91949+21*k*(-2865+2288*k))))))
        ))))))))))))))))))))/(9.47419160688e13L*k20)-(Power(k-1,18)*(588
        +5*k*(-359+k*(98+9*k*(17+4*k))))*LOG20(-ik))/1.62792e6L+((252-k*
        (3755+k*(-25200+k*(98610+k*(-238260+k*(316863+k2*(-969000+k*(
        2209320+k*(-2645370+k*(1478048+k*(755820+k*(-2519400+k*(2781030
        +k*(-1860480+k*(767448+k*(-158916+k*(-113601+k*(34428+5*k*(8090
        +k*(2852+9*k*(55+4*k)))))))))))))))))))))*LOG20(ik))/1.62792e6L;
    case 15:
      return (836215380+k*(-14494399920+k*(117987569700+k*(-598544906960
        +19*k*(111408197901+k*(-290957510844+k*(579466873524+k*(
        -894894888888+17*k*(63322342155+2*k*(-29382700344+k*(20304381348
        +k*(-9277798968+k*(1081085776+9*k*(313643400+k*(-425983720+k*(
        424168776+k*(-408192954+k*(402505272+k*(-376139790+k*(340566032
        +k*(-289604979+k*(227678824+k*(-161427633+19*k*(5304000+k*(
        -2763169+136*k*(8866+k*(-2403+572*k)))))))))))))))))))))))))))
        /(2.36854790172e13L*k19)+(Power(k-1,18)*(6+5*k*(-4+3*k))*(21+k*(
        14+3*k))*LOG20(-ik))/813960.0L-((6+5*k*(4+3*k))*(21+k*(-364+k*(
        2964+k*(-15048+k*(53295+k*(-139536+k*(279072+k*(-434112+k*(
        529074+k*(-503880+k*(369512+k*(-201552+k*(75582+k*(-15504+k3*(
        969+k*(684+k*(228+k*(40+3*k))))))))))))))))))*LOG20(ik))
        /813960.0L;
    }
    break;

  case 5:
    switch (s) {
    case 0:
      return (-1510088580+k*(26630243640+k*(-219848068440+k*(
        1125481036680+k*(-3985229876628+k*(10288325711664+19*k*(
        -1042907442192+k*(1495379239200+k*(-1530445979160+k*(
        937784192688+17*k*(2875904528+k*(-51083602416+k*(64677340200+k
        *(-45766828800+k*(14807607936+k*(10759855392+k*(-26978006363+k
        *(38858418522+k*(-52846712919+2*k*(37257674518+k*(-55008316103
        +k*(82358706758+11*k*(-10703338116+k*(13962044222+k*(
        -15379614593+76*k*(198567502+33*k*(-4267946+17*k*(184777+6*k*(
        -10304+6591*k)))))))))))))))))))))))))))))/(2.17116890991e13L
        *k19*Power(1+k,2))+(Power(k-1,20)*(493+20*k*(9+k))*LOG20(-ik))
        /746130.0L-((493+k*(-9680+k*(90090+k*(-528220+k*(2187185+k*(
        -6794172+k*(16414860+k*(-31550640+k*(48924810+k*(-61680080+k*(
        63371308+k*(-52907400+k*(35565530+k*(-18901960+k*(7674480+k*(
        -2217072+k*(373065+k2*(-14630+k*(-1540+k*(693+20*k*(11+k))))))))
        )))))))))))))*LOG20(ik))/746130.0L;
    case 1: case 2: case 4: case 8:
      return (17215009812+k*(-292606770456+k*(2313819800292+k*(
        -11252648206800+k*(37419575238045+k*(-89192760110454+k
        *(154233732803007+19*k*(-9894705738456+k*(7406895564585+k*(
        -754704561810+k*(-6184945793985+34*k*(262193585290+k*(
        -192329943190+k*(60470421880+k*(35775677300+k*(-70987112768+k
        *(76278807308+k*(-84232011864+k*(106154831244+k*(-147088669728+k
        *(208835049743+k*(-293250231602+k*(374940702137+11*k*(
        -35771141564+7*k*(4116167069+19*k*(-68096078+99*k*(-149533+136*k
        *(12463-4461*k+8112*k2))))))))))))))))))))))))))))
        /(1.0421610767568e15L*k20*Power(1+k,2))-(Power(k-1,20)*(475+194
        *k+24*k2)*LOG20(-ik))/1.790712e6L-((583-11550*k+108570*k2-643720
        *k3+2699235*k4-8505882*k5+20891640*k6-40930560*k7+64913310*k8
        -84063980*k9+89237148*k10-77597520*k11+54964910*k12-31337460*k13
        +14069880*k14-4775232*k15+1119195*k16-131670*k17-73150*k18+5775
        *k20+1430*k21+120*k22)*LOG20(ik))/8.95356e6L;
    case 3: case 6: case 9: case 12:
      return (-1467648+k*(25713324+k*(-210665728+k*(1068749500+k*(
        -3742921000+k*(9530122147+k*(-18022828096+k*(25171379857+k*(
        -24627379530+k*(13357928025+k*(3739852116+k*(-16706659333+2*k*(
        9508627583+k*(-6213525500+k*(1730802320+k*(1549173028+k*(
        -3344121352+k*(4510321452+k*(-5902213044+k*(7982600808+k*( 
        -11418726228+k*(16572696087+k*(-23116438542+k*(29270283829+77*k
        *(-413771504+k*(381951665+57*k*(-4796470+3*k*(996253+68*k*(
        -5123+2216*k)))))))))))))))))))))))))))))/(1.466593128e11L*k20
        *Power(1+k,2))+(Power(k-1,20)*(448+5*k*(43+6*k))*LOG20(-ik))
        /4.47678e6L+(0.00010007192669731369L-53*k/27132.0L+35*k2
        /1938.0L-611*k3/5814.0L+22*k4/51.0L-451*k5/340.0L+19*k6/6.0L-6
        *k7+64*k8/7.0L-203*k9/18.0L+169*k10/15.0L-299*k11/33.0L+52*k12
        /9.0L-17*k13/6.0L+k14-22*k15/105.0L+k17/68.0L+k18/306.0L-k19
        /1938.0L-2*k20/4845.0L-k21/11628.0L-k22/149226.0L)*LOG20(ik);
     case 5: case 10:
       return (-3098285190+k*(54452017620+k*(-447779552220+k*(
         2281944204540+k*(-8036690515854+k*(20611535256312+19*k*(
         -2071911278976+k*(2936841840240+k*(-2950914375780+k*(
         1729467434184+17*k*(12956354824+k*(-104445050248+k*(
         126102873820+k*(-85510189600+k*(24439302448+k*(23967514736+k*(
         -53895585744+k*(76019590416+k*(-101740480842+k*(137841831148+k
         *(-185623028758+k*(240986186968+11*k*(-26403283891+k*(
         28817270482+k*(-27036492658+57*k*(396358746+11*k*(-22546769+136
         *k*(101726+9*k*(-3463+1352*k)))))))))))))))))))))))))))))
         /(2.605402691892e14L*k19*Power(1+k,2))+(Power(k-1,20)*(2023+20
         *k*(43+6*k))*LOG20(-ik))/1.790712e7L+((-2023-k*(-39600+k*(
         367290+k*(-2145220+k*(8843835+k*(-27334692+k*(65659440+k*(
         -125349840+k*(192821310+k*(-240751280+k*(244432188+k*(
         -201048120+k*(132562430+k*(-68643960+k*(26860680+k*(-7333392+k
         *(1119195+k2*(-14630+k*(13860+k*(7623+20*k*(77+6*k)))))))))))))
         ))))))))*LOG20(ik))/1.790712e7L;
    case 7: case 11: case 13: case 14:
      return (1582989408+k*(-24818137344+k*(176119891948+k*(
        -733930036840+k*(1911450721180+k*(-2830694762896+k*(460027134853
        +19*k*(462223303586+k*(-1247983625295+k*(1840250772060+k*(
        -1640808234555+17*k*(34938775690+k*(36940951355+4*k*(
        -18211648400+k*(15382503950+k*(-6769438676+k*(-756057214+k*(
        4950697752+k*(-7236413877+k*(9016261254+k*(-10698255979+k*(
        11190397276+k*(-8947245421+11*k*(376789682+7*k*(21457468+19*k*(
        -2530966+33*k*(104077+68*k*(131+k*(293+624*k))))))))))))))))))))
        )))))))))/(1.0421610767568e15L*k20*Power(1+k,2))-(Power(k-1,20)
        *(48+k*(25+4*k))*LOG20(-ik))/1.193808e6L+((176-k*(3465+k*(-32340
        +k*(190190+k*(-790020+k*(2462229+k*(-5969040+k*(11511720+k*(
        -17907120+k*(22632610+k*(-23279256+k*(19399380+k*(-12932920+3*k
        *(2238390+k*(-852720+k*(198968+5*k2*(-4389+k2*(770+k*(308+k*(55
        +4*k))))))))))))))))))))*LOG20(ik))/1.790712e7L;
    case 15:
      return (-208288080+k*(3626663040+k*(-29501862390+k*(148430662380
        +k*(-514711425228+k*(1294761672204+19*k*(-126887528922+k*(
        173443245360+k*(-163811351760+k*(80697786288+17*k*(2160624508+k
        *(-6954911656+k*(7271715880+k*(-4214986600+k*(556972636+k*(
        1956932912+k*(-3358898928+k*(4464165552+k*(-5878666794+k*(
        7886006836+k*(-10179223156+k*(12256375036+11*k*(-1201558867+k*(
        1135829254+k*(-905100406+57*k*(11090142+11*k*(-525163+136*k*(
        2012+3*k*(-163+52*k)))))))))))))))))))))))))))))
        /(1.302701345946e14L*k19*Power(1+k,2))+(Power(k-1,20)*(136+5*k
        *(16+3*k))*LOG20(-ik))/8.95356e6L+((-136-k*(-2640+k*(24255+k*(
        -140140+k*(570570+k*(-1738044+k*(4103715+k*(-7674480+k*(11511720
        +k*(-13927760+k*(13579566+k*(-10581480+k*(6466460+k*(-2984520+k
        *(959310+k*(-170544+k3*(7315+k*(4620+k*(1386+5*k*(44+3*k))))))))
        ))))))))))))*LOG20(ik))/8.95356e6L;
    }  
    break;

  case 6:
    switch (s) {
    case 0:
      return (-10336+182648*k-1427524*k2+6276264*k3-15489584*k4+11900272
        *k5+57584599*k6-228225588*k7+403488036*k8+95596004*k9-1043924851
        *k10+4814188880*k11-4811068194*k12+10339640248*k13-3813464768
        *k14+7659893320*k15-420065376*k16+1988713584*k17)/(3.6038079e9L
        *k5*Power(1+k,4));
    case 1: case 2: case 4: case 8:
      return (103360-1788536*k+13877304*k2-62339272*k3+170543048*k4
        -245983608*k5-57866199*k6+1121847948*k7-2471506344*k8+2532485932
        *k9+1831605598*k10-5217050532*k11+16713822896*k12-2473823788*k13
        +18505736451*k14+13270392984*k15+8450695512*k16+8566766208*k17)
        /(8.64913896e10L*k6*Power(1 + k,4));
    case 3: case 6: case 9: case 12:
      return (-31620+555560*k-4392392*k2+20196136*k3-56911648*k4
        +85680408*k5+16321224*k6-408172929*k7+981765144*k8-991287252*k9
        -422348908*k10+4040664610*k11-6932409116*k12+10432285076*k13
        -7284854676*k14+7135263021*k15-1820313012*k16+1382863776*k17)
        /(4.32456948e10L*k6*Power(1+k,4));
    case 5: case 10:
      return (56202-994296*k+8195938*k2-41592608*k3+144097338*k4
        -354040504*k5+605932262*k6-603992304*k7-99716517*k8+2003546652
        *k9-4301152743*k10+7614546520*k11-7273215227*k12+8929383564*k13
        -4160754349*k14+4203031040*k15-851006808*k16+611911872*k17)
        /(4.32456948e10L*k5*Power(1+k,4));
    case 7: case 11: case 13: case 14:
      return (4845-103632*k+998138*k2-5741784*k3+21891036*k4-57374456*k5
        +101128972*k6-100763824*k7-14216253*k8+256021464*k9-455964489
        *k10+574003916*k11-183141663*k12+263280264*k13+418764427*k14
        +160852708*k15+240669884*k16+109830336*k17)/(4.32456948e10L*k6
        *Power(1+k,4));
    case 15:
      return (5814-101592*k+821406*k2-4048176*k3+13419766*k4-30800328*k5
        +46870544*k6-34290168*k7-37300899*k8+174210924*k9-308744571*k10
        +414074520*k11-350308929*k12+310516468*k13-129045663*k14
        +93308240*k15-16224936*k16+7845024*k17)/(2.16228474e10L*k5
        *Power(1+k,4));
    }
    break;

  case 7:
    switch (s) {
    case 0:
      return (-94962+1680892*k-13260442*k2+59739496*k3-158187958*k4
        +188656548*k5+241240302*k6-1421582432*k7+2411193893*k8-356916474
        *k9-3989251069*k10+18158657180*k11+11873717815*k12+75835847710
        *k13+122508462585*k14+333133796760*k15+243083407400*k16
        +632878397840*k17+239176672480*k18+438554076880*k19+98319546240
        *k20+86177588640*k21)/(1.56165009e11L*k7*Power(1+k,6));
    case 1: case 2: case 4: case 8:
      return (664734-11590532*k+90298526*k2-404008400*k3+1079444036*k4
        -1411701216*k5-939988140*k6+7770266296*k7-13874847529*k8
        +5513593914*k9+24082479029*k10-44325650260*k11+48439885750*k12
        +102963921520*k13+106306307970*k14+419529266700*k15+743613310385
        *k16+754376463410*k17+1094130960475*k18+748966521400*k19
        +468270374040*k20+185613267840*k21)/(1.873980108e12L*k8
        *Power(1+k,6));
    case 3: case 6: case 9: case 12:
      return (-614992+10830190*k-85414120*k2+388294450*k3-1062159660*k4
        +1463996208*k5+725259032*k6-7598562080*k7+14546372184*k8
        -7434599501*k9-21840865778*k10+54025436807*k11-33349208740*k12
        +5813352360*k13+120280932480*k14+151486906060*k15+29699470360
        *k16+630111939095*k17-86848475410*k18+506496335055*k19
        -30481038420*k20+103714783200*k21)/(2.810970162e12L*k8
        *Power(1+k,6));
    case 5: case 10:
      return (257754-4567220*k+36837504*k2-176271436*k3+540021116*k4
        -1025012076*k5+835801260*k6+1351544404*k7-5504187696*k8
        +8097121404*k9-2877695101*k10-5696482890*k11+23621193050*k12
        +19067541190*k13+10542936420*k14+154383989970*k15-15909205570
        *k16+210004095090*k17-8752264405*k18+93487130320*k19-2917814040
        *k20+13258090560*k21)/(9.36990054e11L*k7*Power(1+k,6));
    case 7: case 11: case 13: case 14: 
      return (94962-1927664*k+17685219*k2-95889010*k3+333553379*k4
        -733618272*k5+807644058*k6+559561324*k7-3831777831*k8+6482382750
        *k9-3332389615*k10-6425822100*k11+17935849460*k12-6215482280*k13
        +10230857700*k14+54693769080*k15+24647340155*k16+77063110230*k17
        +60316498925*k18+40531881100*k19+27934971540*k20+7138971840*k21)
        /(2.810970162e12L*k8*Power(1+k,6));
    case 15:
      return (81396-1426368*k+11327610*k2-52986212*k3+156499314*k4
        -276064224*k5+160852062*k6+539326020*k7-1661300025*k8+1994035854
        *k9-49815576*k10-3059499150*k11+5382999955*k12+1514020080*k13
        -4482407325*k14+22094976450*k15-11339922990*k16+21668037570*k17
        -4105905605*k18+6368733120*k19-510818040*k20+509926560*k21)
        /(1.405485081e12L*k7*Power(1+k,6));
    }
    break;
  }

  fprintf(stderr, "***ERROR in K2\n");
  abort();
  
  return 0.0L;
}


/*----------------------------------------------------------------------
  K3
  
  k1 = k2 = k, k3 = k4 = k+1
----------------------------------------------------------------------*/
static long double K3(int ik, int l, int s)
{
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;

  switch (l) {
  case 0:
    switch (s) {
    case 0:
      return 192277/1074427200.0L+106998337*k/4.88864376e10L+11570527
        *k2/9.585576e8L+68748571*k3/1.7459442e9L+1653011*k4/1.979208e7L
        +1412979*k5/1.19119e7L+2795237*k6/2.52252e7L+98752*k7
        /1.576575e6L+81*k8/4900.0L;
    case 1: case 2:
      return -4141/128931264.0L-19842989*k/4.88864376e10L-11320189*k2
        /4.88864376e9L-286817347*k3/3.66648282e10L-17705377*k4
        /1.02918816e9L-4350551*k5/1.7153136e8L-7963*k6/323400.0L-9131*k7
        /630630.0L-39*k8/9800.0L;
    case 3:
      return 1103/189604800.0L+3707573*k/4.88864376e10L+18722131*k2
        /4.19026608e10L+228865501*k3/1.466593128e11L+457027*k4
        /1.2864852e8L+13969159*k5/2.5729704e9L+3137*k6/573300.0L+4689*k7
        /1.4014e6L+169*k8/176400.0L;
    case 4: case 8:
      return 1667/32232816.0L+391949*k/6.267492e8L+6162249*k2
        /1.8106088e9L+400731467*k3/3.66648282e10L+23416279*k4
        /1.02918816e9L+1806379*k5/5.717712e7L+2173211*k6/7.56756e7L
        +24758*k7/1.576575e6L+39*k8/9800.0L;
    case 5: case 10:
      return -19487/2148854400.0L-425449*k/3.7604952e9L-5920261*k2
        /9.3117024e9L-1442719*k3/6.837264e8L-9352099*k4/2.05837632e9L
        -84869*k5/1.29948e7L-34481*k6/5.6056e6L-88009*k7/2.52252e7L-9*k8
        /9800.0L;
    case 6: case 9:
      return -7459/805820400.0L-19193*k/1.659042e8L-17372353*k2
        /2.66653296e10L-105967387*k3/4.88864376e10L-4817719*k4
        /1.02918816e9L-263209*k5/3.89844e7L-241471*k6/3.78378e7L-137407
        *k7/3.78378e7L-169*k8/176400.0L;
    case 7: case 11:
      return 10559/6446563200.0L+16693*k/7.916832e8L+4787351*k2
        /3.910915008e10L+246516799*k3/5.866372512e11L+1929107*k4
        /2.05837632e9L+1598033*k5/1.1435424e9L+379*k6/277200.0L+2713*k7
        /3.36336e6L+13*k8/58800.0L;
    case 12:
      return 287/19186200.0L+14411*k/8.058204e7L+1969973*k2/2.0511792e9L
        +445256501*k3/1.466593128e11L+1063789*k4/1.7153136e8L+21655709
        *k5/2.5729704e9L+112627*k6/1.513512e7L+1013*k7/257400.0L+169*k8
        /176400.0L;
    case 13: case 14:
      return -31/11850300.0L-1261*k/3.907008e7L-3117*k2/1.74097e7L
        -343275059*k3/5.866372512e11L-254627*k4/2.05837632e8L-5965909*k5
        /3.4306272e9L-6893*k6/4.32432e6L-4903*k7/5.6056e6L-13*k8/58800.0L;
    case 15:
      return 31/67151700.0L+1261*k/2.1488544e8L+9351*k2/2.785552e8L
        +11073389*k3/9.77728752e10L+254627*k4/1.02918816e9L+205721*k5
        /5.717712e8L+6893*k6/2.018016e7L+4903*k7/2.52252e7L+k8/19600.0L; 
    }
    break;

  case 1:
    switch (s) {
    case 0:
      return 227061/1810608800.0L+1404307*k/9.053044e8L+481021*k2
        /5.54268e7L+5595679*k3/1.939938e8L+25625*k4/408408.0L+942727*k5
        /1.02102e7L+15599*k6/171600.0L+29254*k7/525525.0L+81*k8/4900.0L;
    case 1: case 2:
      return -4825/217273056.0L-154289*k/5.4318264e8L-3367*k2
        /2.046528e6L-6597793*k3/1.1639628e9L-69697*k4/5.44544e6L-22849
        *k5/1.16688e6-53*k6/2640.0L-1348*k7/105105.0L-39*k8/9800.0L;
    case 3:
      return 7281/1810608800.0L+22753*k/4.288284e8L+158309*k2
        /4.988412e8L+2631467*k3/2.3279256e9L+3083*k4/1.16688e6L+6667*k5
        /1.5912e6L+487*k6/109200.0L+12457*k7/4.2042e6L+169*k8/176400.0L;
    case 4: case 8:
      return 953/26114550.0L+606443*k/1.3579566e9L+547031*k2/2.217072e8L
        +447823*k3/5.54268e7L+845987*k4/4.900896e7L+3044227*k5
        /1.225224e8L+21431*k6/900900.0L+14767*k7/1.05105e6L+39*k8
        /9800.0L;
    case 5: case 10:
      return -22947/3621217600.0L-72503*k/9.053044e8L-67607*k2
        /1.478048e8L-218669*k3/1.410864e8L-7471*k4/2.178176e6L-34843*k5
        /6.8068e6L-4073*k6/800800.0L-26233*k7/8.4084e6L-9*k8/9800.0L;
    case 6: case 9:
      return -26161/4073869800.0L-382*k/4.700619e6L-84329*k2/1.813968e8L
        -524911*k3/3.325608e8L-34325*k4/9.801792e6L-5099*k5/972400.0L
        -943*k6/180180.0L-40739*k7/1.26126e7L-169*k8/176400.0L;
    case 7: case 11:
      return 4129/3621217600.0L+484493*k/3.25909584e10L+23329*k2
        /2.6604864e8L+2860973*k3/9.3117024e9L+4601*k4/6.534528e6L+8479
        *k5/7.7792e6L+13*k6/11550.0L+2419*k7/3.36336e6L+13*k8/58800.0L;
    case 12:
      return 743/69837768.0L+196817*k/1.527701175e9L+53957*k2/7.67448e7L
        +2263483*k3/9.976824e8L+6251*k4/1.31274e6L+492887*k5/7.351344e7L
        +5623*k6/900900.0L+6401*k7/1.8018e6L+169*k8/176400.0L; 
    case 13: case 14:
      return -1/544635.0L-22679*k/9.876048e8L-2693*k2/2.078505e7L-575639
        *k3/1.3302432e9L-1025*k4/1.089088e6L-44921*k5/3.267264e7L-1369
        *k6/1.0296e6L-4413*k7/5.6056e6L-13*k8/58800.0L;
    case 15:
      return 1/3086265.0L+22679*k/5.4318264e9L+2693*k2/1.108536e8L+18569
        *k3/2.217072e8L+205*k4/1.089088e6L+1549*k5/5.44544e6L+1369*k6
        /4.8048e6L+1471*k7/8.4084e6L+k8/19600.0L;
    }
    break;
 
  case 2:
    switch (s) {
    case 0:
      return (22221108+k*(-5157936+k*(-63281820+k*(-37508408+k*(7740341
        +2*k*(5274108+k*(2453906+k*(9632056+k*(54500403+k*(182147506+19
        *k*(21199499+k*(31794096+17*k*(1904125+4*k*(307088+104247*k)))))
        )))))))))/(1.62954792e10L*k6)-(k14*Power(1+k,2)*(-1+2*k)*(3+2*k)
        *(-91+k*(-14+k*(18+k*(8+k))))*LOG20(ik))/10010.0L;
    case 1: case 2:
      return (-44442216-k*(-1199520+k*(-133880712+k*(-96589136+k*(
        11963868+k*(25937856+k*(10445400+k*(21502824+k*(123068319+k*(
        426138556+19*k*(51529413+2*k*(40189425+68*k*(627328+3*k*(140900
        +50193*k))))))))))))))/(9.77728752e10L*k6)+(k14*Power(1+k,2)*(-1
        +2*k)*(-1092+k*(-1120+k*(48+k*(288+k*(103+12*k)))))*LOG20(ik))
        /120120.0L;
    case 3:
      return (7407036+k*(1319472+k*(-23343040+k*(-19893496+k*(994833+k*(
        5169906+k*(2086917+k*(2153340+k*(11867867+k*(42480851+19*k*(
        5317293+k*(8597111+34*k*(278630+3*k*(65082+24167*k))))))))))))))
        /(4.88864376e10L*k6)-(k14*Power(1+k,3)*(-1+2*k)*(-182+k*(-42+k*(
        36+k*(20+3*k))))*LOG20(ik))/60060.0L;
    case 4: case 8:
      return (-82047168+k*(-123636240+k*(-4900896+k*(71780982+k*(
        42217056+k*(13227872+k*(33956664+k*(186708729+k*(615952964+57*k
        *(23475721+2*k*(17240277+34*k*(502235+2*k*(156344+50193*k)))))))
        ))))))/(9.77728752e10L*k5)-(k15*Power(1+k,2)*(4+k)*(3+2*k)*(-168
        +k*(60+k*(53+12*k)))*LOG20(ik))/120120.0L;
    case 5: case 10:
      return (18232704-k*(-32738328+k*(-6618304+k*(17066042+k*(11764620
        +k*(3464412+k*(4232816+k*(22943859+k*(78375374+19*k*(9289645+4*k
        *(3540258+17*k*(214309+2*k*(69389+23166*k)))))))))))))
        /(6.51819168e10*k5)+(k15*Power(1+k,2)*(-448+k*(-380+k*(200+k*(
        229+8*k*(9+k)))))*LOG20(ik))/80080.0L;
    case 6: case 9:
      return (27349056-k*(-44630712+k*(-3895584+k*(27545546+k*(17838660
        +k*(5226494+k*(6379184+k*(34845635+k*(119045802+19*k*(14147265+4
        *k*(5404533+17*k*(328626+k*(214054+72501*k))))))))))))) 
        /(9.77728752e10L*k5)+(k15*Power(1+k,2)*(-672+k*(-460+k*(344+k*(
        349+12*k*(9+k)))))*LOG20(ik))/120120.0L;
    case 7: case 11:
      return (-18232704+k*(-35017416+k*(-8614648+k*(19160442+k*(14649390
        +k*(4438896+k*(2638692+k*(13168947+k*(46470989+19*k*(5704035+k*(
        9014553+68*k*(141778+9*k*(10625+3718*k)))))))))))))
        /(1.955457504e11L*k5)-(k15*Power(1+k,3)*(-448+k*(12+k*(192+k*(85
        +12*k))))*LOG20(ik))/240240.0L;
    case 12:
      return (15383844+k*(36551088+k*(32142852+k*(13309506+k*(3187737+k
        *(5043556+k*(26728520+k*(87093079+57*k*(3262176+k*(4694983+34*k
        *(133087+k*(79954+24167*k))))))))))))/(4.88864376e10L*k4)-(k16
        *Power(1+k,2)*(3+2*k)*(126+k*(96+k*(28+3*k)))*LOG20(ik))/60060.0L;
    case 13: case 14:
      return (-20511792-k*(53130168+k*(50928192+k*(22788990+k*(5394508+k
        *(3860556+k*(19549824+k*(65935729+57*k*(2561150+k*(3826869+68*k
        *(56413+3*k*(11769+3718*k))))))))))))/(1.955457504e11L*k4)+(k16
        *Power(1+k,2)*(504+k*(828+k*(464+k*(119+12*k))))*LOG20(ik))
        /240240.0L;
    case 15:
      return (1139544+k*(3195864+k*(3313912+k*(1604610+k*(398062+k*(
        141652+k*(610932+k*(2126959+19*k*(256115+k*(395883+17*k*(24177+4
        *k*(3923+1287*k))))))))))))/(3.25909584e10L*k4)-(k16
        *Power(1+k,3)*(28+k*(24+k*(8+k)))*LOG20(ik))/40040.0L;
    }
    break; 

  case 3:
    switch (s) {
    case 0:
      return (80582040+k*(87085152+k*(-1532496042+k*(-1832464788+k*(
        90380598+k*(257629680+k*(177379745+2*k*(14971055+k*(300185325+2
        *k*(422947870+19*k*(52995860+3*k*(25930074+17*k*(1628453+4*k*(
        263128+104247*k))))))))))))))/(9.77728752e10L*k6)-(k14*(-1+2*k*(
        1+k)*(-1+10*k*(1+k)))*(-3060+k*(-408+k*(459+20*k*(9+k))))
        *LOG20(ik))/185640.0L;
    case 1: case 2:
      return (-80582040-k*(130062240+k*(-1566990810+k*(-2253917820+k*(
        -34008282+k*(350713440+k*(186849625+k*(16664010+k*(345621870+k*(
        975811600+57*k*(42846145+6*k*(10850615+34*k*(356525+2*k*(120040
        +50193*k))))))))))))))/(2.933186256e11L*k6)+(k14*(3060+k*(8160+k
        *(-54315+k*(-148230+k*(-91198+25*k*(480+k*(837+242*k+24*k2))))))
        )*LOG20(ik))/556920.0L;
    case 3:
      return (26860680+k*(57679776+k*(-532037814+k*(-893675244+k*(
        -69111504+k*(147096180+k*(70820665+k*(7743450+k*(66739035+2*k*(
        97215425+57*k*(4421470+k*(6961019+102*k*(79167+k*(55422+24167*k)
        )))))))))))))/(2.933186256e11L*k6)-(k14*(1+k)*(-1+2*k*(-1+5*k*(2
        +3*k)))*(-1020+k*(-204+k*(153+5*k*(15+2*k))))*LOG20(ik))
        /556920.0L;
    case 4: case 8:
      return (-386793792+k*(-1639632456+k*(-900340056+k*(919123128+k*(
        614053440+k*(216337380+k*(48669060+k*(499741440+k*(1451851280+57
        *k*(58658035+2*k*(42534687+68*k*(646357+3*k*(135484+50193*k)))))
        ))))))))/(2.933186256e11L*k5)+(k15*(14688+k*(77724-k*(-99684+k*(
        -948+25*k*(2400+k*(1299+286*k+24*k2))))))*LOG20(ik))/556920.0L;
    case 5: case 10:
      return (85954176-k*(-418531806+k*(-330448272+k*(198683562+k*(
        176784300+k*(56451815+k*(8847740+k*(61248165+k*(183967850+19*k*(
        23132275+12*k*(2903298+17*k*(183459+2*k*(60079+23166*k))))))))))
        )))/(1.955457504e11L*k5)+(k15*(-3264+k*(-19329+k*(-29268+k*(
        -5627+5*k*(2930+k*(1833+40*k*(11+k)))))))*LOG20(ik))/371280.0L;
    case 6: case 9:
      return (128931264-k*(-562660560+k*(-331497936+k*(357663306+k*(
        267844500+k*(88264190+3*k*(3257320+k*(32058345+k*(91696550+19*k
        *(11755475+4*k*(4399317+17*k*(280640+k*(183674+72501*k))))))))))
        )))/(2.933186256e11L*k5)+(k15*(-4896+k*(-26520+k*(-35064+k*(379
        +15*k*(1630+k*(933+20*k*(11+k)))))))*LOG20(ik))/556920.0L;
    case 7: case 11:
      return (-85954176+k*(-429276078+k*(-355400430+k*(224915922+k*(
        223527150+k*(72601235+k*(8970555+k*(35780085+k*(108187475+57*k*(
        4737275+9*k*(818173+68*k*(13468+k*(9155+3718*k)))))))))))))/(
        5.866372512e11L*k5)-(k15*(1+k)*(-3264+k*(-16473+k*(-14172+5*k*(
        1709+k*(1786+545*k+60*k2)))))*LOG20(ik))/1.11384e6L; 
    case 12: 
      return (145047672+k*(717604272+k*(913832766+k*(410436180+k*(
        99228990+k*(22656480+k*(138320325+2*k*(208838905+57*k*(8146435+k
        *(11699565+34*k*(343922+3*k*(70294+24167*k))))))))))))/(
        2.933186256e11L*k4)-(k16*(612+k*(408+5*k*(21+2*k)))*(9+2*k*(24+5
        *k*(7+3*k)))*LOG20(ik))/556920.0L;
    case 13: case 14:
      return (-96698448-k*(514664766+k*(723841272+k*(355197150+k*(
        86535540+k*(12154065+k*(50936160+k*(156612775+57*k*(6371750+k*(
        9479955+68*k*(145243+9*k*(10299+3718*k))))))))))))/(
        5.866372512e11L*k4)+(k16*(3672+k*(23409+k*(48048+5*k*(8475+k*(
        3486+715*k+60*k2)))))*LOG20(ik))/1.11384e6L;
    case 15:
      return (5372136+k*(30607038+k*(46879602+k*(25121250+k*(6417810+k*(
        767305+k*(1591755+k*(5052025+19*k*(637175+3*k*(326895+17*k*(
        20749+4*k*(3433+1287*k))))))))))))/(9.77728752e10L*k4)-(k16*(1+k
        )*(1+5*k*(1+k))*(204+k*(153+5*k*(9+k)))*LOG20(ik))/185640.0L;
    }  
    break;

  case 4:
    switch (s) {
    case 0:
      return (-52378326+k*(-404972568+k*(-721062342+k*(-363447084+k*(
        59572149+k*(98552496+k*(215017813+2*k*(405482735+k*(1118173815+2
        *k*(1099030768+19*k*(82534738+9*k*(9468366+17*k*(404041+4*k*(
        47518+11583*k))))))))))))))/(3.25909584e10L*k4*Power(1+k,2))-(
        k16*(3+20*k*(1+k))*(-969+2*k*(-57+k*(57+2*k*(10+k))))*LOG20(ik))
        /90440.0L;
    case 1: case 2:
      return (52378326-k*(-450161712+k*(-866575710+k*(-476143668+k*(
        57894291+k*(104205780+k*(132466085+k*(464379450+k*(1317390285+k
        *(2677212020+57*k*(69299387+2*k*(37041867+68*k*(409220+9*k*(
        22174+5577*k))))))))))))))/(9.77728752e10L*k4*Power(1+k,2))+(k16
        *(-2907+k*(-22230+k*(-27018+5*k*(-144+k*(773+242*k+24*k2)))))
        *LOG20(ik))/271320.0L;
    case 3: 
      return (-52378326+k*(-442972530+k*(-574251678+k*(-29723694+k*(
        77572677+k*(41797707+k*(52246645+k*(223989685+6*k*(95955205+19*k
        *(9647617+k*(12693943+102*k*(117130+k*(69929+24167*k))))))))))))
        )/(2.933186256e11L*k4*(1+k))-(k16*(3+5*k*(5+6*k))*(-969+k*(-171
        +2*k*(57+k*(25+3*k))))*LOG20(ik))/813960.0L;
    case 4: case 8:
      return (-69837768+k*(-366053688+k*(-473489016+k*(-23183160+k*(
        262293486+k*(157508260+k*(198354160+k*(694376100+k*(1888512025
        +k*(3657447100+57*k*(89958761+2*k*(45461593+34*k*(944959+18*k*(
        23890+5577*k))))))))))))))/(9.77728752e10L*k4*Power(1+k,2))+(k16
        *(3876-k*(-16644+k*(-10488+5*k*(1744+k*(1235+286*k+24*k2)))))
        *LOG20(ik))/271320.0L;
    case 5: case 10:
      return (17459442-k*(-102810708+k*(-154750596+k*(-33969936+k*(
        65153907+k*(40790386+k*(28950608+k*(86947520+k*(242388670+k*(
        484634012+19*k*(36908873+12*k*(3211658+17*k*(137819+6*k*(10789
        +2574*k))))))))))))))/(6.51819168e10L*k4*Power(1+k,2))+(k16*(
        -969+k*(-4788+k*(-4047+k*(1970+k*(1737+40*k*(11+k))))))
        *LOG20(ik))/180880.0L;
    case 6: case 9:
      return (69837768-k*(-370161792+k*(-472245774+k*(26882856+k*(
        338164008+k*(187527704+k*(131036637+k*(393708480+k*(1099576130+9
        *k*(244498092+19*k*(18658173+4*k*(4881849+17*k*(210397+k*(99432
        +24167*k))))))))))))))/(2.933186256e11L*k4*Power(1+k,2))+(k16*(
        -3876+k*(-16872+k*(-10203+k*(11230+k*(8053+180*k*(11+k))))))
        *LOG20(ik))/813960.0L;
    case 7: case 11:
      return (-17459442+k*(-86378292+k*(-69345276+k*(44108064+k*(
        37921611+k*(12751466+k*(10402260+k*(40989080+k*(104487610+57*k*(
        3424912+k*(4414651+68*k*(59273+3*k*(11403+3718*k)))))))))))))/(
        1.955457504e11L*k4*(1+k))+(k16*(969-k*(-4845+k*(-4047+k*(2455+k
        *(2219+605*k+60*k2)))))*LOG20(ik))/542640.0L;
    case 12:
      return (73945872+k*(437386950+k*(827530704+k*(666441958+k*(
        263483584+k*(191836827+k*(597673780+k*(1602440455+6*k*(510008468
        +19*k*(37005929+k*(36651125+102*k*(247540+k*(108968+24167*k)))))
        ))))))))/(2.933186256e11L*k3*Power(1+k,2))-(k17*(8+5*k*(7+6*k))
        *(513+2*k*(152+k*(35+3*k)))*LOG20(ik))/813960.0L;
    case 13: case 14:
      return (-17459442-k*(109657548+k*(219399180+k*(185578120+k*(
        73375757+k*(31129826+k*(74251540+k*(203376140+k*(400771634+57*k
        *(10006904+k*(10240313+68*k*(107213+3*k*(16265+3718*k)))))))))))
        ))/(1.955457504e11L*k3*Power(1+k,2))+(k17*(969+k*(5168+k*(7275+k
        *(3374+715*k+60*k2))))*LOG20(ik))/542640.0L;
    case 15:
      return (1027026+k*(5765760+k*(8492484+k*(4160884+k*(1042769+k*(
        540540+k*(1877920+k*(4776090+19*k*(458644+3*k*(193588+17*k*(
        10079+12*k*(470+143*k))))))))))))/(3.25909584e10L*k3*(1+k))-(k17
        *(1+5*k*(1+k))*(57+k*(38+k*(10+k)))*LOG20(ik))/90440.0L;
    }
    break;

  case 5:
    switch (s) {
    case 0:
      return (-263603340+k*(-804683880+k*(-792783589+k*(223734990+k*(
        3724483651+2*k*(8733281274+k*(28629114336+11*k*(6230920914+k*(
        11181756103+k*(15215294154+19*k*(824276336+33*k*(19095958+17*k*(
        611103+4*k*(53836+9477*k))))))))))))))/(5.377508136e11L*k2
        *Power(1+k,4))+(k18*(7315-k*(-770+k*(693+20*k*(11+k))))
        *LOG20(ik))/746130.0L;
    case 1: case 2:
      return (659008350-k*(-2025583560+k*(-2064676120+k*(-157990408+k*(
        4169995593+k*(19726718920+k*(66246497305+11*k*(14826986250+k*(
        27382131343+k*(38360014372+57*k*(713352265+22*k*(25528015+34*k*(
        420415+6*k*(25372+4563*k))))))))))))))/(3.2265048816e12L*k2
        *Power(1+k,4))+(k18*(-7315+k*(-924+k*(693+242*k+24*k2)))
        *LOG20(ik))/1.790712e6L;
    case 3:
      return (-131801670+k*(-277477200+k*(-150982637+k*(53699901+k*(
        357870108+k*(1577772941+11*k*(458097225+k*(1060428264+k*(
        1817466455+57*k*(40716184+11*k*(3479313 + 34*k*(68260+3*k*(9873
        +2197*k)))))))))))))/(1.6132524408e12L*k2*Power(1+k,3))+(k18*(
        7315-k*(-1155+k*(693+5*k*(55+6*k))))*LOG20(ik))/4.47678e6L;
    case 4: case 8:
      return (-395405010+k*(-1082161080+k*(-747038370+k*(1026310116+k*(
        6802364673+k*(30139536280+k*(97496060155+11*k*(20970417390+k*(
        37134293609+k*(49759987772+57*k*(882715195+66*k*(10014575+68*k*(
        78190+k*(26776+4563*k))))))))))))))/(3.2265048816e12L*k2
        *Power(1+k,4))-(k18*(-4389+k*(924+k*(1155+286*k+24*k2)))
        *LOG20(ik))/1.790712e6L;
    case 5: case 10:
      return (131801670-k*(-374594220+k*(-306536828+k*(105212432+k*(
        880848241+k*(3733282526+k*(12306525029+11*k*(2718429996+k*(
        4947350135+11*k*(619518762+19*k*(33888569+12*k*(2172478+17*k*(
        69613+24386*k+4212*k2)))))))))))))/(2.1510032544e12L*k2
        *Power(1+k,4))+(k18*(-7315+k*(770+k*(1617+40*k*(11+k))))
        *LOG20(ik))/5.96904e6L;
    case 6: case 9:
      return (131801670-k*(-346846500+k*(-205202894+k*(287676376+k*(
        1352414973+k*(5630269258+k*(18566272807+33*k*(1368143456+k*(
        2492492905+k*(3438010622+19*k*(188426849+44*k*(3303699+17*k*(
        106334+k*(37538+6591*k))))))))))))))/(3.2265048816e12L*k2
        *Power(1+k,4))+(k18*(-7315+3*k*(770+k*(847+20*k*(11+k))))
        *LOG20(ik))/8.95356e6L;
    case 7: case 11:
      return (-131801670+k*(-235855620+k*(-44567276+k*(153898056+k*(
        435269289+k*(1759217498+11*k*(502647855+k*(1146716592+k*(
        1932701687+57*k*(42460562+11*k*(3545859+68*k*(33845+9*k*(1579
        +338*k)))))))))))))/(6.4530097632e12L*k2*Power(1+k,3))-(k18
        *(-7315+k*(1155+k*(2079+605*k+60*k2)))*LOG20(ik))/1.790712e7L;
    case 12:
      return (62432370+k*(224743428+k*(407446962+k*(1098018519+k*(
        4366060049+k*(13903129061+11*k*(2957006379+k*(5171031796+k*(
        6829879483+57*k*(119137912+11*k*(7950143+102*k*(40398+k*(13422
        +2197*k)))))))))))))/(1.6132524408e12L*k*Power(1+k,4))-(k19*(
        3465+k*(1848+5*k*(77+6*k)))*LOG20(ik))/4.47678e6L;
    case 13: case 14:
      return (-104053950-k*(369972564+k*(570350196+k*(960747840+k*(
        3236907947+k*(10405383158+11*k* (2271227997+k*(4081083328+k*(
        5540346181+57*k*(99348206+99*k*(757161+68*k*(5929+k*(2021+338*k)
        ))))))))))))/(6.4530097632e12L*k*Power(1+k,4))+(k19*(5775+k*(
        3234+715*k+60*k2))*LOG20(ik))/1.790712e7L;
    case 15:
      return (6936930+k*(17818476+k*(18085808+k*(23730876+k*(82401137+11
        *k*(23118090+k*(52017063+k*(86304904+19*k*(5584589+33*k*(152106
        +17*k*(5655+4*k*(574+117*k))))))))))))/(1.0755016272e12L*k
        *Power(1+k,3))-(k19*(385+k*(231+5*k*(11+k)))*LOG20(ik))
        /2.98452e6L;
    } 
    break;

  case 6:
    switch (s) {
    case 0:
      return (1331897+24719620*k+215784435*k2+1175671610*k3+4472587933
        *k4+12587691168*k5+27077215635*k6+45321215490*k7+59464765800*k8
        +61052582392*k9+48469381292*k10+28997038860*k11+12441569170*k12
        +3459536720*k13+476585208*k14)/(2.88304632e10L*Power(1+k,6));
    case 1: case 2:
      return (-460250-8742640*k-78140925*k2-436092602*k3-1700026192*k4
        -4904613480*k5-10818415010*k6-18572554740*k7-24997993578*k8
        -26328268728*k9-21436499675*k10-13144518630*k11-5773521640*k12
        -1639491152*k13-229466952*k14)/(5.76609264e10L*Power(1+k,6));
    case 3:
      return (41657+766687*k+6616432*k2+35503687*k3+132394925*k4
        +363044903*k5+754930403*k6+1208749843*k7+1495351348*k8
        +1417750360*k9+1006953317*k10+511323212*k11+168326282*k12
        +27621022*k13)/(2.88304632e10L*Power(1+k,5));
    case 4: case 8:
      return (2349676+43330824*k+375590691*k2+2030493574*k3+7658242626
        *k4+21347584344*k5+45430890576*k6+75132751464*k7+97256051586*k8
        +98342577704*k9+76741982091*k10+45028751814*k11+18903268796*k12
        +5130289104*k13+688400856*k14)/(1.729827792e11L*Power(1+k,6));
    case 5: case 10:
      return (-266658-5029876*k-44611149*k2-246857040*k3-953286334*k4
        -2721510612*k5-5932926084*k6-10051996176*k7-13330125180*k8
        -13805554776*k9-11028179633*k10-6617057496*k11-2835173844*k12
        -782422280*k13-105907824*k14)/(1.153218528e11L*Power(1+k,6));
    case 6: case 9:
      return (-401212-7569832*k-67158543*k2-371754940*k3-1436215792*k4
        -4102356648*k5-8949033798*k6-15174775872*k7-20145598800*k8
        -20895093728*k9-16726341497*k10-10066974972*k11-4334167748*k12
        -1205875280*k13-165726132*k14)/(1.729827792e11L*Power(1+k,6));
    case 7: case 11:
      return (47726+871588*k+7457495*k2+39638268*k3+146257050*k4
        +396338484*k5+813279172*k6+1282796480*k7+1560256612*k8
        +1451102200*k9+1008413771*k10+499615108*k11+160020660*k12
        +25496328*k13)/(1.153218528e11L*Power(1+k,5));
    case 12:
      return (346500+6351212*k+54686929*k2+293486333*k3+1097990811*k4
        +3033309420*k5+6391048893*k6+10451645556*k7+13359765027*k8
        +13317736644*k9+10224891040*k10+5887963637*k11+2417812624*k12
        +638834568*k13+82863066*k14)/(8.64913896e10L*Power(1+k,6));
    case 13: case 14:
      return (-233240-4371486*k-38501628*k2-211422311*k3-809587404*k4
        -2289835530*k5-4940626044*k6-8275362372*k7-10834833744*k8
        -11062280836*k9-8696597160*k10-5124988731*k11-2151527548*k12
        -580086036*k13-76488984*k14)/(3.459655584e11L*Power(1+k,6));
    case 15:
      return (6860+124362*k+1055496*k2+5560289*k3+20313900*k4+54443070
        *k5+110341488*k6+171633084*k7+205483656*k8+187691500*k9
        +127743840*k10+61757277*k11+19196536*k12+2941884*k13)
        /(5.76609264e10L*Power(1+k,5));
    }
    break;
 
  case 7:
    switch (s) {
    case 0:
      return (50902047+1047075594*k+10198615802*k2+62484512936*k3
        +269765250066*k4+871091019084*k5+2178675364080*k6+4310767069200
        *k7+6829235283945*k8+8706128907630*k9+8919537468930*k10
        +7284506731080*k11+4663585794480*k12+2271881875320*k13
        +798726706300*k14+182508540160*k15+20652025680*k16)
        /(1.249320072e12L*Power(1+k,8));
    case 1: case 2:
      return (-26337675-553262046*k-5504793412*k2-34462814712*k3
        -152078259111*k4-502065768006*k5-1284111794760*k6-2598681253440
        *k7-4211167549095*k8-5491471692630*k9-5754166913220*k10 
        -4804868165880*k11-3143399861715*k12-1563386129310*k13
        -560322667300*k14-130198871040*k15-14915351880*k16)
        /(3.747960216e12L*Power(1+k,8));
    case 3:
      return (14302575+292003275*k+2815877708*k2+17030333738*k3
        +72319163604*k4+228679091770*k5+556999012300*k6+1065801412320*k7
        +1618168009205*k8+1953451492685*k9+1864441390680*k10
        +1386032484850*k11+780170484550*k12+315769347840*k13+82885798740
        *k14+10772198580*k15)/(1.1243880648e13L*Power(1+k,7));
    case 4: case 8:
      return (44927736+918900956*k+8894551850*k2+54125400888*k3
        +231944189127*k4+742878177606*k5+1841414349360*k6+3607638290880
        *k7+5653359437430*k8+7120847132520*k9+7199049718410*k10
        +5793794795880*k11+3649812325545*k12+1746797887290*k13
        +602361446760*k14+134788210080*k15+14915351880*k16)
        /(3.747960216e12L*Power(1+k,8));
    case 5: case 10:
      return (-5090841-106271852*k-1050182365*k2-6525925666*k3
        -28564299543*k4-93463231632*k5-236712695010*k6-473888118780*k7
        -758822925645*k8-976547479380*k9-1008419790675*k10-828539057670
        *k11-532406598195*k12-259575535320*k13-90989351940*k14
        -20619635240*k15-2294669520*k16)/(2.498640144e12L*Power(1+k,8));
    case 6: case 9:
      return (-22965600-479499448*k-4739462913*k2-29458944078*k3
        -128982140781*k4-422185926084*k5-1069732934490*k6-2142720727500
        *k7-3433423773270*k8-4422450757560*k9-4572106738155*k10
        -3762490606530*k11-2423115721695*k12-1185266054100*k13
        -417574580880*k14-95412429960*k15-10772198580*k16)
        /(1.1243880648e13L*Power(1+k,8));
    case 7: case 11:
      return (2732625+55403188*k+530234433*k2+3180350770*k3+13382993983
        *k4+41896609380*k5+100928201570*k6+190781773100*k7+285774638325
        *k8+339875003380*k9+319083173695*k10+232940905350*k11
        +128538801995*k12+50916137300*k13+13060627740*k14+1657261320*k15
        )/(7.495920432e12L*Power(1+k,7));
    case 12:
      return (39779768+809217920*k+7787003145*k2+47084019174*k3
        +200369459250*k4+636883249014*k5+1565565502470*k6+3039227846700
        *k7+4714811891280*k8+5872869443100*k9+5864541395295*k10
        +4655431322730*k11+2887949394450*k12+1358295413130*k13
        +459054253800*k14+100287042600*k15+10772198580*k16)
        /(1.1243880648e13L*Power(1+k,8));
    case 13: case 14:
      return (-4454408-92454103*k-907978552*k2-5604390591*k3-24351852726
        *k4-79047749373*k5-198471021900*k6-393577570830*k7-623719056780
        *k8-793618267035*k9-809399603760*k10-656031783345*k11
        -415313197650*k12-199190894985*k13-68566846980*k14-15227815620
        *k15-1657261320*k16)/(7.495920432e12L*Power(1+k,8));
    case 15:
      return (131012+2638867*k+25075746*k2+149239887*k3+622691234*k4
        +1931316335*k5+4605129780*k6+8607334390*k7+12733449520*k8
        +14936577315*k9+13809547070*k10+9910253225*k11+5363890950*k12
        +2078002735*k13+519138520*k14+63740820*k15)/(1.249320072e12L
        *Power(1+k,7));
    }
    break;
  }

  fprintf(stderr, "***ERROR in K3\n");
  abort();

  return 0.0L;
}


/*---------------------------------------------------------------------- 
  K4
  
  k1 = k2 = k3 = k, k4 = k+1
---------------------------------------------------------------------- */
static long double K4(int ik, int l, int s)
{
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;

  switch (l) {
  case 0:
    switch (s) {
    case 0:
      return 655531/9311702400.0L+32521123*k/3.25909584e10L+52842707*k2
        /8.1477396e9L+88992073*k3/3.4918884e9L+2656133*k4/3.958416e7L
        +17668067*k5/1.429428e8L+54941*k6/350350.0L+262391*k7/2.1021e6L
        +117*k8/2450.0L;
    case 1:
      return -218473/16761064320.0L-12380917*k/6.51819168e10L-62079323
        *k2/4.88864376e10L-15369773*k3/2.9930472e9L-31483*k4/2.261952e6L
        -2264011*k5/8.576568e7L-72697*k6/2.1021e6L-60169*k7/2.1021e6L
        -169*k8/14700.0L;
    case 2: case 8:
      return 5413549/293318625600.0L+49885181*k/1.955457504e11L+2181629
        *k2/1.3579566e9L+68359129*k3/1.12814856e10L+15527321*k4
        /1.02918816e9L+10983013*k5/4.288284e8L+45904*k6/1.576575e6L+7127
        *k7/350350.0L+33*k8/4900.0L;
    case 3: case 9:
      return -332989/97772875200.0L-2184887*k/4.51259424e10L-184108877
        *k2/5.866372512e11L-119216131*k3/9.77728752e10L-3221513*k4
        /1.02918816e9L-14114467*k5/2.5729704e9L-253*k6/39200.0L-118211
        *k7/2.52252e7L-143*k8/88200.0L;
    case 4:
      return 84377/4190266080.0L+55217759*k/1.955457504e11L+656553*k2
        /3.6212176e8L+1030433377*k3/1.466593128e11L+18737669*k4
        /1.02918816e9L+70607*k5/2.144142e6L+257599*k6/6.3063e6L+1104
        *k7/35035.0L+169*k8/14700.0L;
    case 5:
      return -353/96996900.0L-1138261*k/2.17273056e10L-1608163*k2 
        /4.6558512e9L-2100247*k3/1.527701175e9L-7531961*k4/2.05837632e9L
        -277621*k5/4.08408e7L-36511*k6/4.2042e6L-29209*k7/4.2042e6L-13
        *k8/4900.0L;
    case 6: case 12:
      return 3727/705092850.0L+1927427*k/2.66653296e10L+52754447*k2
        /1.1732745024e11L+490404053*k3/2.933186256e11L+421357*k4
        /1.02918816e8L+5843879*k5/8.576568e8L+381349*k6/5.04504e7L+1229
        *k7/240240.0L+143*k8/88200.0L;
    case 7: case 13:
      return -21487/22562971200.0L-1307609*k/9.77728752e10L-1857313*k2
        /2.17273056e10L-6846607*k3/2.09513304e10L-847241*k4
        /1.02918816e9L-4840039*k5/3.4306272e9L-5821*k6/3.6036e6L-433*k7
        /382200.0L-11*k8/29400.0L;
    case 10:
      return 31667/6446563200.0L+21107*k/3.174444e8L+79831117*k2
        /1.955457504e11L+437759423*k3/2.933186256e11L+668587*k4
        /1.8712512e8L+4948693*k5/8.576568e8L+14821*k6/2.4024e6L+14459
        *k7/3.6036e6L+3*k8/2450.0L;
    case 11:
      return -233/257862528.0L-2460779*k/1.955457504e11L-4243549*k2
        /5.33306592e10L-175911607*k3/5.866372512e11L-217873*k4
        /2.9405376e8L-363707*k5/2.9405376e8L-34561*k6/2.52252e7L-5191
        *k7/5.6056e6L-13*k8/44100.0L;
    case 14:
      return 227/161164080.0L+851131*k/4.51259424e10L+6710047*k2
        /5.866372512e10L+80594489*k3/1.955457504e11L+499157*k4
        /5.1459408e8L+150419*k5/9.801792e7L+242119*k6/1.513512e8L
        +30497*k7/3.027024e7L+13*k8/44100.0L;
    case 15:
      return -233/920937600.0L-52357*k/1.50419808e10L-4243549*k2
        /1.955457504e11L-4290527*k3/5.33306592e10L-11467*k4/5.8810752e7L
        -363707*k5/1.1435424e9L-34561*k6/1.009008e8L-179*k7/800800.0L-k8
        /14700.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0:
      return 157461/3621217600.0L+1120073*k/1.8106088e9L+3134249*k2
        /7.759752e8L+6194117*k3/3.879876e8L+4325*k4/102102.0L+268427*k5
        /3.4034e6L+1489*k6/14300.0L+66767*k7/700700.0L+117*k8/2450.0L;
    case 1:
      return -79/9529520.0L-439239*k/3.6212176e9L-948599*k2/1.1639628e9L
        -1102897*k3/3.325608e8L-29609*k4/3.267264e6L-141943*k5
        /8.16816e6L-1013*k6/42900.0L-15553*k7/700700.0L-169*k8/14700.0L;
    case 2: case 8:
      return 1931/164600800.0L+5330863*k/3.25909584e10L+4850753*k2
        /4.6558512e9L+9290819*k3/2.3279256e9L+55337*k4/5.44544e6L+104239
        *k5/5.8344e6L+25967*k6/1.2012e6L+5937*k7/350350.0L+33*k8/4900.0L;
    case 3: case 9:
      return -1447/651819168.0L-1556143*k/4.88864376e10L-1457173*k2
        /6.9837768e9L-182471*k3/2.217072e8L-21163*k4/9.801792e6L-47939
        *k5/1.225224e7L-41*k6/8400.0L-99521*k7/2.52252e7L-143*k8/88200.0L;
    case 4:
      return 44743/3621217600.0L+5662831*k/3.25909584e10L+521219*k2
        /4.6558512e8L+1449077*k3/3.325608e8L+185527*k4/1.633632e7L+40277
        *k5/1.9448e6L+668*k6/25025.0L+500*k7/21021.0L+169*k8/14700.0L;
    case 5:
      return -16707/7242435200.0L-361943*k/1.08636528e10L-6711*k2
        /3.04304e7L-21389*k3/2.4249225e7L-77179*k4/3.267264e7L-2151*k5
        /486200.0L-1171*k6/200200.0L-7473*k7/1.4014e6L-13*k8/4900.0L;
    case 6: case 12:
      return 449/134303400.0L+265093*k/5.7513456e9L+8897*k2/3.069792e7L
        +462941*k3/4.232592e8L+11177*k4/4.08408e6L+192499*k5/4.08408e7L
        +10027*k6/1.8018e6L+3053*k7/720720.0L+143*k8/88200.0L;
    case 7: case 13:
      return -10099/16295479200.0L-190667*k/2.17273056e10L-131983*k2
        /2.3279256e9L-511943*k3/2.3279256e9L-9241*k4/1.633632e7L-163607
        *k5/1.633632e8L-139*k6/114400.0L-121*k7/127400.0L-11*k8/29400.0L;
    case 10:
      return 451/141086400.0L+59459*k/1.3579566e9L+363161*k2
        /1.3302432e9L+677563*k3/6.651216e8L+27329*k4/1.089088e7L+172541
        *k5/4.08408e7L+11573*k6/2.4024e6L+4163*k7/1.2012e6L+3*k8/2450.0L;
    case 11: 
      return -87/144848704.0L-275749*k/3.25909584e10L-19697*k2
        /3.627936e8L-1942949*k3/9.3117024e9L-17309*k4/3.267264e7L-991*k5
        /1.07712e6L-1301*k6/1.2012e6L-13543*k7/1.68168e7L-13*k8/44100.0L;
    case 14:
      return 3721/4073869800.0L+1209937*k/9.77728752e10L+2339*k2
        /3.069792e7L+372727*k3/1.3302432e9L+475*k4/700128.0L+60939*k5
        /5.44544e7L+8971*k6/7.2072e6L+8761*k7/1.009008e7L+13*k8/44100.0L;
    case 15:
      return -87/517316800.0L-5867*k/2.5069968e9L-19697*k2/1.3302432e9L
        -47389*k3/8.465184e8L-911*k4/6.534528e6L-991*k5/4.1888e6L-1301
        *k6/4.8048e6L-467*k7/2.4024e6L-k8/14700.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0:
      return (496271412+k*(-195341832+k*(-1349792248+k*(-610613248+k*(
        193629513+2*k*(53904389+k*(-2051070+k*(-14814324+k*(2512953+k*(
        209976+k*(55191171+k*(172709930+19*k*(26693133+2*k*(23546667+68
        *k*(486236+414633*k+301158*k2)))))))))))))))/(3.25909584e10L*k8)
        +(k12*Power(1+k,2)*(-1+2*k)*(6097+k*(4018+k*(-837+k*(-680+k*(23
        +2*k*(60+k*(19+2*k)))))))*LOG20(ik))/20020.0L;
    case 1: 
      return (-496271412-k*(-31817268+k*(-1521923368+k*(-936179160+k*(
        204696975+k*(148508087+k*(-2490516+k*(-29461068+k*(-3603505+k*(
        1780992+k*(65400528+k*(222819478+19*k*(34708221+16*k*(3991305+17
        *k*(336433+9*k*(33149+24167*k))))))))))))))))/(9.77728752e10L*k8
        )-(k12*Power(1+k,3)*(-1+2*k)*(12194+k*(-140+k*(-1534+k*(-236+k*(
        197+k*(91+12*k))))))*LOG20(ik))/120120.0L;
    case 2: case 8:
      return (-155547756+k*(-17453016+k*(440500872+k*(417369288+k*(
        75665233+k*(-63713664+k*(-43838676+k*(-7679364+k*(8388522+k*(
        76269435+k*(279823582+19*k*(38653677+8*k*(8523924+17*k*(630607
        +54*k*(9494+4719*k)))))))))))))))/(9.77728752e10L*k7)+(k13
        *Power(1+k,2)*(-1+2*k)*(-3822+k*(-4452+k*(-1082+k*(636+k*(486+k
        *(125+12*k))))))*LOG20(ik))/120120.0L;
    case 3: case 9:
      return (51849252-k*(-21201516+k*(156949576+k*(175311472+k*(
        37884861+k*(-23978493+k*(-18084906+k*(-4280724+k*(1491661+k*(
        15371479+k*(58813179+19*k*(8332447+2*k*(7595729+51*k*(192621+2*k
        *(80831+40898*k)))))))))))))))/(9.77728752e10L*k7)-(k13
        *Power(1+k,3)*(-1+2*k)*(-1274+k*(-588+k*(78+k*(152+51*k+6*k2))))
        *LOG20(ik))/120120.0L;
    case 4:
      return (-1471721076+k*(-2365641936+k*(-773476200+k*(254414034+k*(
        60933131+k*(-58751784+k*(-49964796+k*(28727+k*(-954576+k*(
        92528286+k*(279358462+19*k*(42746937+4*k*(18326979+68*k*(370495
        +9*k*(33760+24167*k)))))))))))))))/(9.77728752e10L*k7)+(k13
        *Power(1+k,2)*(36162+k*(23868+k*(-3818+k*(-2525+k*(786+k*(847
        +238*k+24*k2))))))*LOG20(ik))/120120.0L;
    case 5:
      return (327049128-k*(-634443264+k*(-274760528+k*(66370374+k*(
        25408110+k*(-14285040+k*(-12449304+k*(-2150171+k*(101616+k*(
        11886618+k*(39189520+57*k*(2005741+4*k*(895461+68*k*(18391+k*(
        15629+11154*k)))))))))))))))/(6.51819168e10L*k7)-(k13
        *Power(1+k,3)*(8036+k*(-60+k*(-792+k*(-23+k*(165+8*k*(8+k))))))
        *LOG20(ik))/80080.0L;
    case 6: case 12:
      return (138454596+k*(268401168+k*(135948456+k*(-35705082+k*(
        -60534411+k*(-27081474+k*(-4776408+k*(2107981+k*(21194378+k*(
        76270303+114*k*(1727188+k*(2979055+17*k*(215055+2*k*(84665+40898
        *k))))))))))))))/(9.77728752e10L*k6)+(k14*Power(1+k,2)*(-3402+k
        *(-3372+k*(-126+k*(993+2*k*(276+k*(65+6*k))))))*LOG20(ik))
        /120120.0L;
    case 7: case 13:
      return (-92303064-k*(208887840+k*(129102624+k*(-17281474+k*(
        -47456262+k*(-22915116+k*(-4873008+k*(602601+k*(8357034+k*(
        31356716+19*k*(4358262+k*(7751493+272*k*(35833+99*k*(293+143*k))
        ))))))))))))/(1.955457504e11L*k6)-(k14*Power(1+k,3)*(-2268+k*(
        -716+k*(384+k*(357+k*(107+12*k)))))*LOG20(ik))/240240.0L;
    case 10:
      return (88884432+k*(161695296+k*(26883528+k*(-133613536+k*(
        -129666642+k*(-54005868+k*(-10288908+k*(4910304+k*(39641907+k*(
        146898746+19*k*(19452021+4*k*(8327103+17*k*(580531+6*k*(73633
        +30888*k))))))))))))))/(1.955457504e11L*k6)+(k14*Power(1+k,2)*(
        -2184+k*(-1904+k*(1156+k*(2184+k*(1115+4*k*(65+6*k))))))
        *LOG20(ik))/240240.0L;
    case 11:
      return (-29628144-k*(59976000+k*(11680088+k*(-53730040+k*(
        -53501994+k*(-23268546+k*(-4969020+k*(697284+k*(7964627+k*(
        30510437+19*k*(4157067+k*(7344785+68*k*(131926+3*k*(34539+14872
        *k))))))))))))))/(1.955457504e11L*k6)-(k14*Power(1+k,3)*(-728+k
        *(-56+k*(516+k*(368+k*(107+12*k)))))*LOG20(ik))/240240.0L;
    case 14:
      return (-54698112+k*(-150796800+k*(-167181672+k*(-97693484+k*(
        -32689650+k*(-5730536+k*(1176700+k*(11041798+k*(40284417+19*k*(
        5244180+k*(8794443+68*k*(149535+k*(110345+44616*k)))))))))))))
        /(1.955457504e11L*k5)+(k15*Power(1+k,2)*(28+3*k*(6+k))*(48+k*(56
        +k*(23+4*k)))*LOG20(ik))/240240.0L;
    case 15:
      return (18232704-k*(-55529208+k*(-66140200+k*(-40489974+k*(
        -14233590+k*(-2726220+k*(88956+k*(2172171+k*(8185727+19*k*(
        1093965+k*(1888659+34*k*(65963+18*k*(2779+1144*k)))))))))))))
        /(1.955457504e11L*k5)-(k15*Power(1+k,3)*(448+k*(492+k*(240+k*(59
        +6*k))))*LOG20(ik))/240240.0L;
    }
    break;

  case 3:
    switch (s) {
    case 0:
      return (2444321880+k*(2315673360+k*(-46750149600+k*(-49251595008+k
        *(7117417692+k*(2860506264+k*(2709173544+k*(-3330890640+k*(
        2620372725+2*k*(-1456720605+4*k*(441407880+k*(-194704765+19*k*(
        51975145+36*k*(908791+17*k*(159233+k*(76121+100386*k))))))))))))
        ))))/(1.955457504e11L*k8)+(k12*(-92820+k*(-185640+k*(1683000+k*(
        3738912+k*(1693047+2*k*(-196263+4*k*(-27147+5*k*(465+k*(837+20*k
        *(11+k))))))))))*LOG20(ik))/371280.0L;
    case 1:
      return (-2444321880-k*(4759995240+k*(-49323120000+k*(-71150808960
        +k*(5666209164+k*(6063506988+k*(1440114984+k*(-1864086840+k*(
        911273725+k*(-1294391745+2*k*(852592230+k*(-250330990+57*k*(
        41619035+6*k*(5660465+136*k*(107791+58917*k+72501*k2))))))))))))
        )))/(5.866372512e11L*k8)-(k12*(1+k)*(-92820+k*(-185640+k*(
        1868640+k*(2805000+k*(-182121+2*k*(-137913+25*k*(-385+k*(619+218
        *k+24*k2))))))))*LOG20(ik))/1.11384e6L;
    case 2: case 8:
      return (-564074280+k*(-856714320+k*(11225313792+k*(17700497892+k*(
        4717150284+k*(-1019411316+k*(-1763241480+k*(145906430+k*( 
        -473629035+2*k*(431639130+k*(384801010+57*k*(34999385+6*k*(
        8102923+306*k*(39973+4*k*(7114+4719*k)))))))))))))))/(
        5.866372512e11L*k7)+(k13*(21420+k*(55080+k*(-392088+k*(-1120878
        +k*(-885321+2*k*(-74678+25*k*(2154+k*(1299+286*k+24*k2))))))))
        *LOG20(ik))/1.11384e6L;
    case 3: case 9:
      return (62674920-k*(-157865400+k*(1294129760+k*(2551101476+k*(
        747680472+k*(-162001224+k*(-205650060+k*(-28443010+k*(-25413465
        +k*(49796385+2*k*(31901015+19*k*(7389245+4*k*(2795950+51*k*(
        82368+k*(62141+40898*k)))))))))))))))/(1.955457504e11L*k7)-(k13
        *(1+k)*(2380+k*(6120+k*(-48960+k*(-99654+k*(-30549+10*k*(692+k*(
        733+20*k*(10+k))))))))*LOG20(ik))/371280.0L;
    case 4:
      return (-21998896920+k*(-94170716640+k*(-70495996032+k*(8854023024
        +k*(2189293260+k*(922166784+k*(-3784861080+k*(2249752680+k*(
        -2658128385+2*k*(1569438930+k*(-755479430+57*k*(57405725+6*k*(
        5222647+136*k*(123086+52560*k+72501*k2))))))))))))))
        /(5.866372512e11L*k7)+(k13*(835380+k*(4455360+k*(6438648+k*(
        2473704+k*(-435915+2*k*(-61458+25*k*(2154+k*(1299+286*k+24*k2)))
        )))))*LOG20(ik))/1.11384e6L;
    case 5:
      return (2444321880-k*(-12092960880+k*(-11549401248+k*(925351812+k
        *(674036748+k*(-43361472+k*(-289757160+k*(67581710+k*(-129884865
        +k*(164099910+k*(-56747840+19*k*(22293425+72*k*(225799+34*k*(
        17911+k*(8839+11154*k)))))))))))))))/(1.955457504e11L*k7)-(k13
        *(1+k)*(92820+k*(464100+k*(457572+k*(-32130+k*(-30537+5*k*(1087
        +k*(1433+40*k*(10+k))))))))*LOG20(ik))/371280.0L;
    case 6: case 12:
      return (1692222840+k*(7759626336+k*(7405630848+k*(178613820+k*(
        -1527097572+k*(-973484820+k*(-34344800+k*(-149702190+k*(
        249954705+2*k*(98523025+114*k*(4731575+3*k*(2095399+34*k*(92123
        +k*(62475+40898*k))))))))))))))/(5.866372512e11L*k6)+(k14*(
        -64260+k*(-362304+k*(-591192+k*(-302130+k*(50923+10*k*(9795+2*k
        *(2007+5*k*(77+6*k))))))))*LOG20(ik))/1.11384e6L;
    case 7: case 13:
      return (-564074280-k*(2962591632+k*(3434538492+k*(318476928+k*(
        -627486552+k*(-376007940+k*(-63115990+k*(-25524840+k*(41844165+k
        *(49469300+57*k*(3889450+3*k*(1887587+204*k*(13643+44*k*(223+143
        *k))))))))))))))/(5.866372512e11L*k6)-(k14*(1+k)*(-21420+k*(
        -113628+k*(-135150+k*(-13962+5*k*(5029+k*(2831+655*k+60*k2))))))
        *LOG20(ik))/1.11384e6L;
    case 10:
      return (161164080+k*(1033712064+k*(632804634+k*(-1533801192+k*(
        -1856168622+k*(-889631820+k*(-134367085+k*(-46301640+k*(
        140010255+k*(295912550+57*k*(16295465+12*k*(2140721+17*k*(162673
        +2*k*(59843+30888*k))))))))))))))/(5.866372512e11L*k6)+(k14*(
        -6120+k*(-45696+k*(-65331+k*(33048+k*(131753+5*k*(21570+k*(8127
        +20*k*(77+6*k))))))))*LOG20(ik))/1.11384e6L;
    case 11: 
      return (-53721360-k*(373222080+k*(223752606+k*(-721149198+k*(
        -812970774+k*(-376389090+k*(-76162375+k*(-11424765+k*(26137815+k
        *(64299275+57*k*(3474055+k*(5753405+204*k*(37122+k*(28449+14872
        *k))))))))))))))/(5.866372512e11L*k6)-(k14*(1+k)*(-2040+k*(
        -14280+k*(-9129+k*(27606+25*k*(1283+k*(586+k*(131+12*k)))))))
        *LOG20(ik))/1.11384e6L;
    case 14:
      return (-257862528+k*(-1415416464+k*(-2206990170+k*(-1456197204+k
        *(-529902450+k*(-83078240+3*k*(-5639425+k*(13180530+k*(26780325
        +19*k*(4395700+k*(6754671+68*k*(125551+k*(89275+44616*k)))))))))
        ))))/(5.866372512e11L*k5)+(k15*(9792+k*(64056+k*(140355+k*(
        143386+75*k*(1043+k*(324+k*(55+4*k)))))))*LOG20(ik))/1.11384e6L;
    case 15:
      return (85954176+k*(525974526-k*(-906327114+k*(-618916914+k*(
        -228655350+k*(-42919345+k*(-4840935+k*(7128495+k*(17251025+57*k
        *(914225+9*k*(164383+34*k*(6187+4578*k+2288*k2))))))))))))/(
        5.866372512e11L*k5)-(k15*(1+k)*(3264+k*(20145+k*(35286+5*k*(4879 
        +k*(1796+5*k*(71+6*k))))))*LOG20(ik))/1.11384e6L;
    }
    break;
  
  case 4:
    switch (s) {
    case 0:
      return (-2095133040+k*(-15952416480+k*(-27055090062+k*(
        -12084768696+k*(2166646482+k*(1141944804+k*(-76772793+k*(
        -116598912+k*(117581539+2*k*(347641325+2*k*(560321145+k*(
        1285309354+19*k*(116138455+6*k*(24436155+68*k*(338794+9*k*(23867
        +11154*k))))))))))))))))/(6.51819168e10L*k6*Power(1+k,2))+(k14*(
        3+20*k*(1+k))*(38760+k2*(-2907+2*k*(-171+2*k*(57+2*k*(10+k)))))
        *LOG20(ik))/180880.0L;
    case 1:
      return (2095133040-k*(-17349171840+k*(-19997799822+k*(1454076624+k
        *(1254347094+k*(328696368+k*(-407449545+k*(270731580+k*(
        -253473535+2*k*(366543440+k*(349812025+57*k*(24823335+4*k*(
        7116333+68*k*(154238+90888*k+72501*k2))))))))))))))/(
        1.955457504e11L*k6*(1+k))-(k14*(116280+k*(969000+k*(1154079+k*(
        -69027+2*k*(-42978+5*k*(-486+k*(773+242*k+24*k2)))))))*LOG20(ik)
        )/542640.0L;
    case 2: case 8:
      return (419026608+k*(3434915484+k*(6981416442+k*(5037904872+k*(
        831945114+k*(-461669502+k*(-246610325+k*(67961980+k*(550987155+2
        *k*(866053445+k*(1941702344+57*k*(56192551+4*k*(17024987+68*k*(
        218524+27*k*(4724+1573*k)))))))))))))))/(1.955457504e11L*k5
        *Power(1+k,2))+(k15*(-23256+k*(-168606+k*(-227601+2*k*(-31578+5
        *k*(1402+k*(1235+286*k+24*k2))))))*LOG20(ik))/542640.0L;
    case 3: case 9:
      return (-419026608-k*(3714266556+k*(5508336834+k*(1384407024+k*(
        -271993722+k*(-317362794+k*(-16426609+k*(12559960+k*(337387685+6
        *k*(130268030+19*k*(15920743+2*k*(11239660+51*k*(249493+2*k*(
        84349+40898*k))))))))))))))/(5.866372512e11L*k5*(1+k))+(k15*(
        23256+k*(207366-k*(-315267+k*(-82821+2*k*(8665+k*(8053+180*k*(11
        +k)))))))*LOG20(ik))/1.62792e6L;
    case 4:
      return (-5587021440+k*(-29736186480+k*(-43861385568+k*(
        -18461188746+k*(2137955820+k*(674239566+k*(-459316272+k*(
        -250311275+k*(77210230+k*(555475905+2*k*(902332805+k*(2008954820
        +57*k*(59585683+8*k*(9090417+34*k*(247280+3*k*(49614+24167*k))))
        ))))))))))))/(1.955457504e11L*k6*Power(1+k,2))+(k14*(310080+k*(
        1356600+k*(1147296+k*(-68115+2*k*(-27018+5*k*(1402+k*(1235+286*k
        +24*k2)))))))*LOG20(ik))/542640.0L;
    case 5:
      return (698377680-k*(-3455131680+k*(-3314026716+k*(199008810+k*(
        135957822+k*(3459456+k*(-61121538+k*(27317747+k*(-30206905+k*(
        72682730+k*(51828765+19*k*(13102147+12*k*(1123699+204*k*(8400+k
        *(4401+3718*k)))))))))))))))/(6.51819168e10L*k6*(1+k))-(k14*(
        38760+k*(193800+k*(191862+k*(-9747+k*(-9462+k*(1400+k*(1737+40*k
        *(11+k))))))))*LOG20(ik))/180880.0L;
    case 6: case 12:
      return (1117404288+k*(6226588368+k*(10191713532+k*(5454078630+k*(
        -462666204+k*(-1240661786+k*(-468256538+k*(15873021+k*(447731590
        +k*(1399411225+6*k*(512841257+114*k*(7279364+k*(8606179+17*k*(
        431265+244162*k+81796*k2))))))))))))))/(5.866372512e11L*k5
        *Power(1+k,2))+(k15*(8+5*k*(7+6*k))*(-7752+k*(-1938+k*(855+4*k*(
        152+k*(35+3*k)))))*LOG20(ik))/1.62792e6L;
    case 7: case 13:
      return (-139675536-k*(725945220+k*(828167340+k*(85447362+k*(
        -108222114+k*(-62821486+k*(-6607966+k*(-386485+k*(31222995+k*(
        67064530+57*k*(2749341+k*(3721583+136*k*(30475+66*k*(296+143*k))
        ))))))))))))/(1.955457504e11L*k5*(1+k))+(k15*(7752+k*(40698-k*(
        -47652+k*(-5377+k*(6420+k*(3374+715*k+60*k2))))))*LOG20(ik))
        /542640.0L;
    case 10:
      return (52378326+k*(24972948+k*(-340396056+k*(-646557912+k*(
        -470596581+k*(-158701914+k*(4376698+k*(146231220+k*(453282540+k
        *(992717932+57*k*(27872155+4*k*(8136945+17*k*(396509+18*k*(11981
        +3432*k))))))))))))))/(1.955457504e11L*k4*Power(1+k,2))+(k16*(
        -2907+k*(1368+k*(17613+k*(19170+k*(7903+20*k*(77+6*k))))))
        *LOG20(ik))/542640.0L;
    case 11:
      return (-52378326-k*(-108540432+k*(-524395872+k*(-456419964+k*(
        -186149145+k*(-31167780+k*(10872950+k*(79077650+3*k*(70357990+19
        *k*(7766260+k*(11043913+204*k*(56441+k*(37231+14872*k)))))))))))
        ))/(5.866372512e11L*k4*(1+k))-(k16*(-2907+k*(5871+k*(29583+5*k*(
        5271+k*(2119+429*k+36*k2)))))*LOG20(ik))/1.62792e6L;
    case 14:
      return (-139675536+k*(-897458562+k*(-1931409480+k*(-1920826908+k*(
        -993875064+k*(-285326405+k*(-18095160+k*(119153100+k*(368767000
        +9*k*(88224550+19*k*(7295634+k*(8333853+68*k*(99155+k*(52479
        +14872*k))))))))))))))/(5.866372512e11L*k4*Power(1+k,2))+(k16*(
        7752+k*(42465+k*(66918+5*k*(8603+k*(2852+9*k*(55+4*k))))))
        *LOG20(ik))/1.62792e6L;
    case 15:
      return (17459442-k*(-104864760+k*(-173128956+k*(-108756648+k*(
        -36974301+k*(-6018376+k*(672540+k*(7186520+k*(18517990+57*k*(
        673048+k*(930067+34*k*(27835+6*k*(2943+1144*k)))))))))))))/(
        1.955457504e11L*k4*(1+k))-(k16*(969+k*(5871+k*(9861+k*(6275+k*(
        2119+5*k*(77+6*k))))))*LOG20(ik))/542640.0L;
    }
    break;

  case 5:
    switch (s) {
    case 0:
      return (-6721885170+k*(-19811872080+k*(-19101992910+k*(-5193448260
        +k*(1251283554+k*(571688562+k*(1257536469+k*(6981914582+k*(
        26041189068+11*k*(6462335822+k*(13364127637+2*k*(10634684931+38
        *k*(343154237+33*k*(9722881+34*k*(197863+3*k*(32135+9126*k))))))
        ))))))))))/(5.377508136e11L*k4*Power(1+k,4))+(k16*(373065+k2*(
        -21945+2*k*(-1155+k*(693+20*k*(11+k)))))*LOG20(ik))/1.49226e6L;
    case 1:
      return (20165655510-k*(-39269960730+k*(-18167819670+k*(2199006810
        +k*(1169817925+k*(184843497+k*(1321386465+k*(7540101185+11*k*(
        2410416975+k*(6254747205+k*(12099549413+57*k*(312747245+22*k*(
        15725715+136*k*(94745+9*k*(5920+2197*k)))))))))))))))/(
        3.2265048816e12L*k4*Power(1+k,3))-(k16*(223839+k2*(-11704+k*(
        -1386+k*(693+242*k+24*k2))))*LOG20(ik))/1.790712e6L;
    case 2: case 8:
      return (2372430060+k*(7915037130+k*(9559089540+k*(4722408692+k*(
        679345220+k*(1723478037+k*(11079083470+k*(40978445245+11*k*(
        9957513798+k*(20081283355+k*(30977043778+57*k*(640768025+22*k*(
        25839865+34*k*(488047+108*k*(1930+429*k)))))))))))))))/(
        3.2265048816e12L*k3*Power(1+k,4))+(k17*(-26334+k*(-10241+k*(462
        +k*(1155+286*k+24*k2))))*LOG20(ik))/1.790712e6L;
    case 3: case 9:
      return (-395405010-k*(901800900+k*(626636010+k*(97375346+k*(
        -23887979+k*(107556928+k*(653853081+11*k* (207314820+k*(
        523272281+k*(984924290+19*k*(73465807+132*k*(587717+17*k*(26370
        +k*(13405+3718*k))))))))))))))/(1.0755016272e12L*k3*Power(1+k,3)
        )-(k17*(-21945+k*(-7315+k*(385+k*(847+20*k*(11+k)))))*LOG20(ik))
        /2.98452e6L;
    case 4:
      return (-20165655510+k*(-59435616240+k*(-57701383740+k*(
        -16912235340+k*(2076088785+k*(563758518+k*(1749637632+k*(
        11303792570+k*(41880699290+11*k*(10230748200+k*(20761068832+k*(
        32350465618+57*k*(679259765+22*k*(28134185+136*k*(139370+9*k*(
        7428+2197*k))))))))))))))))/(3.2265048816e12L*k4*Power(1+k,4))+(
        k16*(223839+k2*(-8778+k*(462+k*(1155+286*k+24*k2))))*LOG20(ik))
        /1.790712e6L;
    case 5:
      return (6721885170-k*(-13089986910+k*(-6143807670+k*(520269750+k*(
        215316853+k*(-21556149+k*(206754813+k*(1343607951+11*k*(
        419967665+k*(1072842519+k*(2023935025+19*k*(153185897+1584*k*(
        103526+17*k*(4850+k*(2611+1014*k)))))))))))))))/(
        2.1510032544e12L*k4*Power(1+k,3))-(k16*(373065+k2*(-14630+k2*(
        1617+40*k*(11+k))))*LOG20(ik))/5.96904e6L;
    case 6: case 12:
      return (1186215030+k*(3759816060+k*(4113599490+k*(1524789448+k*(
        -189640757+k*(259966512+k*(2979700407+k*(11077504488+11*k*(
        2654144763+k*(5263812736+k*(7968955729+114*k*(80717828+33*k*(
        2122249+34*k*(39220+k*(16517+3718*k)))))))))))))))/(
        3.2265048816e12*k3*Power(1+k,4))+(k17*(-65835+k*(-14630+k*(5775
        +2*k*(1848+5*k*(77+6*k)))))*LOG20(ik))/8.95356e6L;
    case 7: case 13:
      return (-1186215030-k*(2573601030+k*(1560809250+k*(31206214+k*(
        -177926511+k*(130287021+k*(1047889937+11*k* (328894395+k*(
        816156705+k*(1504824293+57*k*(36567973+11*k*(3418191+68*k*(37295
        +36*k*(512+143*k))))))))))))))/(6.4530097632e12L*k3*Power(1+k,3)
        )-(k17*(-65835+k*(-14630+k*(4620+k*(3234+715*k+60*k2))))
        *LOG20(ik))/1.790712e7L;
    case 10:
      return (-131801670+k*(-679819140+k*(-1370520788+k*(-1262892704+k*(
        412357017+k*(5853895106+k*(21734269349+11*k*(5188102536+k*(
        10238111471+k*(15387671306+57*k*(308381753+44*k*(5975713+17*k*(
        214313+6*k*(14165+2808*k))))))))))))))/(6.4530097632e12L*k2
        *Power(1+k,4))+(k18*(7315+k*(16170+k*(7623+20*k*(77+6*k))))
        *LOG20(ik))/1.790712e7L;
    case 11:
      return (131801670-k*(-402341940+k*(-470228720+k*(-246894564+k*(
        125672565+k*(1033439270+11*k*(324445335+k*(801563280+k*(
        1471760759+57*k*(35483050+11*k*(3279675+68*k*(34985+3*k*(5535
        +1352*k)))))))))))))/(6.4530097632e12L*k2*Power(1+k,3))-(k18*(
        1463+k*(1617+k*(693+k*(143+12*k))))*LOG20(ik))/3.581424e6L;
    case 14:
      return (-263603340+k*(-1005854850+k*(-1520323480+k*(-1135697836+k
        *(-177413574+k*(1549814675+k*(5916396140+33*k*(465058915+k*(
        903902482+k*(1335675063+19*k*(78809080+11*k*(5986245+68*k*(52585
        +k*(20509+4056*k))))))))))))))/(6.4530097632e12L*k2*Power(1+k,4)
        )+(k18*(2926+3*k*(847+k*(308+k*(55+4*k))))*LOG20(ik))
        /3.581424e6L;
    case 15:
      return (131801670-k*(-360720360+k*(-365299844+k*(-171646488+k*(
        8113521+k*(275997032+11*k*(86522235+k*(210409458+k*(379213415+57
        *k*(8953028+11*k*(807951+34*k*(16795+18*k*(431+104*k))))))))))))
        )/(6.4530097632e12L*k2*Power(1+k,3))-(k18*(7315+k*(5775+k*(2079
        +5*k*(77+6*k))))*LOG20(ik))/1.790712e7L;
    }
    break;

  case 6:
    switch (s) {
    case 0:
      return (257229+5215128*k+49961157*k2+300291466*k3+1267660373*k4
        +3985359096*k5+9649873307*k6+18343204858*k7+27622192884*k8
        +32967560692*k9+30927863054*k10+22378701828*k11+12063269632*k12
        +4500072176*k13+917867808*k14)/(1.92203088e10L*Power(1+k,6));
    case 1:
      return (-153790-3029543*k-28115821*k2-163121108*k3-661787065*k4
        -1988673380*k5-4570634497*k6-8171790704*k7-11430840817*k8
        -12450827015*k9-10382158600*k10-6408195968*k11-2741479296*k12
        -662904528*k13)/(5.76609264e10L*Power(1+k,5));
    case 2: case 8:
      return (215014+4314513*k+40858107*k2+242394391*k3+1008180768*k4
        +3116079621*k5+7397571512*k6+13739719323*k7+20126164674*k8
        +23226569872*k9+20888804079*k10+14290644648*k11+7099806152*k12
        +2327832576*k13+388328688*k14)/(5.76609264e10L*Power(1+k,6));
    case 3: case 9:
      return (-42330-824446*k-7553952*k2-43193361*k3-172340330*k4
        -507985935*k5-1141411524*k6-1986538543*k7-2689589364*k8
        -2812888155*k9-2224679600*k10-1273434606*k11-481310052*k12
        -93486536*k13)/(5.76609264e10L*Power(1+k,5));
    case 4:
      return (217091+4361139*k+41357253*k2+245783073*k3+1024544808*k4
        +3175829619*k5+7568710326*k6+14133326097*k7+20863842087*k8
        +24362463972*k9+22329346611*k10+15788955084*k11+8359942848*k12
        +3121606368*k13+662904528*k14)/(5.76609264e10L*Power(1+k,6));
    case 5:
      return (-28413-553923*k-5081442*k2-29100692*k3-116347460*k4
        -343885282*k5-775650942*k6-1357485668*k7-1853584513*k8
        -1965558265*k9-1593059028*k10-957092928*k11-404494192*k12
        -101985312*k13)/(3.84406176e10L*Power(1+k,5));
    case 6: case 12:
      return (60792+1209503*k+11346851*k2+66621856*k3+273934701*k4
        +835961998*k5+1956718077*k6+3577804054*k7+5151308459*k8
        +5834977709*k9+5146210762*k10+3454491118*k11+1690350086*k12
        +551558676*k13+93486536*k14)/(5.76609264e10L*Power(1+k,6));
    case 7: case 13:
      return (-23588-454993*k-4124417*k2-23303862*k3-91752240*k4
        -266444052*k5-588749722*k6-1005636038*k7-1333400068*k8
        -1363048885*k9-1052583983*k10-589208848*k11-219660672*k12
        -43147632*k13)/(1.153218528e11L*Power(1+k,5));
    case 10:
      return (120526+2395652*k+22448403*k2+131614084*k3+540193382*k4
        +1644695124*k5+3838051768*k6+6989120952*k7+10005611616*k8
        +11239873128*k9+9788460451*k10+6437864772*k11+3041521748*k12
        +931715464*k13+141210432*k14)/(1.153218528e11L*Power(1+k,6));
    case 11:
      return (-23450-451984*k-4093243*k2-23099944*k3-90808550*k4
        -263164160*k5-579884096*k6-986609932*k7-1300635556*k8-1317694840
        *k9-1002571955*k10-546447524*k11-193243148*k12-33995104*k13)
        /(1.153218528e11L*Power(1+k,5));
    case 14:
      return (102704+2025206*k+18812202*k2+109242947*k3+443669372*k4
        +1335206706*k5+3076189464*k6+5523464628*k7+7786696368*k8
        +8603344948*k9+7362891514*k10+4758983751*k11+2213752852*k12
        +671492452*k13+101985312*k14)/(3.459655584e11L*Power(1+k,6));
    case 15:
      return (-6566-125416*k-1124519*k2-6276534*k3-24373650*k4-69680064
        *k5-151230244*k6-252997976*k7-327354676*k8-324953200*k9
        -241950731*k10-129071446*k11-44841444*k12-7845024*k13)
        /(1.153218528e11L*Power(1+k,5));
    }
    break;

  case 7:
    switch (s) {
    case 0:
      return (29192499+650282802*k+6883842524*k2+46021317368*k3
        +217759955598*k4+774413949612*k5+2144709070140*k6+4727393663040
        *k7+8400101931645*k8+12104353247910*k9+14143572839880*k10
        +13318143548040*k11+9973049181080*k12+5803896006720*k13
        +2520312617120*k14+749852355360*k15+119322815040*k16)
        /(2.498640144e12L*Power(1+k,8));
    case 1:
      return (-17497935-379745082*k-3907035526*k2-25314605330*k3
        -115696700098*k4-395801222660*k5-1049197266400*k6-2199703025540
        *k7-3687821813665*k8-4960593656150*k9-5332722723790*k10
        -4526892260030*k11-2967354711060*k12-1443735219040*k13
        -480571906560*k14-86177588640*k15)/(7.495920432e12L
        *Power(1+k,7));
    case 2: case 8:
      return (24436098+539328659*k+5651260712*k2+37354148664*k3
        +174517581006*k4+611821887018*k5+1667196141000*k6+3607554326580
        *k7+6275315210250*k8+8821587004725*k9+10011462291600*k10
        +9101848383660*k11+6523149844950*k12+3582213518780*k13
        +1433362254920*k14+378123024960*k15+50482729440*k16)
        /(7.495920432e12L*Power(1+k,8));
    case 3: case 9:
      return (-14466450-310780459*k-3161569808*k2-20227499626*k3
        -91142943452*k4-306826089910*k5-798542450960*k6-1639172395000*k7
        -2681364078470*k8-3503925540445*k9-3638556154640*k10
        -2959575180250*k11-1834945152920*k12-825352688040*k13
        -243968657160*k14-36459749040*k15)/(2.2487761296e13L
        *Power(1+k,7));
    case 4:
      return (24623172+543904049*k+5704846412*k2+37753980864*k3
        +176650709106*k4+620477659638*k5+1694936376600*k6+3679512232380
        *k7+6428877949800*k8+9094037910375*k9+10415682056340*k10
        +9604063641660*k11+7043742963150*k12+4026173530380*k13
        +1732036285120*k14+520253464320*k15+86177588640*k16)
        /(7.495920432e12L*Power(1+k,8));
    case 5:
      return (-1616076-34741957*k-353731105*k2-2265556939*k3-10222200428
        *k4-34473003670*k5-89930220770*k6-185197136890*k7-304340470280
        *k8-400422046255*k9-420256153295*k10-347935731385*k11
        -222681862760*k12-106489068960*k13-35514780240*k14-6629045280
        *k15)/(2.498640144e12L*Power(1+k,7));
    case 6: case 12:
      return (20717024+453694118*k+4713708453*k2+30868865502*k3
        +142758406680*k4+494927187552*k5+1332252812370*k6+2844383326860
        *k7+4875898997340*k8+6746666621970*k9+7528706837355*k10
        +6726477658050*k11+4739283432660*k12+2564086496340*k13
        +1016134756920*k14+268245440520*k15+36459749040*k16)
        /(2.2487761296e13L*Power(1+k,8));
    case 7: case 13:
      return (-1343422-28608124*k-288241627*k2-1824731657*k3-8126676937
        *k4-27007614800*k5-69295536730*k6-140022444050*k7-225117588520
        *k8-288678566000*k9-293787214105*k10-234057862055*k11
        -142273878055*k12-62997474260*k13-18514256640*k14-2804596080
        *k15)/(7.495920432e12L*Power(1+k,7));
    case 10:
      return (6857781+150077424*k+1557945539*k2+10192298430*k3
        +47077822335*k4+162962451324*k5+437809334430*k6+932369437860*k7
        +1592928767385*k8+2194016070840*k9+2432555710365*k10
        +2152884812730*k11+1495140339675*k12+790592440740*k13
        +301741167140*k14+74821024920*k15+9178678080*k16)
        /(7.495920432e12L*Power(1+k,8));
    case 11:
      return (-4011525-85377858*k-859645519*k2-5437510546*k3-24191381361
        *k4-80288699120*k5-205644210470*k6-414566866740*k7-664364273485
        *k8-848027954950*k9-857134218105*k10-675596821550*k11
        -403491751145*k12-173293459980*k13-48225947820*k14-6629045280
        *k15)/(2.2487761296e13L*Power(1+k,7));
    case 14:
      return (5841528+126911969*k+1307118858*k2+8478372837*k3
        +38797607028*k4+132942107811*k5+353226671880*k6+743243317530*k7
        +1253376035520*k8+1702337111805*k9+1859617042470*k10
        +1620695162955*k11+1108412664540*k12+577814053335*k13
        +218047538940*k14+53774519340*k15+6629045280*k16)
        /(2.2487761296e13L*Power(1+k,8));
    case 15:
      return (-374409-7903582*k-78871005*k2-494032804*k3-2174551771*k4
        -7132915350*k5-18035865530*k6-35849754080*k7-56572546965*k8
        -71016611710*k9-70511326435*k10-54556347300*k11-31988764895*k12
        -13512088070*k13-3716134380*k14-509926560*k15)/(7.495920432e12L
        *Power(1+k,7));
    }
    break; 
  }

  fprintf(stderr, "***ERROR in K4\n");
  abort();

  return 0.0L;
}


/*----------------------------------------------------------------------
  K5

  k1 = k, k2 = k3 = k4 = k+1
----------------------------------------------------------------------*/
static long double K5(int ik, int l, int s)
{
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
  const long double k19 = k9*k10;

  switch (l) {
  case 0:
    switch (s) {
    case 0:
      return 4051331/2148854400.0L+184867*k/9.0288e6L+31178509*k2
        /3.133746e8L+89955713*k3/3.174444e8L+10274767*k4/1.979208e7L
        +1069333*k5/1.7017e6L+1550713*k6/3.15315e6L+478607*k7/2.1021e6L
        +117*k8/2450.0L;
    case 1: case 4:
      return -233687/1289312640.0L-399967*k/1.953504e8L-110708*k2
        /1.0683225e7L-4515949367*k3/1.466593128e11L-60683663*k4
        /1.02918816e9L-2565685*k5/3.4306272e7L-4658579*k6/7.56756e7L
        -10559*k7/350350.0L-33*k8/4900.0L;
    case 2:
      return -140317/429770880.0L-7920527*k/2.1488544e9L-70109969*k2
        /3.7604952e9L-8117431667*k3/1.466593128e11L-108903143*k4
        /1.02918816e9L-22944377*k5/1.7153136e8L-1033853*k6/9.45945e6L
        -9238*k7/175175.0L-169*k8/14700.0L;
    case 3: case 6:
      return 12049/379209600.0L+14811*k/3.97936e7L+88271221*k2
        /4.51259424e10L+1772776183*k3/2.933186256e11L+1550231*k4
        /1.2864852e8L+13663607*k5/8.576568e8L+67*k6/4900.0L+175949*k7
        /2.52252e7L+143*k8/88200.0L;
    case 5:
      return 158017/6446563200.0L+17069*k/5.96904e7L+2500397*k2
        /1.6713312e9L+79357199*k3/1.72540368e10L+383189*k4/4.200768e7L
        +3430541*k5/2.858856e8L+17071*k6/1.6632e6L+132101*k7/2.52252e7L
        +3*k8/2450.0L;
    case 7:
      return -509/117210240.0L-337999*k/6.4465632e9L-12827767*k2
        /4.51259424e10L-2815123*k3/3.1039008e9L-6041*k4/3.23136e6L
        -1756019*k5/6.8612544e8L-172541*k6/7.56756e7L-183373*k7
        /1.513512e8L-13*k8/44100.0L;
    case 8:
      return 985/1790712.0L+2536427*k/4.2977088e8L+562447*k2/1.98968e7L
        +81159569*k3/1.0255896e9L+29223013*k4/2.05837632e8L+28711*k5
        /171360.0L+120368*k6/945945.0L+119657*k7/2.1021e6L+169*k8
        /14700.0L;
    case 9: case 12:
      return -7081/134303400.0L-68993*k/1.1721024e8L-2705843*k2
        /9.209376e8L-64578637*k3/7.5209904e9L-5528291*k4/3.4306272e8L
        -2332763*k5/1.169532e8L-26813*k6/1.68168e6L-190283*k7/2.52252e7L
        -143*k8/88200.0L;
    case 10:
      return -24991/268606800.0L-445241*k/4.2977088e8L-41963*k2
        /8.1396e6L-56442983*k3/3.7604952e9L-11538337*k4/4.11675264e8L
        -447571*k5/1.29948e7L-11447*k6/420420.0L-53233*k7/4.2042e6L-13
        *k8/4900.0L;
    case 11: case 14:
      return 383/42411600.0L+787*k/7.53984e6L+9751*k2/1.80576e7L
        +36943561*k3/2.25629712e10L+1313825*k4/4.11675264e8L+1280431*k5
        /3.118752e8L+14753*k6/4.32432e6L+641*k7/382200.0L+11*k8/29400.0L;
    case 13:
      return 191/26860680.0L+105737*k/1.28931264e9L+1364401*k2
        /3.2232816e9L+57784679*k3/4.51259424e10L+854443*k4/3.4306272e8L
        +598721*k5/1.8712512e8L+80501*k6/3.027024e7L+66167*k7/5.04504e7L
        +13*k8/44100.0L;
    case 15:
      return -62/50363775.0L-1261*k/8.5954176e7L-3117*k2/3.97936e7L
        -11073389*k3/4.51259424e10L-254627*k4/5.1459408e8L-205721*k5
        /3.118752e8L-6893*k6/1.2108096e7L-4903*k7/1.68168e7L-k8/14700.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0:
      return 899271/724243520.0L+6744077*k/4.938024e8L+2492453*k2
        /3.69512e7L+1952633*k3/9.9484e6L+1008009*k4/2.72272e6L+1907869
        *k5/4.08408e6L+16691*k6/42900.0L+138839*k7/700700.0L+117*k8
        /2450.0L;
    case 1: case 4:
      return -275311/2172730560.0L-15724543*k/1.08636528e10L-149771*k2
        /2.01552e7L-447567*k3/1.98968e7L-55507*k4/1.25664e6L-355639*k5
        /6.12612e6L-12967*k6/257400.0L-9369*k7/350350.0L-33*k8/4900.0L;
    case 2:
      return -25583/127807680.0L-74822999*k/3.25909584e10L-33673*k2
        /2.8424e6L-587593*k3/1.62792e7L-3502333*k4/4.900896e7L-1162883
        *k5/1.225224e7L-10741*k6/128700.0L-23654*k7/525525.0L-169*k8
        /14700.0L;
    case 3: case 6:
      return 9193/434546112.0L+215251*k/8.576568e8L+1337351*k2
        /9.976824e8L+1514897*k3/3.581424e8L+5063*k4/583440.0L+19517*k5
        /1.633632e6L+61*k6/5600.0L+51253*k7/8.4084e6L+143*k8/88200.0L;
    case 5:
      return 7597/434546112.0L+560359*k/2.7159132e9L+485599*k2
        /4.434144e8L+409459*k3/1.193808e8L+97627*k4/1.400256e7L+232541
        *k5/2.450448e7L+4729*k6/554400.0L+39437*k7/8.4084e6L+3*k8
        /2450.0L;
    case 7:
      return -19561/6518191680.0L-3585163*k/9.77728752e10L-805433*k2
        /3.9907296e9L-470353*k3/7.162848e8L-27241*k4/1.9603584e7L-3071
        *k5/1.55584e6L-239*k6/128700.0L-54101*k7/5.04504e7L-13*k8
        /44100.0L;
    case 8:
      return 33443/90530440.0L+26214451*k/6.51819168e9L+250957*k2
        /1.27908e7L+159861*k3/2.8424e6L+5102789*k4/4.900896e7L+71437*k5
        /556920.0L+1865*k6/18018.0L+106147*k7/2.1021e6L+169*k8/14700.0L;
    case 9: case 12:
      return -30493/814773960.0L-752281*k/1.77768864e9L-4282789*k2
        /1.9953648e9L-3309*k3/516800.0L-1815349*k4/1.4702688e8L-105979
        *k5/6.68304e6L-1283*k6/96096.0L-171593*k7/2.52252e7L-143*k8
        /88200.0L;
    case 10:
      return -961/16460080.0L-1438139*k/2.17273056e9L-11669*k2
        /3.464175e6L-7823*k3/775200.0L-213089*k4/1.089088e7L-18787*k5
        /742560.0L-43*k6/2002.0L-15481*k7/1.4014e6L-13*k8/4900.0L;
    case 11: case 14:
      return 3329/543182640.0L+468511*k/6.51819168e9L+251881*k2
        /6.651216e8L+20047*k3/1.70544e7L+231613*k4/9.801792e7L+2173*k5
        /685440.0L+1007*k6/360360.0L+571*k7/382200.0L+11*k8/29400.0L;
    case 13:
      return 12541/2444321880.0L+234425*k/3.910915008e9L+313649*k2
        /9.976824e8L+297937*k3/3.069792e8L+284987*k4/1.4702688e8L+9853
        *k5/3.81888e6L+3247*k6/1.44144e6L+60077*k7/5.04504e7L+13*k8
        /44100.0L;
    case 15:
      return -8/9258795.0L-22679*k/2.17273056e9L-18851*k2/3.325608e8L
        -18569*k3/1.023264e8L-205*k4/544544.0L-1549*k5/2.97024e6L-1369
        *k6/2.88288e6L-1471*k7/5.6056e6L-k8/14700.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0:
      return (35529354+k*(14612724+k*(-6034014+k*(-14814324+k*(15905176
        +k*(156984296+k*(812432181+2*k*(1185065951+38*k*(59985602+3*k*(
        25749801+34*k*(655868+k*(354427+100386*k))))))))))))/(
        1.62954792e10L*k4)+(k16*(3+2*k)*(-291+k*(-232+k*(79+2*k*(60+k*(19
        +2*k)))))*LOG20(ik))/20020.0L;
    case 1: case 4:
      return (-71058708-k*(51934932+k*(-17292366+k*(-27559476+k*(2426060
        +k*(104935488+k*(559071543+k*(1701171098+57*k*(59575469+2*k*(
        39891555+34*k*(1052927+36*k*(16358+4719*k))))))))))))/(
        9.77728752e10L*k4)-(k16*(1+k)*(3+2*k)*(-582+k*(-68+k*(197+k*(91
        +12*k))))*LOG20(ik))/120120.0L;
    case 2:
      return (-71058708-k*(35818524+k*(-11519676+k*(-33685596+k*(
        13018180+k*(145633488+k*(819355602+k*(2490247178+57*k*(88688309
        +6*k*(20077805+136*k*(135936+k*(78376+24167*k))))))))))))/(
        9.77728752e10L*k4)+(k16*(1746-k*(-2718+k*(-638+k*(1122+k*(885
        +242*k+24*k2)))))*LOG20(ik))/120120.0L;
    case 3: case 6:
      return (23686236+k*(19509336+k*(-4957302+k*(-10802274+k*(-1832376
        +k*(17050992+k*(97921589+k*(309680021+114*k*(5680403+k*(7975627
        +34*k*(222079+k*(131569+40898*k))))))))))))/(9.77728752e10L*k4)
        +(k16*(1+k)*(-582+k*(-510+k*(129+2*k*(132+k*(49+6*k)))))
        *LOG20(ik))/120120.0L;
    case 5:
      return (47372472+k*(49762944+k*(-10032414+k*(-21447468+k*(-4723964
        +k*(29817984+k*(166922727+k*(526527562+57*k*(19082507+4*k*(
        6619503+17*k*(361543+209042*k+61776*k2)))))))))))/(
        1.955457504e11L*k4)+(k16*Power(1+k,2)*(-1164+k*(-120+k*(367+4*k
        *(43+6*k))))*LOG20(ik))/240240.0L;
    case 7:
      return (-15790824-k*(18052776+k*(-2390234+k*(-8270850+k*(-2863700
        +k*(4842292+k*(30184429+k*(98678637+19*k*(11182425+k*(16192605
        +68*k*(231918+k*(141233+44616*k))))))))))))/(1.955457504e11L*k4) 
        -(k16*Power(1+k,2)*(-388+k*(-76+k*(133+3*k*(25+4*k))))*LOG20(ik)
        )/240240.0L;
    case 8:
      return (-59337684+k*(-75261312+k*(-45740268+k*(19743612+k*(
        285323136+k*(1433488056+k*(4147331266+57*k*(137378687+2*k*(
        86707005+136*k*(535460+3*k*(92637+24167*k)))))))))))/(
        9.77728752e10L*k3)+(k17*(3+2*k)*(486+k*(804+k*(486+k*(125+12*k))
        ))*LOG20(ik))/120120.0L;
    case 9: case 12:
      return (19779228-k*(-30703428+k*(-18084906+k*(-1677396+k*(31088512
        +k*(162967336+3*k*(163607483+19*k*(16890429+4*k*(5543186+17*k*(
        284250+k*(152903+40898*k)))))))))))/(9.77728752e10L*k3)-(k17*(1
        +k)*(3+2*k)*(162+k*(152+51*k+6*k2))*LOG20(ik))/120120.0L;
    case 10:
      return (13186152-k*(-18678240+k*(-12449304+k*(237860+k*(29247008+k
        *(156994824+k*(474328660+19*k*(49579043+48*k*(1373163+17*k*(
        71850+k*(39653+11154*k)))))))))))/(6.51819168e10L*k3)-(k17*(324
        +k*(800+k*(788+k*(381+8*k*(11+k)))))*LOG20(ik))/80080.0L;
    case 11: case 14:
      return  (-13186152+k*(-22422456+k*(-14775516+k*(-3238844+k*(
        9870756+k*(55898808+k*(175084078+57*k*(6294559+k*(8647221+68*k*(
        116515+132*k*(501+143*k)))))))))))/(1.955457504e11L*k3)+(k17*(1
        +k)*(324+k*(568+k*(374+k*(109+12*k))))*LOG20(ik))/240240.0L;
    case 13:
      return (-13186152+k*(-24213168+k*(-15128946+k*(-3432576+k*(8592556
        +k*(48264608+k*(150565129+57*k*(5365254+k*(7299325+68*k*(96935+k
        *(53987+14872*k)))))))))))/(1.955457504e11L*k3)+(k17
        *Power(1+k,2)*(324+k*(288+k*(97+12*k)))*LOG20(ik))/240240.0L;
    case 15:
      return (4395384-k*(-8722224+k*(-6093990+k*(-1770608+k*(1310820+k*(
        8553048+k*(27650467+57*k*(1024460+k*(1451571+34*k*(40295+6*k*(
        3923+1144*k)))))))))))/(1.955457504e11L*k3)-(k17*Power(1+k,2)*(
        108+k*(112+k*(43+6*k)))*LOG20(ik))/240240.0L;
    }
    break;

  case 3:
    switch (s) {
    case 0:
      return (1462564026+k*(4831180200+k*(4889074806+k*(619605756+k*(
        -1621825961+k*(704922400+k*(10429007075+2*k*(19421554845+2*k*(
        23813422105+k*(41245467170+19*k*(2707673951+6*k*(403869995+68*k
        *(3698687+3*k*(493109+100386*k))))))))))))))/(1.955457504e11L*k4
        *Power(1+k,2))+(k16*(-55539+2*k*(-65421+4*k*(-7672+5*k*(875+k*(
        837+20*k*(11+k))))))*LOG20(ik))/371280.0L;
    case 1: case 4:
      return (-1462564026-k*(4175621604+k*(2393881182+k*(-936098856+k*(
        -1029323505+k*(269969490+k*(3155457855+2*k*(5482809860+k*(
        12420204995+57*k*(343591925+4*k*(94554991+68*k*(1051913+27*k*(
        18697+4719*k)))))))))))))/(5.866372512e11L*k4*(1+k))+(k16*(55539 
        -k*(-161487+2*k*(-48078+25*k*(726+k*(837+242*k+24*k2)))))
        *LOG20(ik))/1.11384e6L;
    case 2:
      return (-1462564026-k*(4939965954+k*(5211760554+k*(794074050+k*(
        -1949698961+k*(-552751465+k*(4749881245+k*(19258245325+2*k*(
        24665084455+k*(44681693420+57*k*(1025422369+8*k*(120716807+34*k
        *(2337148+9*k*(110470+24167*k))))))))))))))/(5.866372512e11L*k4
        *Power(1+k,2))+(k16*(55539-k*(-134973+2*k*(-34858+25*k*(726+k*(
        837+242*k+24*k2)))))*LOG20(ik))/1.11384e6L;
    case 3: case 6:
      return (487521342+k*(1428135786+k*(887562522+k*(-308145222+k*(
        -429705535+k*(-42635005+k*(522818505+k*(1934877505+2*k*(
        2278090925+114*k*(32931020+3*k*(12645208+17*k*(592417+300554*k
        +81796*k2))))))))))))/(5.866372512e11L*k4*(1+k))+(k16*(-18513+k
        *(-55206+k*(-35527+10*k*(1185+2*k*(852+5*k*(55+6*k))))))
        *LOG20(ik))/1.11384e6L;
    case 5:
      return (487521342+k*(1173354336+k*(50094198+k*(-362951820+k*(
        -75391575+k*(17708460+k*(448249395+k*(1247991910+57*k*(47919245
        +4*k*(16405005+17*k*(934297+6*k*(90731+30888*k))))))))))))/(
        5.866372512e11L*k4)+(k16*(-18513+k*(-64044+k*(-48747+5*k*(2370+k
        *(3507+20*k*(55+6*k))))))*LOG20(ik))/1.11384e6L;
    case 7:
      return (-162507114-k*(403205418+k*(40578846+k*(-134892450+k*(
        -41275325+3*k*(-2455915+k*(28336605+k*(75205625+19*k*(9299125+k
        *(13103865+68*k*(197644+k*(120163+44616*k))))))))))))/(
        5.866372512e11L*k4)+(k16*(6171-k*(-21807+k*(-17639+75*k*(47+k*(
        93+k*(33+4*k))))))*LOG20(ik))/1.11384e6L;
    case 8:
      return (-979071786+k*(-4064553108+k*(-6461619318+k*(-4899516216+k
        *(-400380015+k*(9174756150+k*(34479133205+2*k*(41823749495+k*(
        71465540080+57*k*(1537712725+4*k*(336637609+68*k*(2993539+9*k*(
        127461+24167*k)))))))))))))/(5.866372512e11L*k3*Power(1+k,2))+(
        k17*(37179+2*k*(59562+25*k*(2646+k*(1299+286*k+24*k2))))
        *LOG20(ik))/1.11384e6L;
    case 9: case 12:
      return (108785754-k*(-397737648+k*(-468355734+k*(-232275540+k*(
        731465+k*(308826980+k*(1064599525+2*k*(1192288680+19*k*(97311775
        +6*k*(17510540+17*k*(757591+350222*k+81796*k2)))))))))))/(
        1.955457504e11L*k3*(1+k))-(k17*(4131+k*(15321+10*k*(1835+k*(933
        +20*k*(11+k)))))*LOG20(ik))/371280.0L;
    case 10:
      return (108785754-k*(-467733420+k*(-785893878+k*(-653423652+k*(
        -207256735+k*(436953790+k*(1861783605+k*(4715497435+k*(
        8414074790+19*k*(568514345+36*k*(14516769+68*k*(136000+k*(55171
        +11154*k)))))))))))))/(1.955457504e11L*k3*Power(1+k,2))-(k17*(
        4131+k*(13848+5*k*(3340+k*(1833+40*k*(11+k)))))*LOG20(ik))
        /371280.0L;
    case 11: case 14:
      return (-108785754+k*(-413854056+k*(-525442302+k*(-293784120+k*(
        -48911415+k*(149055060+k*(552551065+k*(1286166475+57*k*(36482380
        +k*(41145875+136*k*(233423+198*k*(574+143*k))))))))))))/(
        5.866372512e11L*k3*(1+k))+(k17*(4131+k*(15933+5*k*(4110+k*(2331
        +605*k+60*k2))))*LOG20(ik))/1.11384e6L;
    case 13:
      return (-108785754+k*(-343858284+k*(-244719090+k*(-61148640+k*(
        8850765+k*(123832170+k*(365517295+57*k*(13470490+k*(18358395+68
        *k*(252115+3*k*(47897+14872*k)))))))))))/(5.866372512e11L*k3)+(
        k17*(4131+k*(17406+25*k*(921+k*(486+k*(121+12*k)))))*LOG20(ik))
        /1.11384e6L;
    case 15:
      return (36261918-k*(-119991564+k*(-96985350+k*(-28985040+k*(
        -1657425+k*(22284570+k*(65676325+57*k*(2548700+k*(3595845+34*k*(
        103745+18*k*(3433+1144*k)))))))))))/(5.866372512e11L*k3)-(k17*(
        1377+k*(6006+5*k*(1695+k*(996+5*k*(55+6*k)))))*LOG20(ik))
        /1.11384e6L;
    }
    break;
  
  case 4:
    switch (s) {
    case 0:
      return (52378326+k*(509729220+k*(1577223648+k*(2149138992+k*(
        1225734237+k*(418382874+k*(3407205022+k*(16541991936+k*(
        49740615729+2*k*(53654885263+4*k*(21575953035+k*(26197719087+19
        *k*(1261299491+12*k*(71196042+17*k*(2017066+9*k*(70199+11154*k))
        ))))))))))))))/(6.51819168e10L*k4*Power(1+k,4))+(k16*(3+20*k*(1
        +k))*(-969+2*k*(-57+2*k*(57+2*k*(10+k))))*LOG20(ik))/180880.0L;
    case 1: case 4:
      return (-52378326-k*(544648104+k*(1474737264+k*(1482208728+k*(
        339673425+k*(-169293852+k*(986805950+k*(4919772080+k*(
        13907512575+2*k*(13962451405+k*(20593132543+57*k*(393667891+2*k
        *(156487663+34*k*(2572249+108*k*(8585+1573*k)))))))))))))))/(
        1.955457504e11L*k4*Power(1+k,3))+(k16*(2907-k*(-24567+2*k*(
        -15618+5*k*(198+k*(773+242*k+24*k2)))))*LOG20(ik))/542640.0L;
    case 2:
      return (-52378326-k*(512810298+k*(1606845240+k*(2228346120+k*(
        1268321145+k*(65893191+k*(1203878078+k*(7681393930+k*(
        24456866855+k*(54870322355+2*k*(45910732928+k*(58111369474+57*k
        *(974592865+2*k*(346015665+136*k*(1289840+427692*k+72501*k2)))))
        ))))))))))/(1.955457504e11L*k4*Power(1+k,4))+(k16*(2907-k*(
        -19893+2*k*(-11058+5*k*(198+k*(773+242*k+24*k2)))))*LOG20(ik))
        /542640.0L;
    case 3: case 6:
      return (52378326+k*(547729182+k*(1506412908+k*(1555145592+k*(
        366839109+k*(-317325099+k*(315196648+k*(2472502086+k*(7354085805
        +k*(15305200429+6*k*(3904104403+114*k*(38803441+k*(32164813+34*k
        *(553287+k*(209883+40898*k)))))))))))))))/(5.866372512e11L*k4
        *Power(1+k,3))+(k16*(3+5*k*(5+6*k))*(-969+k*(-171+4*k*(57+k*(25
        +3*k))))*LOG20(ik))/1.62792e6L;
    case 5:
      return (17459442+k*(193188996+k*(445777332+k*(259363104+k*(
        -61798191+k*(-52096954+k*(136152848+k*(606224520+k*(1595526020+k
        *(2953188092+57*k*(69190903+4*k*(16556671+17*k*(649547+18*k*(
        15413+3432*k))))))))))))))/(1.955457504e11L*k4*Power(1+k,2))+(
        k16*(-969+k*(-9804+k*(-15447+k*(290+k*(3283+20*k*(55+6*k))))))
        *LOG20(ik))/542640.0L;
    case 7:
      return (-17459442-k*(194216022+k*(457020564+k*(280900620+k*(
        -61939059+k*(-81002285+k*(53689740+k*(324894900+k*(889862050+9*k
        *(189288590+19*k*(13797657+k*(13719147+68*k*(140348+k*(62775
        +14872*k))))))))))))))/(5.866372512e11L*k4*Power(1+k,2))+(k16*(
        969+k*(9861+k*(16017-5*k*(-73+k*(773+9*k*(33+4*k))))))*LOG20(ik)
        )/1.62792e6L;
    case 8:
      return (-27729702+k*(-299459160+k*(-1098016920+k*(-1955023980+k*(
        -1504631765+k*(2337037248+k*(14651314470+k*(44231815730+k*(
        94574643075+2*k*(75207663780+k*(90079624154+57*k*(1421409895+2*k
        *(471089935+136*k*(1620530+486855*k+72501*k2))))))))))))))/(
        1.955457504e11L*k3*Power(1+k,4))+(k17*(1539+2*k*(6042+5*k*(2086
        +k*(1235+286*k+24*k2))))*LOG20(ik))/542640.0L;
    case 9: case 12:
      return (27729702-k*(-317945628+k*(-1037548512+k*(-1498705572+k*(
        -949026169+k*(654182438+k*(4307904666+k*(12163184330+k*(
        24159710069+6*k*(5861123802+19*k*(330764005+4*k*(64416421+51*k*(
        687048+k*(238217+40898*k))))))))))))))/(5.866372512e11L*k3
        *Power(1+k,3))-(k17*(1539+k*(14649+2*k*(13795+k*(8053+180*k*(11
        +k)))))*LOG20(ik))/1.62792e6L;
    case 10:
      return (3081078-k*(-33729696+k*(-127351224+k*(-239054634+k*(
        -232140181+k*(16713788+k*(723839558+k*(2368508357+k*(5266899363
        +k*(8703646142+k*(10850029320+19*k*(535702579+24*k*(15471079+102
        *k*(74442+k*(23563+3718*k)))))))))))))))/(6.51819168e10L*k3
        *Power(1+k,4))-(k17*(171+k*(1368+k*(2540+k*(1737+40*k*(11+k)))))
        *LOG20(ik))/180880.0L;
    case 11: case 14:
      return (-3081078+k*(-35783748+k*(-120576456+k*(-183942486+k*(
        -137132177+k*(-2178358+k*(232579334+k*(701219695+k*(1442500961+k
        *(2175638018+57*k*(42479190+k*(34419809+68*k*(287189+132*k*(790
        +143*k))))))))))))))/(1.955457504e11L*k3*Power(1+k,3))+(k17*(171
        +k*(1653+k*(3310+k*(2219+605*k+60*k2))))*LOG20(ik))/542640.0L;
    case 13:
      return (-9243234+k*(-112144032+k*(-319531212+k*(-337319164+k*(
        -133414645+k*(103505220+k*(527083700+k*(1374495650+3*k*(
        837718970+19*k*(57961580+k*(54441683+204*k*(173441+k*(71551
        +14872*k)))))))))))))/(5.866372512e11L*k3*Power(1+k,2))+(k17*(
        513+k*(5738+5*k*(2459+k*(1426+363*k+36*k2))))*LOG20(ik))
        /1.62792e6L;
    case 15:
      return (1027026-k*(-12612600+k*(-37249212+k*(-42182504+k*(
        -20289997+k*(2630264+k*(30719260+k*(83603960+k*(157945766+57*k*(
        3770456+k*(3671555+34*k*(72955+6*k*(5231+1144*k)))))))))))))/(
        1.955457504e11L*k3*Power(1+k,2))-(k17*(57+k*(646+k*(1455+k*(964
        +5*k*(55+6*k)))))*LOG20(ik))/542640.0L;
    }
    break;

  case 5:
    switch (s) {
    case 0:
      return (263603340+k*(1331890560+k*(3168863893+k*(11133474092+k*(
        69898380903+2*k*(172447462113+k*(601392111996+k*(1546589949390+k
        *(3030007627728+11*k*(418677925446+k*(497563294122+k*(
        461817525394+19*k*(17432331641+198*k*(47653093+68*k*(273489+k*(
        70043+9126*k))))))))))))))))/(1.0755016272e12L*k2*Power(1+k,6))
        +(k18*(-7315+2*k*(-385+k*(693+20*k*(11+k))))*LOG20(ik))
        /1.49226e6L;
    case 1: case 4:
      return (-395405010-k*(1602430830+k*(2587072995+k*(4365455484+k*(
        22287458544+k*(106150640617+k*(353203906775+k*(858227031685+11*k
        *(143013775031+k*(200842077001+k*(216850653323+57*k*(3133576805
        +264*k*(7296835+51*k*(62773+6*k*(2983+429*k)))))))))))))))/(
        3.2265048816e12L*k2*Power(1+k,5))+(k18*(4389-k*(-462+k*(693+242
        *k+24*k2)))*LOG20(ik))/1.790712e6L;
    case 2:
      return (-263603340-k*(1345764420+k*(2935380955+k*(6338324889+k*(
        31646569695+k*(158089532973+k*(569551733337+k*(1513818781140+k*(
        3067645003611+11*k*(438992543010+k*(541144761537+k*(521952151628
        +57*k*(6840528725+132*k*(29302121+68*k*(176380+47634*k+6591*k2))
        )))))))))))))/(3.2265048816e12L*k2*Power(1+k,6))+(k18*(2926-k*(
        -462+k*(693+242*k+24*k2)))*LOG20(ik))/1.790712e6L;
    case 3: case 6:
      return (131801670+k*(541080540+k*(854432189+k*(972083593+k*(
        3610031672+k*(17573673969+k*(60438171285+k*(151416575910+11*k*(
        26027629037+k*(37747876828+k*(42149770631+114*k*(315507845+33*k
        *(6103462+17*k*(164117+48998*k+7436*k2))))))))))))))/(
        3.2265048816e12L*k2*Power(1+k,5))+(k18*(-7315+k*(-1155+2*k*(693
        +5*k*(55+6*k))))*LOG20(ik))/8.95356e6L;
    case 5:
      return (395405010+k*(1207025820+k*(1264754712+k*(1105708344+k*(
        5677095177+k*(26679529826+k*(84233395049+11*k*(17474964816+k*(
        29742436003+k*(38187618106+57*k*(647107613+132*k*(3497141+17*k*(
        103921+33946*k+5616*k2)))))))))))))/(6.4530097632e12L*k2
        *Power(1+k,4))+(k18*(-21945+k*(-2310+k*(3003+20*k*(55+6*k))))
        *LOG20(ik))/1.790712e7L;
    case 7:
      return (-131801670-k*(409278870+k*(430862900+k*(247360412+k*(
        937327911+k*(4677078155+k*(15267513905+33*k*(1086922465+k*(
        1905567531+k*(2523067523+19*k*(132458665+11*k*(8887635+68*k*(
        68500+k*(23317+4056*k))))))))))))))/(6.4530097632e12L*k2
        *Power(1+k,4))+(k18*(1463-3*k*(-77+k*(77+k*(33+4*k))))*LOG20(ik)
        )/3.581424e6L;
    case 8:
      return (-124864740+k*(-263651310+k*(6059403675+k*(60261922032+k*(
        309912216029+k*(1079013706467+k*(2756114615040+k*(5355250643425
        +11*k*(732895252218+k*(861188811111+k*(788548190368+57*k*(
        9757524125+528*k*(9791920+17*k*(218501+53757*k+6591*k2))))))))))
        ))))/(3.2265048816e12L*k*Power(1+k,6))+(k19*(1386+k*(1155+286*k
        +24*k2))*LOG20(ik))/1.790712e6L;
    case 9: case 12:
      return (20810790-k*(-80030652+k*(100729265+k*(2021526208+k*(
        10411943958+k*(34604153660+k*(83418621915+11*k*(13768202020+k*(
        19120667317+k*(20373656842+19*k*(869189725+22*k*(23805296+51*k*(
        199455+54786*k+7436*k2)))))))))))))/(1.0755016272e12L*k
        *Power(1+k,5))-(k19*(1155+k*(847+20*k*(11+k)))*LOG20(ik))
        /2.98452e6L;
    case 10:
      return (27747720-k*(-124283562+k*(316063618+k*(5502527706+k*(
        30776706194+k*(111445303953+k*(294188374350+k*(590720097895+11*k
        *(83632986084+k*(101801746101+k*(96716322802+19*k*(3732373817
        +396*k*(5204177+68*k*(30348+k*(7837+1014*k)))))))))))))))/(
        2.1510032544e12L*k*Power(1+k,6))-(k19*(1540+k*(1617+40*k*(11+k))
        )*LOG20(ik))/5.96904e6L;
    case 11: case 14:
      return (-41621580+k*(-186013074+k*(-135660330+k*(1684093086+k*(
        10069469228+k*(34886997505+k*(86697573320+11*k*(14745732565+k*(
        21121889834+k*(23240392537+57*k*(341759860+33*k*(6464171+816*k*(
        3515+k*(1006+143*k))))))))))))))/(6.4530097632e12L*k
        *Power(1+k,5))+(k19*(2310+k*(2079+605*k+60*k2))*LOG20(ik))
        /1.790712e7L;
    case 13:
      return (-62432370+k*(-210943980+k*(-81262740+k*(1509830556+k*(
        7753516355+k*(24396366770+11*k*(5014341525+k*(8440877200+k*(
        10698479389+57*k*(178507030+11*k*(11357225+204*k*(27435+k*(8655
        +1352*k)))))))))))))/(6.4530097632e12L*k*Power(1+k,4))+(k19*(
        693+k*(462+k*(121+12*k)))*LOG20(ik))/3.581424e6L;
    case 15:
      return (20810790-k*(-75624744+k*(-75926916+k*(208047528+k*(
        1326531713+k*(4344780872+11*k*(918795243+k*(1591287922+k*(
        2076928471+57*k*(35723684+99*k*(260679+34*k*(3907+2*k*(639+104*k
        )))))))))))))/(6.4530097632e12L*k*Power(1+k,4))-(k19*(1155+k*(
        924+5*k*(55+6*k)))*LOG20(ik))/1.790712e7L;
    }
    break;

  case 6:
    switch (s) {
    case 0:
      return (22864135+438151546*k+3967130546*k2+22543952796*k3
        +90069777782*k4+268533171940*k5+618711512208*k6+1125325459224*k7
        +1635636686949*k8+1910146973798*k9+1791159417890*k10
        +1339148225332*k11+786335985576*k12+352971545896*k13
        +115394122112*k14+24938261520*k15+2753603424*k16)
        /(5.76609264e10L*Power(1+k,8));
    case 1: case 4:
      return (-8040235-149502136*k-1309185585*k2-7168042261*k3
        -27467985646*k4-78118674534*k5-170555818164*k6-291540810660*k7
        -394121880459*k8-422364141314*k9-356978343581*k10-234652679361
        *k11-116805148790*k12-41932871636*k13-9843008976*k14-1164986064
        *k15)/(1.729827792e11L*Power(1+k,7));
    case 2:
      return (-9450035-185806201*k-1727200070*k2-10083611318*k3
        -41418471319*k4-127049972864*k5-301431573228*k6-565070814600*k7
        -847396340427*k8-1022260491751*k9-991611591626*k10-768256983742
        *k11-468511351135*k12-219054428018*k13-74886221824*k14
        -17013717504*k15-1988713584*k16)/(1.729827792e11L*Power(1+k,8));
    case 3: case 6:
      return (406735+7744499*k+69483411*k2+389993456*k3+1532913100*k4
        +4474600286*k5+10033949866*k6+17629447284*k7+24517437029*k8
        +27056521385*k9+23577658409*k10+16004693244*k11+8244948046*k12
        +3072812696*k13+752260540*k14+93486536*k15)/(5.76609264e10L
        *Power(1+k,7));
    case 5:
      return (2299430+41392764*k+349681041*k2+1839273788*k3+6736908246
        *k4+18201573180*k5+37465446396*k6+59797537344*k7+74540976552*k8
        +72435162424*k9+54230347485*k10+30503986164*k11+12289972996*k12
        +3218777688*k13+423631296*k14)/(3.459655584e11L*Power(1+k,6));
    case 7:
      return (-372890-6865838*k-59355339*k2-319648013*k3-1199373722*k4
        -3321353490*k5-7011508692*k6-11484870876*k7-14703719772*k8
        -14687958148*k9-11317039255*k10-6561690033*k11-2731613584*k12
        -742097668*k13-101985312*k14)/(3.459655584e11L*Power(1+k,6));
    case 8:
      return (21023940+401308425*k+3617966604*k2+20462930006*k3
        +81330368179*k4+241077778196*k5+551869567764*k6+996457552176*k7
        +1436359991154*k8+1661468982711*k9+1540650987144*k10
        +1136588932806*k11+656567371189*k12+288663950006*k13+91802946304
        *k14+19085873616*k15+1988713584*k16)/(1.729827792e11L
        *Power(1+k,8));
    case 9: case 12:
      return (-2428860-44957255*k-391729618*k2-2133069839*k3-8124625314
        *k4-22951653048*k5-49735350270*k6-84299389872*k7-112868874786*k8
        -119621101971*k9-99793022932*k10-64575924815*k11-31522647148*k12
        -11031440964*k13-2499085884*k14-280459608*k15)/(1.729827792e11L
        *Power(1+k,7));
    case 10:
      return (-1902460-37240770*k-344509032*k2-2000647132*k3-8169765429
        *k4-24898644780*k5-58647407196*k6-109050726888*k7-162029911158
        *k8-193395625166*k9-185275790220*k10-141424422804*k11
        -84685684187*k12-38686708032*k13-12823761504*k14-2790332400*k15
        -305955936*k16)/(1.153218528e11L*Power(1+k,8));
    case 11: case 14:
      return (726740+13767530*k+122838492*k2+685277090*k3+2675527049*k4
        +7751921796*k5+17238945360*k6+30005596128*k7+41285419290*k8
        +45001964266*k9+38649912724*k10+25779333570*k11+12991831267*k12
        +4704191380*k13+1106148384*k14+129442896*k15)/(3.459655584e11L
        *Power(1+k,7));
    case 13:
      return (686280+12290390*k+103242832*k2+539678399*k3+1963201032*k4
        +5263731090*k5+10742317440*k6+16980495108*k7+20934558876*k8
        +20084104068*k9+14809772680*k10+8176888355*k11+3216667852*k12
        +815080164*k13+101985312*k14)/(3.459655584e11L*Power(1+k,6));
    case 15:
      return (-109760-2009490*k-17264184*k2-92338181*k3-343852002*k4
        -944210070*k5-1974504000*k6-3199820076*k7-4046732424*k8
        -3985128028*k9-3018761040*k10-1713968625*k11-694358914*k12
        -181594476*k13-23535072*k14)/(3.459655584e11L*Power(1+k,6));
    }
    break;

  case 7:
    switch (s) {
    case 0:
      return (862626681+18261257784*k+183711793145*k2+1167607479126*k3
        +5256713924622*k4+17814408610976*k5+47131137342702*k6
        +99625332659772*k7+170744417332755*k8+239347720704240*k9
        +275506401584475*k10+260305217300730*k11+200889493966260*k12
        +125343802431960*k13+62128975312140*k14+23764207478760*k15
        +6671162361120*k16+1245500971680*k17+119322815040*k18)
        /(2.498640144e12L*Power(1+k,10));
    case 1: case 4:
      return (-307185585-6329158038*k-61810716673*k2-380219937940*k3
        -1651011744864*k4-5374372352406*k5-13591797477738*k6
        -27303788884680*k7-44158151769615*k8-57904519356990*k9
        -61668857296035*k10-53153966668380*k11-36728723819490*k12
        -19996506552510*k13-8329668868920*k14-2522225138880*k15
        -502035179040*k16-50482729440*k17)/(7.495920432e12L
        *Power(1+k,9));
    case 2:
      return (-352097025-7628413233*k-78580740727*k2-511650190703*k3
        -2361137097904*k4-8206468300770*k5-22281041720544*k6
        -48364475783598*k7-85181869507695*k8-122808031681755*k9
        -145522835546625*k10-141697735653165*k11-112849942770930*k12
        -72785701684200*k13-37375874175030*k14-14854117692000*k15
        -4349837589120*k16-851705728320*k17-86177588640*k18)
        /(7.495920432e12L*Power(1+k,10));
    case 3: case 6:
      return (138850425+2922539655*k+29168925869*k2+183448437249*k3
        +814777813584*k4+2714080442478*k5+7027238403882*k6
        +14459886426870*k7+23967817361175*k8+32231190180285*k9
        +35228199089655*k10+31189065272715*k11+22161583668690*k12
        +12425909659560*k13+5341844125200*k14+1674369266940*k15
        +346519167480*k16+36459749040*k17)/(2.2487761296e13L
        *Power(1+k,9));
    case 5:
      return (44037903+881355244*k+8337462415*k2+49517470782*k3
        +206809682367*k4+644587270524*k5+1552414443630*k6+2950283512740
        *k7+4477312544055*k8+5452766950620*k9+5322321320745*k10
        +4130738072730*k11+2507965267335*k12+1157415956700*k13
        +385818726060*k14+83999703000*k15+9178678080*k16)
        /(7.495920432e12L*Power(1+k,8));
    case 7:
      return (-21337575-435776579*k-4208187237*k2-25522621329*k3
        -108894178785*k4-346855421901*k5-854042570430*k6-1660051682970
        *k7-2577839843715*k8-3214054739775*k9-3213589882635*k10
        -2556764501655*k11-1592914077285*k12-755437654545*k13
        -259358276400*k14-58363858380*k15-6629045280*k16)
        /(2.2487761296e13L*Power(1+k,8));
    case 8:
      return (794584728+16763289599*k+168017491026*k2+1063568060127*k3
        +4767352125782*k4+16078663508574*k5+42315161104902*k6
        +88925711281758*k7+151422472602900*k8+210726821417565*k9
        +240576716556570*k10+225172505357685*k11+171880654231530*k12
        +105853921884480*k13+51637945263390*k14+19355944597800*k15
        +5289881853600*k16+951234121440*k17+86177588640*k18)
        /(7.495920432e12L*Power(1+k,10));
    case 9: case 12:
      return (-278768280-5720672839*k-55626479860*k2-340572436875*k3
        -1471298189832*k4-4762643464440*k5-11970989816472*k6
        -23885465468610*k7-38340274434600*k8-49854211426665*k9
        -52592484169500*k10-44839305121785*k11-30591545781240*k12
        -16403186704950*k13-6704799369840*k14-1980884160900*k15
        -381138901800*k16-36459749040*k17)/(2.2487761296e13L
        *Power(1+k,9));
    case 10: 
      return (-35506884-766305631*k-7860804826*k2-50951346463*k3
        -233971912383*k4-808840399796*k5-2183133185343*k6-4708106800698
        *k7-8232490597230*k8-11773391549385*k9-13824143814840*k10
        -13320627963765*k11-10480137879525*k12-6661977642090*k13
        -3360605515785*k14-1305731851380*k15-371069635200*k16
        -69679859760*k17-6629045280*k18)/(2.498640144e12L*Power(1+k,10));
    case 11: case 14:
      return (13804140+289256807*k+2873086132*k2+17975117023*k3
        +79383032319*k4+262793944314*k5+675803842827*k6+1380195376350*k7
        +2268742503270*k8+3022610281245*k9+3268978245390*k10
        +2859300465765*k11+2003050703805*k12+1104073726500*k13
        +464621828265*k14+141628536120*k15+28202861280*k16+2804596080
        *k17)/(7.495920432e12L*Power(1+k,9));
    case 13:
      return (39469864+786348913*k+7402199322*k2+43728523029*k3
        +181572243966*k4+562336753671*k5+1344891807840*k6+2536257555930
        *k7+3816152271540*k8+4603174482885*k9+4444486356750*k10
        +3406576595355*k11+2038080010050*k12+923872771395*k13
        +301006884780*k14+63523744620*k15+6629045280*k16)
        /(2.2487761296e13*Power(1+k,8));
    case 15:
      return (-2096192-42596281*k-409115518*k2-2466709197*k3-10457092548
        *k4-33075613131*k5-80814991830*k6-155753215770*k7-239584946400
        *k8-295557784005*k9-291969364830*k10-229075378035*k11
        -140378602500*k12-65236808655*k13-21818304390*k14-4735987500*k15
        -509926560*k16)/(7.495920432e12L*Power(1+k,8));
    }
    break;
  }

  fprintf(stderr, "***ERROR in K5\n");
  abort();

  return 0.0L;
}


/*----------------------------------------------------------------------
  K6

  k1 = k3 = k, k2 = k4 = k+1
----------------------------------------------------------------------*/
static long double K6(int ik, int l, int s)
{
  const long double k = (long double)ik;
  const long double k2  = k*k;
  const long double k3  = k*k2;
  const long double k4  = k2*k2;
  const long double k5  = k2*k3;
  const long double k6  = k3*k3;
  const long double k7  = k3*k4;
  const long double k8  = k4*k4;
  const long double k9  = k4*k5;
  const long double k10 = k5*k5;
  const long double k11 = k5*k6;
  const long double k12 = k6*k6;
  const long double k13 = k6*k7;
  const long double k14 = k7*k7;
  const long double k15 = k7*k8;
  const long double k16 = k8*k8;
  const long double k17 = k8*k9;
  const long double k18 = k9*k9;
 
  switch (l) {
  case 0:
    switch (s) {
    case 0:
      return 715537/1074427200.0L+59599231*k/6.9837768e9L+47838487*k2
        /9.585576e8L+43787443*k3/2.494206e8L+8144699*k4/1.979208e7L
        +11985461*k5/1.786785e7L+1542161*k6/2.1021e6L+1503809*k7
        /3.15315e6L+169*k8/1225.0L;
    case 1: case 4:
      return -42631/644656320.0L-12220639*k/1.39675536e10L-258417199*k2
        /4.88864376e10L-1409595839*k3/7.33296564e10L-47973469*k4
        /1.02918816e9L-136313*k5/1.73264e6L-7267*k6/80850.0L-71177*k7
        /1.1466e6L-143*k8/7350.0L;
    case 2: case 8:
      return 63149/358142400.0L+324391*k/1.474704e8L+50642567*k2
        /4.0738698e9L+3064373479*k3/7.33296564e10L+94626559*k4
        /1.02918816e9L+116833109*k5/8.576568e8L+2498233*k6/1.89189e7L
        +86777*k7/1.1466e6L+143*k8/7350.0L;
    case 3: case 6: case 9: case 12:
      return -56303/3223281600.0L-2540143*k/1.12814856e10L-55144319*k2
        /4.19026608e10L-20358467*k3/4.4442216e9L-398477*k4/3.811808e7L
        -287887*k5/1.786785e7L-1237969*k6/7.56756e7L-375827*k7
        /3.78378e7L-121*k8/44100.0L;
    case 5:
      return 14857/1611640800.0L+62461*k/4.988412e8L+8439449*k2
        /1.08636528e10L+1925621*k3/6.636168e8L+2975891*k4/4.11675264e8L
        +827129*k5/6.59736e7L+26657*k6/1.8018e6L+14951*k7/1.4014e6L+13
        *k8/3675.0L;
    case 7: case 13:
      return 7829/3223281600.0L+11338*k/3.52546425e8+56664941*k2
        /2.933186256e11L+202846559*k3/2.933186256e11L+1671847*k4
        /1.02918816e9L+1265057*k5/4.900896e8L+136867*k6/5.04504e7L+14437
        *k7/8.4084e6L+11*k8/22050.0L;
    case 10:
      return 14543/306979200.0L+22177*k/3.83724e7L+11949251*k2
        /3.7604952e9L+47300083*k3/4.583103525e9L+44731231*k4
        /2.05837632e9L+26102827*k5/8.576568e8L+74623*k6/2.7027e6L+2953
        *k7/200200.0L+13*k8/3675.0L;
    case 11: case 14:
      return -6037/1289312640.0L-14101*k/2.387616e8L-7583987*k2
        /2.25629712e10L-41497163*k3/3.66648282e10L-5096321*k4
        /2.05837632e9L-2480041*k5/6.8612544e8L-173477*k6/5.04504e7L
        -49031*k7/2.52252e7L-11*k8/22050.0L;
    case 15: 
      return 4187/6446563200.0L+2467*k/2.930256e8L+2223733*k2
        /4.51259424e10L+8349673*k3/4.88864376e10L+792593*k4
        /2.05837632e9L+995651*k5/1.7153136e9L+173189*k6/3.027024e8L
        +3187*k7/9.45945e6L+k8/11025.0L;
    }
    break;

  case 1:
    switch (s) {
    case 0:
      return 282/870485.0L+11316983*k/2.7159132e9L+182741*k2
        /7.4613e6L+1530259*k3/1.76358e7L+279823*k4/1.36136e6L+57971*k5
        /170170.0L+5721*k6/14300.0L+343793*k7/1.05105e6L+169*k8/1225.0L;
    case 1: case 4:
      return -78973/2172730560.0L-2613407*k/5.4318264e9L-6774269*k2
        /2.3279256e9L-12343963*k3/1.1639628e9L-10817*k4/418880.0L-17091
        *k5/388960.0L-349*k6/6600.0L-16889*k7/382200.0L-143*k8/7350.0L;
    case 2: case 8:
      return 26557/296281440.0L+354283*k/3.133746e8L+7547257*k2
        /1.1639628e9L+785041*k3/3.52716e7L+25013*k4/495040.0L+643661*k5
        /8.16816e6L+151607*k6/1.8018e6L+7363*k7/127400.0L+143*k8/7350.0L;
    case 3: case 6: case 9: case 12:
      return -3617/362121760.0L-12678937*k/9.77728752e10L-2673487*k2
        /3.4918884e9L-81781*k3/3.02328e7L-61933*k4/9.801792e6L-739*k5
        /72930.0L-40063*k6/3.6036e6L-32993*k7/4.2042e6L-121*k8/44100.0L;
    case 5: 
      return 7873/1448487040.0L+75167*k/1.01846745e9L+2131841*k2
        /4.6558512e9L+362137*k3/2.116296e8L+19931*k4/4.66752e6L+263*k5
        /35360.0L+1831*k6/200200.0L+32743*k7/4.2042e6L+13*k8/3675.0L;
    case 7: case 13:
      return 1319/888844320.0L+3865153*k/1.955457504e11L+1672441*k2
        /1.39675536e10L+2020693*k3/4.6558512e9L+3397*k4/3.267264e6L+859
        *k5/502656.0L+13837*k6/7.2072e6L+34981*k7/2.52252e7L+11*k8
        /22050.0L;
    case 10:
      return 1553/62078016.0L+1263343*k/4.0738698e9L+64341*k2
        /3.69512e7L+11843*k3/2.0349e6L+96251*k4/7.53984e6L+465853*k5
        /2.450448e7L+17137*k6/900900.0L+7129*k7/600600.0L+13*k8/3675.0L;
    case 11: case 14:
      return -547/197520960.0L-1151873*k/3.25909584e10L-407773*k2
        /1.9953648e9L-7873*k3/1.119195e7L-3461*k4/2.178176e6L-79601*k5
        /3.267264e7L-6029*k6/2.4024e6L-13567*k7/8.4084e6L-11*k8/22050.0L;
    case 15:
      return 41/100279872.0L+15397*k/2.8756728e9L+11527*k2/3.627936e8L
        +20101*k3/1.790712e8L+5101*k4/1.9603584e7L+2229*k5/5.44544e6L
        +6233*k6/1.44144e7L+899*k7/3.15315e6L+k8/11025.0L;
    }
    break;

  case 2:
    switch (s) {
    case 0:
      return (734029128+k*(205392096+k*(40948614+k*(-149167648+k*(
        127933106+k*(-119541744+k*(148175431+2*k*(-57678560+k*(204123693
        +k*(321334666+19*k*(57909809+4*k*(20024406+17*k*(1622927+2*k*(
        558949+435006*k))))))))))))))/(1.62954792e10L*k6)-(k14*(18036+k
        *(24032+k*(6263+2*k*(-1278+k*(-389+k*(5+2*k)*(46+k*(17+2*k))))))
        )*LOG20(ik))/20020.0L;
    case 1: case 4:
      return (-1468058256-k*(899811360+k*(-106917930+k*(-143496304+k*(
        53155788+k*(-21269304+k*(65752025+k*(-32823192+k*(260793258+k*(
        579511876+19*k*(91808193+2*k*(69284925+68*k*(1363154+33*k*(30157
        +22308*k))))))))))))))/(9.77728752e10L*k6)+(k14*(1+k)*(36072+k*(
        24008+k*(-3473+k*(-2755+k*(281+k*(667+218*k+24*k2))))))
        *LOG20(ik))/120120.0L;
    case 2: case 8:
      return (-404375328+k*(-251719272+k*(-56407428+k*(115585512+k*(
        -4990104+k*(66983812+k*(11915790+k*(449132340+k*(1270725644+57*k
        *(55316041+2*k*(42130317+68*k*(711168+11*k*(45757+22308*k)))))))
        ))))))/(9.77728752e10L*k5)-(k15*(-9936+k*(-16644+k*(-7866+k*(
        1364+k*(2892+k*(1347+286*k+24*k2))))))*LOG20(ik))/120120.0L;
    case 3: case 6: case 9: case 12:
      return (134791776-k*(-128511432+k*(-16658334+k*(29920996+k*(
        15135960+k*(8504944+k*(5211539+k*(54337070+k*(172302042+19*k*(
        22317465+4*k*(8837088+17*k*(600261+436094*k+207636*k2)))))))))))
        )/(9.77728752e10L*k5)+(k15*(1+k)*(-3312+k*(-3332+k*(-391+k*(693
        +k*(457+12*k*(10+k))))))*LOG20(ik))/120120.0L;
    case 5:
      return (978705504+k*(925892352+k*(-88505298+k*(-106806728+k*(
        11061792+k*(12114144+k*(19408993+k*(-2823696+k*(78077010+k*(
        206011244+19*k*(31338867+4*k*(12476511+68*k*(244433+9*k*(20633
        +14872*k))))))))))))))/(1.955457504e11L*k6)-(k14*Power(1+k,2)*(
        24048+k*(-32+k*(-2285+k*(-294+k*(367+4*k*(43+6*k))))))*LOG20(ik)
        )/240240.0L;
    case 7: case 13:
      return (-89861184+k*(-115410960+k*(-19565742+k*(21299964+k*(
        14582400+k*(5066208+k*(2408159+k*(17317654+k*(58470458+57*k*(
        2550490+k*(4164727+68*k*(71723+53302*k+25168*k2))))))))))))/(
        1.955457504e11L*k5)-(k15*Power(1+k,2)*(-2208+k*(-744+k*(237+k*(
        288+k*(97+12*k)))))*LOG20(ik))/240240.0L;
    case 10:
      return (142117416+k*(190312416+k*(145064808+k*(44672544+k*(
        23895620+k*(28043904+k*(225757788+k*(711536416+57*k*(28639489
        +12*k*(3599849+68*k*(56322+k*(37793+14872*k))))))))))))/(
        1.955457504e11L*k4)-(k16*(3492+k*(8352+k*(8476+k*(4848+k*(1659+4
        *k*(77+6*k))))))*LOG20(ik))/240240.0L;
    case 11: case 14:
      return (-47372472-k*(78577128+k*(57352764+k*(22722000+k*(6051780+k
        *(4489716+k*(28146754+k*(94107808+57*k*(3837935+k*(5960845+68*k
        *(94613+64742*k+25168*k2)))))))))))/(1.955457504e11L*k4)+(k16*(1
        +k)*(1164+k*(1992+k*(1446+k*(594+k*(131+12*k)))))*LOG20(ik))
        /240240.0L;
    case 15:
      return (15790824+k*(31238928+k*(23776438+k*(10011120+k*(2448124+k
        *(909832+k*(4525285+k*(15727236+19*k*(1962285+2*k*(1564233+17*k
        *(101181+32*k*(2207+858*k))))))))))))/(1.955457504e11L*k4)-(k16
        *Power(1+k,2)*(388+k*(400+k*(203+6*k*(9+k))))*LOG20(ik))
        /240240.0L;
    }
    break;

  case 3:
    switch (s) {
    case 0:
      return (21998896920+k*(69727497840+k*(69422770494+k*(18460087800
        +k*(-5120309964+k*(-1378984992+k*(771050273+k*(740895680+k*(
        1381243615+2*k*(2717061445+k*(7766140265+2*k*(8082442760+19*k*(
        653495459+12*k*(62536610+17*k*(3142469+2*k*(956531+435006*k)))))
        )))))))))))/(9.77728752e10L*k6*Power(1+k,2))-(k14*(835380+k*(
        1856400+k*(872661+2*k*(-65421+k*(-34819+10*k*(670+k*(837+20*k
        *(11+k))))))))*LOG20(ik))/185640.0L;
    case 1: case 4:
      return (-21998896920-k*(59950210320+k*(33272536374+k*(-3969615804
        +k*(-2843734740+k*(571237128+k*(247761465+k*(555798810+k*(
        93443955+k*(2296053670+k*(4410336625+57*k*(176637655+6*k*(
        35862307+68*k*(602106+121*k*(2905+2028*k)))))))))))))))/(
        2.933186256e11L*k6*(1+k))+(k14*(835380+k*(2320500+k*(1336761+k*(
        -161487+k*(-104418+25*k*(480+k*(837+242*k+24*k2)))))))*LOG20(ik)
        )/556920.0L;
    case 2: case 8:
      return (-5076668520+k*(-17553525066+k*(-20645405550+k*(-8505506394
        +k*(1222823448+k*(2365167833+k*(1164088135+k*(1270097615+k*(
        4488204545+k*(12609886495+k*(25533336680+57*k*(665478027+2*k*(
        364383515+68*k*(4248203+33*k*(69863+22308*k)))))))))))))))/(
        2.933186256e11L*k5*Power(1+k,2))+(k15*(192780+k*(483939-k*(
        -324387+k*(-14168+25*k*(2400+k*(1299+286*k+24*k2))))))
        *LOG20(ik))/556920.0L;
    case 3: case 6: case 9: case 12:
      return (1692222840-k*(-5099075982+k*(-3882495078+k*(-317447130+k*(
        639102156+k*(340980115+k*(120111775+3*k*(38701465+k*(174207495+k
        *(427616075+19*k*(43053995+4*k*(13921391+17*k*(793097+485950*k
        +207636*k2)))))))))))))/(2.933186256e11L*k5*(1+k))+(k15*(-64260
        +k*(-197013+k*(-154044+k*(-12841+15*k*(1630+k*(933+20*k*(11+k)))
        ))))*LOG20(ik))/556920.0L;
    case 5:
      return (14665931280+k*(33448615200+k*(525856716+k*(-3877037472+k*(
        1312420032+k*(-1177046640+k*(1641517290+k*(-1489164600+k*(
        1661978670+k*(-1000228300+57*k*(47931175+12*k*(1158265+68*k*(
        84973+3*k*(8523+14872*k))))))))))))))/(5.866372512e11L*k6)-(k14
        *(556920+k*(1856400+k*(1355274+k*(-128088+k*(-103002+5*k*(1140+k
        *(3507+20*k*(55+6*k))))))))*LOG20(ik))/1.11384e6L;
    case 7: case 13:
      return (-1128148560+k*(-2897984628+k*(-776816964+k*(355362084+k*(
        154539000+k*(155003170+k*(-73127670+k*(124422210+k*(48638150+57
        *k*(7257950+3*k*(2718835+68*k*(61459+36642*k+25168*k2)))))))))))
        )/(5.866372512e11L*k5)+(k15*(42840+k*(155142-k*(-145206+k*(
        -17374+25*k*(798+k*(486+k*(121+12*k)))))))*LOG20(ik))/1.11384e6L;
    case 10:
      return (1462564026+k*(5919037740+k*(9421361334+k*(7867118952+k*(
        3931113907+k*(1312559150+k*(817297135+k*(2502192880+k*(
        6890691995+k*(13636713730+57*k*(345301703+4*k*(90914465+68*k*(
        1002296+9*k*(55427+14872*k))))))))))))))/(5.866372512e11L*k4
        *Power(1+k,2))-(k16*(55539+k*(172152+k*(194348+5*k*(22800+k*(
        8127+20*k*(77+6*k))))))*LOG20(ik))/1.11384e6L;
    case 11: case 14:
      return (-487521342-k*(1754493048+k*(2129124690+k*(1268181684+k*(
        451646825+k*(104099240+k*(75653325+k*(277224310+k*(697682425+57
        *k*(22276940+3*k*(9434957+68*k*(125571+73250*k+25168*k2)))))))))
        )))/(5.866372512e11L*k4*(1+k))+(k16*(18513+k*(67599+k*(83326+25
        *k*(1992+k*(717+k*(143+12*k))))))*LOG20(ik))/1.11384e6L;
    case 15:
      return (162507114+k*(511991172+k*(400553538+k*(156063600+k*(
        45679795+3*k*(-798490+k*(6052035+k*(9529300+19*k*(1653025+2*k*(
        1158165+17*k*(84053+32*k*(1717+858*k))))))))))))/(5.866372512e11
        *k4)-(k16*(6171+k*(25938+k*(35657+15*k*(1460+k*(531+10*k*(11+k))
        ))))*LOG20(ik))/1.11384e6L;
    }
    break;

  case 4:
    switch (s) {
    case 0:
      return (1047566520+k*(10071341280+k*(30553717194+k*(41328547260+k
        *(25323764466+k*(4440740304+k*(-1387466301+k*(-127684770+k*(
        617390944+k*(2083335192+k*(6713782623+2*k*(8312012701+k*(
        15625488621+2*k*(11313399798+19*k*(667360369+12*k*(48106898+17*k
        *(1855369+6*k*(150457+48334*k))))))))))))))))))/(3.25909584e10L
        *k6*Power(1+k,4))-(k14*(3+20*k*(1+k))*(19380+k2*(-969+2*k*(-57+k
        *(57+2*k*(10+k)))))*LOG20(ik))/90440.0L;
    case 1: case 4:
      return (-1047566520-k*(10769718960+k*(28421827434+k*(28175563416+k
        *(8497018530+k*(-1544118576+k*(-494844753+k*(164406276+k*(
        270090140+k*(752830310+k*(2333349900+k*(5401008905+k*(9371330039
        +57*k*(215834333+2*k*(107065275+68*k*(1169437+33*k*(18687+7436*k
        )))))))))))))))))/(9.77728752e10L*k6*Power(1+k,3))+(k14*(58140+k
        *(484500+k*(578493+k*(-24567+k*(-31578+5*k*(-144+k*(773+242*k+24
        *k2)))))))*LOG20(ik))/271320.0L;
    case 2: case 8:
      return (-209513304+k*(-2066646582+k*(-6615542934+k*(-9792662880+k
        *(-7059362310+k*(-1863395079+k*(595728523+k*(896812446+k*(
        1856170530+k*(5595213260+k*(13579588701+k*(24992733918+k*(
        35221911196+57*k*(668368145+2*k*(274726925+68*k*(2454238+33*k*(
        31323+7436*k)))))))))))))))))/(9.77728752e10L*k5*Power(1+k,4))+(
        k15*(11628+k*(80427-k*(-96387+k*(-15048+5*k*(1744+k*(1235+286*k
        +24*k2))))))*LOG20(ik))/271320.0L;
    case 3: case 6: case 9: case 12:
      return (209513304-k*(-2206322118+k*(-6224083866+k*(-7056293244+k*(
        -3038441406+k*(130677729+k*(528568157+k*(343297016+k*(644423742
        +k*(1900965710+k*(4313927534+9*k*(810833239+19*k*(54215073+4*k*(
        12860051+17*k*(523133+247814*k+69212*k2)))))))))))))))/(
        2.933186256e11L*k5*Power(1+k,3))+(k15*(-11628+k*(-99807+k*(
        -139992+k*(-23883+k*(11230+k*(8053+180*k*(11+k)))))))*LOG20(ik))
        /813960.0L;
    case 5:
      return (698377680+k*(7645397760+k*(17061039996+k*(9573972408+k*(
        -876743868+k*(-580515936+k*(59642466+k*(53922844+k*(103212652+k
        *(173193480+k*(684062155+k*(1310983618+57*k*(39411761+4*k*(
        10825977+68*k*(148933+78471*k+44616*k2)))))))))))))))/(
        1.955457504e11L*k6*Power(1+k,2))-(k14*(38760+k*(387600+k*(579462
        +k*(-19608+k*(-31122+k*(-1420+k*(3283+20*k*(55+6*k))))))))
        *LOG20(ik))/542640.0L;
    case 7: case 13:
      return (-139675536+k*(-1563998436+k*(-3795299508+k*(-2768453688+k
        *(-257693436+k*(298730614+k*(156608530+k*(78823080+k*(170203000 
        +k*(514823500+3*k*(354351613+57*k*(9704684+k*(10826167+68*k*(
        131199+70318*k+25168*k2))))))))))))))/(5.866372512e11L*k5
        *Power(1+k,2))+(k15*(7752+k*(79458-k*(-135318+k*(-25042+5*k*(
        1946+k*(1426+363*k+36*k2))))))*LOG20(ik))/1.62792e6L;
    case 10:
      return (52378326+k*(540540000+k*(1910412504+k*(3374507136+k*(
        3419812851+k*(2155154820+k*(1053453986+k*(1166837908+k*(
        3140291642+k*(7468228088+k*(13488650578+k*(18572475068+57*k*(
        342077993+4*k*(67536243+68*k*(569746+3*k*(73061+14872*k)))))))))
        )))))))/(1.955457504e11L*k4*Power(1+k,4))-(k16*(2907+k*(21432+k
        *(34428+k*(20880+k*(7903+20*k*(77+6*k))))))*LOG20(ik))
        /542640.0L;
    case 11: case 14:
      return (-52378326-k*(575458884+k*(1828466640+k*(2640333696+k*(
        2022321483+k*(916864494+k*(334801870+k*(379288080+k*(1043108550
        +k*(2322691780+3*k*(1281294752+57*k*(27776194+k*(25399825+68*k*(
        245007+2*k*(53463+12584*k)))))))))))))))/(5.866372512e11L*k4
        *Power(1+k,3))+(k16*(2907+k*(26277+k*(47538+5*k*(5784+k*(2119
        +429*k+36*k2)))))*LOG20(ik))/1.62792e6L;
    case 15:
      return (17459442+k*(203459256+k*(570533964+k*(616143528+k*(
        317703477+k*(101607688+k*(30017364+k*(48421560+k*(137426410+9*k
        *(31342824+19*k*(2486289+2*k*(1352241+17*k*(61831+32*k*(981+286
        *k))))))))))))))/(5.866372512e11L*k4*Power(1+k,2))-(k16*(969+k*(
        10374+k*(21831+k*(13460+k*(4811+90*k*(11+k))))))*LOG20(ik))
        /1.62792e6L;
    }
    break;

  case 5:
    switch (s) {
    case 0:
      return (6721885170+k*(33255642420+k*(65579423910+k*(63875251440+k
        *(29645173961+k*(4445597902+k*(5312172267+k*(36298828640+k*(
        141319822969+k*(407247653250+k*(903350470249+22*k*(71480568794+k
        *(98714842020+k*(108547794988+19*k*(4988342147+132*k*(25943622
        +17*k*(808531+318490*k+79092*k2)))))))))))))))))
        /(5.377508136e11L*k4*Power(1+k,6))-(k16*(373065+2*k2*(-7315+k*(
        -770+k*(693+20*k*(11+k)))))*LOG20(ik))/1.49226e6L;
    case 1: case 4:
      return (-20165655510-k*(79601271750+k*(117136999980+k*(74488754340
        +k*(14350477635+k*(-2638610133+k*(4736699232+k*(28712540528+k*(
        105238116880+k*(286741378670+11*k*(54220068958+k*(87517278248+k
        *(110475628057+57*k*(1911259215+22*k*(66326995+68*k*(564556+3*k
        *(80167+22308*k)))))))))))))))))/(3.2265048816e12L*k4
        *Power(1+k,5))+(k16*(223839+k2*(-8778+k*(-924+k*(693+242*k+24
        *k2))))*LOG20(ik))/1.790712e6L;
    case 2: case 8:
      return (-2372430060+k*(-12264492240+k*(-25763758020+k*(
        -27593459972+k*(-13870461930+k*(8467992144+k*(63197995830+k*(
        239000249253+k*(678612821562+k*(1482726761295+11*k*(230493553896
        +k*(311399365005+k*(332999154032+57*k*(4914930259+66*k*(48495435
        +68*k*(347789+k*(118075+22308*k)))))))))))))))))
        /(3.2265048816e12L*k3*Power(1+k,6))-(k17*(-26334+k*(-5852+k*(924
        +k*(1155+286*k+24*k2))))*LOG20(ik))/1.790712e6L;
    case 3: case 6: case 9: case 12:
      return (1186215030-k*(-4946031090+k*(-7935847920+k*(-5891717728+k
        *(-1530250631+k*(1736810279+k*(8225051126+k*(29069714215+k*(
        77925174134+33*k*(4826296427+k*(7627689262+11*k*(852423663+19*k
        *(42735655+4*k*(7768409+17*k*(245519+2*k*(45497+9438*k))))))))))
        ))))))/(3.2265048816e12L*k3*Power(1+k,5))+(k17*(-65835+k*(-14630
        +3*k*(770+k*(847+20*k*(11+k)))))*LOG20(ik))/8.95356e6L;
    case 5:
      return (20165655510+k*(59435616240+k*(57701383740+k*(16787370600
        +k*(-2468101857+k*(-616329300+k*(1832487288+k*(8654320264+k*(
        30092128186+11*k*(7026204744+k*(13651139768+k*(20290727624+57*k
        *(407006287+44*k*(8048777+68*k*(77293+9*k*(3981+1352*k))))))))))
        ))))))/(6.4530097632e12L*k4*Power(1+k,4))-(k16*(1119195+k2*(
        -43890+k*(-4620+k*(3003+20*k*(55+6*k)))))*LOG20(ik))
        /1.790712e7L;
    case 7: case 13:
      return (-1186215030+k*(-3759816060+k*(-4176031860+k*(-1726467496+k
        *(60879377+k*(662851494+k*(2424869940+k*(8157824130+11*k*(
        1873136226+k*(3564682478+k*(5170315196+57*k*(100367090+33*k*(
        2534745+68*k*(22713+9454*k+2288*k2))))))))))))))/(
        6.4530097632e12L*k3*Power(1+k,4))-(k17*(-13167+k*(-2926+k*(462+k
        *(462+k*(121+12*k)))))*LOG20(ik))/3.581424e6L;
    case 10:
      return (527206680+k*(2941258320+k*(6989313968+k*(9832077216+k*(
        13770390036+k*(39188710200+k*(136095731097+k*(379942193130+k*(
        818809126167+11*k*(125288211264+k*(166073808165+k*(173457398038
        +57*k*(2483935999+132*k*(11766649+68*k*(79628+3*k*(8245+1352*k))
        ))))))))))))))/(6.4530097632e12L*k2*Power(1+k,6))-(k18*(29260+k
        *(18480+k*(7623+20*k*(77+6*k))))*LOG20(ik))/1.790712e7L;
    case 11: case 14:
      return (-263603340-k*(1207025820+k*(2266903600+k*(2351148176+k*(
        2167784086+k*(4938940254+k*(16215350055+k*(42740431860+11*k*(
        7818060379+k*(12130084154+k*(14578363651+57*k*(236751800+33*k*(
        5021335+68*k*(37577+22*k*(581+104*k)))))))))))))))/(
        6.4530097632e12L*k2*Power(1+k,5))+(k18*(2926+k*(1848+k*(693+k*(
        143+12*k))))*LOG20(ik))/3.581424e6L;
    case 15:
      return (131801670+k*(471711240+k*(657737132+k*(475141160+k*(
        313185327+k*(697483016+k*(2233171289+33*k*(168127222+k*(
        314279609+k*(446139052+19*k*(25287613+22*k*(924651+17*k*(31511
        +32*k*(379+78*k))))))))))))))/(6.4530097632e12L*k2*Power(1+k,4)) 
        -(k18*(7315+3*k*(1540+k*(539+10*k*(11+k))))*LOG20(ik))
        /1.790712e7L;
    }
    break;

  case 6:
    switch (s) {
    case 0:
      return (547/9152528.0L+364943*k/2.941884e8L+44002153*k2
        /3.6038079e9L+545854279*k3/7.2076158e9L+4782515267*k4
        /1.44152316e10L+714355009*k5/6.552378e8L+20076096541*k6
        /7.2076158e9L+3398348237*k7/6.0063465e8L+113847059*k8/1.22892e7L
        +233978267*k9/1.89924e7L+1176589*k10/88200.0L+3678354773*k11
        /3.133746e8L+274958261*k12/3.29868e7L+19450663*k13/4.12335e6L
        +1009667*k14/485100.0L+165329*k15/242550.0L+169*k16/1225.0L)
        /Power(1+k,8);
    case 1: case 4:
      return (-97/10767680.0L-1493893*k/8.2372752e9L-19903931*k2
        /1.153218528e10L-594663469*k3/5.76609264e10L-1249391581*k4
        /2.88304632e10L-1958862313*k5/1.44152316e10L-13876413*k6
        /4.21498e7L-2629273*k7/4.178328e6L-72952009*k8/7.59696e7L
        -268001291*k9/2.279088e8L-15478163*k10/1.34064e7L-119120233*k11
        /1.319472e8L-7319779*k12/1.319472e7L-254147*k13/970200.0L-8021
        *k14/88200.0L-143*k15/7350.0L)/Power(1+k,7);
    case 2: case 8:
      return (43229/2471182560.0L+ (15601417*k)/4.32456948e10+152469437
        *k2/4.32456948e10L+938252213*k3/4.32456948e10L+16296241873*k4
        /1.729827792e11L+26509866421*k5/8.64913896e10L+11160910829*k6
        /1.44152316e10L+2797047877*k7/1.80190395e9L+1568757803*k8
        /6.267492e8L+139663639*k9/4.27329e7L+295754467*k10/8.54658e7L
        +327289783*k11/1.106028e8L+800582999*k12/3.958416e8L+214419683
        *k13/1.979208e8L+319061*k14/727650.0L+10937*k15/88200.0L+143*k16
        /7350.0L)/Power(1+k,8);
    case 3: case 6: case 9: case 12:
      return (-917/353026080.0L-8979811*k/1.729827792e11L-2421943*k2
        /4.94236512e9L-501954689*k3/1.729827792e11L-1043951291*k4
        /8.64913896e10L-42579991*k5/1.1380446e9L-80598407*k6
        /9.00951975e8L-2110853*k7/1.2534984e7L-12013*k8/47600.0L
        -1737149*k9/5.7456e6L-197452361*k10/6.837264e8L-86199053*k11
        /3.958416e8L-2515357*k12/1.979208e7L-322079*k13/5.8212e6L-48401
        *k14/2.9106e6L-121*k15/44100.0L)/Power(1+k,7);
    case 5:
      return (5099/3294910080.0L+(3947*k)/1.307504e8L+10665407*k2
        /3.84406176e10L+3065449*k3/1.92203088e9L+10278559*k4
        /1.6016924e9L+1315217*k5/6.864396e7L+873291*k6/1.98968e7L+746741
        *k7/9.4962e6L+482893*k8/4.34112e6L+9454969*k9/7.59696e7L+175177
        *k10/1.59936e6L+1645031*k11/2.19912e7L+1789*k12/46200.0L+9*k13
        /616.0L+13*k14/3675.0L)/Power(1+k,6);
    case 7: case 13:
      return (11/24961440.0L+57839*k/6.7836384e9L+559717*k2/7.2076158e9L
        +185305*k3/4.19352192e8L+2815293*k4/1.6016924e9L+439217*k5
        /8.479548e7L+21419*k6/1.8326e6L+2344319*k7/1.139544e8L+64931*k8
        /2.279088e6L+14152603*k9/4.558176e8L+15839*k10/599760.0L+412649
        *k11/2.39904e7L+32107*k12/3.8808e6L+353*k13/129360.0L+11*k14
        /22050.0L)/Power(1+k,6);
    case 10:
      return (2309/449305920.0L+162703*k/1.5444891e9L+88454713*k2
        /8.64913896e10L+12280673*k3/1.9657134e9L+9309053777*k4
        /3.459655584e11L+7503035677*k5/8.64913896e10L+6253040141*k6
        /2.88304632e10L+885030281*k7/2.0593188e9L+214242283*k8
        /3.133746e8L+1532429*k9/1.7442e6L+155998553*k10/1.709316e8L
        +130283677*k11/1.709316e8L+400117951*k12/7.916832e8L+51446501
        *k13/1.979208e8L+1877*k14/18900.0L+8333*k15/323400.0L+13*k16
        /3675.0L)/Power(1+k,8);
    case 11: case 14:
      return (-59/78450240.0L-8783*k/5.883768e8L-537611*k2/3.84406176e9L
        -2785583*k3/3.3918192e9L-11827247*k4/3.4946016e9L-149659847*k5
        /1.44152316e10L-707238593*k6/2.88304632e10L-1903369*k7
        /4.178328e7L-1914391*k8/2.84886e7L-500169*k9/6.3308e6L-5603033
        *k10/7.59696e7L-792023*k11/1.46608e7L-320189*k12/1.0555776e7L
        -16169*k13/1.2936e6L-6703*k14/1.9404e6L-11*k15/22050.0L)
        /Power(1+k,7);
    case 15:
      return (89/706052160.0L+787*k/3.251556e8L+2520929*k2
        /1.153218528e11L+387577*k3/3.14514144e9L+3814039*k4/7.8628536e9L
        +2036347*k5/1.44152316e9L+3942017*k6/1.2534984e9L+34501*k7
        /6.3308e6L+4825*k8/651168.0L+675043*k9/8.54658e7L+5533*k10
        /846720.0L+541489*k11/1.319472e8L+43649*k12/2.32848e7L+83*k13
        /145530.0L+k14/11025.0L)/Power(1+k,6);
    }
    break;

  case 7:
    switch (s) {
    case 0:
      return (2989803/59491432000.0L+44721*k/3.91391e7L+15466465337*k2
        /1.249320072e12L+17698243889*k3/2.08220012e11L+64935653*k4
        /1.565564e8L+238815619729*k5/1.56165009e11L+7734431541*k6
        /1.749748e9L+1067908215637*k7/1.04110006e11L+17015140657*k8
        /8.7671584e8L+314778973671*k9/1.04110006e10L+141211174931*k10
        /3.6212176e9L+5809471097*k11/1.392776e8L+66890794473*k12
        /1.8106088e9L+217898101*k13/8.083075e6L+765203801*k14/4.76476e7L
        +137556011*k15/1.786785e7L+6087199*k16/2.1021e6L+121279*k17
        /150150.0L+169*k18/1225.0L)/Power(1+k,10);
    case 1: case 4:
      return (-12031/1551950400.0L-637573*k/3.7182145e9L-6757287853*k2
        /3.747960216e12L-975905869*k3/8.1477396e10L-10059723763*k4
        /1.78474296e11L-22671025499*k5/1.13574552e11L-172287290563*k6
        /3.12330018e11L-230355341*k7/1.8929092e8L-543259899979*k8
        /2.498640144e11L-202711591*k9/6.390384e7L-249703707*k10
        /6.584032e7L-2881160521*k11/7.759752e8L-32175714611*k12
        /1.08636528e10L-14534957*k13/7.623616e6L-278782481*k14
        /2.858856e8L-2435527*k15/6.3063e6L-42793*k16/382200.0L-143*k17
        /7350.0L)/Power(1+k,9);
    case 2: case 8:
      return (877077/59491432000.0L+19218613*k/5.76609264e10+286053*k2
        /7.970144e7L+30578360213*k3/1.249320072e12L+444341912153*k4
        /3.747960216e12L+270812352983*k5/6.24660036e11L+67449265361*k6
        /5.4318264e10L+55675321798*k7/1.9520626125e10L+88853404347*k8
        /1.665760096e10L+292697429407*k9/3.56948592e10L+113121182089*k10
        /1.08636528e10L+808025771*k11/7.39024e7L+102837837313*k12
        /1.08636528e10L+276095377*k13/4.11502e7L+56348133*k14/1.46608e7L
        +497785219*k15/2.858856e8L+2522791*k16/4.2042e6L+55429*k17
        /382200.0L+143*k18/7350.0L)/Power(1+k,10);
    case 3: case 6: case 9: case 12:
      return (-10429/4655851200.0L-553517953*k/1.1243880648e13L
        -5779014721*k2/1.1243880648e13L-4232974549*k3/1.249320072e12L
        -234196351*k4/1.4814072e10L-41571705341*k5/7.495920432e11L
        -94726185293*k6/6.24660036e11L-7288209349*k7/2.20468248e10L
        -3670154927*k8/6.2990928e9L-478465157*k9/5.717712e8L-2909809597
        *k10/2.9628144e9L-30639845723*k11/3.25909584e10L-8122593*k12
        /1.1142208e7L-775150109*k13/1.7153136e9L-93892349*k14
        /4.288284e8L-669817*k15/8.4084e6L-252653*k16/1.26126e7L-121*k17
        /44100.0L)/Power(1+k,9);
    case 5:
      return (19997/14872858000.0L+3869317*k/1.33855722e11L+1104294427
        *k2/3.747960216e12L+1181699039*k3/6.24660036e11L+931317287*k4
        /1.08636528e11L+6061689521*k5/2.08220012e11L+63081947*k6
        /8.219211e8L+626753132*k7/3.904125225e9L+977313761*k8
        /3.6212176e9L+90618421*k9/2.469012e8L+133051307*k10/3.292016e8L
        +1949519333*k11/5.4318264e9L+291103679*k12/1.1435424e9L+8113159
        *k13/5.717712e7L+763901*k14/1.26126e7L+26353*k15/1.4014e6L+13
        *k16/3675.0L)/Power(1+k,8);
    case 7: case 13:
      return (616159/1606268664000.0L+46032551*k/5.621940324e12L+1404269
        *k2/1.6959096e10L+1975350991*k3/3.747960216e12L+805995899*k4
        /3.40723656e11L+348379263*k5/4.3835792e10L+775874081*k6
        /3.747960216e10L+142934747*k7/3.34639305e9L+767941403*k8
        /1.08636528e10L+279644431*k9/2.9628144e9L+301253599*k10
        /2.9628144e9L+318175841*k11/3.6212176e9L+103150097*k12
        /1.7153136e9L+109327891*k13/3.4306272e9L+10061*k14/800800.0L
        +12239*k15/3.6036e6L+11*k16/22050.0L)/Power(1+k,8);
    case 10:
      return (221327/50992656000.0L+17420537*k/1.78474296e11L+7833497293
        *k2/7.495920432e12L+6636009067*k3/9.36990054e11L+255389886457*k4
        /7.495920432e12L+154484117101*k5/1.249320072e12L+877614350819*k6
        /2.498640144e12L+166430960481*k7/2.08220012e11L+2688352517*k8
        /1.817192832e9L+561583591243*k9/2.498640144e11L+20391827259*k10
        /7.2424352e9L+564494927*k11/1.939938e8L+53706382189*k12
        /2.17273056e10L+18591163631*k13/1.08636528e10L+1087784063*k14
        /1.1435424e9L+118326889*k15/2.858856e8L+1816*k16/13475.0L+125963
        *k17/4.2042e6L+13*k18/3675.0L)/Power(1+k,10);
    case 11: case 14:
      return (-2017/3103900800.0L-1900169*k/1.33855722e11L-3311150681*k2
        /2.2487761296e13L-3612064007*k3/3.747960216e12L-586728203*k4
        /1.31507376e11L-29077911971*k5/1.873980108e12L-105000424271*k6
        /2.498640144e12L-5661466777*k7/6.24660036e10L-236608988783*k8
        /1.4991840864e12L-202403287*k9/9.053044e8L-39237829*k10
        /1.519392e8L-465077671*k11/1.9171152e9L-3994325927*k12
        /2.17273056e10L-10548547*k13/9.52952e7L-177001933*k14
        /3.4306272e9L-298913*k15/1.68168e7L-11553*k16/2.8028e6L-11*k17
        /22050.0L)/Power(1+k,9);
    case 15:
      return (117707/1070845776000.0L+13613*k/5.84097696e9L+7015301*k2
        /2.9983681728e11L+3042343*k3/2.0593188e10L+3184733*k4
        /4.845456e9L+1170831041*k5/5.35422888e11L+1055897033*k6
        /1.873980108e11L+56638759*k7/4.9315266e9L+71879851*k8
        /3.8342304e9L+72995347*k9/2.9628144e9L+128971*k10/4.9504e6L
        +51251149*k11/2.3279256e9L+652831*k12/4.45536e7L+2567603*k13
        /3.4306272e8L+282647*k14/1.009008e8L+2207*k15/3.15315e6L+k16
        /11025.0L)/Power(1+k,8);
    }
    break;
  }

  fprintf(stderr, "***ERROR in K6\n");
  abort();
 
  return 0.0L;
}


/*----------------------------------------------------------------------
  Log10

  Log10(k) = Log(1+1/k) - Sum_{n=1}^{n=10} -(-k)^(-n)/n
----------------------------------------------------------------------*/
static long double Log10_p(int ik)
{
  const long double k = (long double)ik;
  
  switch (ik) {
  case  1: return  0.047512259925024674497L;
  case  2: return  0.000030460290704064517696L;
  case  3: return  3.9323298846371458031e-7L;
  case  4: return  1.7637820091827803027e-8L;
  case  5: return  1.5736371658942577077e-9L;
  case  6: return  2.1739635751582420319e-10L;
  case  7: return  4.0656117927529053413e-11L;
  case  8: return  9.495928859220929287e-12L;
  case  9: return  2.6293159473253921615e-12L;
  case -2: return -0.000082324409151658623581L;
  case -3: return -7.3993612969569245044e-7L;
  case -4: return -2.8134929023917393609e-8L;
  case -5: return -2.2805589621155014395e-9L;
  case -6: return -2.9583038879468987501e-10L;
  case -7: return -5.2911852516572046791e-11L;
  case -8: return -1.1954180362685154289e-11L;
  case -9: return -3.2257464786154677734e-12L;
  }
  
  return (1/11.0L-(1/12.0L-(1/13.0L-(1/14.0L-(1/15.0L-(1/16.0L-(1/17.0L
    -(1/18.0L-(1/19.0L-(1/20.0L-(1/21.0L-(1/22.0L-(1/23.0L-(1/24.0L
    -(1/25.0L-(1/26.0L-(1/27.0L-(1/28.0L-(1/29.0L-(1/30.0L/k)/k)/k)/k)
    /k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/k)/Power(k,11);
}


/*----------------------------------------------------------------------
  Log20

  Log20(k) = Log(1+1/k) - Sum_{n=1}^{n=20} -(-k)^(-n)/n
----------------------------------------------------------------------*/
static long double Log20_p(int ik)
{
  const long double k = (long double)ik;

  switch (ik) {
  case  1: return 0.24375777384517366187e-1L;
  case  2: return 0.15373987349227200010e-7L;
  case  3: return 0.34539215997191893165e-11L;
  case  4: return 0.87420032720118016847e-14L;
  case  5: return 0.83860107068515694217e-16L;
  case  6: return 0.18728546595548083514e-17L;
  case  7: return 0.75027005129753430217e-19L;
  case  8: return 0.46126205086322442056e-20L;
  case  9: return 0.39347524476167441189e-21L;
  case -2: return -0.43508916373144593619e-7L;
  case -3: return -0.66796680783346570899e-11L;
  case -4: return -0.14223779326883030407e-13l;
  case -5: return -0.12344185873684561726e-15L;
  case -6: return -0.25815890775266002772e-17L;
  case -7: return -0.98721805194977473078e-19L;
  case -8: return -0.58625732161813776697e-20L;
  case -9: return -0.48684696965060483753e-21L;
  }

  return -(133573286426580-(136998242488800-(140603459396400
    -(144403552893600-(148414762696200-(152655184487520
    -(157145042854800-(161907013850400-(166966608033225
    -(172352627647200-(178097715235440-(184239015760800
    -(190818980609400-(197886350261600-(205497363733200
    -(213717258282528-(222622144044300-(232301367698400
    -(242860520775600-254425307479200*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)*k)
    *k)*k)*k)*k)*k)*k)*k)*k)/(5.3429314570632e15L*Power(k,40));
}


