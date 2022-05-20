/**********************************************************************
  Calc_Vlocal.c:

     Calc_Vlocal.c is a subroutine to calculate the local potential
     in the KB form.

  Log of Calc_Vlocal.c:

     12/Feb/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

static void Spline_Func(double R, double fv[3], double *af);

void Calc_Vlocal(double Vlocal[ASIZE1], double *V0, int local_switch)
{
  static int i,j,ip,jp,po,m,n,l;
  static double Vi,dVi,d2Vi,d3Vi,ri,rj,r;
  static double Vj,dVj,d2Vj,d3Vj;
  static double c0,c1,c2,c3,c4,c5,c6,c7;
  static double r2,r3,r4,r5,r6,r7,r8;
  static double dis_r,sum,p,p2,MinD,tmp1,tmp2;
  static double a[ASIZE6][ASIZE6],c[ASIZE6];
  static double tmp0,Z,dr,fv[3];
  static double vm1,vm2,vp1,vp2;

  Z = (double)AtomNum - (total_electron - valence_electron);

  /****************************************************
    find ip where r_i is close to specified r
  ****************************************************/

  MinD = 10000.0;
  for (i=0; i<Grid_Num; i++){
    tmp0 = fabs(MRV[i]-Vlocal_cutoff);
    if (tmp0<MinD){
      MinD = tmp0;
      ip = i;
    }
  }

  /****************************************************
    Value, first, second and third derivatives at ri
  ****************************************************/

  ri = MRV[ip];

  if (local_switch==0){

    Vi   = -Z/ri;
    dVi  = Z/ri/ri;
    d2Vi = -2.0*Z/ri/ri/ri;
    d3Vi =  6.0*Z/ri/ri/ri/ri;
  }

  else if (local_switch==4){

    Spline_Func(ri,fv,V0);

    Vi   = fv[0];
    dVi  = fv[1];
    d2Vi = fv[2];

    Spline_Func(ri-0.001,fv,V0);
    tmp1 = fv[2];
    Spline_Func(ri+0.001,fv,V0);
    tmp2 = fv[2];
    d3Vi = (tmp2-tmp1)/0.002;
  }

  /****************************************************
    Value, first, second and third derivatives at rj
  ****************************************************/

  jp = 0;
  rj = MRV[jp];
  Vj = Vlocal_origin*Vi;
  dVj  =  0.0;
  d2Vj =  0.0;
  d3Vj =  0.0;

  /****************************************************
                           At ri
  ****************************************************/

  r2 = ri*ri;
  r3 = ri*r2;
  r4 = ri*r3;
  r5 = ri*r4;
  r6 = ri*r5;
  r7 = ri*r6;
 
  a[0][0] = 1.0;
  a[0][1] = ri;
  a[0][2] = r2;
  a[0][3] = r3;
  a[0][4] = r4;
  a[0][5] = r5;
  a[0][6] = r6;
  a[0][7] = r7;

  a[1][0] = 0.0;
  a[1][1] = 1.0;
  a[1][2] = 2.0*ri;
  a[1][3] = 3.0*r2;
  a[1][4] = 4.0*r3;
  a[1][5] = 5.0*r4;
  a[1][6] = 6.0*r5;
  a[1][7] = 7.0*r6;

  a[2][0] = 0.0;
  a[2][1] = 0.0;
  a[2][2] = 2.0;
  a[2][3] = 6.0*ri;
  a[2][4] = 12.0*r2;
  a[2][5] = 20.0*r3;
  a[2][6] = 30.0*r4;
  a[2][7] = 42.0*r5;

  a[3][0] = 0.0;
  a[3][1] = 0.0;
  a[3][2] = 0.0;
  a[3][3] = 6.0;
  a[3][4] = 24.0*ri;
  a[3][5] = 60.0*r2;
  a[3][6] = 120.0*r3;
  a[3][7] = 210.0*r4;

  a[0][8] = Vi;
  a[1][8] = dVi;
  a[2][8] = d2Vi;
  a[3][8] = d3Vi;

  /****************************************************
                         At rj
  ****************************************************/

  r2 = rj*rj;
  r3 = rj*r2;
  r4 = rj*r3;
  r5 = rj*r4;
  r6 = rj*r5;
  r7 = rj*r6;
 
  a[4][0] = 1.0;
  a[4][1] = rj;
  a[4][2] = r2;
  a[4][3] = r3;
  a[4][4] = r4;
  a[4][5] = r5;
  a[4][6] = r6;
  a[4][7] = r7;

  a[5][0] = 0.0;
  a[5][1] = 1.0;
  a[5][2] = 2.0*rj;
  a[5][3] = 3.0*r2;
  a[5][4] = 4.0*r3;
  a[5][5] = 5.0*r4;
  a[5][6] = 6.0*r5;
  a[5][7] = 7.0*r6;

  a[6][0] = 0.0;
  a[6][1] = 0.0;
  a[6][2] = 2.0;
  a[6][3] = 6.0*rj;
  a[6][4] = 12.0*r2;
  a[6][5] = 20.0*r3;
  a[6][6] = 30.0*r4;
  a[6][7] = 42.0*r5;

  a[7][0] = 0.0;
  a[7][1] = 0.0;
  a[7][2] = 0.0;
  a[7][3] = 6.0;
  a[7][4] = 24.0*rj;
  a[7][5] = 60.0*r2;
  a[7][6] = 120.0*r3;
  a[7][7] = 210.0*r4;

  a[4][8] = Vj;
  a[5][8] = dVj;
  a[6][8] = d2Vj;
  a[7][8] = d3Vj;

  /****************************************************
              c0,c1,c2,c3,c4,c5,c6, and c7
  ****************************************************/

  Gauss_LEQ(7,a,c);
  c0 = c[0];
  c1 = c[1];
  c2 = c[2];
  c3 = c[3];
  c4 = c[4];
  c5 = c[5];
  c6 = c[6];
  c7 = c[7];

  /****************************************************
                         Vlocal
  ****************************************************/

  for (i=0; i<=(ip-1); i++){
    r = MRV[i];
    r2 = r*r;
    r3 = r2*r;
    r4 = r2*r2;
    r5 = r4*r;
    r6 = r2*r4;
    r7 = r6*r;
    Vlocal[i] = c0 + c1*r  + c2*r2 + c3*r3 + c4*r4
      + c5*r5 + c6*r6 + c7*r7;
  }

  if (local_switch==0){

    for (i=ip; i<Grid_Num; i++){
      r = MRV[i];
      Vlocal[i] = -Z/r;
    }
  }

  else if (local_switch==4){
    for (i=ip; i<Grid_Num; i++){
      Vlocal[i] = V0[i];
    }
  }

}



static void Spline_Func(double R, double fv[3], double *af)
{
  int mp_min,mp_max,m;
  double h1,h2,h3,f1,f2,f3,f4;
  double g1,g2,x1,x2,y1,y2,f,Df,D2f;
  double result;

  mp_min = 4;
  mp_max = Grid_Num - 1;
 
  if (R<MRV[4]){
    m = 5;
  }
  else if (MRV[Grid_Num-1]<R){
    m = Grid_Num - 3;
  } 
  else{
    do{ 
      m = (mp_min + mp_max)/2;
      if (MRV[m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;
  }

  /****************************************************
                 Spline like interpolation
  ****************************************************/

  h1 = MRV[m-1] - MRV[m-2];
  h2 = MRV[m]   - MRV[m-1];
  h3 = MRV[m+1] - MRV[m];

  f1 = af[m-2];
  f2 = af[m-1];
  f3 = af[m];
  f4 = af[m+1];

  /****************************************************
                 Treatment of edge points
  ****************************************************/

  if (m==1){
    h1 = -(h2+h3);
    f1 = f4;
  }
  if (m==(Grid_Num-1)){
    h3 = -(h1+h2);
    f4 = f1;
  }

  /****************************************************
              calculate the values at R
  ****************************************************/

  g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
  g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

  x1 = R - MRV[m-1];
  x2 = R - MRV[m];
  y1 = x1/h2;
  y2 = x2/h2;

  f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

  Df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + y2*y2*(2.0*f2 + h2*g1)/h2
    + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
    - y1*y1*(2.0*f3 - h2*g2)/h2;
    
  D2f =  2.0*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + 4.0*y2*(2.0*f2 + h2*g1)
    + 2.0*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
    - 4.0*y1*(2.0*f3 - h2*g2);
  D2f = D2f/h2/h2;

  fv[0] = f;
  fv[1] = Df;
  fv[2] = D2f;
}
