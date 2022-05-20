/**********************************************************************
  Hamming_O.c:

     Hamming_O.c is a subroutine to solve a second order differential
     equation from origin using Hamming method.

  Log of Hamming_O.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

#define MSIZE2            10          /* Matrix size of in RIH */

static void RIH(int n,double a[MSIZE2][MSIZE2],
                double b[MSIZE2],double xx[MSIZE2]);
static double Poly_Func(int n,double r,double b[MSIZE2]);
static double DPoly_Func(int n,double r,double b[MSIZE2]);
static double DMF_Func(int eq_type,
                       int cp, double ep, int NL, double kappa,
                       double M, double L, 
                       double VPL[ASIZE1],
                       double VNLp[ASIZE1]);
static double DMF_Func0(int cp, double ep, int NL,
                        double M, double L,
                        double VPL[ASIZE1], double VNLp[ASIZE1]);
static double DMF_Func1(int cp, double ep, int NL,
                        double M, double L,
                        double VPL[ASIZE1], double VNLp[ASIZE1]);
static double DMF_Func2(int cp, double ep, int NL,
                        double M, double L,
                        double VPL[ASIZE1], double VNLp[ASIZE1]);
static double DMF_Func3(int cp, double ep, int NL, double kappa,
                        double M, double L,
                        double VPL[ASIZE1],
                        double VNLp[ASIZE1]);

static double dVF(double R, double VPL[ASIZE1]);

void Hamming_O(int eq_type,
               int NL, double ep, double kappa,
               double Mo[ASIZE1],  double Lo[ASIZE1], double DM[ASIZE1], 
               double VPL[ASIZE1], double VNLp[ASIZE1],
               double Reduce_Num_grid)
{
  /*****************************************************
    eq_type
     0: the usual radial differential equation without 
        relativistic corrections 

        dM/dx = -(2NL + 1)M + 2r^2(V - ep)L
        dL/dx = M

     1: the radial differential equation for evaluating
        logarithmic derivatives of wave functions
        for separable pseudo potentials without 
        relativistic corrections 

        dM/dx = -(2NL+1)M + 2r^2*(VL-ep)L + 2r^(2-NL)*VNLp
        dL/dx = M

     2: the radial differential equation with scalar 
        relativistic corrections

        dM/dx = - (2NL + 1 + r*alpha^2/2/MQ*dV/dr)M 
                + 2MQ*r^2(V - ep)L
                - r*alpha^2/2/MQ*dV/dr*NL*L
        dL/dx = M

     3: the radial differential equation with a full
        relativistic treatment

        dM/dx = - (2NL + 1 + r*alpha^2/2/MQ*dV/dr)M 
                + 2MQ*r^2(V - ep)L
                - r*alpha^2/2/MQ*dV/dr*(NL+1+kappa)*L
        dL/dx = M
  *****************************************************/

  int i,j,k,po,po0,Grid_Num1;
  double a[MSIZE2],b[MSIZE2],c[MSIZE2];
  double tm[MSIZE2][MSIZE2],vv[MSIZE2],tmp0;
  double DMp;
  double Lp,Lm,Lc,Mp,Mm,Lp2,Lc2,Mp2,al,dum;
  double k1M,k2M,k3M,k4M;
  double k1L,k2L,k3L,k4L;
  double r0,f0,Mi,Li,Mtmp,r1,f1;

  Grid_Num1 = (int)((double)Grid_Num*Reduce_Num_grid);

  for (i=0; i<Grid_Num; i++){
    Lo[i] = 0.0;
    Mo[i] = 0.0;
    DM[i] = 0.0;
  }

  /****************************************************
                       From origin
  ****************************************************/

  tm[0][0] = 1.0; tm[0][1] = MRV[0]; tm[0][2] = MRV[0]*MRV[0];
  tm[1][0] = 1.0; tm[1][1] = MRV[1]; tm[1][2] = MRV[1]*MRV[1];
  tm[2][0] = 1.0; tm[2][1] = MRV[2]; tm[2][2] = MRV[2]*MRV[2];
  vv[0] = VPL[0];
  vv[1] = VPL[1];
  vv[2] = VPL[2];

  RIH(2,tm,vv,a);

  if (eq_type==0 || eq_type==2 || eq_type==3){
    b[0] = 1.0;
    b[1] = 0.0;
    b[2] = (a[0]-ep)/(2.0*(double)NL+3.0)*b[0];
    b[3] = a[1]/(3.0*(double)NL+6.0)*b[0];
    b[4] = (a[0]*b[2]+a[2]*b[0]-ep*b[2])/(4.0*(double)NL+10.0);
  }
  else if (eq_type==1){
    tm[0][0] = 1.0; tm[0][1] = MRV[0]; tm[0][2] = MRV[0]*MRV[0];
    tm[1][0] = 1.0; tm[1][1] = MRV[1]; tm[1][2] = MRV[1]*MRV[1];
    tm[2][0] = 1.0; tm[2][1] = MRV[2]; tm[2][2] = MRV[2]*MRV[2];
    vv[0] = VNLp[0];
    vv[1] = VNLp[1];
    vv[2] = VNLp[2];
    RIH(2,tm,vv,c);
    b[0] = 1.0;
    b[1] = 0.0;
    b[2] = (a[0]*b[0]-ep*b[0]+c[0])/(2.0*(double)NL+3.0);
    b[3] = (a[1]*b[0]+c[1])/(3.0*(double)NL+6.0);
    b[4] = (a[0]*b[2]+a[2]*b[0]-ep*b[2]+c[2])/(4.0*(double)NL+10.0);
  }

  Lo[0] = Poly_Func(4,MRV[0],b);
  Lo[1] = Poly_Func(4,MRV[1],b);
  Mo[0] = DPoly_Func(4,MRV[0],b)*MRV[0];
  Mo[1] = DPoly_Func(4,MRV[1],b)*MRV[1];
  DM[0] = DMF_Func(eq_type,0,ep,NL,kappa,Mo[0],Lo[0],VPL,VNLp);
  DM[1] = DMF_Func(eq_type,1,ep,NL,kappa,Mo[1],Lo[1],VPL,VNLp);

  /****************************************************
              a predictor-corrector method
  ****************************************************/

  po = 0; 
  po0 = 0; 
  node = 0;

  /* Schrodinger */

  if (eq_type==0){

    for (i=1; i<=(Grid_Num1-5); i++){

      /* check the classical turning point */

      if (0.0<(VPL[i]-ep) && po==0){
	CTP = i;
	po = 1;
      }

      /* evolve L, M, and DM */

      Lp = 32.0*Lo[i]-31.0*Lo[i-1] - dx*(16.0*Mo[i]+14.0*Mo[i-1]) + dx*dx*(4.0*DM[i]-2.0*DM[i-1]);
      Mp = -4.0*Mo[i]+5.0*Mo[i-1] + dx*(4.0*DM[i]+2.0*DM[i-1]);
      DMp = DMF_Func0(i+1,ep,NL,Mp,Lp,VPL,VNLp);
      Mo[i+1] = Mo[i] + dx/12.0*(8.0*DM[i]-DM[i-1]+5.0*DMp);
      Lo[i+1] = Lo[i] + dx/12.0*(8.0*Mo[i]-Mo[i-1]+5.0*Mo[i+1]);
      DM[i+1] = DMF_Func0(i+1,ep,NL,Mo[i+1],Lo[i+1],VPL,VNLp);

      /* count the number of nodes */

      if (Lo[i]*Lo[i+1]<0.0){
	node++;
      }

      /* check divergence */

      if (LimitE<fabs(Mo[i+1]) ||
	  LimitE<fabs(Lo[i+1]) ||
	  LimitE<fabs(DM[i+1])){

	po0 = 1;
	UL_grid = i;

	for (j=i+2; j<=(Grid_Num1-5); j++){
	  Lo[j] = Lo[i+1];           
	}
      }

      if (po0==1) break;
    }
  }

  /* Schrodinger eq. for a non-local potential */

  else if (eq_type==1){

    for (i=1; i<=(Grid_Num1-5); i++){

      /* check the classical turning point */

      if (0.0<(VPL[i]-ep) && po==0){
	CTP = i;
	po = 1;
      }

      /* evolve L, M, and DM */

      Lp = 32.0*Lo[i]-31.0*Lo[i-1] - dx*(16.0*Mo[i]+14.0*Mo[i-1]) + dx*dx*(4.0*DM[i]-2.0*DM[i-1]);
      Mp = -4.0*Mo[i]+5.0*Mo[i-1] + dx*(4.0*DM[i]+2.0*DM[i-1]);
      DMp = DMF_Func1(i+1,ep,NL,Mp,Lp,VPL,VNLp);
      Mo[i+1] = Mo[i] + dx/12.0*(8.0*DM[i]-DM[i-1]+5.0*DMp);
      Lo[i+1] = Lo[i] + dx/12.0*(8.0*Mo[i]-Mo[i-1]+5.0*Mo[i+1]);
      DM[i+1] = DMF_Func1(i+1,ep,NL,Mo[i+1],Lo[i+1],VPL,VNLp);

      /* count the number of nodes */

      if (Lo[i]*Lo[i+1]<0.0){
	node++;
      }

      /* check divergence */

      if (LimitE<fabs(Mo[i+1]) ||
	  LimitE<fabs(Lo[i+1]) ||
	  LimitE<fabs(DM[i+1])){

	po0 = 1;
	UL_grid = i;

	for (j=i+2; j<=(Grid_Num1-5); j++){
	  Lo[j] = Lo[i+1];           
	}
      }

      if (po0==1) break;
    }
  }

  /* scalar relativistic */

  else if (eq_type==2){

    for (i=1; i<=(Grid_Num1-5); i++){

      /* check the classical turning point */

      if (0.0<(VPL[i]-ep) && po==0){
	CTP = i;
	po = 1;
      }

      /* evolve L, M, and DM */

      Lp = 32.0*Lo[i]-31.0*Lo[i-1] - dx*(16.0*Mo[i]+14.0*Mo[i-1]) + dx*dx*(4.0*DM[i]-2.0*DM[i-1]);
      Mp = -4.0*Mo[i]+5.0*Mo[i-1] + dx*(4.0*DM[i]+2.0*DM[i-1]);
      DMp = DMF_Func2(i+1,ep,NL,Mp,Lp,VPL,VNLp);
      Mo[i+1] = Mo[i] + dx/12.0*(8.0*DM[i]-DM[i-1]+5.0*DMp);
      Lo[i+1] = Lo[i] + dx/12.0*(8.0*Mo[i]-Mo[i-1]+5.0*Mo[i+1]);
      DM[i+1] = DMF_Func2(i+1,ep,NL,Mo[i+1],Lo[i+1],VPL,VNLp);

      /* count the number of nodes */

      if (Lo[i]*Lo[i+1]<0.0){
	node++;
      }

      /* check divergence */

      if (LimitE<fabs(Mo[i+1]) ||
	  LimitE<fabs(Lo[i+1]) ||
	  LimitE<fabs(DM[i+1])){

	po0 = 1;
	UL_grid = i;

	for (j=i+2; j<=(Grid_Num1-5); j++){
	  Lo[j] = Lo[i+1];           
	}
      }

      if (po0==1) break;
    }
  }

  /* full relativistic */

  else if (eq_type==3){

    for (i=1; i<=(Grid_Num1-5); i++){

      /* check the classical turning point */

      if (0.0<(VPL[i]-ep) && po==0){
	CTP = i;
	po = 1;
      }

      /* evolve L, M, and DM */

      Lp = 32.0*Lo[i]-31.0*Lo[i-1] - dx*(16.0*Mo[i]+14.0*Mo[i-1]) + dx*dx*(4.0*DM[i]-2.0*DM[i-1]);
      Mp = -4.0*Mo[i]+5.0*Mo[i-1] + dx*(4.0*DM[i]+2.0*DM[i-1]);
      DMp = DMF_Func3(i+1,ep,NL,kappa,Mp,Lp,VPL,VNLp);
      Mo[i+1] = Mo[i] + dx/12.0*(8.0*DM[i]-DM[i-1]+5.0*DMp);
      Lo[i+1] = Lo[i] + dx/12.0*(8.0*Mo[i]-Mo[i-1]+5.0*Mo[i+1]);
      DM[i+1] = DMF_Func3(i+1,ep,NL,kappa,Mo[i+1],Lo[i+1],VPL,VNLp);

      /* count the number of nodes */

      if (Lo[i]*Lo[i+1]<0.0){
	node++;
      }

      /* check divergence */

      if (LimitE<fabs(Mo[i+1]) ||
	  LimitE<fabs(Lo[i+1]) ||
	  LimitE<fabs(DM[i+1])){

	po0 = 1;
	UL_grid = i;

	for (j=i+2; j<=(Grid_Num1-5); j++){
	  Lo[j] = Lo[i+1];           
	}
      }

      if (po0==1) break;
    }
  }


  if (po0==0) UL_grid = Grid_Num1-6;

}


double Poly_Func(int n,double r,double b[MSIZE2])
{
  static int i,j;
  static double sum,di;

  sum = 0.0;
  di = -1.0;
  for (i=0; i<=n; i++){
    di = di + 1.0;
    sum = sum + b[i]*pow(r,di);
  }

  return sum;
}


double DPoly_Func(int n,double r,double b[MSIZE2])
{
  static int i,j;
  static double sum,di;

  sum = 0.0;
  di = -1.0;
  for (i=1; i<=n; i++){
    di = di + 1.0;
    sum = sum + (di+1.0)*b[i]*pow(r,di);
  }

  return sum;
}

double DMF_Func(int eq_type,
                int cp, double ep, int NL, double kappa,
                double M, double L,
                double VPL[ASIZE1],
                double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Schrodinger equation */

  if (eq_type==0){
    d1 = -(2.0*(double)NL+1.0)*M;
    d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L;
  }

  /* Schrodinger equation for non-local potentials */

  else if (eq_type==1){
    d1 = -(2.0*(double)NL+1.0)*M;
    d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L 
       + 2.0*pow(MRV[cp],2.0-(double)NL)*VNLp[cp];
  }

  /* scalar relativistic equation */

  else if (eq_type==2){
    c = 137.0359895;
    al = 0.0072971395213076740;         /* 1.0/c */
    tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
    r = MRV[cp];
    MQ = 1.0 + tmp0*(ep-VPL[cp]);
    dV = dVF(r,VPL);
    d3 = tmp0*r/MQ*dV;
    d1 = -(2.0*(double)NL + 1.0 + d3)*M;
    /* scalar treatment */
    d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*(double)NL)*L;
  }

  /* full relativistic equation */

  else if (eq_type==3){
    c = 137.0359895;
    al = 0.0072971395213076740;         /* 1.0/c */
    tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
    r = MRV[cp];
    MQ = 1.0 + tmp0*(ep-VPL[cp]);
    dV = dVF(r,VPL);
    d3 = tmp0*r/MQ*dV;
    d1 = -(2.0*(double)NL + 1.0 + d3)*M;
    d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*((double)NL + 1.0 + kappa))*L;
  }

  result = d1 + d2;
  return result;
}

double DMF_Func0(int cp, double ep, int NL,
                 double M, double L,
                 double VPL[ASIZE1],
                 double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Schrodinger */
  d1 = -(2.0*(double)NL+1.0)*M;
  d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L;
  return d1 + d2;
}

double DMF_Func1(int cp, double ep, int NL,
                 double M, double L,
                 double VPL[ASIZE1],
                 double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Schrodinger for non-local potentials */

  d1 = -(2.0*(double)NL+1.0)*M;
  d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L 
     + 2.0*pow(MRV[cp],2.0-(double)NL)*VNLp[cp];

  return d1 + d2;
}

double DMF_Func2(int cp, double ep, int NL,
                 double M, double L,
                 double VPL[ASIZE1],
                 double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* scalar relativistic */

  c = 137.0359895;
  al = 0.0072971395213076740;         /* 1.0/c */
  tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
  r = MRV[cp];
  MQ = 1.0 + tmp0*(ep-VPL[cp]);
  dV = dVF(r,VPL);
  d3 = tmp0*r/MQ*dV;
  d1 = -(2.0*(double)NL + 1.0 + d3)*M;
  /* scalar treatment */
  d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*(double)NL)*L;

  return d1 + d2;
}

double DMF_Func3(int cp, double ep, int NL, double kappa,
                 double M, double L,
                 double VPL[ASIZE1],
                 double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* full relativistic */

  c = 137.0359895;
  al = 0.0072971395213076740;         /* 1.0/c */
  tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
  r = MRV[cp];
  MQ = 1.0 + tmp0*(ep-VPL[cp]);
  dV = dVF(r,VPL);
  d3 = tmp0*r/MQ*dV;
  d1 = -(2.0*(double)NL + 1.0 + d3)*M;
  /* dependence on all angular momentum */
  d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*((double)NL + 1.0 + kappa))*L;

  return d1 + d2;
}


void RIH(int n,double a[MSIZE2][MSIZE2],
         double b[MSIZE2],double xx[MSIZE2])
{

  /****************************************************
                     From 0 to n
  ****************************************************/

  static int i,j,k,l,m;
  static double w,sum;
  static double y[MSIZE2],da[MSIZE2][MSIZE2];

  if (n==-1){
    for (i=0; i<MSIZE2; i++){
      for (j=0; j<MSIZE2; j++){
	a[i][j] = 0.0;
      }
    }
  }
  else{

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	da[i][j] = a[i][j];
      }
    }

    /****************************************************
                      LU factorization
    ****************************************************/

    for (k=0; k<=n-1; k++){
      w = 1.0/a[k][k];
      for (i=k+1; i<=n; i++){
	a[i][k] = w*a[i][k];
	for (j=k+1; j<=n; j++){
	  a[i][j] = a[i][j] - a[i][k]*a[k][j];
	}
      }
    }

    /****************************************************
                           Ly = b
    ****************************************************/

    for (i=0; i<=n; i++){
      y[i] = b[i];
      for (j=0; j<=i-1; j++){
	y[i] = y[i] - a[i][j]*y[j];
      }
    }

    /****************************************************
                          Ux = y
    ****************************************************/

    for (i=n; 0<=i; i--){
      xx[i] = y[i];
      for (j=n; (i+1)<=j; j--){
	xx[i] = xx[i] - a[i][j]*xx[j];
      }
      xx[i] = xx[i]/a[i][i];
    }

    for (i=0; i<=n; i++){
      for (j=0; j<=n; j++){
	a[i][j] = da[i][j];
      }
    }

  }

}



double dVF(double R, double VPL[ASIZE1])
{
  static int mp_min,mp_max,m,po;
  static double h1,h2,h3,f1,f2,f3,f4;
  static double g1,g2,x1,x2,y1,y2,f,df;
  static double rm,a,b,y12,y22,result;

  mp_min = 0;
  mp_max = Grid_Num-1;

  if (MRV[Grid_Num-1]<R){
    result = 0.0;
  }
  else if (R<MRV[0]){
    po = 1;
    m = 4;
    rm = MRV[m];

    h1 = MRV[m-1] - MRV[m-2];
    h2 = MRV[m]   - MRV[m-1];
    h3 = MRV[m+1] - MRV[m];

    f1 = VPL[m-2];
    f2 = VPL[m-1];
    f3 = VPL[m];
    f4 = VPL[m+1];

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = rm - MRV[m-1];
    x2 = rm - MRV[m];
    y1 = x1/h2;
    y2 = x2/h2;
    y12 = y1*y1;
    y22 = y2*y2;

    f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

    df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y22*(2.0*f2 + h2*g1)/h2
       + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
       - y12*(2.0*f3 - h2*g2)/h2;

    a = 0.5*df/rm;
    b = f - a*rm*rm;      
    result = 2.0*a*R;
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

    /****************************************************
                   Spline like interpolation
    ****************************************************/

    h1 = MRV[m-1] - MRV[m-2];
    h2 = MRV[m]   - MRV[m-1];
    h3 = MRV[m+1] - MRV[m];

    f1 = VPL[m-2];
    f2 = VPL[m-1];
    f3 = VPL[m];
    f4 = VPL[m+1];

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
                Calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - MRV[m-1];
    x2 = R - MRV[m];
    y1 = x1/h2;
    y2 = x2/h2;

    f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

    df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y2*y2*(2.0*f2 + h2*g1)/h2
       + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
       - y1*y1*(2.0*f3 - h2*g2)/h2;

    result = df;

  }
  return result;
}
 
