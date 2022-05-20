/**********************************************************************
  Hamming_I.c:

     Hamming_I.c is a subroutine to solve a second order differential
     equation from infinity using the Hamming method.

  Log of Hamming_I.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

static double DMF_Func(int eq_type, int cp, double ep, int NL, double kappa,
                       double M, double L,
                       double VPL[ASIZE1], double VNLp[ASIZE1]);
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
                        double VPL[ASIZE1], double VNLp[ASIZE1]);

static double dVF(double R, double VPL[ASIZE1]);

void Hamming_I(int eq_type, int NL, double ep, double kappa,
               double Mi[ASIZE1],  double Li[ASIZE1], double DM[ASIZE1],
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

     2: the scaler relativistic radial differential equation
        ref: D.D.Koelling and B.N.Harmon,
             J.Phys.C: Solid State Phys. 10, 3107 (1977).

        dM/dx = -(2*NL + 1 + W*dV/dr)M
                +(2*MQ*r^2(V - ep) - W*dV/dr*NL)L
        dL/dx = M

     3: the radial differential equation with a full
        relativistic treatment

        dM/dx = - (2NL + 1 + r*alpha^2/2/MQ*dV/dr)M 
                + 2MQ*r^2(V - ep)L
                - r*alpha^2/2/MQ*dV/dr*(NL+1+kappa)*L
        dL/dx = M
  *****************************************************/

  int i,j,k,po,po0,Grid_Num1;
  double DMp,Lp,Mp,al,MinV;
  double r0,r1,f0,f1,Mtmp; 

  Grid_Num1 = (int)((double)Grid_Num*Reduce_Num_grid);

  /****************************************************
                    From an infinity
  ****************************************************/

  al = sqrt(2.0*(Vinf-ep));
  MinV = 1.0e-30;

  for (i=0; i<Grid_Num; i++){
    Li[i] = 0.0;
    Mi[i] = 0.0;
    DM[i] = 0.0;
  }

  Li[Grid_Num1-1] = exp(-al*MRV[Grid_Num1-1])/pow(MRV[Grid_Num1-1],NL+1.0);
  Li[Grid_Num1-2] = exp(-al*MRV[Grid_Num1-2])/pow(MRV[Grid_Num1-2],NL+1.0);

  if (fabs(Li[Grid_Num1-1])<MinV){
    Li[Grid_Num1-1] = MinV;
    Li[Grid_Num1-2] = MinV;
  }

  Mi[Grid_Num1-1] = -Li[Grid_Num1-1]*(al*MRV[Grid_Num1-1]+NL+1.0);
  Mi[Grid_Num1-2] = -Li[Grid_Num1-2]*(al*MRV[Grid_Num1-2]+NL+1.0);

  if (fabs(Mi[Grid_Num1-1])<MinV){
    Mi[Grid_Num1-1] = MinV;
    Mi[Grid_Num1-2] = MinV;
  }
  
  DM[Grid_Num1-1] = DMF_Func(eq_type,Grid_Num1-1,ep,NL,kappa,
                             Mi[Grid_Num1-1],Li[Grid_Num1-1],VPL,VNLp);
  DM[Grid_Num1-2] = DMF_Func(eq_type,Grid_Num1-2,ep,NL,kappa,
                             Mi[Grid_Num1-2],Li[Grid_Num1-2],VPL,VNLp);

  if (fabs(DM[Grid_Num1-1])<MinV){
    DM[Grid_Num1-1] = MinV;
    DM[Grid_Num1-2] = MinV;
  }

  /****************************************************
              a predictor-corrector method
  ****************************************************/

  po = 0; 
  po0 = 0; 

  /****************************************************
     0: the usual radial differential equation without 
        relativistic corrections 

        dM/dx = -(2NL + 1)M + 2r^2(V - ep)L
        dL/dx = M
  *****************************************************/

  if (eq_type==0){

    for (i=(Grid_Num1-2); 2<=i; i--){

      Lp = 32.0*Li[i]-31.0*Li[i+1] + dx*(16.0*Mi[i]+14.0*Mi[i+1]) + dx*dx*(4.0*DM[i]-2.0*DM[i+1]);
      Mp = -4.0*Mi[i]+5.0*Mi[i+1] - dx*(4.0*DM[i]+2.0*DM[i+1]);
      DMp = DMF_Func0(i-1,ep,NL,Mp,Lp,VPL,VNLp);
      Mi[i-1] = Mi[i] - dx/12.0*(8.0*DM[i]-DM[i+1]+5.0*DMp);
      Li[i-1] = Li[i] - dx/12.0*(8.0*Mi[i]-Mi[i+1]+5.0*Mi[i-1]);
      DM[i-1] = DMF_Func0(i-1,ep,NL,Mi[i-1],Li[i-1],VPL,VNLp);

      if (LimitE<fabs(Mi[i-1]) ||
	  LimitE<fabs(Li[i-1]) ||
	  LimitE<fabs(DM[i-1])){
    
	po = 1;
	LL_grid = i - 2;
	for (j=i-2; 2<=j; j--){
	  Li[j] = Li[i-1];
	}                      
      }
      if (po==1) break;
    }
  }

  /****************************************************
     1: the radial differential equation for evaluating
        logarithmic derivatives of wave functions
        for separable pseudo potentials without 
        relativistic corrections 

        dM/dx = -(2NL+1)M + 2r^2*(VL-ep)L + 2r^(2-NL)*VNLp 
        dL/dx = M
  *****************************************************/

  else if (eq_type==1){

    for (i=(Grid_Num1-2); 2<=i; i--){

      Lp = 32.0*Li[i]-31.0*Li[i+1] + dx*(16.0*Mi[i]+14.0*Mi[i+1]) + dx*dx*(4.0*DM[i]-2.0*DM[i+1]);
      Mp = -4.0*Mi[i]+5.0*Mi[i+1] - dx*(4.0*DM[i]+2.0*DM[i+1]);
      DMp = DMF_Func1(i-1,ep,NL,Mp,Lp,VPL,VNLp);
      Mi[i-1] = Mi[i] - dx/12.0*(8.0*DM[i]-DM[i+1]+5.0*DMp);
      Li[i-1] = Li[i] - dx/12.0*(8.0*Mi[i]-Mi[i+1]+5.0*Mi[i-1]);
      DM[i-1] = DMF_Func1(i-1,ep,NL,Mi[i-1],Li[i-1],VPL,VNLp);

      if (LimitE<fabs(Mi[i-1]) ||
	  LimitE<fabs(Li[i-1]) ||
	  LimitE<fabs(DM[i-1])){
    
	po = 1;
	LL_grid = i - 2;
	for (j=i-2; 2<=j; j--){
	  Li[j] = Li[i-1];
	}                      
      }
      if (po==1) break;
    }
  }

  /****************************************************
     2: the scalar relativistic radial differential equation
        ref: D.D.Koelling and B.N.Harmon,
             J.Phys.C: Solid State Phys. 10, 3107 (1977).

        dM/dx = -(2*NL + 1 + W*dV/dr)M
                +(2*MQ*r^2(V - ep) - W*dV/dr*NL)L
        dL/dx = M
  *****************************************************/

  else if (eq_type==2){

    for (i=(Grid_Num1-2); 2<=i; i--){

      Lp = 32.0*Li[i]-31.0*Li[i+1] + dx*(16.0*Mi[i]+14.0*Mi[i+1]) + dx*dx*(4.0*DM[i]-2.0*DM[i+1]);
      Mp = -4.0*Mi[i]+5.0*Mi[i+1] - dx*(4.0*DM[i]+2.0*DM[i+1]);
      DMp = DMF_Func2(i-1,ep,NL,Mp,Lp,VPL,VNLp);
      Mi[i-1] = Mi[i] - dx/12.0*(8.0*DM[i]-DM[i+1]+5.0*DMp);
      Li[i-1] = Li[i] - dx/12.0*(8.0*Mi[i]-Mi[i+1]+5.0*Mi[i-1]);
      DM[i-1] = DMF_Func2(i-1,ep,NL,Mi[i-1],Li[i-1],VPL,VNLp);

      if (LimitE<fabs(Mi[i-1]) ||
	  LimitE<fabs(Li[i-1]) ||
	  LimitE<fabs(DM[i-1])){
    
	po = 1;
	LL_grid = i - 2;
	for (j=i-2; 2<=j; j--){
	  Li[j] = Li[i-1];
	}                      
      }
      if (po==1) break;
    }
  }

  /****************************************************
     3: the radial differential equation with a full
        relativistic treatment

        dM/dx = - (2NL + 1 + r*alpha^2/2/MQ*dV/dr)M 
                + 2MQ*r^2(V - ep)L
                - r*alpha^2/2/MQ*dV/dr*(NL+1+kappa)*L
        dL/dx = M
  *****************************************************/

  else if (eq_type==3){

    for (i=(Grid_Num1-2); 2<=i; i--){

      Lp = 32.0*Li[i]-31.0*Li[i+1] + dx*(16.0*Mi[i]+14.0*Mi[i+1]) + dx*dx*(4.0*DM[i]-2.0*DM[i+1]);
      Mp = -4.0*Mi[i]+5.0*Mi[i+1] - dx*(4.0*DM[i]+2.0*DM[i+1]);
      DMp = DMF_Func3(i-1,ep,NL,kappa,Mp,Lp,VPL,VNLp);
      Mi[i-1] = Mi[i] - dx/12.0*(8.0*DM[i]-DM[i+1]+5.0*DMp);
      Li[i-1] = Li[i] - dx/12.0*(8.0*Mi[i]-Mi[i+1]+5.0*Mi[i-1]);
      DM[i-1] = DMF_Func3(i-1,ep,NL,kappa,Mi[i-1],Li[i-1],VPL,VNLp);

      if (LimitE<fabs(Mi[i-1]) ||
	  LimitE<fabs(Li[i-1]) ||
	  LimitE<fabs(DM[i-1])){
    
	po = 1;
	LL_grid = i - 2;
	for (j=i-2; 2<=j; j--){
	  Li[j] = Li[i-1];
	}                      
      }
      if (po==1) break;
    }
  }

  if (po==0) LL_grid = 2; 
}




double DMF_Func(int eq_type,
                int cp, double ep, int NL, double kappa,
                double M, double L,
                double VPL[ASIZE1],
                double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Schrodinger */

  if (eq_type==0){
    d1 = -(2.0*(double)NL + 1.0)*M;
    d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp] - ep)*L;
  }

  /* Schrodinger for non-local potentials */

  else if (eq_type==1){
    d1 = -(2.0*(double)NL+1.0)*M;
    d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L 
       + 2.0*pow(MRV[cp],2.0-(double)NL)*VNLp[cp];
  }

  /* scalar relativistic */

  else if (eq_type==2){
    c = 137.0359895;
    al = 0.0072971395213076740;         /* 1.0/c */
    tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
    r = MRV[cp];
    MQ = 1.0 + tmp0*(ep-VPL[cp]);
    dV = dVF(r,VPL);
    d3 = tmp0*r/MQ*dV;
    d1 = -(2.0*(double)NL + 1.0 + d3)*M;
    /* scaler treatment */
    d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*(double)NL)*L;
  }

  /* full relativistic */

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
                 double VPL[ASIZE1], double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Schrodinger */
  d1 = -(2.0*(double)NL + 1.0)*M;
  d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp] - ep)*L;
  return d1 + d2;
}


double DMF_Func1(int cp, double ep, int NL,
                 double M, double L,
                 double VPL[ASIZE1], double VNLp[ASIZE1])
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
                 double VPL[ASIZE1], double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Dirac */

  c = 137.0359895;
  al = 0.0072971395213076740;         /* 1.0/c */
  tmp0 = 0.0000266241225967150;       /* 0.5*al*al */
  r = MRV[cp];
  MQ = 1.0 + tmp0*(ep-VPL[cp]);
  dV = dVF(r,VPL);
  d3 = tmp0*r/MQ*dV;
  d1 = -(2.0*(double)NL + 1.0 + d3)*M;
  /* scaler treatment */
  d2 = (2.0*r*r*MQ*(VPL[cp]-ep)-d3*(double)NL)*L;

  return d1 + d2;
}

double DMF_Func3(int cp, double ep, int NL, double kappa,
                 double M, double L,
                 double VPL[ASIZE1], double VNLp[ASIZE1])
{
  static double d1,d2,result,MQ,c,al,d3,dV,r,tmp0;

  /* Dirac */

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
 

