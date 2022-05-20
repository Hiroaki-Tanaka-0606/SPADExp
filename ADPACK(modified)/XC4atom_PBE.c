/**********************************************************************
  XC4atom_PBE.c:

     XC4atom_PBE.c is a subroutine to calculate the exchange-correlation
     potential within GGA by PBE.

  Log of XC4atom_PBE.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

static void F_rho(double R,
                  double rh[ASIZE1],
                  double Frho[2]);
static void F_grho(double x, double y, double z,
                   double rh[ASIZE1],
                   double GEdens[4][2]);

void XC4atom_PBE(double rh[ASIZE1], double xc[ASIZE1], int P_switch)
{
  /********************************************************
   GGA exchange-correlation potential constructed by PBE.
   This routine dose not take into account spin-poralized
   states.
          P_switch:
              0  \epsilon_XC (XC energy density)  
              1  \mu_XC      (XC potential)  
  *********************************************************/

  int i,IX1,ID1,Ip;
  double r,delta1,d1[2],CXYZ1[3],GEdens[4][2];
  double dens[2],GDENS[3][2],Exc[2][3][2];
  double DEXDD[2],DECDD[2];
  double DEXDGD[3][2],DECDGD[3][2];
  double dfdgn[2][3][2],gdfdgn[3][2];
  double dExdn[2][3][2],dEcdn[2][3][2];
  double r0,r1,r2,a[ASIZE6][ASIZE6],x[ASIZE6];

  delta1 = 0.01;  /* default=0.01 */
  d1[0] = -delta1;
  d1[1] =  delta1;

  for (i=0; i<Grid_Num; i++){
    r = MRV[i];

    for (IX1=0; IX1<=2; IX1++){
      for (ID1=0; ID1<=1; ID1++){

        CXYZ1[0] = 0.0;
        CXYZ1[1] = 0.0;
        CXYZ1[2] = r;
        CXYZ1[IX1] = CXYZ1[IX1] + d1[ID1];
        F_grho(CXYZ1[0],CXYZ1[1],CXYZ1[2],rh,GEdens);
        dens[0] = GEdens[0][0];
        dens[1] = GEdens[0][1];

        GDENS[0][0] = GEdens[1][0];
        GDENS[1][0] = GEdens[2][0];
        GDENS[2][0] = GEdens[3][0];
        GDENS[0][1] = GEdens[1][1];
        GDENS[1][1] = GEdens[2][1];
        GDENS[2][1] = GEdens[3][1];

        XC_PBE(dens,GDENS,Exc[ID1][IX1],DEXDD,DECDD,DEXDGD,DECDGD);

        dfdgn[ID1][IX1][0] = DEXDGD[IX1][0] + DECDGD[IX1][0];
        dfdgn[ID1][IX1][1] = DEXDGD[IX1][1] + DECDGD[IX1][1];
        dExdn[ID1][IX1][0] = DEXDD[0];
        dExdn[ID1][IX1][1] = DEXDD[1];
        dEcdn[ID1][IX1][0] = DECDD[0];
        dEcdn[ID1][IX1][1] = DECDD[1];

      } 
    }

    if (P_switch==0){
      xc[i] = (Exc[0][0][0] + Exc[1][0][0]
	     + Exc[0][1][0] + Exc[1][1][0] 
             + Exc[0][2][0] + Exc[1][2][0]
             + Exc[0][0][1] + Exc[1][0][1]
	     + Exc[0][1][1] + Exc[1][1][1] 
             + Exc[0][2][1] + Exc[1][2][1])/6.0;
    }
    else if (P_switch==1){
      for (IX1=0; IX1<=2; IX1++){
        gdfdgn[IX1][0] = 0.50*(dfdgn[1][IX1][0] - dfdgn[0][IX1][0])/delta1;
        gdfdgn[IX1][1] = 0.50*(dfdgn[1][IX1][1] - dfdgn[0][IX1][1])/delta1;
      }

      DEXDD[0] = (dExdn[0][0][0] + dExdn[1][0][0]
	        + dExdn[0][1][0] + dExdn[1][1][0] 
	        + dExdn[0][2][0] + dExdn[1][2][0])/6.0;
      DEXDD[1] = (dExdn[0][0][1] + dExdn[1][0][1]
	        + dExdn[0][1][1] + dExdn[1][1][1]
	        + dExdn[0][2][1] + dExdn[1][2][1])/6.0;

      DECDD[0] = (dEcdn[0][0][0] + dEcdn[1][0][0]
	        + dEcdn[0][1][0] + dEcdn[1][1][0]
	        + dEcdn[0][2][0] + dEcdn[1][2][0])/6.0;
      DECDD[1] = (dEcdn[0][0][1] + dEcdn[1][0][1]
	        + dEcdn[0][1][1] + dEcdn[1][1][1]
	        + dEcdn[0][2][1] + dEcdn[1][2][1])/6.0;

      xc[i] = DEXDD[0] + DECDD[0] - (gdfdgn[0][0]
			           + gdfdgn[1][0]
				   + gdfdgn[2][0]);
    }
  }

  /****************************************************
           Correction of around the origin
  ****************************************************/

  Ip = Grid_Num/200;
  r0 = MRV[Ip];    
  r1 = MRV[Ip+1];  
  r2 = MRV[Ip+2];  

  a[0][0] = r0*r0;
  a[0][1] = r0;
  a[0][2] = 1.0;
  a[1][0] = r1*r1;
  a[1][1] = r1;
  a[1][2] = 1.0;
  a[2][0] = r2*r2;
  a[2][1] = r2;
  a[2][2] = 1.0;

  a[0][3] = xc[Ip];
  a[1][3] = xc[Ip+1];
  a[2][3] = xc[Ip+2];

  Gauss_LEQ(2,a,x);

  for (i=0; i<=Ip-1; i++){
    r = MRV[i];
    xc[i] = x[0]*r*r + x[1]*r + x[2];
  }

  /****************************************************
          Correction of around the end point
  ****************************************************/

  Ip = Grid_Num - 20;
  r0 = MRV[Ip];     
  r1 = MRV[Ip-1];   
  r2 = MRV[Ip-2];   
  a[0][0] = r0*r0;
  a[0][1] = r0;
  a[0][2] = 1.0;
  a[1][0] = r1*r1;
  a[1][1] = r1;
  a[1][2] = 1.0;
  a[2][0] = r2*r2;
  a[2][1] = r2;
  a[2][2] = 1.0;

  a[0][3] = xc[Ip];
  a[1][3] = xc[Ip-1];
  a[2][3] = xc[Ip-2];

  Gauss_LEQ(2,a,x);

  for (i=(Ip+1); i<Grid_Num; i++){
    r = MRV[i];
    xc[i] = x[0]*r*r + x[1]*r + x[2];
  }

}

void F_grho(double x, double y, double z,
            double rh[ASIZE1],
            double GEdens[4][2])
{

  static double r,Frho[2];
  
  r = sqrt(x*x + y*y + z*z);
  F_rho(r,rh,Frho);

  GEdens[0][0] = 0.50*Frho[0];
  GEdens[1][0] = 0.50*Frho[1]*x/r;
  GEdens[2][0] = 0.50*Frho[1]*y/r;
  GEdens[3][0] = 0.50*Frho[1]*z/r;

  GEdens[0][1] = GEdens[0][0]; 
  GEdens[1][1] = GEdens[1][0]; 
  GEdens[2][1] = GEdens[2][0]; 
  GEdens[3][1] = GEdens[3][0];
}

void F_rho(double R,
           double rh[ASIZE1],
           double Frho[2])
{

  static int mp_min,mp_max,m;
  static double h1,h2,h3,f1,f2,f3,f4;
  static double g1,g2,x1,x2,y1,y2,f,Df;
  static double result;

  mp_min = 0;
  mp_max = Grid_Num - 1;

  if (MRV[Grid_Num-1]<R){
    m = Grid_Num - 3;
  }
  else if (R<MRV[0]){
    m = 2;
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

  f1 = rh[m-2];
  f2 = rh[m-1];
  f3 = rh[m];
  f4 = rh[m+1];

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

  Df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + y2*y2*(2.0*f2 + h2*g1)/h2
    + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
    - y1*y1*(2.0*f3 - h2*g2)/h2;

  Frho[0] = f;
  Frho[1] = Df;
}
