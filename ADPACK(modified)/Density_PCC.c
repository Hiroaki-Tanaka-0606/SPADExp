/**********************************************************************
  Density_PCC.c:

     Density_PCC.c is a subroutine to calculate densities on grids
     which are used in the partial core correction.

  Log of Density_PCC.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Density_PCC(double rho_PCC[ASIZE1], int state_num, double OcpN0[ASIZE15][2][ASIZE3+1][ASIZE3])
{
  static int i,j,ip,jp,po,m,n,l,so;
  static double mag,Vi,dVi,d2Vi,d3Vi,ri,rj,r;
  static double Vj,dVj,d2Vj,d3Vj;
  static double dVix,d2Vix,d3Vix; 
  static double dVjx,d2Vjx,d3Vjx; 
  static double c0,c1,c2,c3,c4,c5,c6,c7;
  static double r2,r3,r4,r5,r6,r7,r8;
  static double dis_r,sum,p,p2;
  static double a[ASIZE6][ASIZE6],c[ASIZE6];
  static double rho_V_tmp[ASIZE1];

  if (fabs(total_electron-valence_electron)<10e-14){
    for (i=0; i<Grid_Num; i++){
      rho_PCC[i] = 0.0;
    }
  }
  else {
    mag = PCC_Ratio;

    Density(state_num);

    /************************************************
       the valence density calculated by the wave
       functions in the all electron calculations
    ************************************************/

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      sum = 0.0;
      for (so=0; so<SOI_switch; so++){
        for (m=0; m<Number_VPS; m++){
          n = NVPS[m];
          l = LVPS[m];
          p = PF[so][n][l][i];
          p2 = p*p;
          sum = sum + OcpN0[state_num][so][n][l]*p2/4.0/PI/r/r;
        }
      }
      rho_V_tmp[i] = sum;
    }

    /****************************************************
       Searching of ip where rho_core = mag*rho_V_tmp
    ****************************************************/

    for (i=0; i<Grid_Num; i++){
      rho_PCC[i] = rho[0][i] - rho_V_tmp[i];
    }

    po = 0;
    i = 0;
    do {
      if (rho_PCC[i]<(mag*rho_V_tmp[i])){
	ip = i;
	po = 1;
      }
      i++;      
    }
    while(po==0);

    /****************************************************
       Value, first, second and third derivatives at ri
    ****************************************************/

    ri = MRV[ip];
    Vi = rho_PCC[ip];
    dVix  = (rho_PCC[ip-2] - 8.0*rho_PCC[ip-1]
	     +8.0*rho_PCC[ip+1] - rho_PCC[ip+2])/12.0/dx;
    d2Vix = (-rho_PCC[ip-2] + 16.0*rho_PCC[ip-1] - 30.0*rho_PCC[ip]
	     +16.0*rho_PCC[ip+1] - rho_PCC[ip+2])/12.0/dx/dx;
    d3Vix = (-rho_PCC[ip-2] + 2.0*rho_PCC[ip-1]
	     -2.0*rho_PCC[ip+1] + rho_PCC[ip+2])/2.0/dx/dx/dx;

    dVi  = dVix/ri;
    d2Vi = -dVix/ri/ri + d2Vix/ri/ri;
    d3Vi = 2.0*dVix/ri/ri/ri - 3.0*d2Vix/ri/ri/ri + d3Vix/ri/ri/ri;

    /****************************************************
       Value, first, second and third derivatives at rj
    ****************************************************/

    jp = 0;
    rj = MRV[jp];
    Vj = PCC_Ratio_Origin*Vi;

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
                          rho_PCC
    ****************************************************/

    for (i=0; i<=(jp-1); i++){
      r = MRV[i];
      r2 = r*r;
      r4 = r2*r2;
      r6 = r2*r4;
      r8 = r2*r6;
      rho_PCC[i] = -(rho_V_tmp[i] - rho_V_tmp[jp]) + Vj;
    }
    for (i=jp; i<=(ip-1); i++){
      r = MRV[i];
      r2 = r*r;
      r3 = r2*r;
      r4 = r2*r2;
      r5 = r4*r;
      r6 = r2*r4;
      r7 = r6*r;
      rho_PCC[i] = c0 + c1*r  + c2*r2 + c3*r3 + c4*r4
                              + c5*r5 + c6*r6 + c7*r7;
    }
  }

}







