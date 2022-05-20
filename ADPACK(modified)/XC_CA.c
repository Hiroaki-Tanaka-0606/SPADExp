/**********************************************************************
  XC_CA.c:

     XC_CA.c is a subroutine to calculate the exchange-correlation
     potential within LDA constructed by Ceperly and Alder,
     and parametrized by Perdew and Zunger.

  Log of XC_CA.c:

     10/Dec/2002  Released by T.Ozaki,

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void XC_CA(double rh[ASIZE1], double xc[ASIZE1], int P_switch)
{
  /****************************************************
          P_switch:
              0  \epsilon_XC (XC energy density)  
              1  \mu_XC      (XC potential)  
  ****************************************************/

  static int i,j;
  static double dum,dum2,rs,coe;
  static double Ex,Ec,dEx,dEc;
  static double be,fe,fm,d1,d2; 

  /****************************************************
                    Non-relativisic
  ****************************************************/

  coe = pow(3.0/4.0/PI,1.0/3.0);
  for (i=0; i<Grid_Num; i++){
    if (rh[i]<1.0e-14) rh[i] = 1.0e-14;
    rs = coe*pow(rh[i],-1.0/3.0);
    Ex = -0.4582/rs;
    dEx = 0.4582/rs/rs;
    if (1.0<=rs){
      dum = (1.0+1.0529*sqrt(rs)+0.3334*rs);
      Ec = -0.1423/dum;
      dEc = 0.1423/dum/dum*(1.0529*0.5/sqrt(rs)+0.3334);
    }
    else{
      Ec = -0.0480+0.0311*log(rs)-0.0116*rs+0.0020*rs*log(rs);
      dEc = 0.0311/rs + 0.0020*log(rs) - 0.0096;
    }

    if (P_switch==0)
      xc[i] = Ex + Ec;
    else if (P_switch==1)
      xc[i] = Ex + Ec - 0.333333333333333*rs*(dEx + dEc);
  }

  /****************************************************
                    Non-relativisic
  ****************************************************/

  /*
  if (Equation_Type==0){
    coe = pow(3.0/4.0/PI,1.0/3.0);
    for (i=0; i<Grid_Num; i++){
      if (rh[i]<1.0e-14) rh[i] = 1.0e-14;
      rs = coe*pow(rh[i],-1.0/3.0);
      Ex = -0.4582/rs;
      dEx = 0.4582/rs/rs;
      if (1.0<=rs){
	dum = (1.0+1.0529*sqrt(rs)+0.3334*rs);
	Ec = -0.1423/dum;
	dEc = 0.1423/dum/dum*(1.0529*0.5/sqrt(rs)+0.3334);
      }
      else{
	Ec = -0.0480+0.0311*log(rs)-0.0116*rs+0.0020*rs*log(rs);
	dEc = 0.0311/rs + 0.0020*log(rs) - 0.0096;
      }

      if (P_switch==0)
	xc[i] = Ex + Ec;
      else if (P_switch==1)
	xc[i] = Ex + Ec - 0.333333333333333*rs*(dEx + dEc);
    }
  }
  */

  /****************************************************
           Relativistic exchange-correlation
  ****************************************************/

  /*
  else if (Equation_Type==1){
    coe = pow(3.0/4.0/PI,0.333333333333333);
    for (i=0; i<Grid_Num; i++){
      if (rh[i]<1.0e-14) rh[i] = 1.0e-14;
      rs = coe*pow(rh[i],-0.333333333333333);
      be = 0.0140/rs;
      dum2 = log(be+sqrt(1.0+be*be));
      dum = sqrt(1.0+be)/be - dum2/be/be;
      fe = 1.0 - 1.50*dum*dum;
      fm = -0.50 + 1.50*dum2/be/sqrt(1.0+be*be);

      Ex = -0.4582/rs;
      dEx = 0.4582/rs/rs;
      if (1.0<=rs){
	dum = (1.0+1.0529*sqrt(rs)+0.3334*rs);
	Ec = -0.1423/dum;
	dEc = 0.1423/dum/dum*(1.0529*0.5/sqrt(rs)+0.3334);
      }
      else{
	Ec = -0.0480+0.0311*log(rs)-0.0116*rs+0.0020*rs*log(rs);
	dEc = 0.0311/rs + 0.0020*log(rs) - 0.0096;
      }

      if (P_switch==0){
	xc[i] = fe*Ex + Ec;
      }
      else if (P_switch==1){
	xc[i] = fm*(Ex - 0.333333333333333*rs*dEx)
                  + Ec - 0.333333333333333*rs*dEc;
      }
    }
  }
  */

} 


