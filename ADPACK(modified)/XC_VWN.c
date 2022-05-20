/**********************************************************************
  XC_VWN.c:

     XC_VWN.c is a subroutine to calculate the exchange-correlation
     potential within LDA constructed by Ceperly and Alder,
     and parametrized by Vosko, Wilk, and Nusair (VWN).

  Log of XC_VWN.c:

     13/May/2010  Released by T.Ozaki,

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void XC_VWN(double rh[ASIZE1], double xc[ASIZE1], int P_switch)
{
  /****************************************************
          P_switch:
              0  \epsilon_XC (XC energy density)  
              1  \mu_XC      (XC potential)  
  ****************************************************/

  int i,j;
  double dum,dum2,rs,coe;
  double Ex,Ec,dEx,dEc;
  double tmp,X0,x0,X,Q,x;
  double A,b,c;

  /****************************************************
                    Non-relativisic
  ****************************************************/

  coe = pow(3.0/4.0/PI,1.0/3.0);

  for (i=0; i<Grid_Num; i++){
    if (rh[i]<1.0e-15) rh[i] = 1.0e-15;
    rs = coe*pow(rh[i],-1.0/3.0);

    /* the exchange part */

    tmp = 3.0/4.0*pow(9.0/(4.0*PI*PI),1.0/3.0);
    Ex = -tmp/rs; 
    dEx = tmp/rs/rs;

    /* the correlation part */

    A = 0.0310907;
    b = 3.72744;
    c = 12.9352;
    x0 = -0.10498;
    X0 = x0*x0 + b*x0 + c;  

    x = sqrt(rs);  
    X = x*x + b*x + c;
    Q = sqrt(4.0*c-b*b);

    Ec = A*( log(x*x/X)
             + 2.0*b/Q*atan(Q/(2.0*x+b)) 
             - b*x0/X0*(log((x-x0)*(x-x0)/X) 
	     + 2.0*(b+2.0*x0)/Q*atan(Q/(2.0*x+b)))
             );

    dEc = (A*((2.0*c + b*x)/(c + rs + b*x) - (4.0*b*x)/(b*b + Q*Q + 4.0*rs + 4.0*b*x) 
	 - (b*x*x0*((-4.0*(b + 2.0*x0))/(b*b + Q*Q + 4.0*rs + 4.0*b*x) 
	 + (2.0*c + 2.0*x*x0 + b*(x + x0))/((c + rs + b*x)*(x - x0))))/X0))/(2.0*rs);

    /* exchange-correlation energy density */

    if (P_switch==0)
      xc[i] = Ex + Ec;

    /* exchange-correlation potential */

    else if (P_switch==1)
      xc[i] = Ex + Ec - 0.333333333333333*rs*(dEx + dEc);

  }

} 


