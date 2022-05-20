/**********************************************************************
  XC_PBE.c:

     XC_PBE.c is a subroutine to calculate the exchange-correlation
     potential developed by Perdew, Burke and Ernzerhof within
     generalized gradient approximation.

     Ref: J.P.Perdew, K.Burke & M.Ernzerhof, PRL 77, 3865 (1996)
     for given densities

  Log of XC_PBE.c:

     10/Dec/2002  Released by T.Ozaki (but this code is originally 
                  wrriten by L.C.Balbas and J.M.Soler for SIESTA code)

***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "adpack.h"

#define FOUTHD  (4.0/3.0) 
#define HALF    0.50
#define THD     (1.0/3.0)
#define THRHLF  1.50
#define TWO     2.0
#define TWOTHD  (2.0/3.0)
#define beta    0.0667250
#define kappa   0.8040

void XC_PBE(double dens[2], double GDENS[3][2], double Exc[2],
            double DEXDD[2], double DECDD[2],
            double DEXDGD[3][2], double DECDGD[3][2])
{
  int IS,IX;
  double Ec_unif[1],Vc_unif[2];
  double dt,rs,zeta;
  double den_min,gd_min,phi,t,ks,kF,f1,f2,f3,f4;
  double A,H,Fc,Fx;
  double GDMT,GDT[3];
  double DRSDD,DKFDD,DKSDD,DZDD[2],DPDZ;
  double DECUDD,DPDD,DTDD,DF1DD,DF2DD,DF3DD,DF4DD,DADD;
  double DHDD,DFCDD[2],DTDGD,DF3DGD,DF4DGD,DHDGD,DFCDGD[3][2];
  double DS[2],GDMS,KFS,s,f,DFDD,DFXDD[2],Vx_unif[2],Ex_unif[1];
  double GDS,DSDGD,DSDD,DF1DGD,DFDGD,DFXDGD[3][2];
  double D[2],GD[3][2],GDM[2];
  double gamma,mu;

  gamma = (1.0 - log(TWO))/(PI*PI);
  mu = beta*PI*PI/3.0;

  /****************************************************
         Lower bounds of density and its gradient
              to avoid divisions by zero
  ****************************************************/

  den_min = 1.0e-40;
  gd_min  = 1.0e-40;

  /****************************************************
   Translate density and its gradient to new variables
  ****************************************************/

  D[0] = dens[0];
  D[1] = dens[1];
  dt = largest(den_min,dens[0] + dens[1]);

  for (IX=0; IX<=2; IX++){
    GD[IX][0] = GDENS[IX][0];
    GD[IX][1] = GDENS[IX][1];
    GDT[IX] = GDENS[IX][0] + GDENS[IX][1];
  } 
  GDM[0] = sqrt(GD[0][0]*GD[0][0] + GD[1][0]*GD[1][0] + GD[2][0]*GD[2][0]);
  GDM[1] = sqrt(GD[0][1]*GD[0][1] + GD[1][1]*GD[1][1] + GD[2][1]*GD[2][1]);
  GDMT   = sqrt(GDT[0]*GDT[0] + GDT[1]*GDT[1] + GDT[2]*GDT[2]);
  GDMT = largest(gd_min, GDMT);

  /****************************************************
          Local correlation energy and potential 
  ****************************************************/

  XC_PW91C(dens,Ec_unif,Vc_unif);

  /****************************************************
                Total correlation energy
  ****************************************************/

  rs = pow(3.0/(4.0*PI*dt),THD);
  kF = pow(3.0*PI*PI*dt,THD);
  ks = sqrt(4.0*kF/PI);
  zeta = (dens[0] - dens[1])/dt;
  zeta = largest(-1.0 + den_min,zeta);
  zeta = smallest( 1.0 - den_min,zeta);
  phi = 0.50*(pow(1.0 + zeta,TWOTHD)
            + pow(1.0 - zeta,TWOTHD));
  t = GDMT/(2.0*phi*ks*dt);
  f1 = Ec_unif[0]/(gamma*pow(phi,3.0));
  f2 = exp(-f1);
  A = beta/gamma/(f2 - 1.0);
  f3 = t*t + A*t*t*t*t;
  f4 = beta/gamma * f3/(1.0 + A*f3);
  H = gamma*pow(phi,3.0)*log(1.0 + f4);
  Fc = Ec_unif[0] + H;

  /****************************************************
              Correlation energy derivatives
  ****************************************************/

  DRSDD = -(THD*rs/dt);
  DKFDD =   THD*kF/dt;
  DKSDD = HALF*ks*DKFDD/kF;
  DZDD[0] = 1.0/dt - zeta/dt;
  DZDD[1] = -(1.0/dt) - zeta/dt;
  DPDZ = HALF*TWOTHD*(1.0/pow(1.0 + zeta,THD) - 1.0/pow(1.0 - zeta,THD));
  for (IS=0; IS<=1; IS++){
    DECUDD = (Vc_unif[IS] - Ec_unif[0])/dt;
    DPDD = DPDZ*DZDD[IS];
    DTDD = (-t)*(DPDD/phi + DKSDD/ks + 1.0/dt);
    DF1DD = f1*(DECUDD/Ec_unif[0] - 3.0*DPDD/phi);
    DF2DD = (-f2)*DF1DD;
    DADD = (-A)*DF2DD/(f2 - 1.0);
    DF3DD = (2.0*t + 4.0*A*t*t*t) * DTDD + DADD*t*t*t*t;
    DF4DD = f4*(DF3DD/f3 - (DADD*f3+A*DF3DD)/(1.0 + A*f3));
    DHDD = 3.0*H*DPDD/phi;
    DHDD = DHDD + gamma*phi*phi*phi*DF4DD/(1.0 + f4);
    DFCDD[IS] = Vc_unif[IS] + H + dt * DHDD;
    for (IX=0; IX<=2; IX++){
      DTDGD = (t/GDMT)*GDT[IX]/GDMT;
      DF3DGD = DTDGD*(2.0*t + 4.0*A*t*t*t);
      DF4DGD = f4*DF3DGD*(1.0/f3 - A/(1.0 + A*f3));
      DHDGD = gamma*phi*phi*phi*DF4DGD/(1.0 + f4);
      DFCDGD[IX][IS] = dt*DHDGD;
    }
  }

  /****************************************************
              Exchange energy and potential
  ****************************************************/

  Fx = 0.0;
  for (IS=0; IS<=1; IS++){
    DS[IS] = largest(den_min, 2.0*D[IS]);
    GDMS = largest(gd_min, 2.0*GDM[IS]);
    KFS = pow(3.0*PI*PI*DS[IS],THD);
    s = GDMS/(2.0*KFS*DS[IS]);
    f1 = 1.0 + mu*s*s/kappa;
    f = 1.0 + kappa - kappa/f1;

    /****************************************************
                Note nspin=1 in call to XC_EX
    ****************************************************/

    XC_EX(1,DS[IS],DS,Ex_unif,Vx_unif);

    Fx = Fx + DS[IS]*Ex_unif[0]*f;
    DKFDD = THD * KFS/DS[IS];
    DSDD = s*(-(DKFDD/KFS) - 1.0/DS[IS]);
    DF1DD = 2.0*(f1 - 1.0)*DSDD/s;
    DFDD = kappa*DF1DD/(f1*f1);
    DFXDD[IS] = Vx_unif[0]*f + DS[IS]*Ex_unif[0]*DFDD;
    for (IX=0; IX<=2; IX++){
      GDS = 2.0*GD[IX][IS];
      DSDGD = (s/GDMS)*GDS/GDMS;
      DF1DGD = 2.0*mu*s*DSDGD/kappa;
      DFDGD = kappa*DF1DGD/(f1*f1);
      DFXDGD[IX][IS] = DS[IS]*Ex_unif[0]*DFDGD;
    }
  }
  Fx = HALF*Fx/dt;

  /****************************************************
                   Set output arguments
  ****************************************************/

  Exc[0] = Fx;
  Exc[1] = Fc;
  for (IS=0; IS<=1; IS++){
    DEXDD[IS] = DFXDD[IS];
    DECDD[IS] = DFCDD[IS];
    for (IX=0; IX<=2; IX++){
      DEXDGD[IX][IS] = DFXDGD[IX][IS];
      DECDGD[IX][IS] = DFCDGD[IX][IS];
    } 
  }
}



