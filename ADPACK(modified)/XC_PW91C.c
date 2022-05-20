/**********************************************************************
  XC_PW91C.c:

     XC_PW91C.c is a subroutine to calculate the correlation
     potential developed by Perdew and Wang (Ref: J.P.Perdew
     and Y.Wang, PRB, 45, 13244 (1992)) for given up (dens[0]) 
     and down (dens[1]) densities.

  Log of PW91C.c:

     10/Dec/2002  Released by T.Ozaki (but this code is originally 
                  wrriten by L.C.Balbas and J.M.Soler for SIESTA code)

***********************************************************************/

#include <math.h>
#include <stdlib.h>
#include "adpack.h"

void XC_PW91C(double dens[2], double Ec[1], double Vc[2])
{

  int i;
  double dtot,rs,zeta,coe,den_min;
  double dum,dum1,dum2,b,c,dbdrs,dcdrs;
  double fpp0,f,dfdz;
  double G[3],dGdrs[3];
  double dEcdd[2];
  double dEcdrs,dEcdz;
  double drsdd,dzdd[2];
  double result;

  /****************************************************
                parameters from Table I of
            Perdew & Wang, PRB, 45, 13244 (92)
  ****************************************************/

  double p[3]      = {1.0000000, 1.0000000, 1.0000000};
  double A[3]      = {0.0310910, 0.0155450, 0.0168870};
  double alpha1[3] = {0.2137000, 0.2054800, 0.1112500};
  double beta1[3]  = {7.5957000,14.1189000,10.3570000};
  double beta2[3]  = {3.5876000, 6.1977000, 3.6231000};
  double beta3[3]  = {1.6382000, 3.3662000, 0.8802600};
  double beta4[3]  = {0.4929400, 0.6251700, 0.4967100};

  den_min = 1.0e-14;

  /****************************************************
                      zeta and rs
  ****************************************************/

  dtot = largest(den_min,dens[0] + dens[1]);
  zeta = (dens[0] - dens[1])/dtot;
  rs = pow(3.0/(4.0*PI*dtot),0.333333333333333);
  drsdd = -rs/dtot/3.0;
  dzdd[0] = 1.0/dtot - zeta/dtot;
  dzdd[1] =-1.0/dtot - zeta/dtot;

  /****************************************************
           eps_c(rs,0)=G(0), eps_c(rs,1)=G(1) and
                  -alpha_c(rs)=G(2)
           using eq.(10) of cited reference
         (Perdew & Wang, PRB, 45, 13244 (1992))
  ****************************************************/

  for (i=0; i<=2; i++){
    dum = sqrt(rs);
    b =   beta1[i]*dum
        + beta2[i]*rs
        + beta3[i]*rs*dum
        + beta4[i]*pow(rs,p[i]+1.0);

    dbdrs =  beta1[i]*0.50/dum
           + beta2[i]
           + beta3[i]*1.50*dum
           + beta4[i]*(p[i] + 1.0)*pow(rs,p[i]);

    c = 1.0 + 1.0/(2.0*A[i]*b);
    dcdrs = -(c - 1.0)*dbdrs/b;
    dum = log(c);
    dum1 = 1.0 + alpha1[i]*rs;
    G[i] = -2.0*A[i]*dum1*dum;
    dGdrs[i] = -2.0*A[i]*(alpha1[i]*dum + dum1*dcdrs/c);
  }

  /****************************************************
            f''(0) and f(zeta) from eq.(9)
  ****************************************************/

  c = 1.92366105093154;
  fpp0 = 1.70992093416137;
  dum1 = 1.0 + zeta;
  dum2 = 1.0 - zeta;
  f = (pow(dum1,1.333333333333333)
     + pow(dum2,1.333333333333333)
     - 2.0)*c;
  dfdz = 1.333333333333333*(pow(dum1,0.333333333333333)
                          - pow(dum2,0.333333333333333))*c;

  /****************************************************
               eps_c(rs,zeta) from eq.(8)
  ****************************************************/

  dum = pow(zeta,4.0);
  dum1 = pow(zeta,3.0);

  Ec[0] =  G[0] - G[2]*f/fpp0*(1.0 - dum) + (G[1] - G[0])*f*dum;

  dEcdrs =   dGdrs[0] - dGdrs[2]*f/fpp0*(1.0 - dum)
          + (dGdrs[1] - dGdrs[0])*f*dum;
  dEcdz = -G[2]/fpp0*(dfdz*(1.0 - dum) - f*4.0*dum1)
          + (G[1] - G[0])*(dfdz*dum + f*4.0*dum1);

  /****************************************************
                Find correlation potential
  ****************************************************/

  dum = dEcdrs*drsdd;
  dEcdd[0] = dum + dEcdz*dzdd[0];
  dEcdd[1] = dum + dEcdz*dzdd[1];
  Vc[0] = Ec[0] + dtot*dEcdd[0];
  Vc[1] = Ec[0] + dtot*dEcdd[1];
} 
