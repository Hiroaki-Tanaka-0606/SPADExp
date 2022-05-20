/*************************************************************************
  Simple_Mixing.c:

  Simple_Mixing.c is a subroutine to achieve self-consistent field using
  a simple mixing method. In this simple mixing scheme, a mixing weight
  is changed by monitoring the historical behavior of the residual norm
  of the charge density and the sum of the Kohn-Sham eigenvalues to
  accelerate the convergence. The idea was proposed by Dr. Dam Hieu Chi
  when he stayed in AIST at summer 2002.

  Log of Simple_Mixing.c:

     10/Dec/2002  Released by T.Ozaki

*************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Simple_Mixing(int state_num,
                   double rho0[ASIZE1],
                   double rho1[ASIZE1],
                   double rho2[ASIZE1],
                   double Rrho2[ASIZE1])
{
  static int i;
  static double Mix_wgt,Mix_wgt2;
  static double Nrm,s0,s1,s2,tmp0;
  static double Min_Weight,Max_Weight;

  Min_Weight = Min_Mixing_weight;
  Max_Weight = Max_Mixing_weight;
  
  /****************************************************
                        NormRD
  ****************************************************/

  for (i=0; i<Grid_Num; i++){
    Rrho2[i] = rho0[i] - rho1[i];
  }

  Nrm = 0.0;

  s0 = Rrho2[0]*Rrho2[0]*MRV[0]*MRV[0]*MRV[0]
     + Rrho2[Grid_Num-1]*Rrho2[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + Rrho2[i]*Rrho2[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + Rrho2[i]*Rrho2[i]*MRV[i]*MRV[i]*MRV[i];
  }
  Nrm = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  for (i=4; 1<=i; i--){
    NormRD[i] = NormRD[i-1];
    History_Uele[i] = History_Uele[i-1];
  }
  NormRD[0] = Nrm;
  History_Uele[0] = Eeigen[state_num];

  /****************************************************
                Change Mixing_weight
  ****************************************************/

  if (sgn(History_Uele[0]-History_Uele[1])
      ==sgn(History_Uele[1]-History_Uele[2])
      && NormRD[0]<NormRD[1]){

    tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

    if (tmp0<Max_Weight){
      if (Min_Weight<tmp0) Mixing_weight = tmp0;
      else                 Mixing_weight = Min_Weight;
    }
    else                   Mixing_weight = Max_Weight;
  }
   
  else if (sgn(History_Uele[0]-History_Uele[1])
	   ==sgn(History_Uele[1]-History_Uele[2])
	   && NormRD[1]<NormRD[0]){

    tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

    if (tmp0<Max_Weight){
      if (Min_Weight<tmp0) Mixing_weight = tmp0;
      else                 Mixing_weight = Min_Weight;
    }
    else                   Mixing_weight = Max_Weight;

  }

  else if (sgn(History_Uele[0]-History_Uele[1])
	   !=sgn(History_Uele[1]-History_Uele[2])
	   && NormRD[0]<NormRD[1]){

    tmp0 = NormRD[1]/(largest(NormRD[1]-NormRD[0],10e-10))*Mixing_weight;

    if (tmp0<Max_Weight){
      if (Min_Weight<tmp0) Mixing_weight = tmp0;
      else                 Mixing_weight = Min_Weight;
    }
    else                   Mixing_weight = Max_Weight;
  }

  else if (sgn(History_Uele[0]-History_Uele[1])
	   !=sgn(History_Uele[1]-History_Uele[2])
	   && NormRD[1]<NormRD[0]){

    tmp0 = NormRD[1]/(largest(NormRD[1]+NormRD[0],10e-10))*Mixing_weight;

    if (tmp0<Max_Weight){
      if (Min_Weight<tmp0) Mixing_weight = tmp0;
      else                 Mixing_weight = Min_Weight;
    }
    else                   Mixing_weight = Max_Weight;
  }

  Mix_wgt = Mixing_weight;

  /****************************************************
                        Mixing
  ****************************************************/
 
  Mix_wgt2 = 1.0 - Mix_wgt;
  for (i=0; i<Grid_Num; i++){
    rho0[i] = Mix_wgt*rho0[i] + Mix_wgt2*rho1[i];
    rho2[i] = rho1[i];
    rho1[i] = rho0[i];
  }
}
