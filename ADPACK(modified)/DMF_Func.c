/**********************************************************************
  DMF_Func.c:

     DMF_Func.c is a subroutine to calculate the derivative of MF.

  Log of DMF_Func.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vps_pao.h"

double DMF_Func(int cp, double ep, int NL,
                double M, double L, double VPL[BHS_YOUSO1])
{
  static double d1,d2,result,MQ,c,al,d3,dV;

  /*-- Non-relativistic --*/

  d1 = -(2.0*NL+1.0)*M;
  d2 = 2.0*MRV[cp]*MRV[cp]*(VPL[cp]-ep)*L;

  /*-- Relativistic --*/ 

  /*
    c = 137.0359895;
    al = 1.0/c;
    MQ = 1.0 + (ep-VPL[cp])/2.0/c/c; 
    dV = (VPL[cp+1]-VPL[cp-1])/dx/2.0/MRV[cp];
    d3 = al*al*MRV[cp]/2.0/MQ*dV;
    d1 = -(2.0*NL + 1.0 + d3)*M;
    d2 = (2.0*MRV[cp]*MRV[cp]*MQ*(VPL[cp]-ep)-d3*NL)*L;
  */

  result = d1 + d2;

  return result;
}

