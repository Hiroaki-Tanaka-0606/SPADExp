/**********************************************************************
  XC_EX.c:

     XC_EX.c is a subroutine to calculate exchange energy density
     and potential by a generalized gradient correction.

  Log of XC_EX.c:

     10/Dec/2002  Released by T.Ozaki, (but this is originally adapted
                  by J.M.Soler for SIESTA code from routine "velect" of
                  Froyen's pseudopotential generation program) 

***********************************************************************/

#include <math.h>
#include "adpack.h"

#define ZERO   0.00
#define ONE    1.00
#define PFIVE  0.50
#define OPF    1.50
#define C014   0.0140
#define TRD    0.333333333333333
#define FTRD   1.333333333333333
#define TFTM   0.519842099789746   /* 2**FTRD-2       */
#define A0     0.521061761197848   /* (4/(9*PI))**TRD */

void XC_EX(int NSP, double DS0, double DS[2], double EX[1], double VX[2])
{
  double den_min,ALP,D0,D1,D,Z,FZ,FZP;
  double RS,VXP,EXP,VXF,EXF;

  den_min = 1.0e-40;
  ALP = 2.0*TRD;

  if (NSP==2){
    D0 = largest(DS[0],den_min);
    D1 = largest(DS[1],den_min);
    D = D0 + D1;
    Z = (D0 - D1)/D;
    FZ = (pow(1.0 + Z,FTRD) + pow(1 - Z,FTRD) - 2.0)/TFTM;
    FZP = FTRD*(pow(1.0 + Z,TRD) - pow(1.0 - Z,TRD))/TFTM;
  }    
  else{ 
    D = largest(DS0,den_min);
    Z = ZERO;
    FZ = ZERO;
    FZP = ZERO;
  }
  RS = pow(3.0/(4.0*PI*D),TRD);
  VXP = -(3.0*ALP/(2.0*PI*A0*RS));
  EXP = 3.0*VXP/4.0;
  VXF = pow(2.0,TRD)*VXP;
  EXF = pow(2.0,TRD)*EXP;

  if (NSP==2){
    VX[0] = VXP + FZ*(VXF - VXP) + (1.0 - Z)*FZP*(EXF - EXP);
    VX[1] = VXP + FZ*(VXF - VXP) - (1.0 + Z)*FZP*(EXF - EXP);
    EX[0] = EXP + FZ*(EXF - EXP);
  }
  else{
    VX[0] = VXP;
    EX[0] = EXP;
  } 
}


