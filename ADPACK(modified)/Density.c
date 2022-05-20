/**********************************************************************
  Density.c:

     Density.c is a subroutine to calculate charge densities.

  Log of Density.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

#define kB            0.00008617251324000000   /* eV/K           */   

void Density(int state_num)
{
  static int i,j,k,l,n,so,po,loopN;
  static double sum,r,p,p2,f,f2;
  static double Beta,ChemP,ChemP_MIN,ChemP_MAX;
  static double E_Temp,degeneracy;
  static double e,Num_State,max_e=50.0;
  static double FermiF,Dnum;

  /**************************************
    If EDPP, update OcpN
  **************************************/
  
  if (VPP_switch==3 && Calc_Type==1) {

    total_electron = (double)AtomNum - charge_states[state_num];

    po = 0;
    loopN = 0;

    E_Temp = 3000.0/27.211386; /* 3000(K) */
    Beta = 1.0/kB/E_Temp;
    ChemP_MAX =  15.0;
    ChemP_MIN = -15.0;
 
    if (SOI_switch==1) degeneracy = 2.0;
    else               degeneracy = 1.0;

    do { 
      ChemP = 0.50*(ChemP_MAX + ChemP_MIN);
      Num_State = 0.0;

      for (so=0; so<SOI_switch; so++){
        for (n=1; n<=max_ocupied_N; n++){
          for (l=0; l<n; l++){

	    if (0.0<OcpN[state_num][so][n][l]){
              e = (E[state_num][so][n][l] - ChemP)*Beta;
              if (e<=-max_e) e = -max_e;
              if (max_e<=e)  e =  max_e;
              FermiF = 1.0/(1.0 + exp(e));
              OcpN[state_num][so][n][l] = degeneracy*(2.0*(double)l+1.0)*FermiF;
              Num_State = Num_State + OcpN[state_num][so][n][l];

	      /*
              printf("so=%2d n=%2d l=%2d E=%15.12f e=%15.12f FermiF=%15.12f\n",
                      so,n,l,E[state_num][so][n][l],e,FermiF); 
	      */
            }

          }
        }
      }
 
      Dnum = ((double)AtomNum - Num_State) - charge_states[state_num];
      if (0.0<=Dnum) ChemP_MIN = ChemP;
      else           ChemP_MAX = ChemP;
      if (fabs(Dnum)<1.0e-11) po = 1;
      loopN++;

    } while (po==0 && loopN<1000);
  }

  /**************************************
        calculate electron density
  **************************************/

  for (i=0; i<Grid_Num; i++){
    r = MRV[i];
    sum = 0.0;
    for (so=0; so<SOI_switch; so++){
      for (n=1; n<=max_ocupied_N; n++){
	for (l=0; l<n; l++){

	  if (0.0<OcpN[state_num][so][n][l]){ 

	    if (TwoComp_frag==0){      /* one-component representation */
	      p = PF[so][n][l][i];
	      p2 = p*p;
	      sum = sum + OcpN[state_num][so][n][l]*p2/4.0/PI/r/r;
	    }

	    else{                      /* two-component representation */
	      p = PF[so][n][l][i];
	      p2 = p*p;
	      f = FF[so][n][l][i];
	      f2 = f*f;
	      sum = sum + OcpN[state_num][so][n][l]*(p2+f2)/4.0/PI/r/r;
	    }

	  }
	}  
      }
    }
    rho[0][i] = sum;
  }

  /*
  for(i=0;i<Grid_Num;i+=Grid_Num/100) {
    printf("%f ", PF[0][1][0][i]);
  }
  printf("\n");
  exit(0);
  */

  /*
  for(i=0;i<Grid_Num;i+=Grid_Num/10) {
    printf("%f ", rho[0][i]);
  }
  printf("\n");
  exit(0);
  */


}
