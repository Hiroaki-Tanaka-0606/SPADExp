/**********************************************************************
  Total_Energy.c:

     Total_Energy.c is a subroutine to calculate the total energy
     of atom.

  Log of All_Electron.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Total_Energy(int state_num)
{
  static int i,n,l,so;
  static double s0,s1,s2,XC,tmp1,sum;
    
  /****************************************************
                    Hartree energy
             0.5*4pi\int rho(r)*VH(r)*r*r dr
            = 2pi\int rho(r)*VH(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho[0][0]*Vh[0]*MRV[0]*MRV[0]*MRV[0]
     + rho[0][Grid_Num-1]*Vh[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[0][i]*Vh[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[0][i]*Vh[i]*MRV[i]*MRV[i]*MRV[i];
  }
  EHart[state_num] = 2.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
      Exchange-correlation part (this is not energy)
              4pi\int rho(r)*Vxc(r)*r*r dr
            = 4pi\int rho(r)*Vxc(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho[0][0]*Vxc[0]*MRV[0]*MRV[0]*MRV[0]
     + rho[0][Grid_Num-1]*Vxc[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[0][i]*Vxc[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[0][i]*Vxc[i]*MRV[i]*MRV[i]*MRV[i];
  }
  XC = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
              Exchange-correlation energy
              4pi\int rho(r)*Exc(r)*r*r dr
            = 4pi\int rho(r)*Exc(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  if      (XC_switch==0) XC_CA(rho[0],Vxc,0);
  else if (XC_switch==1) XC4atom_PBE(rho[0],Vxc,0);
  else if (XC_switch==2) XC_VWN(rho[0],Vxc,0);

  s0 = rho[0][0]*Vxc[0]*MRV[0]*MRV[0]*MRV[0]
     + rho[0][Grid_Num-1]*Vxc[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[0][i]*Vxc[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[0][i]*Vxc[i]*MRV[i]*MRV[i]*MRV[i];
  }
  Exc[state_num] = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
        Calculate exchange-correlation potential
  ****************************************************/

  if      (XC_switch==0) XC_CA(rho[0],Vxc,1);
  else if (XC_switch==1) XC4atom_PBE(rho[0],Vxc,1);
  else if (XC_switch==2) XC_VWN(rho[0],Vxc,1);

  /****************************************************
            Electron-core coulomb energy
              4pi*Z*\int rho(r)*r dr
            = 4pi*Z*\int rho(r)*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho[0][0]*MRV[0]*MRV[0]
     + rho[0][Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[0][i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[0][i]*MRV[i]*MRV[i];
  }
  Eec[state_num] = -4.0*PI*(double)AtomNum*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
               Kohn-Sham eigenvalue energy
                        \sum ep
  ****************************************************/

  Eeigen[state_num] = 0.0;
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
        if (0.0<OcpN[state_num][so][n][l]){ 
          Eeigen[state_num] += OcpN[state_num][so][n][l]*E[state_num][so][n][l];
        }
      }
    }
  }

  /****************************************************
     Kinetic energy = Eeigen - (2*EHart + XC + Eec) 
  ****************************************************/

  /* frozen core approximation */
  if (Calc_Type==3 && state_num==1){

    for (so=0; so<SOI_switch; so++){
      for (n=1; n<=max_ocupied_N; n++){
	for (l=0; l<n; l++){

	  if (Relaxed_Switch[n][l]==0){ 
	    Energy_Kin[1][so][n][l] = Energy_Kin[0][so][n][l]; 
	  }
          else if (0.0<OcpN[state_num][so][n][l] && Relaxed_Switch[n][l]==1){

	    if (TwoComp_frag==0){      /* one-component representation */
	      s0 = PF[so][n][l][0]*PF[so][n][l][0]*V[0]*MRV[0]
	         + PF[so][n][l][Grid_Num-1]*PF[so][n][l][Grid_Num-1]*V[Grid_Num-1]*MRV[Grid_Num-1];
	    }
	    else{                      /* two-component representation */
	      s0 = (PF[so][n][l][0]*PF[so][n][l][0]+FF[so][n][l][0]*FF[so][n][l][0])*V[0]*MRV[0]
	         + (PF[so][n][l][Grid_Num-1]*PF[so][n][l][Grid_Num-1]
                   +FF[so][n][l][Grid_Num-1]*FF[so][n][l][Grid_Num-1])*V[Grid_Num-1]*MRV[Grid_Num-1];
	    }

	    s1 = 0.0;
	    for (i=1; i<=(Grid_Num-2); i=i+2){
 	      if (TwoComp_frag==0){      /* one-component representation */
 	        s1 += PF[so][n][l][i]*PF[so][n][l][i]*V[i]*MRV[i];
	      }
	      else{                      /* two-component representation */
 	        s1 += (PF[so][n][l][i]*PF[so][n][l][i]+FF[so][n][l][i]*FF[so][n][l][i])*V[i]*MRV[i];
	      }
	    }

	    s2 = 0.0; 
	    for (i=2; i<=(Grid_Num-3); i=i+2){
 	      if (TwoComp_frag==0){      /* one-component representation */
	        s2 += PF[so][n][l][i]*PF[so][n][l][i]*V[i]*MRV[i];
	      }
	      else{                      /* two-component representation */
	        s2 += (PF[so][n][l][i]*PF[so][n][l][i]+FF[so][n][l][i]*FF[so][n][l][i])*V[i]*MRV[i];
	      }
	    }

	    tmp1 = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
	    Energy_Kin[state_num][so][n][l] = E[state_num][so][n][l] - tmp1;
          }
	}
      }
    }
  }

  else {
    for (so=0; so<SOI_switch; so++){
      for (n=1; n<=max_ocupied_N; n++){
	for (l=0; l<n; l++){
	  if (0.0<OcpN[state_num][so][n][l]){ 

	    if (TwoComp_frag==0){      /* one-component representation */
	      s0 = PF[so][n][l][0]*PF[so][n][l][0]*V[0]*MRV[0]
	         + PF[so][n][l][Grid_Num-1]*PF[so][n][l][Grid_Num-1]*V[Grid_Num-1]*MRV[Grid_Num-1];
	    }
	    else{                      /* two-component representation */
	      s0 = (PF[so][n][l][0]*PF[so][n][l][0]+FF[so][n][l][0]*FF[so][n][l][0])*V[0]*MRV[0]
	         + (PF[so][n][l][Grid_Num-1]*PF[so][n][l][Grid_Num-1]
                   +FF[so][n][l][Grid_Num-1]*FF[so][n][l][Grid_Num-1])*V[Grid_Num-1]*MRV[Grid_Num-1];
	    }

	    s1 = 0.0;
	    for (i=1; i<=(Grid_Num-2); i=i+2){
 	      if (TwoComp_frag==0){      /* one-component representation */
 	        s1 += PF[so][n][l][i]*PF[so][n][l][i]*V[i]*MRV[i];
	      }
	      else{                      /* two-component representation */
 	        s1 += (PF[so][n][l][i]*PF[so][n][l][i]+FF[so][n][l][i]*FF[so][n][l][i])*V[i]*MRV[i];
	      }
	    }

	    s2 = 0.0; 
	    for (i=2; i<=(Grid_Num-3); i=i+2){
 	      if (TwoComp_frag==0){      /* one-component representation */
	        s2 += PF[so][n][l][i]*PF[so][n][l][i]*V[i]*MRV[i];
	      }
	      else{                      /* two-component representation */
	        s2 += (PF[so][n][l][i]*PF[so][n][l][i]+FF[so][n][l][i]*FF[so][n][l][i])*V[i]*MRV[i];
	      }
	    }

	    tmp1 = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
	    Energy_Kin[state_num][so][n][l] = E[state_num][so][n][l] - tmp1;
	  }
	}
      }
    }
  }

  Ekin[state_num] = 0.0;
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
        if (0.0<OcpN[state_num][so][n][l]){ 
          Ekin[state_num] += OcpN[state_num][so][n][l]*Energy_Kin[state_num][so][n][l];
	}
      }
    }
  }

  /*
  Ekin[state_num] = Eeigen[state_num] - (2.0*EHart[state_num] + XC + Eec[state_num]);
  printf("sum=%15.12f  Ekin=%15.12f\n",sum,Ekin[state_num]);
  */

  /****************************************************
        Total energy = Ekin + EHart + Exc + Eec
  ****************************************************/
  
  Etot[state_num] = Ekin[state_num] + EHart[state_num] + Exc[state_num] + Eec[state_num];

  printf("<ALL>  **** Energies of atom ****\n");
  printf("<ALL>  Ekin   = %19.12f (Hartree)\n",Ekin[state_num]);
  printf("<ALL>  EHart  = %19.12f (Hartree)\n",EHart[state_num]);
  printf("<ALL>  Exc    = %19.12f (Hartree)\n",Exc[state_num]);
  printf("<ALL>  Eec    = %19.12f (Hartree)\n",Eec[state_num]);
  printf("<ALL>  Etot   = %19.12f (Hartree)\n",Etot[state_num]);
  printf("<ALL>  Eeigen = %19.12f (Hartree)\n\n",Eeigen[state_num]);
  /* printf("%18.12f\n",Eeigen[state_num]-EHart[state_num]+Exc[state_num]-XC); */

}
