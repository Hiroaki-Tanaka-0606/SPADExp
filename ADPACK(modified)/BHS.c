/**********************************************************************
  BHS.c:

     BHS.c is a subroutine to calculate norm conseving pseudo
     potentials using BHS's scheme.

  Log of BHS.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

static void Simpson(int so, int m);

void BHS(int state_num)
{
  static int i,j,k,l,n,m,po,poi,po_node;
  static int loop_max,loop_num,so;
  static double a,b,c,sum;
  static double ld,rp,f,g,cl[ASIZE2]; 
  static double ep,dum,kappa;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double ratio,Deriv_O,Deriv_I,DD_new,DD_old,kei;
  static double W1[ASIZE2][ASIZE1];
  static double Gamma[ASIZE2],Delta[ASIZE2];
  static double dum2; 

  ld = 3.50;
  loop_max = 10000;

  for (so=0; so<SOI_switch; so++){
    for (m=0; m<Number_VPS; m++){

      printf("<VPS>  Generating pseudo potentials %2d by the BHS scheme....\n",m);

      n = NVPS[m];
      l = LVPS[m];

      ep = E[state_num][so][n][l];
      cl[m] = -10.0;
      kei = 0.1;
      po = 0;
      poi= 0;
      po_node = 0;

      do{

	for (i=0; i<Grid_Num; i++){
	  rp = MRV[i]/VPS_Rcut[m];
	  dum = pow(rp,ld);
	  f = exp(-dum);
	  V1[m][i] = (1.0-f)*V[i] + cl[m]*f;
	}

	Hamming_O(0,l,ep,kappa,Mo,Lo,DMo,V1[m],V1[m],RNG[n][l]);

	if (node==0 || po_node==1){

	  po_node = 1;

	  Hamming_I(0,l,ep,kappa,Mi,Li,DMi,V1[m],V1[m],RNG[n][l]);
	  ratio = Lo[MCP[so][n][l]]/Li[MCP[so][n][l]];
	  Deriv_O = (Lo[MCP[so][n][l]+1]-Lo[MCP[so][n][l]-1])/2.0/dx;
	  Deriv_I = ratio*(Li[MCP[so][n][l]+1]-Li[MCP[so][n][l]-1])/2.0/dx;

	  if (Check_switch==1){ 
	    printf("BHS n=%2d l=%2d Deriv_O-Deriv_I=%20.17f\n",
		   n,l,Deriv_O-Deriv_I);
	  }

	  if (fabs(Deriv_O-Deriv_I)<=10e-12){
	    po = 1;
	  }
	  else{

	    DD_old = DD_new;
	    DD_new = Deriv_O-Deriv_I;
	    if (poi==0) DD_old = DD_new;

	    cl[m] = cl[m] + kei*(Deriv_I-Deriv_O);

	    if ((fabs(DD_old)-fabs(DD_new))<0.0){
	      kei = kei - 0.02;
	    } 
	    else if ((fabs(DD_new-DD_old)/DD_new)<0.02){
	      kei = kei + 0.02;
	    }

	    poi++;
	  }
	}
	else{
	  cl[m] = cl[m] + 0.2;
	}

	loop_num++;      

	if (loop_max<loop_num){
	  printf("********************************************************\n");
	  printf("Failure in searching of an appropriate cl[m]\n");
	  printf("in the pseudo potential generation by the BHS scheme.\n");
	  printf("Please check cutoff radii of both pseudo potentials and\n");
	  printf("pseudo atomic oribtals.\n");
	  printf("********************************************************\n");
	  exit(0); 
	}
      } while(po==0);


      for (i=0; i<=MCP[so][n][l]; i++){
	W1[m][i] = pow(MRV[i],l+1.0)*Lo[i];
      }
      for (i=(MCP[so][n][l]+1); i<Grid_Num; i++){
	W1[m][i] = pow(MRV[i],l+1.0)*ratio*Li[i];
      }

      /****************************************************
	   Gamma*W1l = Pl  ->   Gamma = Pl/W1l
      ****************************************************/

      Gamma[m] = PF[so][n][l][Grid_Num-1]/W1[m][Grid_Num-1];

      /****************************************************
	   Delta
      ****************************************************/

      sum = 0.0;
      for (i=0; i<Grid_Num; i++){
	sum = sum + W1[m][i]*W1[m][i]*MRV[i];
      }
      a = sum*dx;

      sum = 0.0;
      for (i=0; i<Grid_Num; i++){
	rp = MRV[i]/VPS_Rcut[m];
	dum = pow(rp,ld);
	g = exp(-dum)*pow(MRV[i],l+1.0);
	sum = sum + W1[m][i]*g*MRV[i];
      }
      b = sum*dx;

      sum = 0.0;
      for (i=0; i<Grid_Num; i++){
	rp = MRV[i]/VPS_Rcut[m];
	dum = pow(rp,ld);
	g = exp(-dum)*pow(MRV[i],l+1.0);
	sum = sum + g*g*MRV[i];
      }
      c = sum*dx;

      Delta[m] = (-b+sqrt(b*b-c*(a-1.0/Gamma[m]/Gamma[m])))/c;

      /****************************************************
	   W2l
      ****************************************************/

      for (i=0; i<Grid_Num; i++){
	rp = MRV[i]/VPS_Rcut[m];
	dum = pow(rp,ld);
	g = exp(-dum)*pow(MRV[i],l+1.0);
	W2[so][m][i] = Gamma[m]*(W1[m][i]+Delta[m]*g);
      }
      Simpson(so, m);

      /****************************************************
	   V2
      ****************************************************/

      for (i=0; i<Grid_Num; i++){
	rp = MRV[i]/VPS_Rcut[m];
	a = pow(rp,ld);
	g = exp(-a)*pow(MRV[i],l+1.0);
	dum = Gamma[m]*Delta[m]*g/W2[so][m][i];
	b = pow(rp,2.0*ld); 
	dum2 = (ld*ld*b - (2.0*ld*l+ld*(ld+1.0))*a)/2.0/MRV[i]/MRV[i];
	V2[so][m][i] = V1[m][i] + dum*(E[state_num][so][n][l] - V1[m][i] + dum2);
      }

    } /* m */
  } /* so */ 

  Density_V(state_num,OcpN);
  Hartree(rho_V[state_num],Vh_V);
  
  if (PCC_switch==0){ 
    if      (XC_switch==0) XC_CA(rho_V[state_num],Vxc_V,1);
    else if (XC_switch==1) XC4atom_PBE(rho_V[state_num],Vxc_V,1);
  }
  else {
    Density_PCC(rho_PCC,state_num,OcpN);

    for (i=0; i<Grid_Num; i++){
      rho[1][i] = rho_V[state_num][i] + rho_PCC[i];
    }
    if      (XC_switch==0) XC_CA(rho[1],Vxc_V,1);
    else if (XC_switch==1) XC4atom_PBE(rho[1],Vxc_V,1);
  }

  for (so=0; so<SOI_switch; so++){
    for (m=0; m<Number_VPS; m++){
      for (i=0; i<Grid_Num; i++){
        VPS[0][so][m][i] = V2[so][m][i] - Vh_V[i] - Vxc_V[i];
      }
    }
  }

  /* for std output */
  printf("\n");
}


void Simpson(int so, int m)
{

  /****************************************************
           Normalization by Simpson Integration
  ****************************************************/

  static int i;
  static double s0,s1,s2,sum,sum2;

  sum = 0.0;
  s0 = W2[so][m][0]*W2[so][m][0]*MRV[0]
     + W2[so][m][Grid_Num-1]*W2[so][m][Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + W2[so][m][i]*W2[so][m][i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + W2[so][m][i]*W2[so][m][i]*MRV[i];
  }
  sum = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  if (0.0<=W2[so][m][Grid_Num-200]){
    sum2 = sqrt(sum);
    for (i=0; i<=(Grid_Num-1); i++){
      W2[so][m][i] = W2[so][m][i]/sum2;
    }
  }
  else {
    sum2 = sqrt(sum);
    for (i=0; i<=(Grid_Num-1); i++){
      W2[so][m][i] = -W2[so][m][i]/sum2;
    }
  }

}
