/**********************************************************************
  Make_EDPP.c:

     Make_EDPP.c is a subroutine to generate environment dependent
     pseudo potentials. 

  Log of Make_EDPP.c:

     6/Sep/2004  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

static void SCF_NL_Equation(int cs,double dE[ASIZE15], double EZ); 
static double Solve_NL_Equation(int MCP_on, int cs, int so, int n, int l, int mul, 
                                double ene, int *node,
                                double U[ASIZE1], 
                                double Vh_NL[ASIZE1], double Vxc_NL[ASIZE1],
                                double Vextra[ASIZE1]);
static double Solve_NL_Equation2(int MCP_on, int cs, int so, int n, int l, int mul, 
                                 double ene, int *node,
                                 double U[ASIZE1], 
                                 double Vh_NL[ASIZE1], double Vxc_NL[ASIZE1]);
static double phi0_phi1(double phi0[ASIZE1], double phi1[ASIZE1], double rmax);
static double Int_RadF(double R, double RadF[ASIZE1]);
static void Simpson_Normalization(double U[ASIZE1]);
static double Func_dEdQ(double Q);
static double Func_Eextra(double Q);
static double CalcQ(int cs);
static void CubicSpline1(int cs0, int cs1, double Eextra[ASIZE15],
                         double GradEext[ASIZE15]);

static double Int_EDPP_Proj[ASIZE15];
static double Evps[ASIZE15][ASIZE13][ASIZE4];
static double CoesQ[ASIZE15+1];

void Make_EDPP()
{
  static int i,cs,cs0,po,SCF,SCF_MAX;
  static double scaleA;
  static double GradEext[ASIZE15];

  static double dE[ASIZE15];
  static double dE1;
  static double dE2;
  static double InpE1;
  static double InpE2;

  static int num_div,div,po_div0,po_div;
  static int Min_cs,num_plus,num_minus;
  static int found_switch[ASIZE15];
  static int Plus_Order[ASIZE15];
  static int Minus_Order[ASIZE15];
  static double MinC,scale1,StepE,StartE;

  /* find plus charge states */
  for (cs=0; cs<charge_state_num; cs++){
    found_switch[cs] = 0;
  }

  Plus_Order[0] = 0;
  num_plus = 1;  
  do{
    MinC = 1.0e+100;
    po = 1;
    for (cs=1; cs<charge_state_num; cs++){
      if (  found_switch[cs]==0
         && 0.0<charge_states[cs]
         && fabs(charge_states[cs])<MinC){

	MinC = fabs(charge_states[cs]);
        Min_cs = cs;
        po = 0;
        found_switch[cs] = 1;
      }
    }

    if (po==0){
      num_plus++;
      Plus_Order[num_plus-1] = Min_cs;
    }
  }
  while(po==0);

  /* find minus charge states */
  for (cs=0; cs<charge_state_num; cs++){
    found_switch[cs] = 0;
  }

  Minus_Order[0] = 0;
  num_minus = 1;  
  do{
    MinC = 1.0e+100;
    po = 1;
    for (cs=1; cs<charge_state_num; cs++){
      if (  found_switch[cs]==0
         && charge_states[cs]<0.0
         && fabs(charge_states[cs])<MinC){

	MinC = fabs(charge_states[cs]);
        Min_cs = cs;
        po = 0;
        found_switch[cs] = 1;
      }
    }

    if (po==0){
      num_minus++;
      Minus_Order[num_minus-1] = Min_cs;
    }
  }
  while(po==0);

  /*
  for (i=0; i<num_plus; i++){
    printf("plus  %2d %2d\n",i,Plus_Order[i]);  
  }
  for (i=0; i<num_minus; i++){
    printf("minus %2d %2d\n",i,Minus_Order[i]);  
  }
  exit(0);
  */

  /***************************************************
     find Eextra and GradEext for plus charge states
  ***************************************************/

  Eextra[0]   = 0.0;
  GradEext[0] = 0.0;
  SCF_NL_Equation(0, dE, GradEext[0]);
  Int_EDPP_Proj[0] = CalcQ(0);
  printf("cs=%2d Int_EDPP_Proj=%15.12f\n",0,Int_EDPP_Proj[0]);

  for (i=1; i<num_plus; i++){

    cs0 = Plus_Order[i-1]; 

    /* first trial */
    cs = Plus_Order[i];
    Eextra[cs]   = 0.0;
    GradEext[cs] = 0.0;
    SCF_NL_Equation(cs, dE, GradEext[cs]);
    Int_EDPP_Proj[cs] = CalcQ(cs);

    InpE1 = Eextra[cs];
    dE1 = dE[cs];

    /* second trial */

    scale1 = 0.04;
    po_div0 = 0;

    do{

      po_div = 0; 
      num_div = 6; 
      StepE = 2.0*scale1*fabs(dE[cs])/(double)num_div;
      StartE = -scale1*fabs(dE[cs]);
      div = -1;

      do{ 
	div++;

	InpE1 = InpE2;
	dE1 = dE2;

	Eextra[cs] = StartE + (double)div*StepE;
	CubicSpline1(cs0,cs,Eextra,GradEext);
	SCF_NL_Equation(cs, dE, GradEext[cs]);
	Int_EDPP_Proj[cs] = CalcQ(cs);

	InpE2 = Eextra[cs];
	dE2 = dE[cs];

	if (div!=0 && dE1*dE2<0.0){
	  po_div = 1;
	}

	printf("div=%2d Eextra=%15.12f dE=%15.12f\n",div,Eextra[cs],dE[cs]);


      }
      while (div<num_div && po_div==0);

      if (po_div==0) scale1 = scale1*2.0;
      else           po_div0 = 1; 
    }
    while (po_div0==0);


    InpE2 = Eextra[cs];
    dE2 = dE[cs];

    printf("InpE1=%15.12f InpE2=%15.12f\n",InpE1,InpE2); 
    printf("dE1=%15.12f dE2=%15.12f\n",dE1,dE2); 






 
    SCF = 0;
    SCF_MAX = 10;
    po = 0;
    scaleA = 1.2;

    /* Self consistent */
    do{
      SCF++;

      if (dE1*dE2<0.0){
        Eextra[cs] = (InpE1*dE2-InpE2*dE1)/(dE2-dE1);
      }
      else{
        if (fabs(dE1)<fabs(dE2)){
          Eextra[cs] = InpE1
                      - sgn((dE2-dE1)/(InpE2-InpE1))*fabs((InpE2-InpE1))*scaleA;
	}
        else{
          Eextra[cs] = InpE2
                      - sgn((dE2-dE1)/(InpE2-InpE1))*fabs((InpE2-InpE1))*scaleA;
        }

	/*
        scaleA = 1.4*scaleA;
	*/

      }
      
      CubicSpline1(cs0,cs,Eextra,GradEext);
      SCF_NL_Equation(cs, dE, GradEext[cs]);
      Int_EDPP_Proj[cs] = CalcQ(cs);

      /*
      if (0.0<dE1*dE[cs] && 0.0<dE2*dE[cs]){

        if (fabs(dE1)<fabs(dE2)){
          InpE2 = Eextra[cs];
          dE2 = dE[cs]; 
        }
        else{ 
          InpE1 = Eextra[cs];
          dE1 = dE[cs]; 
	}
      }
      else if (0.0<dE1*dE[cs]) {
        InpE1 = Eextra[cs];
        dE1 = dE[cs]; 
      }
      else {
        InpE2 = Eextra[cs];
        dE2 = dE[cs]; 
      }
      */

  
    } while (po==0 && SCF<SCF_MAX);  

  }

}

void CubicSpline1(int cs0, int cs1, double Eextra[ASIZE15], double GradEext[ASIZE15])
{
  static int i,j;
  static double x[ASIZE6],y[ASIZE6],y0p;
  static double Amat[ASIZE6][ASIZE6];
  static double Cvec[ASIZE6];

  x[0] = Int_EDPP_Proj[cs0];
  x[1] = Int_EDPP_Proj[cs1];

  y[0] = Eextra[cs0];
  y[1] = Eextra[cs1];

  y0p = GradEext[cs0];

  for (i=0; i<=1; i++){
    for (j=0; j<=2; j++){
      Amat[i][j] = pow(x[i],(double)j);
    }
  }

  Amat[2][0] = 0.0;
  Amat[2][1] = 1.0;
  Amat[2][2] = 2.0*x[0];

  Amat[0][3] = y[0];
  Amat[1][3] = y[1];
  Amat[2][3] = y0p;

  Gauss_LEQ(2,Amat,Cvec);
  GradEext[cs1] = Cvec[1] + 2.0*Cvec[2]*x[1];    

  /*
  printf("c %15.12f %15.12f %15.12f\n",Cvec[0],Cvec[1],Cvec[2]);
  printf("x0=%15.12f x1=%15.12f\n",x[0],x[1]);
  printf("E0=%15.12f E1=%15.12f\n",Eextra[cs0],Eextra[cs1]);
  printf("GradE0=%15.12f GradE1=%15.12f\n",GradEext[cs0],GradEext[cs1]);
  */

}


double CalcQ(int cs)
{
  static int i;
  static double s0,s1,s2,tmp;

  /* calculate Int_EDPP_Proj */

  s0 = Vedpp[cs][0]*MRV[0]*MRV[0]*MRV[0]
    + Vedpp[cs][Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 += Vedpp[cs][i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 += Vedpp[cs][i]*MRV[i]*MRV[i]*MRV[i];
  }
  tmp = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  return tmp;  
}




void SCF_NL_Equation(int cs, double dE[ASIZE15], double EZ)
{
  static int so,L,mul,m,n,i,po,refine_loop;
  static int SCF,SCF_OK,num_ene,l,node,p;
  static double sum,mixw,tmp,alp;
  static double dene,ene,Trial_DD;
  static double DD_min,DD_max,mc;
  static double rc,r,ep,sum_p,sum_c;
  static double MinE,MaxE,tmp0;
  static double drho_p,drho_c;
  static double s0,s1,s2,XC;
  static double Vsl[ASIZE1];
  static double VNLp[ASIZE1];
  static double rho_NL[ASIZE1];
  static double rho_p[ASIZE1];
  static double Vh_NL[ASIZE1];
  static double Vxc_NL[ASIZE1];
  static double Vextra[ASIZE1];
  static double NLW[ASIZE13][ASIZE2][ASIZE1];

  /* set initial valence charge and potential */

  i = Restart_load(cs);
  Density_V(cs);
  Hartree(rho_V[cs],Vh_NL);

  if (PCC_switch==0){
    if      (XC_switch==0) XC_CA(rho_V[cs],Vxc_NL,1);
    else if (XC_switch==1) XC4atom_PBE(rho_V[cs],Vxc_NL,1);
  }
  else {
    Density_PCC(rho_PCC,cs);
    for (i=0; i<Grid_Num; i++){
      rho[1][i] = rho_V[cs][i] + rho_PCC[i];
    }
    if      (XC_switch==0) XC_CA(rho[1],Vxc_NL,1);
    else if (XC_switch==1) XC4atom_PBE(rho[1],Vxc_NL,1);
  }

  for (i=0; i<Grid_Num; i++){
    Vh_V[i]  = Vh_NL[i];
    Vxc_V[i] = Vxc_NL[i];
    rho_NL[i] = rho_V[cs][i];
  }  

  /*
  Q = Int_EDPP_Proj[cs];
  EZ = Func_dEdQ(Q);
  printf("cs=%2d EZ=%15.12f\n",cs,EZ);
  */

  for (i=0; i<Grid_Num; i++){
    r = MRV[i];
    Vextra[i] = EZ*erf(r/Core_R[0])/r;
  }  

  /* start SCF calculation */

  SCF = 0;
  SCF_MAX = 27;
  SCF_OK = 0;
  mixw = 0.3;

  do{
    SCF++;

    for (so=0; so<SOI_switch; so++){
      for (L=0; L<ASIZE2; L++){

	for (mul=0; mul<NumVPS_L[L]; mul++){
          
	  m = GI_VPS[L][mul];
	  n = NVPS[m];
	  rc = VPS_Rcut[m];

          if (mul==0){
            DD_min = Solve_NL_Equation(1,cs,so,n,L,mul,E[cs][so][n][L],&node,NLW[so][m],
                                       Vh_NL,Vxc_NL,Vextra);
	  }

          if (SCF==1){   
  	    MinE = E[cs][so][n][L] - 0.20;
	    MaxE = smallest(E[cs][so][n][L]+0.20,Vinf-1.0e-14);
	  }
          else{
  	    MinE = Evps[cs][so][m] - 0.20;
	    MaxE = smallest(Evps[cs][so][m]+0.20,Vinf-1.0e-14);
	  }

          DD_min = Solve_NL_Equation(0,cs,so,n,L,mul,MinE,&node,NLW[so][m],Vh_NL,Vxc_NL,Vextra);
          DD_max = Solve_NL_Equation(0,cs,so,n,L,mul,MaxE,&node,NLW[so][m],Vh_NL,Vxc_NL,Vextra); 
         

          /****************************************************
 	             Refinement of the eigenvalue by 
                          a bisection method
	  ****************************************************/

          po = 0; 
          refine_loop = 0;

          do {

	    if (5<refine_loop && DD_max*DD_min<0.0){
	      ene = (MinE*DD_max - MaxE*DD_min)/(DD_max - DD_min);
	    }
	    else{ 
	      ene = 0.5*(MinE + MaxE);
	    }

            Trial_DD = Solve_NL_Equation(0,cs,so,n,L,mul,ene,&node,NLW[so][m],Vh_NL,Vxc_NL,Vextra);

            if (0.0<DD_min*Trial_DD){
              DD_min = Trial_DD;
              MinE = ene;
	    }
            else{
              DD_max = Trial_DD;
              MaxE = ene;
            }

	    /*
            printf("cs=%2d SCF=%2d so=%2d L=%2d mul=%2d ene=%15.12f Trial_DD=%15.12f\n",
                    cs,SCF,so,L,mul,ene,Trial_DD);
	    */

            if (fabs(Trial_DD)<1.0e-8) po = 1;

            refine_loop++;
	  } while(po==0 && refine_loop<50);          

          Evps[cs][so][m] = ene;
          
	} /* mul */
      } /* L */
    } /* so */  

    /****************************************************
                 Reconstruction of potential
    ****************************************************/

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      sum = 0.0;
      for (so=0; so<SOI_switch; so++){
	for (m=0; m<Number_VPS; m++){
	  n = NVPS[m];
	  l = LVPS[m];
	  sum = sum + OcpN[cs][so][n][l]*NLW[so][m][i]*NLW[so][m][i]/4.0/PI/r/r;
	}
      }

      rho_p[i] = rho_NL[i];
      rho_NL[i] = (1.0-mixw)*rho_NL[i] + mixw*sum;
    }

     sum = 0.0;
     for (i=0; i<Grid_Num; i++){
       tmp = rho_NL[i] - rho_p[i];
       sum += tmp*tmp*MRV[i]*MRV[i]*MRV[i]*dx;
     }      
     drho_p = drho_c;
     drho_c = sqrt(sum);

     if (drho_c<1.0e-5 ) SCF_OK = 1; 

     if (drho_p<drho_c) mixw = mixw/2.0;
     else               mixw = 1.1*mixw;
     if (0.7<mixw)      mixw = 0.7;

     /*
     printf("cs=%2d SCF =%2d drho=%15.12f\n",cs,SCF,drho_c);
     */

     Hartree(rho_NL,Vh_NL);

     if (PCC_switch==0){
       for (i=0; i<Grid_Num; i++){
	 rho[1][i] = rho_NL[i];
       }
     }
     else {
       for (i=0; i<Grid_Num; i++){
	 rho[1][i] = rho_NL[i] + rho_PCC[i];
       }
     }
     if      (XC_switch==0) XC_CA(rho[1],Vxc_NL,1);
     else if (XC_switch==1) XC4atom_PBE(rho[1],Vxc_NL,1);

  } while(SCF<SCF_MAX && SCF_OK==0);

  /****************************************************
     store Vedpp
  ****************************************************/

  Gaussian_Alp = 2.0;

  /*
  printf("ABC Av_CoreR=%15.12f\n",Av_CoreR);
  */

  alp = 1.0/Core_R[0]/Core_R[0];

  for (i=0; i<Grid_Num; i++){
    r = MRV[i];
    mc = pow(alp/PI,1.5)*exp(-alp*r*r);
    Vedpp[cs][i] = mc*(Vh_NL[i] + Vlocal[i]);
  }

  /****************************************************
     calculation of total energy for valence electrons
  ****************************************************/

  /****************************************************
                    Hartree energy
             0.5*4pi\int rho(r)*VH(r)*r*r dr
            = 2pi\int rho(r)*VH(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho_NL[0]*Vh_NL[0]*MRV[0]*MRV[0]*MRV[0]
     + rho_NL[Grid_Num-1]*Vh_NL[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho_NL[i]*Vh_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho_NL[i]*Vh_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }
  EHart_VPS[cs] = 2.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
  
  /****************************************************
      Exchange-correlation part (this is not energy)
              4pi\int rho(r)*Vxc(r)*r*r dr
            = 4pi\int rho(r)*Vxc(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho_NL[0]*Vxc_NL[0]*MRV[0]*MRV[0]*MRV[0]
     + rho_NL[Grid_Num-1]*Vxc_NL[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho_NL[i]*Vxc_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho_NL[i]*Vxc_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }
  XC = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
              Exchange-correlation energy
              4pi\int rho(r)*Exc(r)*r*r dr
            = 4pi\int rho(r)*Exc(r)*r*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  if      (XC_switch==0) XC_CA(rho[1],Vxc_NL,0);
  else if (XC_switch==1) XC4atom_PBE(rho[1],Vxc_NL,0);

  s0 = rho[1][0]*Vxc_NL[0]*MRV[0]*MRV[0]*MRV[0]
     + rho[1][Grid_Num-1]*Vxc_NL[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[1][i]*Vxc_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[1][i]*Vxc_NL[i]*MRV[i]*MRV[i]*MRV[i];
  }
  Exc_VPS[cs] = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
            Electron-local part coulomb energy
              4pi*Z*\int rho(r)*r dr
            = 4pi*Z*\int rho(r)*r*r dx
        Radial integrations by the Simpson method
  ****************************************************/

  s0 = rho_NL[0]*Vlocal[0]*MRV[0]*MRV[0]*MRV[0]
     + rho_NL[Grid_Num-1]*Vlocal[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho_NL[i]*Vlocal[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho_NL[i]*Vlocal[i]*MRV[i]*MRV[i]*MRV[i];
  }
  Eec_VPS[cs] = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
               energy fo non-local part
  ****************************************************/

  Enl_VPS[cs] = 0.0;
  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){
      for (mul=0; mul<NumVPS_L[L]; mul++){

        m = GI_VPS[L][mul];
        n = NVPS[m];

        sum = 0.0;
        for (p=0; p<projector_num[L]; p++){
          tmp = phi0_phi1(VNL_W2[so][L][p],NLW[so][m],MRV[Grid_Num-1]);  
          sum += OcpN[cs][so][n][L]*proj_ene[0][so][L][p]*tmp*tmp;
	}  
        Enl_VPS[cs] += sum;
      }
    }
  }

  /****************************************************
   Eeigen
  ****************************************************/

  Eeigen_VPS[cs] = 0.0;
  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){
      for (mul=0; mul<NumVPS_L[L]; mul++){

        m = GI_VPS[L][mul];
        n = NVPS[m];
        Eeigen_VPS[cs] += OcpN[cs][so][n][L]*Evps[cs][so][m];

        printf("so=%2d L=%2d mul=%2d Evps=%15.12f\n",so,L,mul,Evps[cs][so][m]);
      }
    }
  }

  /****************************************************
   Eextra2
  ****************************************************/

  s0 = rho_NL[0]*Vextra[0]*MRV[0]*MRV[0]*MRV[0]
     + rho_NL[Grid_Num-1]*Vextra[Grid_Num-1]*
       MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 += rho_NL[i]*Vextra[i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 += rho_NL[i]*Vextra[i]*MRV[i]*MRV[i]*MRV[i];
  }
  Eextra2[cs] = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /****************************************************
   Kinetic energy = Eeigen - (2*EHart + XC + Eec + Enl + Eextra) 
  ****************************************************/

  Ekin_VPS[cs] = Eeigen_VPS[cs] - (2.0*EHart_VPS[cs]+XC+Eec_VPS[cs]+Enl_VPS[cs]+Eextra2[cs]);

  /****************************************************
     Total energy0 = Ekin + EHart + Exc + Eec + Enl
  ****************************************************/

  Etot0_VPS[cs] = Ekin_VPS[cs] + EHart_VPS[cs] + Exc_VPS[cs] + Eec_VPS[cs] + Enl_VPS[cs];

  /****************************************************
     Total energy = Total energy0 + Eextra
  ****************************************************/

  Etot_VPS[cs] = Etot0_VPS[cs] + Eextra[cs];

  /****************************************************
     difference of Energy
  ****************************************************/

  dE[cs] = (Etot[cs]-Etot[0]) - (Etot_VPS[cs]-Etot_VPS[0]);

  /****************************************************
      standard output 
  ****************************************************/

  /*
  printf("<EDPP>  **** Energies of pseudo atom charge state=%2d ****\n",cs);
  printf("<EDPP>  Eeigen  = %19.12f (Hartree)\n",Eeigen_VPS[cs]);
  printf("<EDPP>  Ekin    = %19.12f (Hartree)\n",Ekin_VPS[cs]);
  printf("<EDPP>  EHart   = %19.12f (Hartree)\n",EHart_VPS[cs]);
  printf("<EDPP>  Exc     = %19.12f (Hartree)\n",Exc_VPS[cs]);
  printf("<EDPP>  Eec     = %19.12f (Hartree)\n",Eec_VPS[cs]);
  printf("<EDPP>  Enl     = %19.12f (Hartree)\n",Enl_VPS[cs]);
  printf("<EDPP>  Eextra  = %19.12f (Hartree)\n",Eextra[cs]);
  printf("<EDPP>  Etot0   = %19.12f (Hartree)\n",Etot0_VPS[cs]);
  printf("<EDPP>  Etot    = %19.12f (Hartree)\n",Etot_VPS[cs]);
  */

  printf("<EDPP>  **** Energies of pseudo atom charge state=%2d ****\n",cs);
  printf("<EDPP>  Eextra2 = %19.12f (Hartree)\n",Eextra2[cs]);
  printf("<EDPP>  Eextra  = %19.12f (Hartree)\n",Eextra[cs]);
  printf("<EDPP>  Etot0   = %19.12f (Hartree)\n",Etot0_VPS[cs]);
  printf("<EDPP>  Etot    = %19.12f (Hartree)\n",Etot_VPS[cs]);

} 

double Func_Eextra(double Q)
{
  static int cs;
  static double sum;

  sum = 0.0;
  for (cs=0; cs<=charge_state_num; cs++){
    sum += CoesQ[cs]*pow(Q,(double)cs); 
  }

  return sum;
}

double Func_dEdQ(double Q)
{
  static int cs;
  static double sum;

  sum = 0.0;
  for (cs=1; cs<=charge_state_num; cs++){
    sum += (double)cs*CoesQ[cs]*pow(Q,(double)(cs-1.0)); 
  }

  return sum;
}


double Solve_NL_Equation(int MCP_on, int cs, int so, int n, int l, int mul,
                         double ene, int *node, 
                         double U[ASIZE1], 
                         double Vh_NL[ASIZE1], double Vxc_NL[ASIZE1],
                         double Vextra[ASIZE1])
{
  static int SCF,SCF_Max,SCF_po,m,i,j,p;
  static int is,ie;
  static double Trial_D,DD_min,DD_max;
  static double Deriv_O,Deriv_I;
  static double Reduce_Num_grid,kappa,mixing;
  static double tmp0,sum_p,sum_c,ratio,r;
  static double dsum,dsum_p;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double Vsl[ASIZE1],Up[ASIZE1];
  static double VNLp[ASIZE1];

  kappa= 0.0;
  Reduce_Num_grid = RNG[n][l];

  /****************************************
        semi-local pseudo potentials
  ****************************************/

  m = GI_VPS[l][0]; 
  for (i=0; i<Grid_Num; i++){
    Vsl[i] = VPS[0][so][m][i] + Vh_V[i] + Vxc_V[i];
  }

  Hamming_O(0,l,ene,kappa,Mo,Lo,DMo,Vsl,Vsl,Reduce_Num_grid);
  Hamming_I(0,l,ene,kappa,Mi,Li,DMi,Vsl,Vsl,Reduce_Num_grid);

  /* Matching point */
  if (LL_grid<=UL_grid){
    if (LL_grid<=MCP[so][n][l] && MCP[so][n][l]<=UL_grid){
      j = MCP[so][n][l];
    }
    else{
      j = (LL_grid + UL_grid)/2;
    }
  }
  else{
    printf("error(2)\n");
    exit(0);
  }

  ratio = Lo[j]/Li[j];
  for (i=0; i<=j; i++){
    U[i] = pow(MRV[i],(double)l+1.0)*Lo[i];
  }
  for (i=(j+1); i<Grid_Num; i++){
    U[i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
  }
  Simpson_Normalization(U);

  /*
  for (i=0; i<Grid_Num; i++){
    U[i] = W2[so][m][i];
  }
  */

  ratio = Lo[j]/Li[j];
  Deriv_O = Mo[j];
  Deriv_I = ratio*Mi[j];
  Trial_D = Deriv_O - Deriv_I;

  /*
  if (l==0){
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f\n",MRV[i],U[i]/MRV[i],W2[so][m][i]/MRV[i]);
  }
  exit(0);
  }
  */

  /*
  printf("Semi Trial_D=%15.12f\n",Trial_D);
  */

  /****************************************
      fully separable pseudo potentials
  ****************************************/

  SCF_po = 0;
  SCF_Max = 1000;
  SCF = 0;
  mixing = 0.0001;

  /* screened Vlocal */
  for (i=0; i<Grid_Num; i++){
    Vsl[i] = Vlocal[i] + Vh_NL[i] + Vxc_NL[i] + Vextra[i];
  }
 
  /* SCF iteration loop */
  do{
    SCF++; 

    /* sum over projectors */
    for (i=0; i<Grid_Num; i++) VNLp[i] = 0.0;

    sum_p = sum_c;
    sum_c = 0.0; 
    for (p=0; p<projector_num[l]; p++){
    
      /* <v|phi_i> */
      tmp0 = phi0_phi1(VNL_W2[so][l][p],U,MRV[Grid_Num-1]);
    
      sum_c = sum_c + tmp0;
      /* |v>*c*<v|phi_i> */
      for (i=0; i<Grid_Num; i++){
	r = MRV[i];
	VNLp[i] += VNL_W2[so][l][p][i]/r*proj_ene[0][so][l][p]*tmp0;
      }
    }
     
    /*
    printf("so=%2d n=%2d l=%2d SCF=%2d sum_c=%15.12f\n",so,n,l,SCF,sum_c); 
    */

    Hamming_O(1,l,ene,kappa,Mo,Lo,DMo,Vsl,VNLp,Reduce_Num_grid);
    Hamming_I(1,l,ene,kappa,Mi,Li,DMi,Vsl,VNLp,Reduce_Num_grid);

    /* Matching point */
    if (LL_grid<=UL_grid){
      if (LL_grid<=MCP[so][n][l] && MCP[so][n][l]<=UL_grid){
	j = MCP[so][n][l];
      }
      else{
        j = (LL_grid + UL_grid)/2;
      }
    }
    else{
      printf("error(2)\n");
      exit(0);
    }

    /* calculate Up */
    ratio = Lo[j]/Li[j];
    for (i=0; i<=j; i++){
      Up[i] = pow(MRV[i],(double)l+1.0)*Lo[i];
    }
    for (i=(j+1); i<Grid_Num; i++){
      Up[i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
    }

    /* simple mixing */
    for (i=0; i<Grid_Num; i++){
      U[i] = mixing*Up[i] + (1.0-mixing)*U[i];
    }

    dsum_p = dsum; 
    dsum = fabs(sum_p-sum_c);

    if ( dsum<1.0e-12 ) SCF_po = 1;

    if (dsum_p<dsum) mixing = mixing/2.0;
    else             mixing = 1.1*mixing;
    if (0.8<mixing)  mixing = 0.8;

  } while(SCF_po==0 && SCF<SCF_Max);

  if (SCF_po==0){
    printf("dsum is large: dsum=%15.12f\n",dsum);
  }

  /* find the difference in the derivative of wave functions */

  Simpson_Normalization(U);

  if (MCP_on==1){

    Trial_D = 1.0e+20;

    if (0<(j-500)) is = j - 500;
    else           is = 1;   
    if ((j+500)<(Grid_Num-1)) ie = j + 500;
    else                      ie = Grid_Num;   

    for (i=is; i<ie; i++){

      ratio = Lo[i]/Li[i];
      Deriv_O = Mo[i];
      Deriv_I = ratio*Mi[i];
      if (fabs(Deriv_O-Deriv_I)<fabs(Trial_D)){
        Trial_D = Deriv_O - Deriv_I;
        j = i;
      }
    }
    MCP[so][n][l] = j;
  }  
  else{ 
    ratio = Lo[j]/Li[j];
    Deriv_O = Mo[j];
    Deriv_I = ratio*Mi[j];
    Trial_D = Deriv_O - Deriv_I;
  }

  /* Recalculation of the number of nodes */
  *node = 0;
  j = 0;
  for (i=1; i<=(Grid_Num-2); i++){
    if (U[i]*U[i+1]<0.0) {
      (*node)++;
      j++;
    }
  }


  /*
  if (l==0){
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f\n",MRV[i],U[i]/MRV[i],W2[so][m][i]/MRV[i]);
  }
  exit(0);
  }
  */


  /*
  printf("node=%2d %2d   Trial_D=%15.12f\n",*node,j,Trial_D);
  */

  /*
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f\n",MRV[i],W2[so][0][i]/MRV[i]);
  }
  exit(0);

  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f\n",MRV[i],U[i]/MRV[i]);
  }
  exit(0);
  */

  return Trial_D;
}

void Simpson_Normalization(double U[ASIZE1])
{
  static int i;
  static double tmp0,sumG,s0,s1,s2;

  s0 = U[0]*U[0]*MRV[0] + U[Grid_Num-1]*U[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + U[i]*U[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + U[i]*U[i]*MRV[i];
  }
  sumG = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  tmp0 = 1.0/sqrt(sumG);

  for (i=0; i<=(Grid_Num-1); i++){
    U[i] = tmp0*U[i];
  }
}

double phi0_phi1(double phi0[ASIZE1], double phi1[ASIZE1], double rmax)
{
  static int i,n,j,l,nf,fg;
  static double r,rmin,Sr,Dr,sum,dum;

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  Sr = rmax + rmin;
  Dr = rmax - rmin;
  
  sum = 0.0;
  for (i=0; i<=(n-1); i++){
    r = 0.50*(Dr*x[i] + Sr);
    sum += 0.5*Dr*w[i]*Int_RadF(r,phi0)*Int_RadF(r,phi1);
  }

  return sum;
}


double Int_RadF(double R, double RadF[ASIZE1])
{
  static int mp_min,mp_max,m,po;
  static double h1,h2,h3,f1,f2,f3,f4;
  static double g1,g2,x1,x2,y1,y2,f;
  static double rm,y12,y22,df,a,b;
  static double result;

  mp_min = 0;
  mp_max = Grid_Num - 1;
 
  if (MRV[Grid_Num-1]<R){
    f = 0.0;
  }
  else if (R<MRV[0]){
    po = 1;
    m = 4;
    rm = MRV[m];

    h1 = MRV[m-1] - MRV[m-2];
    h2 = MRV[m]   - MRV[m-1];
    h3 = MRV[m+1] - MRV[m];

    f1 = RadF[m-2];
    f2 = RadF[m-1];
    f3 = RadF[m];
    f4 = RadF[m+1];

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = rm - MRV[m-1];
    x2 = rm - MRV[m];
    y1 = x1/h2;
    y2 = x2/h2;
    y12 = y1*y1;
    y22 = y2*y2;

    f =  y22*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y12*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

    df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y22*(2.0*f2 + h2*g1)/h2
       + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
       - y12*(2.0*f3 - h2*g2)/h2;

    a = 0.5*df/rm;
    b = f - a*rm*rm;      
    f = a*R*R + b;
  } 
  else{
    do{ 
      m = (mp_min + mp_max)/2;
      if (MRV[m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    /****************************************************
                   Spline like interpolation
    ****************************************************/

    if (m!=1){
      h1 = MRV[m-1] - MRV[m-2];
    }

    h2 = MRV[m]   - MRV[m-1];

    if (m!=(Grid_Num-1)){
      h3 = MRV[m+1] - MRV[m];
    }

    if (m!=1){
      f1 = RadF[m-2];
    }

    f2 = RadF[m-1];
    f3 = RadF[m];

    if (m!=(Grid_Num-1)){
      f4 = RadF[m+1];
    }

    /****************************************************
                   Treatment of edge points
    ****************************************************/

    if (m==1){
      h1 = -(h2+h3);
      f1 = f4;
    }
    if (m==(Grid_Num-1)){
      h3 = -(h1+h2);
      f4 = f1;
    }

    /****************************************************
                Calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - MRV[m-1];
    x2 = R - MRV[m];

    y1 = x1/h2;
    y2 = x2/h2;

    f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
      + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }

  result = f;
  return result; 
}


double Solve_NL_Equation2(int MCP_on, int cs, int so, int n, int l, int mul,
                          double ene, int *node, 
                          double U[ASIZE1], 
                          double Vh_NL[ASIZE1], double Vxc_NL[ASIZE1])
{
  static int SCF,SCF_Max,SCF_po,m,i,j,p,q;
  static int is,ie;
  static double Trial_D,DD_min,DD_max;
  static double Deriv_O,Deriv_I,sum;
  static double Reduce_Num_grid,kappa,mixing;
  static double tmp0,sum_p,sum_c,ratio,r;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double Vsl[ASIZE1],Up[ASIZE1];
  static double U0[ASIZE1];
  static double VNLp[ASIZE1];
  static double W[ASIZE4][ASIZE1];
  static double Amat[ASIZE6][ASIZE6];
  static double Cvec[ASIZE6];

  kappa= 0.0;
  Reduce_Num_grid = RNG[n][l];

  /****************************************
                homo equation
  ****************************************/

  m = GI_VPS[l][0]; 
  for (i=0; i<Grid_Num; i++){
    Vsl[i] = Vlocal[i] + Vh_NL[i] + Vxc_NL[i];
  }

  Hamming_O(0,l,ene,kappa,Mo,Lo,DMo,Vsl,Vsl,Reduce_Num_grid);
  Hamming_I(0,l,ene,kappa,Mi,Li,DMi,Vsl,Vsl,Reduce_Num_grid);

  /* Matching point */
  if (LL_grid<=UL_grid){
    if (LL_grid<=MCP[so][n][l] && MCP[so][n][l]<=UL_grid){
      j = MCP[so][n][l];
    }
    else{
      j = (LL_grid + UL_grid)/2;
    }
  }
  else{
    printf("error(2)\n");
    exit(0);
  }

  ratio = Lo[j]/Li[j];
  for (i=0; i<=j; i++){
    U0[i] = pow(MRV[i],(double)l+1.0)*Lo[i];
  }
  for (i=(j+1); i<Grid_Num; i++){
    U0[i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
  }
  Simpson_Normalization(U0);

  /*
  for (i=0; i<Grid_Num; i++){
    U[i] = W2[so][m][i];
  }
  */

  ratio = Lo[j]/Li[j];
  Deriv_O = Mo[j];
  Deriv_I = ratio*Mi[j];
  Trial_D = Deriv_O - Deriv_I;
  printf("Semi Trial_D=%15.12f\n",Trial_D);

  /****************************************
               inhomo equation
  ****************************************/

  for (p=0; p<projector_num[l]; p++){

    for (i=0; i<Grid_Num; i++){
      VNLp[i] = VNL_W2[so][l][p][i]/MRV[i];
    }

    Hamming_O(1,l,ene,kappa,Mo,Lo,DMo,Vsl,VNLp,Reduce_Num_grid);

    for (i=0; i<Grid_Num; i++){
      W[p][i] = pow(MRV[i],(double)l+1.0)*Lo[i];
    }
  }

  for (p=0; p<projector_num[l]; p++){
    for (q=0; q<projector_num[l]; q++){
      Amat[p][q] = proj_ene[0][so][l][p]*phi0_phi1(VNL_W2[so][l][p],W[q],MRV[Grid_Num-1]);
      if (p==q) Amat[p][q] += 1.0;
    }
  }
  
  for (p=0; p<projector_num[l]; p++){
    Amat[p][projector_num[l]] = proj_ene[0][so][l][p]*phi0_phi1(VNL_W2[so][l][p],U0,MRV[Grid_Num-1]);
  }

  for (p=0; p<=projector_num[l]; p++){
    for (q=0; q<=projector_num[l]; q++){
      printf("%15.12f ",Amat[p][q]); 
    }
    printf("\n");
  }

  Gauss_LEQ(projector_num[l]-1,Amat,Cvec);     

  for (p=0; p<projector_num[l]; p++){
    printf("p=%2d %15.12f\n",p,Cvec[p]);
  }

  for (i=0; i<Grid_Num; i++){
    sum = 0.0;  
    for (p=0; p<projector_num[l]; p++){
      sum += Cvec[p]*W[p][i];
    }
    U[i] = U0[i] + sum;
  }
  Simpson_Normalization(U);


  /* Recalculation of the number of nodes */
  *node = 0;
  j = 0;
  for (i=1; i<=(Grid_Num-2); i++){
    if (U[i]*U[i+1]<0.0) {
      (*node)++;
      j++;
    }
  }

  if (l==0){
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f\n",MRV[i],U[i]/MRV[i],W2[so][m][i]/MRV[i]);
  }
  exit(0);
  }


  /*
  printf("node=%2d %2d   Trial_D=%15.12f\n",*node,j,Trial_D);
  */

  /*
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f\n",MRV[i],W2[so][0][i]/MRV[i]);
  }
  exit(0);

  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f\n",MRV[i],U[i]/MRV[i]);
  }
  exit(0);
  */

  return Trial_D;
}


