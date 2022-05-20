/**********************************************************************
  All_Electron.c:

   All_Electron.c is a subroutine to perform the self-consistent
   calculation of an atomic Kohn-Sham equation including all electrons.

  Log of All_Electron.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

static void Check_Sum_Num_Electrons();
static void Simpson_Normalization(int so, int n, int l);
static double Quadratic_Interpolate(double x1, double x2, double x3,
	                            double y1, double y2, double y3);

static double Check_Sum;

void All_Electron(int state_num)
{
  static int i0,i1,imin,p_node,so;
  static int i,j,k,n,l,po,SCF,SCF_OK,po_node,CTP_P;
  static int cNode,pNode,Node_min,Node_max,MatP_fix_switch;
  static int Enum,Divided_Num,Reduce_po;
  static double minMapD,tmp0,tmp1,tmp2,cDD,pDD,pDD2,pOEE,Fugou;
  static double DD_max,DD_min,di,Reduce_Num_grid;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double ratio,Deriv_O,Deriv_I,ep,sum,dG,r,te;
  static double kei,DD_old,DD_new,s0,s1,s2,OEE,MaxE,MinE;
  static double Min_ep,Max_ep,Min_D,Max_D,Trial_ep,Trial_D;
  static double p_ep,p_D,cri,scale,al,c,kappa;
  static double StepE,E_old[ASIZE13][ASIZE3][ASIZE2];

	/* beginning of the addition */
	int CalcType0=0;
	if(Calc_Type==0){
		CalcType0=1;
		Calc_Type=2;
	}
	/* end of the addition */
	
  MatP_fix_switch = 0;
  
  if (Calc_Type==0 || Calc_Type==1 || Calc_Type==3)
    CTP_P = 40;
  else 
    CTP_P = 20;
      
  for (so=0; so<SOI_switch; so++){
		for (n=1; n<=max_ocupied_N; n++){
			for (l=0; l<n; l++){
        tmp0 = (double)AtomNum/(double)n;
        E[state_num][so][n][l] = -0.5*tmp0*tmp0 - 10.0;
      }
    }
  }

  /* calculation of frozen core */
  if (Calc_Type==3 && state_num==1)
		for (so=0; so<SOI_switch; so++){
			for (n=1; n<=max_ocupied_N; n++){
				for (l=0; l<n; l++){
        E[1][so][n][l] = E[0][so][n][l];
      }
    }
  }

  Initial_Density();
  
  SCF = 1;
  SCF_OK = 0;

	for (n=1; n<=max_ocupied_N; n++){
		for (l=0; l<n; l++){
      RNG[n][l] = 1.0; 
    }
  }

  /* reduction of the grid range for atoms with a heavy weight */

  if (83<=AtomNum){  
    if (1.9<=Grid_Xmax){
      RNG[1][0] = (log(4.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
      RNG[2][0] = (log(6.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
    if (2.3<=Grid_Xmax){
      RNG[2][1] = (log(7.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
      RNG[3][0] = (log(7.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
      RNG[3][1] = (log(9.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
  }

  else if (50<=AtomNum){  
    if (2.2<=Grid_Xmax){
      RNG[1][0] = (log(8.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
      RNG[2][0] = (log(9.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
    if (2.3<=Grid_Xmax){
      RNG[2][1] = (log(10.0)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
  }

  /* reduction of the grid range for 1s state */

  else if (15<=AtomNum) {
    if (2.8<=Grid_Xmax){
      RNG[1][0] = (log(16.445)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
  }

  else if (AtomNum!=0) {
    if (3.1<=Grid_Xmax){
      RNG[1][0] = (log(22.198)-Grid_Xmin)/(Grid_Xmax-Grid_Xmin);
    }
  }

  /* start SCF calculation */

	
  do{
    Core(SCF,(double)AtomNum,Vcore);
    Hartree(rho[0],Vh);
    if (SCF==1){
      XC_Xa(rho[0],Vxc);
    }
    else{
      if      (XC_switch==0) XC_CA(rho[0],Vxc,1);
      else if (XC_switch==1) XC4atom_PBE(rho[0],Vxc,1);
      else if (XC_switch==2) XC_VWN(rho[0],Vxc,1);
    }

    VP();

    /****************************************************
             search eigenvalues for n and l  
    ****************************************************/

    for (so=0; so<SOI_switch; so++){
			for (n=1; n<=max_ocupied_N; n++){
				for (l=0; l<n; l++){
          MatP[so][n][l] = -100;
        }
      }
    }
		
    for (so=0; so<SOI_switch; so++){
			for (n=1; n<=max_ocupied_N; n++){
				for (l=0; l<n; l++){

          if (l==0){
            kappa = -1.0;    /* for j=l+1/2 */
          } 
          else{ 
            if      (so==0)  kappa = (double)(-l-1);    /* for j=l+1/2 */
            else if (so==1)  kappa = (double)l;         /* for j=l-1/2 */
	  }

	  if ( (Calc_Type==3
                &&
                state_num==1
                &&
                0.0<OcpN[state_num][so][n][l]
                &&
                Relaxed_Switch[n][l]==1)
              ||
               (Calc_Type==3   
                &&
                state_num==0
                &&
                0.0<OcpN[state_num][so][n][l])
              ||
               (Calc_Type!=3   
                &&
                0.0<OcpN[state_num][so][n][l])
	       ){ 

	    /****************************************************
 	      a lower bound and a upper bound of the solution 
	    ****************************************************/

	    po = 0;
	    scale = 1.0;
	    Reduce_Num_grid = RNG[n][l];
	    Reduce_po = 0;

	    i = 0;
	    do{

	      i++;

	      if (SCF==1){
		tmp0 = (double)AtomNum/(double)n;
		MinE = -0.5*tmp0*tmp0 - 3.0;
		MaxE = smallest(MinE+100.0,Vinf-1.0e-14);
                if (50<=AtomNum && n<=3){  
  	  	  MaxE = smallest(MinE+300.0,Vinf-1.0e-14);
		}
	      }
	      else{
		if (SCF==2){
		  MinE = E[state_num][so][n][l] - 5.0;
		  MaxE = smallest(MinE+100.0,Vinf-1.0e-14);
		}
		else{
		  tmp2 = 10.0*fabs(E_old[so][n][l]-E[state_num][so][n][l]);
		  if (tmp2<1.0e-5) tmp2 = 1.0e-5;
		  tmp2 = scale*tmp2;  
		  MinE = E[state_num][so][n][l] - tmp2;
		  MaxE = E[state_num][so][n][l] + tmp2;
		  if (Vinf<MaxE) MaxE = Vinf - 1.0e-14;
		}
	      }

	      /* MinE */

	      if (Equation_Type==0){      /* Schrodinger equation */
		Hamming_O(0,l,MinE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_min = node;
		Hamming_I(0,l,MinE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==1){ /* scalar relativistic equation */
		Hamming_O(2,l,MinE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_min = node;
		Hamming_I(2,l,MinE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==2){ /* full relativistic equation */
		Hamming_O(3,l,MinE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_min = node;
		Hamming_I(3,l,MinE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }

	      /* Matching point */

	      if (LL_grid<=UL_grid){
		if (LL_grid<=MatP[so][n][l] && MatP[so][n][l]<=UL_grid
		    && MatP_fix_switch==1){
		  j = MatP[so][n][l];
		}
		else{ 
		  if (LL_grid<=(CTP+CTP_P) && (CTP+CTP_P)<=UL_grid){
		    j = CTP + CTP_P;
		    if (UL_grid<j) j = (LL_grid + UL_grid)/2;
		  }
		  else{
		    j = (LL_grid + UL_grid)/2;
		  }
		  MatP[so][n][l] = j;
		}
	      }
	      else{
		Reduce_po = 1;
		Reduce_Num_grid = 0.70*Reduce_Num_grid;
		if (Reduce_Num_grid<RNG[n][l]){
		  RNG[n][l] = Reduce_Num_grid;
		  if (Check_switch==1){ 
		    printf("A n=%i  l=%i  Reduce_Num_grid=%5.2f\n",n,l,RNG[n][l]);
		  }
		}
	      }

	      /* Deviation between the derivatives */

	      if (Reduce_po==0){
		ratio = Lo[j]/Li[j];
		Deriv_O = Mo[j];
		Deriv_I = ratio*Mi[j];
		DD_min = Deriv_O - Deriv_I;
	      }

	      /* MaxE */

	      if (Equation_Type==0){      /* Schrodinger equation */ 
		Hamming_O(0,l,MaxE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_max = node;
		Hamming_I(0,l,MaxE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==1){ /* scalar relativistic equation */
		Hamming_O(2,l,MaxE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_max = node;
		Hamming_I(2,l,MaxE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==2){ /* full relativistic equation */
		Hamming_O(3,l,MaxE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Node_max = node;
		Hamming_I(3,l,MaxE,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }

	      /* matching point */

	      if (LL_grid<=UL_grid){
		if (LL_grid<=MatP[so][n][l] && MatP[so][n][l]<=UL_grid
		    && MatP_fix_switch==1){
		  j = MatP[so][n][l];  
		}
		else{ 
		  if (LL_grid<=(CTP+CTP_P) && (CTP+CTP_P)<=UL_grid){
		    j = CTP + CTP_P;
		    if (UL_grid<j) j = (LL_grid + UL_grid)/2;
		  }
		  else{
		    j = (LL_grid + UL_grid)/2;
		  }
		  MatP[so][n][l] = j;
		}
	      }
	      else{
		Reduce_po = 1;
		Reduce_Num_grid = 0.70*Reduce_Num_grid;
		if (Reduce_Num_grid<RNG[n][l]){
		  RNG[n][l] = Reduce_Num_grid;
		  if (Check_switch==1){ 
		    printf("B n=%i  l=%i  Reduce_Num_grid=%5.2f\n",n,l,RNG[n][l]);
		  }
		}
	      }

	      if (Reduce_po==0){
		ratio = Lo[j]/Li[j];
		Deriv_O = Mo[j];
		Deriv_I = ratio*Mi[j];
		DD_max = Deriv_O - Deriv_I;
	      }

	      if (Reduce_po==1){
		Reduce_po = 0;
	      }

	      else if (0.0<DD_max*DD_min && fabs(MaxE-MinE)<0.1 && 
		       abs(Node_min-(n-l-1))<=1 && abs(Node_max-(n-l-1))<=1 &&
		       (Node_min==(n-l-1) || Node_min==(n-l-1))
		       ||
		       0.0<DD_max*DD_min && ((n-l-1)<Node_min || Node_max<(n-l-1))
		       ){

		scale = 100.0*scale; 
	      }
	      else{
		po = 1;
		Max_ep = MaxE;
		Min_ep = MinE;
		Max_D = DD_max;
		Min_D = DD_min;
	      }

	    }while (po==0 && i<10);

	    /****************************************************
		 refinement of the eigenvalue by 
		 the linear interpolation
	    ****************************************************/

	    if (Check_switch==1){ 
	      printf("i=%i MinE=%10.6f MaxE=%10.6f Min_D=%10.6f Max_D=%10.6f ",
		     i,MinE,MaxE,Min_D,Max_D);
	      printf("Node_min=%2d Node_max=%2d\n",Node_min,Node_max);
	    }

	    i = 0;
	    po = 0;
	    do{

	      i++; 

	      if (1<i && i<20 && Max_D*Min_D<0.0 &&
		  abs(Node_min-(n-l-1))<=1 && abs(Node_max-(n-l-1))<=1 &&
		  (Node_min==(n-l-1) || Node_min==(n-l-1))
		  ){

 	        Trial_ep = (Min_ep*Max_D - Max_ep*Min_D)/(Max_D - Min_D);
	      }
	      else{ 
		Trial_ep = 0.50*(Min_ep + Max_ep);
	      }

	      if (Equation_Type==0){       /* Schrodinger equation */
		Hamming_O(0,l,Trial_ep,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Hamming_I(0,l,Trial_ep,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==1){  /* scalar relativistic equation */
		Hamming_O(2,l,Trial_ep,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Hamming_I(2,l,Trial_ep,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }
	      else if (Equation_Type==2){  /* full relativistic equation */
		Hamming_O(3,l,Trial_ep,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
		Hamming_I(3,l,Trial_ep,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
	      }

	      /* Matching point */

	      if (LL_grid<=UL_grid){
		if (LL_grid<=MatP[so][n][l] && MatP[so][n][l]<=UL_grid
		    && MatP_fix_switch==1){
		  j = MatP[so][n][l];  
		}
		else{ 
		  if (LL_grid<=(CTP+CTP_P) && (CTP+CTP_P)<=UL_grid){
		    j = CTP + CTP_P;
		    if (UL_grid<j) j = (LL_grid + UL_grid)/2;
		  }
		  else{
		    j = (LL_grid + UL_grid)/2;
		  }
		  MatP[so][n][l] = j;
		}
	      }
	      else{
		Reduce_Num_grid = 0.70*Reduce_Num_grid;
		if (Reduce_Num_grid<RNG[n][l]){
		  RNG[n][l] = Reduce_Num_grid;
		  if (Check_switch==1){ 
		    printf("C n=%i  l=%i  Reduce_Num_grid=%5.2f\n",n,l,RNG[n][l]);
		  }
		}
	      }

	      ratio = Lo[j]/Li[j];
	      Deriv_O = Mo[j];
	      Deriv_I = ratio*Mi[j];
	      Trial_D = Deriv_O - Deriv_I;

	      if (Check_switch==1){
		printf("SCF=%2d i=%2d n=%2d l=%2d T_ep=%16.12f T_D=%18.15f ",
		       SCF,i,n,l,Trial_ep,Trial_D);
		printf("node=%2d Node_min=%2d Node_max=%2d MatP=%2d\n",node,Node_min,Node_max,j);
	      }

	      if (SCF<=3)      cri = 1.0e-6;
	      else if (SCF<=6) cri = 1.0e-9;
	      else             cri = 1.0e-12;

	      if (i!=1 && fabs(Deriv_O-Deriv_I)<=cri && node==(n-l-1)
		  ||
		  200<i 
		  ){

		po = 1;
		ep = Trial_ep; 
		MCP[so][n][l] = j;

		if (Check_switch==1){
		  printf("SCF=%2d  i=%2d  n l ep D %i %i %18.14f %18.14f\n",
			 SCF,i,n,l,Trial_ep,Trial_D);
		}
	      }

	      else{

		if ((n-l-1)<node){
		  p_D = Max_D;
		  p_ep = Max_ep;
		  p_node = Node_max; 
		  Max_D = Trial_D;
		  Max_ep = Trial_ep;
		  Node_max = node;
		}
		else if ((n-l-1)==node){

		  /****************************************************
		     if n-l-1 = an even number,
		     minus to plus as the energy decreases 
		     if n-l-1 = an odd number 
		     plus to minus as the energy decreases 
		  ****************************************************/

		  if ((n-l-1)%2==0){
		    if (0.0<=Trial_D){
		      p_D = Min_D;
		      p_ep = Min_ep;
		      p_node = Node_min; 
		      Min_D = Trial_D;
		      Min_ep = Trial_ep;
		      Node_min = node;
		    }
		    else{
		      p_D = Max_D;
		      p_ep = Max_ep;
		      p_node = Node_max; 
		      Max_D = Trial_D;
		      Max_ep = Trial_ep;
		      Node_max = node;
		    }
		  }
		  else{
		    if (0.0<=Trial_D){
		      p_D = Max_D;
		      p_ep = Max_ep;
		      p_node = Node_max; 
		      Max_D = Trial_D;
		      Max_ep = Trial_ep;
		      Node_max = node;
		    }
		    else{
		      p_D = Min_D;
		      p_ep = Min_ep;
		      p_node = Node_min; 
		      Min_D = Trial_D;
		      Min_ep = Trial_ep;
		      Node_min = node;
		    }
		  }
                  
		}
		else{
		  p_D = Min_D;
		  p_ep = Min_ep;
		  p_node = Node_min; 
		  Min_D = Trial_D;
		  Min_ep = Trial_ep;
		  Node_min = node;
		}

	      }

	    } while(po==0);

	    /****************************************************
                     eigenenergy and radial wave function
	    ****************************************************/

	    E_old[so][n][l] = E[state_num][so][n][l];
	    E[state_num][so][n][l] = ep;

	    for (i=0; i<=MCP[so][n][l]; i++){
	      PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*Lo[i];
	    }

	    for (i=(MCP[so][n][l]+1); i<Grid_Num; i++){
	      PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
	    }

	    /* In case of Dirac equation, compute the small wave function */

	    if (TwoComp_frag){

	      c = 137.0359895;
	      al = 1.0/c;
	      for (i=0; i<=MCP[so][n][l]; i++){
		r = MRV[i];
		tmp0 = 1.0/(al*(2.0/al/al + ep - V[i]));
		dG = pow(r,(double)l)*(((double)l+1.0)*Lo[i] + Mo[i]);
		FF[so][n][l][i] = tmp0*(dG + kappa*PF[so][n][l][i]/r);
	      }
	      for (i=(MCP[so][n][l]+1); i<Grid_Num; i++){
		r = MRV[i];
		tmp0 = 1.0/(al*(2.0/al/al + ep - V[i]));
		dG = ratio*pow(r,(double)l)*(((double)l+1.0)*Li[i] + Mi[i]);
		FF[so][n][l][i] = tmp0*(dG + kappa*PF[so][n][l][i]/r);
	      }
	    }

	    /****************************************************
		    Normalization by Simpson Integration
	    ****************************************************/

	    Simpson_Normalization(so, n, l);
	    if (Check_switch==1){
	      printf("SCF=%2d n=%2d l=%2d so=%2d eigenvalue=%15.12f (Hartree)\n",
		     SCF,n,l,so,E[state_num][so][n][l]);
	    }

	  } 
        }  /* l */
      }  /* n */
    } /* so */


    /****************************************************
                     Sum of eigenvalues
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
    HisEeigen[SCF] = Eeigen[state_num];

    /****************************************************
                      Check SCF criterion
    ****************************************************/

    if (SCF<SCF_MAX){

      if      (Mixing_switch==0){
        Density(state_num);
        Check_Sum_Num_Electrons();
        Simple_Mixing(state_num,rho[0],rho[1],rho[2],Rrho[0]);
      }
      else if (Mixing_switch==1){
        Density(state_num);
        Check_Sum_Num_Electrons();
        GR_Pulay(state_num,SCF);
      }
      else{
        printf("Mixing_switch=%i is not available\n",Mixing_switch);  
        exit(0);
      }

      /* Check the total number of electrons */

      if (Calc_Type==3 && state_num==1) te = total_electron1; /* for frozen core */
      else                              te = total_electron;

      if ( 1.0e-5<fabs(Check_Sum-te) ){
        printf("\n\nThis calculation is invalid.\n");
        printf("Calculated total number of electrons is %15.8f\n",
               Check_Sum);
        printf("Please use a larger value for grid.num\n");
        printf("Please check both grid.xmin and grid.xmax\n\n");
        exit(0);
      }

      /* Converge ? */

      HisNormRD[SCF] = NormRD[0];

      if (NormRD[0]<SCF_criterion) SCF_OK = 1;

      printf("<ALL>  SCF=%4d  Eeigen=%18.12f  NormRD=%18.12f\n",
              SCF,Eeigen[state_num],NormRD[0]);

      if (NormRD[0]<1.0e-6){
        MatP_fix_switch = 1;
      }

      if (Check_switch==1){
        printf("SCF=%2d  ",SCF);
        for(i=0; i<Grid_Num; i+=Grid_Num/10){
         printf("%10.6f ", rho[0][i]);
        }
        printf("\n");
      }

    }

    SCF++; 
  } while(SCF<=SCF_MAX && SCF_OK==0);

  SCF_END = SCF;
  Total_Energy(state_num);

	/* beginning of the addition */
	if(CalcType0==1){
		Calc_Type=0;
	}
	/* end of the addition */
} 


void Simpson_Normalization(int so, int n, int l)
{

  /****************************************************
           Normalization by Simpson Integration
  ****************************************************/

  static int i;
  static double s0,s1,s2,sumG,sumF,sum,mul0;
  static double tmp0;

  sumG = 0.0;
  sumF = 0.0;

  s0 = PF[so][n][l][0]*PF[so][n][l][0]*MRV[0]
     + PF[so][n][l][Grid_Num-1]*PF[so][n][l][Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + PF[so][n][l][i]*PF[so][n][l][i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + PF[so][n][l][i]*PF[so][n][l][i]*MRV[i];
  }
  sumG = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;

  /* In case of a full relativistic equation, the contribution
     of the minor wave function */

  if (TwoComp_frag){
    s0 = FF[so][n][l][0]*FF[so][n][l][0]*MRV[0]
       + FF[so][n][l][Grid_Num-1]*FF[so][n][l][Grid_Num-1]*MRV[Grid_Num-1];

    s1 = 0.0;
    for (i=1; i<=(Grid_Num-2); i=i+2){
      s1 = s1 + FF[so][n][l][i]*FF[so][n][l][i]*MRV[i];
    }

    s2 = 0.0; 
    for (i=2; i<=(Grid_Num-3); i=i+2){
      s2 = s2 + FF[so][n][l][i]*FF[so][n][l][i]*MRV[i];
    }
    sumF = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
  }

  /* normalize */
  /* printf("sumG=%15.12f sumF=%15.12f\n",sumG,sumF); */

  tmp0 = 1.0/sqrt(sumG + sumF);

  if (0.0<=PF[so][n][l][Grid_Num-200]){
    for (i=0; i<=(Grid_Num-1); i++){
      PF[so][n][l][i] = tmp0*PF[so][n][l][i];
    }
    if (TwoComp_frag){
      for (i=0; i<=(Grid_Num-1); i++){
        FF[so][n][l][i] = tmp0*FF[so][n][l][i];
      }
    }
  }
  else{
    for (i=0; i<=(Grid_Num-1); i++){
      PF[so][n][l][i] = -tmp0*PF[so][n][l][i];
    }
    if (TwoComp_frag){
      for (i=0; i<=(Grid_Num-1); i++){
        FF[so][n][l][i] = -tmp0*FF[so][n][l][i];
      }
    }
  }

}


double Quadratic_Interpolate(double x1, double x2, double x3,
                             double y1, double y2, double y3)
{

  static double D0,D1,D2,D3,a1,a2,a3;
  static double tmp0,tmp1,tmp2,tmp3,result;

  D0 = x1*x1*x2 + x2*x2*x3 + x3*x3*x1
    -x2*x3*x3 - x3*x1*x1 - x1*x2*x2;

  D1 = y1*x2 + y2*x3 + y3*x1
    -x2*y3 - x3*y1 - x1*y2;

  D2 = x1*x1*y2 + x2*x2*y3 + x3*x3*y1
    -y2*x3*x3 - y3*x1*x1 - y1*x2*x2;

  D3 = x1*x1*x2*y3 + x2*x2*x3*y1 + x3*x3*x1*y2
    -y1*x2*x3*x3 - y2*x3*x1*x1 - y3*x1*x2*x2;

  a1 = D1/D0; 
  a2 = D2/D0; 
  a3 = D3/D0; 

  tmp0 = a2*a2 - 4.0*a1*a3;
  if (0.0<=tmp0){
    tmp1 = sqrt(tmp0);
    tmp2 = (-a2 + tmp1)/(2.0*a1);
    tmp3 = (-a2 - tmp1)/(2.0*a1);

    if (x1<tmp2 && tmp2<x2) result = tmp2;
    else                    result = tmp3;
  }
  else{
    result = (x1*y2 - x2*y1)/(y2 - y1);
  }

  return result;
}


void Check_Sum_Num_Electrons()
{

  /****************************************************
       check the sum of the number of electrons
  ****************************************************/

  static int i;
  static double s0,s1,s2,sum;
  static double tmp0;

  sum = 0.0;

  s0 = rho[0][0]*MRV[0]*MRV[0]*MRV[0]
     + rho[0][Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + rho[0][i]*MRV[i]*MRV[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + rho[0][i]*MRV[i]*MRV[i]*MRV[i];
  }

  Check_Sum = 4.0*PI*(s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
}

