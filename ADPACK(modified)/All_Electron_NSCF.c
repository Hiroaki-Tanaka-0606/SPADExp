/**********************************************************************
  All_Electron_NSCF.c:

   All_Electron.c is a subroutine to calculate occupied and unoccupied 
   electronic states of atomic Kohn-Sham equation including all electrons.

  Log of All_Electron_NSCF.c:

     10/Dec/2002  All_Electron.c is released by T.Ozaki
     21/Sep/2021  All_Electron.c is copied to All_Electron_NSCF.c
     DD/MMM/2021  All_Electron_NSCF.c is modified by H.Tanaka

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

static int count_node(double L[ASIZE1], int Mat_p);

void All_Electron_NSCF(int state_num)
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

	/* new variable */
	static int ie;
	int search_node[ASIZE11];
	int search_node2[ASIZE11];
	int target_node,trial_node;
	int Mat_p;
	static double Trial_dE;
	static double log_diff_outward, log_diff_inward, int_outward, int_inward;
	static double interval_outward, interval_inward;
	static int Grid_Num1;
	static int trial_count;
	static double E_diff;
	
  MatP_fix_switch = 0;
      
  for (so=0; so<SOI_switch; so++){
		for(n=1; n<=max_N; n++){
			for(l=0; (l<n && l<=MaxL_PAO); l++){
				// for (n=1; n<=max_ocupied_N; n++){
				// for (l=0; l<n; l++){
        tmp0 = (double)AtomNum/(double)n;
        E[state_num][so][n][l] = -0.5*tmp0*tmp0 - 10.0;
      }
    }
  }

	Reduce_Num_grid=RNG[NVPS[0]][LVPS[0]];

	/****************************************************
             search eigenvalues for n and l  
	****************************************************/

	/* Matching point is set to the point where V[i]>0 */
	Mat_p=Grid_Num-1;
	for(i=Grid_Num-1; i>=0; i--){
		if(V[i]<0){
			Mat_p=i;
			break;
		}
	}
	printf("Matching point= %d x= %17.14f r= %17.14f\n", Mat_p, MXV[Mat_p], MRV[Mat_p]);
	Grid_Num1=(int)((double)Grid_Num*Reduce_Num_grid);
	printf("Grid_Num= %d RNG= %.1f Reduced Grid_Num= %d\n", Grid_Num, Reduce_Num_grid, Grid_Num1);
		
	/* for (so=0; so<SOI_switch; so++){*/
	for(so=0; so<1; so++){
		for(l=0; l<=MaxL_PAO; l++){
			printf("l= %d\n", l);

			/* rough search, similar to Multiple_PAO.c */
			for(ie=0; ie<Num_Partition; ie++){
				Trial_ep=search_LowerE+(search_UpperE-search_LowerE)*ie/(Num_Partition-1);
				Hamming_O(0,l,Trial_ep,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
			  search_node[ie]=node;
				search_node2[ie]=count_node(Lo, Mat_p);
				Hamming_I(0,l,Trial_ep,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);

				// printf("T_ep= %16.12f node= %d node2= %d\n", Trial_ep, search_node[ie], search_node2[ie]);
			}
			for(n=1; n<=max_N; n++){
				target_node=n-l-1;
				if(target_node<0) continue;
				printf("n= %d node= %d\n", n, target_node);
				int Trial_ep_found=0;
				for(ie=Num_Partition-1; ie>=0; ie--){
					if(search_node[ie]==target_node){
						Trial_ep=search_LowerE+(search_UpperE-search_LowerE)*ie/(Num_Partition-1);
						Trial_ep_found=1;
						break;
					}
				}
				if(Trial_ep_found==0){
					printf("can not find adequate Trial_ep\n");
					continue;
				}
				trial_count=0;
				do{
					/* use PF[so][n][l] as temporary space for wavefunction */

					/* outward */
					Hamming_O(0,l,Trial_ep,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
				  for (i=0; i<=Mat_p; i++){
						PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*Lo[i];
					}
					interval_outward=MRV[Mat_p]-MRV[Mat_p-1];
					log_diff_outward=(PF[so][n][l][Mat_p]-PF[so][n][l][Mat_p-1])/(interval_outward*PF[so][n][l][Mat_p]);
					int_outward=0.0;
					for(i=0; i<Mat_p; i++){
						int_outward+=pow(PF[so][n][l][i], 2)*(MRV[i+1]-MRV[i]);
					}
					int_outward/=pow(PF[so][n][l][Mat_p], 2);

					/* inward */
					Hamming_I(0,l,Trial_ep,kappa,Mi,Li,DMi,V,V,Reduce_Num_grid);
					scale=Lo[Mat_p]/Li[Mat_p];
				  for (i=Mat_p; i<Grid_Num1; i++){
						PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*Li[i]*scale;
						// printf("%8.5e ", PF[so][n][l][i]);
					}
					interval_inward=MRV[Mat_p+1]-MRV[Mat_p];
					log_diff_inward=(PF[so][n][l][Mat_p+1]-PF[so][n][l][Mat_p])/(interval_inward*PF[so][n][l][Mat_p]);
					int_inward=0.0;
					for(i=Mat_p; i<Grid_Num1-1; i++){
						int_inward+=pow(PF[so][n][l][i], 2)*(MRV[i+1]-MRV[i]);
						// printf("%18.12e ",int_inward);
					}
					int_inward/=pow(PF[so][n][l][Mat_p], 2);

					// printf("ldo= %18.12e io= %18.12e ldi= %18.12e ii= %18.12e\n", log_diff_outward, int_outward, log_diff_inward, int_inward);

					Trial_dE=(log_diff_outward-log_diff_inward)/(int_outward+int_inward);
					if(!isfinite(Trial_dE)){
						// printf("Warning: Trial_dE is not finite");
						Trial_dE=log_diff_outward-log_diff_inward;
					}

					Hamming_O(0,l,Trial_ep+Trial_dE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
					trial_node=count_node(Lo,Mat_p);
					while(trial_node!=target_node){
						// printf("Warning: E+dE= %12.5e has %d nodes, but the target is %d\n", Trial_ep+Trial_dE, trial_node, target_node);
						Trial_dE/=2.0;
						Hamming_O(0,l,Trial_ep+Trial_dE,kappa,Mo,Lo,DMo,V,V,Reduce_Num_grid);
						trial_node=count_node(Lo,Mat_p);
					}

					E_diff=fabs(Trial_dE/Trial_ep);
					
					trial_count++;
					// printf("i= %d Trial_ep= %17.14f Trial_dE= %17.14f dE/E= %17.14f\n", trial_count, Trial_ep, Trial_dE, E_diff);
					Trial_ep+=Trial_dE;
					
				}while(trial_count<100 || (trial_count<1000 && E_diff>1e-4));
				printf("i= %d Trial_ep= %17.14f Trial_dE= %17.14f dE/E= %17.14f\n", trial_count, Trial_ep, Trial_dE, E_diff);

				
				/****************************************************
		    Normalization by Simpson Integration
				****************************************************/
					
				Simpson_Normalization(so, n, l);

			}
				

	    /****************************************************
                     eigenenergy and radial wave function
	    ****************************************************/
					/*
	    E_old[so][n][l] = E[state_num][so][n][l];
	    E[state_num][so][n][l] = ep;

	    for (i=0; i<=MCP[so][n][l]; i++){
	      PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*Lo[i];
	    }

	    for (i=(MCP[so][n][l]+1); i<Grid_Num; i++){
	      PF[so][n][l][i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
	    }

	    /* In case of Dirac equation, compute the small wave function */
					/*
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
					/*
	    Simpson_Normalization(so, n, l);
	    if (Check_switch==1){
	      printf("SCF=%2d n=%2d l=%2d so=%2d eigenvalue=%15.12f (Hartree)\n",
		     SCF,n,l,so,E[state_num][so][n][l]);
				 }

				 } */
		}  /* l */
	} /* so */


    /****************************************************
                     Sum of eigenvalues
    ****************************************************/
	/*
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
	/*
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
	/*
		if (Calc_Type==3 && state_num==1) te = total_electron1; /* for frozen core *//*
      else                              te = total_electron;

      if ( 1.0e-5<fabs(Check_Sum-te) ){
        printf("\n\nThis calculation is invalid.\n");
        printf("Calculated total number of electrons is %15.8f\n",
               Check_Sum);
        printf("Please use a larger value for grid.num\n");
        printf("Please check both grid.xmin and grid.xmax\n\n");
        exit(0);
      }

      /* Converge ? *//*

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
  Total_Energy(state_num);*/
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


int count_node(double L[ASIZE1], int Mat_p){
	int i;
	int count=0;
	for(i=2; i<Mat_p; i++){
		if(L[i]*L[i-1]<0.0) {
			count++;
			// printf("%d ",i);
		}
	}
	// printf("\n");
	return count;
}
	
