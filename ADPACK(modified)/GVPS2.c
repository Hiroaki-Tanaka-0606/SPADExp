/**********************************************************************
  GVPS.c:

     GVPS.c is a subroutine to calculate norm conseving pseudo
     potentials using a generalized scheme developed by T.Ozaki.

  Log of GVPS.c:

     29/Apr/2006  Modified by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"



#define iorigin  10 





static double Norm_Core(double af[ASIZE1], double rc);




static double Int_phi0_phi1(double *phi0,double *phi1,int ic);

static void Do_PatchWork(int nloop,
                         int so, int L, int mul,
                         double slope,
                         int icut[20],
                         double ***EAE, double ****WFAE,
                         double *VAE_Tmp,
                         double ***EPS,  
                         double ****WFPS, double ****WFPS0,
                         double ****GVPS);

static double Dlog_Norm(int l,
                        double a[ASIZE6][ASIZE6],
                        double a0[ASIZE6][ASIZE6],
                        double b[ASIZE6], double c[ASIZE6], double c2,
                        double Ncore,
			double rc, double r2, double r3, double r4);

static double Norm_Core2(int l,double rc,
                         double c0,
                         double c2,
                         double c4,
                         double c6,
                         double c8,
                         double c10,
			 double c12);

static void Spline_Func(double R, double fv[3], double *af);


static double Norm_Core0(double *af, int ic);

static void Generate_GVPS(double ***EAE, double ****WFAE, double VAE[], 
                          double ***EPS, double ****WFPS, double ****WFPS0,
                          double ****GVPS, double Rcutoff[20], int icut[20]);

static void Generate_TMVPS(double ***EAE, double ****WFAE, double VAE[], 
                           double ***EPS, double ****WFPS, double ****WFPS0,
                           double ****GVPS, double Rcutoff[20], int icut[20]);

static void Calc_BoundStates(char *mode, double ***EV, double ****WF,
                             double VAE[], double ****GVPS);

static void Simpson(double *PAO);

static int Calc_VPS_PAO( int calc_type,
                         int so, int l, double rc, double ep,
                         double *AEW, double AEV[],
                         double *PW,  double *PV );


static void search_roughly(
        int Num_Part,  /* input, array size of search_?[] */
        double search_e[], int search_node[], double search_D[], int search_MatP[],/* output */
        double min, double max, 
	int so, int GL, int CTP_P, int Grid_Num, int Smallest_MP,double alpha, /* input */
        int asize1, double Mo[], double Lo[], double DMo[], /* work */
        double Mi[],double Li[], double DMi[],
        double Ltmp[], double Vpot[] );

static int search_bisection(
           int NumL, double Min_ep, double Max_ep, double Min_D, double Max_D,
           int Min_MatP, int Max_MatP,
           double Criterion, /* input */
           int retnode, double Trial_ep, double Trial_D, /* output */
	   int so, int GL, int CTP_P, int Grid_Num, int Smallest_MP,double alpha, /* input */
	   int Num_Part, double Mo[], double Lo[], double DMo[], /* work */
	   double Mi[],  double Li[], double DMi[],
	   double Ltmp[], double Vpot[],
           double *E, double *WF);







void GVPS(int state_num)
{ 
  static int m,n,l,i,j,ip,po,po1,numvps,ll,mul,mul1;
  static int m1,n1,loop_num,loop_max,num_nodes,so;
  int icut[20];
  static double max_SF,r,ef;
  static double ep,rc,a2,sum,dnorm,da,dnorm_p;
  double Rcutoff[20];
  char *s_vec[20];

  double ***EAE,****WFAE;
  double ***EPS,****WFPS0,****WFPS;
  double ****GVPS;

  /*****************************************************
                    allocation of arrays
  *****************************************************/

  EAE = (double***)malloc(sizeof(double**)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    EAE[so] = (double**)malloc(sizeof(double*)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      EAE[so][l] = (double*)malloc(sizeof(double)*GVPS_ProNum); 
    }
  }

  WFAE = (double****)malloc(sizeof(double***)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    WFAE[so] = (double***)malloc(sizeof(double**)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      WFAE[so][l] = (double**)malloc(sizeof(double*)*GVPS_ProNum); 
      for (mul=0; mul<GVPS_ProNum; mul++){
        WFAE[so][l][mul] = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
      }
    }  
  }

  EPS = (double***)malloc(sizeof(double**)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    EPS[so] = (double**)malloc(sizeof(double*)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      EPS[so][l] = (double*)malloc(sizeof(double)*GVPS_ProNum); 
    }
  }

  WFPS = (double****)malloc(sizeof(double***)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    WFPS[so] = (double***)malloc(sizeof(double**)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      WFPS[so][l] = (double**)malloc(sizeof(double*)*GVPS_ProNum); 
      for (mul=0; mul<GVPS_ProNum; mul++){
        WFPS[so][l][mul] = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
      }
    }  
  }

  WFPS0 = (double****)malloc(sizeof(double***)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    WFPS0[so] = (double***)malloc(sizeof(double**)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      WFPS0[so][l] = (double**)malloc(sizeof(double*)*GVPS_ProNum); 
      for (mul=0; mul<GVPS_ProNum; mul++){
        WFPS0[so][l][mul] = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
      }
    }  
  }

  GVPS = (double****)malloc(sizeof(double***)*SOI_switch); 
  for (so=0; so<SOI_switch; so++){
    GVPS[so] = (double***)malloc(sizeof(double**)*(GVPS_MaxL+1)); 
    for (l=0; l<=GVPS_MaxL; l++){
      GVPS[so][l] = (double**)malloc(sizeof(double*)*GVPS_ProNum); 
      for (mul=0; mul<GVPS_ProNum; mul++){
        GVPS[so][l][mul] = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
      }
    }  
  }

  /*****************************************************
         calculate all electron eigenfunctions
  *****************************************************/

  Calc_BoundStates("AE", EAE, WFAE, V_all_ele, GVPS);

  /*****************************************************
            generate each TM pseudopotential
         of the lowest L states to be pseudized
  *****************************************************/

  Generate_TMVPS(EAE, WFAE, V_all_ele, EPS, WFPS, WFPS0, GVPS, Rcutoff, icut);

  /*****************************************************
   calculate pseudo-eigenfunctions for the TM potential
  *****************************************************/

  Calc_BoundStates("PS", EPS, WFPS0, V_all_ele, GVPS);

  /*****************************************************
        generate a pseudopotential of each state
  *****************************************************/

  Generate_GVPS(EAE, WFAE, V_all_ele, EPS, WFPS, WFPS0, GVPS, Rcutoff, icut);


  
  for (i=0; i<Grid_Num; i++){
    l = 2; 
    printf("%15.12f %15.12f %15.12f %15.12f %15.12f %15.12f\n",
             MRV[i],V_all_ele[i],GVPS[0][l][0][i],GVPS[0][l][1][i],GVPS[0][l][2][i],GVPS[0][l][3][i]);

  }




  exit(0);





  /*
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f %15.12f %15.12f\n",
             MRV[i],WFAE[0][0][0][i],WFAE[0][0][1][i],WFAE[0][0][2][i],WFAE[0][0][3][i]);
  }

  exit(0); 
  */



  /*
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f %15.12f %15.12f\n",
             MRV[i],V_all_ele[i],GVPS[0][0][0][i],GVPS[0][1][0][i],GVPS[0][2][0][i]);

  }
  */


  /*
  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      for (mul=0; mul<GVPS_ProNum; mul++) {
        printf("so=%2d GL=%2d mul=%2d diff=%15.12f\n",so,l,mul,EPS[so][l][mul]-EAE[so][l][mul]);
      }
    }
  }
  */













  /**************************************
           freeing of arrays
  **************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      free(EAE[so][l]);
    }
    free(EAE[so]);
  }
  free(EAE);

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      for (mul=0; mul<GVPS_ProNum; mul++){
        free(WFAE[so][l][mul]);
      }
      free(WFAE[so][l]);
    }  
    free(WFAE[so]);
  }
  free(WFAE);

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      free(EPS[so][l]);
    }
    free(EPS[so]);
  }
  free(EPS);

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      for (mul=0; mul<GVPS_ProNum; mul++){
        free(WFPS[so][l][mul]);
      }
      free(WFPS[so][l]);
    }  
    free(WFPS[so]);
  }
  free(WFPS);

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      for (mul=0; mul<GVPS_ProNum; mul++){
        free(WFPS0[so][l][mul]);
      }
      free(WFPS0[so][l]);
    }  
    free(WFPS0[so]);
  }
  free(WFPS0);

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<=GVPS_MaxL; l++){
      for (mul=0; mul<GVPS_ProNum; mul++){
        free(GVPS[so][l][mul]);
      }
      free(GVPS[so][l]);
    }  
    free(GVPS[so]);
  }
  free(GVPS);

  /* for std output */
  printf("\n");

}







static void Generate_GVPS(double ***EAE, double ****WFAE, double VAE[], 
                          double ***EPS, double ****WFPS, double ****WFPS0,
                          double ****GVPS, double Rcutoff[20], int icut[20])
{
  int so,L,mul,m,po;
  int i,j,i0,i1,i2;
  int k,k1,k2,k3,k4,k5;
  int k6,k7,k8,k9,k10;
  int knum,knum0,num,ndiv;
  double r0,r1,r2,rcut,dr,fv[3];
  double r,r3,r4,r5,r6,r7,r8,r9;
  double *VAE_Tmp;

  int nloop,mul0,mul1;
  double slope_min,slope_max,slope,dslope,rc; 
  double NormAE,dNormPS,dNormPS_min,dNormPS_max,NormPS;
  double dNormPS0,dNormPS1,opt_slope,opt_dNorm;

  /* allocation of array */

  VAE_Tmp = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
  for (i=0; i<Grid_Num; i++) VAE_Tmp[i] = VAE[i];

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){
      for (mul=1; mul<GVPS_ProNum; mul++){

	/* find the upper and lower bound of slope */

	i = icut[L];
	rc = MRV[i];

	NormAE = Norm_Core0(WFAE[so][L][mul],icut[L]);

        if (mul==1)
  	  slope_min = 0.7;
        else 
  	  slope_min = 0.1;

	slope_max = 3.0;
        ndiv = 700;

        dslope = (slope_max - slope_min)/(double)(ndiv-1); 

        po = 0;
        nloop = 0;
        dNormPS0 = 0.0;
        opt_dNorm = 10000;

        do {

          slope = slope_min + (double)nloop*dslope;

	  Do_PatchWork( nloop, so, L, mul, slope, icut, EAE, WFAE, VAE_Tmp, EPS, WFPS, WFPS0, GVPS );
	  dNormPS1 = Norm_Core0(WFPS[so][L][mul],icut[L]) - NormAE;

          if (fabs(dNormPS1)<fabs(opt_dNorm)){
            opt_dNorm = dNormPS1;
            opt_slope = slope;
          } 

	  /*
          if (so==0 && L==0 && mul==2){
	    for (i=0; i<Grid_Num; i++){
	      printf("%15.12f %15.12f %15.12f %15.12f\n",
		     MRV[i],WFAE[so][L][mul][i],WFPS0[so][L][mul][i],WFPS[so][L][mul][i]);
	    }
            exit(0);
          }
	  */

	  /*
          printf("search the boundary: so=%2d L=%2d mul=%2d NormAE=%15.12f slope=%15.12f dNormPS1=%15.12f\n",
		 so,L,mul,NormAE,slope,dNormPS1);
	  */

	  if (dNormPS1*dNormPS0<0.0){
            po = 1;
            slope_min = slope - dslope;
            slope_max = slope;

	    dNormPS_min = dNormPS0;
	    dNormPS_max = dNormPS1;
	  }

          dNormPS0 = dNormPS1;

          nloop++;

        } while (po==0 && nloop<ndiv);

	/* refine an optimum slope by a bisection method */

	if (po==1){

	  printf("find region so=%2d L=%2d mul=%2d\n",so,L,mul);

	  po = 0;
	  nloop = 0;

	  do {

	    slope = 0.5*(slope_min + slope_max);
	    Do_PatchWork( nloop, so, L, mul, slope, icut, EAE, WFAE, VAE_Tmp, EPS, WFPS, WFPS0, GVPS );

	    dNormPS = Norm_Core0(WFPS[so][L][mul],icut[L]) - NormAE;

	    printf("L=%2d mul=%2d nloop=%4d slope=%15.12f dNormPS=%18.15f\n",
		   L,mul,nloop,slope,dNormPS);

	    if (0.0<dNormPS*dNormPS_min){
	      slope_min = slope;
	      dNormPS_min = dNormPS;
	    }
	    else {
	      slope_max = slope;
	      dNormPS_max = dNormPS;
	    }

	    if (fabs(dNormPS)<1.0e-13){
	      po = 1;
	    }

	    nloop++;

	  } while (po==0 && nloop<200);

	}

        else{
          slope = opt_slope;
	  Do_PatchWork( nloop, so, L, mul, slope, icut, EAE, WFAE, VAE_Tmp, EPS, WFPS, WFPS0, GVPS );
	  dNormPS = Norm_Core0(WFPS[so][L][mul],icut[L]) - NormAE;
          printf("\nso=%2d L=%2d mul=%2d slope=%15.12f diff of Norm=%15.12f\n\n",so,L,mul,slope,dNormPS);
        }

      }
    }
  }

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){

      for (mul0=0; mul0<GVPS_ProNum; mul0++){
        for (mul1=0; mul1<GVPS_ProNum; mul1++){

          NormAE = Int_phi0_phi1(WFAE[so][L][mul0],WFAE[so][L][mul1],icut[L]);
          NormPS = Int_phi0_phi1(WFPS[so][L][mul0],WFPS[so][L][mul1],icut[L]);
          dNormPS = NormPS - NormAE;

          printf("so=%2d L=%2d mul0=%2d mul1=%2d NormAE=%15.12f NormPS=%15.12f dNormPS=%15.12f\n",
                  so,L,mul0,mul1,NormAE,NormPS,dNormPS);
        }
      }
    }
  }

  /* freeing of array */

  free(VAE_Tmp);

}













static void Do_PatchWork(int nloop,
                         int so, int L, int mul,
                         double slope,
                         int icut[20],
                         double ***EAE, double ****WFAE,
                         double *VAE_Tmp,
                         double ***EPS,  
                         double ****WFPS, double ****WFPS0,
                         double ****GVPS)
{
  int i,j,k,i0,i1,k2,po,n;
  double u0,u1,v0,v1,v2,sum,sum0,sum1,q0,q2;
  double a[ASIZE6][ASIZE6];
  double b[ASIZE6];
  double Q_Coes[30];
  double fv[3];
  double r,r0,r1,r2,r3,r4,r5;
  double r6,r7,r8,r9,dr0,dr1;
  double factor,center,rl,rr,width2;
  static double IntAE[20];
  static double IntCore[20];
  static double IntPoly[20][30];

  /*********************************************************
           calculate parameters for the generalized
                norm concerving condition
  *********************************************************/

  if (nloop==0){

    for (k=0; k<mul; k++){ 

      /* calculate Int_0^{MRV[icut[L]]} dr (WFAE_mul*WFAE_k) */

      sum0 = 0.0; 
      for (i=1; i<icut[L]; i++){
        sum0 += WFAE[so][L][mul][i]*WFAE[so][L][k][i]*MRV[i]; 
      }

      sum1 = WFAE[so][L][mul][0]*WFAE[so][L][k][0]*MRV[0]
           + WFAE[so][L][mul][icut[L]]*WFAE[so][L][k][icut[L]]*MRV[icut[L]];

      IntAE[k] = 0.5*dx*(sum1 + 2.0*sum0);

      /* calculate Int_0^{MRV[iorigin]} dr (WFPS0_mul*WFPS_k) */

      sum0 = 0.0; 
      for (i=1; i<iorigin; i++){
        sum0 += WFPS0[so][L][mul][i]*WFPS[so][L][k][i]*MRV[i]; 
      }

      sum1 = WFPS0[so][L][mul][0]*WFPS[so][L][k][0]*MRV[0]
           + WFPS0[so][L][mul][iorigin]*WFPS[so][L][k][iorigin]*MRV[iorigin];

      IntCore[k] = 0.5*dx*(sum1 + 2.0*sum0);

      /* calculate Int_{MRV[iorigin]}^{MRV[icut[L]]} dr (r^n*WFPS_k) */

      for (n=0; n<=(9+mul); n++){

        sum0 = 0.0; 
        for (i=(iorigin+1); i<icut[L]; i++){
          r = MRV[i]; 
          sum0 += pow(r,(double)n)*WFPS[so][L][k][i]*MRV[i]; 
        }

        sum1 = pow(MRV[iorigin],(double)n)*WFPS[so][L][k][iorigin]*MRV[iorigin]
             + pow(MRV[icut[L]],(double)n)*WFPS[so][L][k][icut[L]]*MRV[icut[L]];

        IntPoly[k][n] = 0.5*dx*(sum1 + 2.0*sum0);

      } /* n */
    } /* k */
  } /* if (nloop==0) */

  /*********************************************************
                construct a linear equation
  *********************************************************/

  r0 = MRV[iorigin]; 
  r1 = MRV[icut[L]];

  /* for r0 */

  r = r0;

  /* zeroth derivative */
  for (k=0; k<=(9+mul); k++){
    a[0][k] = pow(r,(double)k);
  }

  /* first derivative */
  for (k=0; k<=(9+mul); k++){
    a[1][k] = (double)k*pow(r,(double)(k-1.0));
  }

  /* second derivative */
  for (k=0; k<=(9+mul); k++){
    a[2][k] = (double)k*((double)k-1.0)*pow(r,(double)(k-2.0));
  }

  /* third derivative */
  for (k=0; k<=(9+mul); k++){
    a[3][k] = (double)k*((double)k-1.0)*((double)k-2.0)*pow(r,(double)(k-3.0));
  }

  /* forth derivative */
  for (k=0; k<=(9+mul); k++){
    a[4][k] = (double)k*((double)k-1.0)*((double)k-2.0)*((double)k-3.0)*pow(r,(double)(k-4.0));
  }

  Spline_Func(r,fv,WFPS0[so][L][mul]);
  u0 = slope*fv[0];
  u1 = slope*fv[1];

  Spline_Func(r,fv,GVPS[so][L][0]);
  v0 = fv[0];
  v1 = fv[1];
  v2 = fv[2];

  a[0][10+mul] = u0;
  a[1][10+mul] = u1;
  a[2][10+mul] = u0*(2.0*v0 - 2.0*EPS[so][L][mul] + (double)L*((double)L+1.0)/r/r);
  a[3][10+mul] = u1*(2.0*v0 - 2.0*EPS[so][L][mul] + (double)L*((double)L+1.0)/r/r)
               + u0*(2.0*v1 - 2.0*(double)L*((double)L+1.0)/r/r/r);
  a[4][10+mul] = a[2][10+mul]*(2.0*v0 - 2.0*EPS[so][L][mul] + (double)L*((double)L+1.0)/r/r)
               + 4.0*u1*(v1 - (double)L*((double)L+1.0)/r/r/r)
               + 2.0*u0*(v2 + 3.0*(double)L*((double)L+1.0)/r/r/r/r);

  /* for r1 */

  r = r1;

  /* zeroth derivative */
  for (k=0; k<=(9+mul); k++){
    a[5][k] = pow(r,(double)k);
  }

  /* first derivative */
  for (k=0; k<=(9+mul); k++){
    a[6][k] = (double)k*pow(r,(double)(k-1.0));
  }

  /* second derivative */
  for (k=0; k<=(9+mul); k++){
    a[7][k] = (double)k*((double)k-1.0)*pow(r,(double)(k-2.0));
  }

  /* third derivative */
  for (k=0; k<=(9+mul); k++){
    a[8][k] = (double)k*((double)k-1.0)*((double)k-2.0)*pow(r,(double)(k-3.0));
  }

  /* forth derivative */
  for (k=0; k<=(9+mul); k++){
    a[9][k] = (double)k*((double)k-1.0)*((double)k-2.0)*((double)k-3.0)*pow(r,(double)(k-4.0));
  }

  Spline_Func(r,fv,WFAE[so][L][mul]);
  u0 = fv[0];
  u1 = fv[1];

  Spline_Func(r,fv,VAE_Tmp);
  v0 = fv[0];
  v1 = fv[1];
  v2 = fv[2];

  a[5][10+mul] = u0;
  a[6][10+mul] = u1;
  a[7][10+mul] = u0*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r/r);
  a[8][10+mul] = u1*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r/r)
               + u0*(2.0*v1 - 2.0*(double)L*((double)L+1.0)/r/r/r);
  a[9][10+mul] = a[7][10+mul]*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r/r)
               + 4.0*u1*(v1 - (double)L*((double)L+1.0)/r/r/r)
               + 2.0*u0*(v2 + 3.0*(double)L*((double)L+1.0)/r/r/r/r);

  /* generalized norm conserving */

  for (i=0; i<mul; i++){ 
    for (k=0; k<=(9+mul); k++){
      a[9+i+1][k] = IntPoly[i][k];
    }

    a[9+i+1][10+mul] = IntAE[i] - slope*IntCore[i];
  }

  /* solve the linear equation */
          
  Gauss_LEQ(9+mul,a,b);

  for (i=0; i<=(9+mul); i++) Q_Coes[i] = b[i];

  /*********************************************************
            construct WFPS by a patchwork method                        
  *********************************************************/
       
  for (i=0; i<iorigin; i++){
    WFPS[so][L][mul][i] = slope*WFPS0[so][L][mul][i];
  }

  for (i=iorigin; i<icut[L]; i++){
    r = MRV[i];
    sum = 0.0;
    for (j=0; j<=(9+mul); j++) sum += Q_Coes[j]*pow(r,(double)j);
    WFPS[so][L][mul][i] = sum;
  }

  for (i=icut[L]; i<Grid_Num; i++){
    WFPS[so][L][mul][i] = WFAE[so][L][mul][i];
  }

  /*********************************************************
              construct GVPS by a patchwork method                         
  *********************************************************/

  /* visinicity of core */ 

  for (i=0; i<iorigin; i++){
    GVPS[so][L][mul][i] = GVPS[so][L][0][i] + (EAE[so][L][mul] - EPS[so][L][mul]);
  }

  /* intermediate region */ 

  for (i=iorigin; i<icut[L]; i++){

    r = MRV[i];

    q0 = 0.0;
    q2 = 0.0;

    for (j=0; j<=(9+mul); j++){
      q0 += Q_Coes[j]*pow(r,(double)j);
      q2 += (double)j*((double)j-1.0)*Q_Coes[j]*pow(r,(double)(j-2));
    }

    GVPS[so][L][mul][i] = EAE[so][L][mul] - 0.5*(double)L*((double)L+1.0)/r/r + 0.5*q2/q0;
  }  

  /* outside of cutoff */ 
  for (i=icut[L]; i<Grid_Num; i++){
    GVPS[so][L][mul][i] = GVPS[so][L][0][i];
  }

}








static void Generate_TMVPS(double ***EAE, double ****WFAE, double VAE[], 
                           double ***EPS, double ****WFPS, double ****WFPS0,
                           double ****GVPS, double Rcutoff[20], int icut[20])
{
  int so,L,m,po,i,i1,node,node0,imin,mul;
  double ep,rc,RMOP[20];

  /* find radius of the most outer peak of the ground state */

  for (L=0; L<=GVPS_MaxL; L++){

    i = Grid_Num - 20;
    po = 0;

    do {

      if (fabs(WFAE[0][L][0][i-1])<fabs(WFAE[0][L][0][i])){
	po = 1;
        i1 = i;
      }

      i--;
    } while (po==0 && 1<i);

    if (po==0){
      printf("could not find radius of most outer peak L=%2d\n",L);
      exit(0);
    }

    RMOP[L] = MRV[i1];
  }

  /* adjust the cutoff radius */
  
  for (L=0; L<=GVPS_MaxL; L++){

    /* find the number of nodes of the lowest ground state */

    node = 0;
    for (i=1; i<=(Grid_Num-2); i++){
      if (WFAE[0][L][0][i]*WFAE[0][L][0][i+1]<0.0) node++;
    }
    node0 = node;

    imin = Grid_Num - 1;

    for (so=0; so<SOI_switch; so++){
      for (mul=1; mul<GVPS_ProNum; mul++){

        i = 1;
        node = 0;
        po = 0;
        do {
          if (WFAE[so][L][mul][i]*WFAE[so][L][mul][i+1]<0.0) node++;

          if (node==(node0+1)){
            i1 = i; 
            po = 1;
          }

          i++;
        } while (po==0 && i<(Grid_Num-1));

        if (po==0){
          printf("could not adjust the cutoff radius at so=%2d L=%2d\n",so,L);
          exit(0);
        }

        if (i1<imin) imin = i1;    

      } /* mul */ 
    } /* so */

    icut[L] = i1 - 5;
    Rcutoff[L] = MRV[icut[L]];

    printf("L=%2d Rcutoff=%15.12f RMOP=%15.12f\n",L,Rcutoff[L],RMOP[L]);

  } /* L */

  /* generate TMVPS */

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){

      ep = EAE[so][L][0];
      m = GI_VPS[L][0]; 
      rc = Rcutoff[L];

      po = Calc_VPS_PAO(0,so,L,rc,ep,WFAE[so][L][0],VAE,WFPS0[so][L][0],GVPS[so][L][0]);

      for (i=0; i<Grid_Num; i++){
        WFPS[so][L][0][i] = WFPS0[so][L][0][i];     
      }  
    }
  }
}






static int Calc_VPS_PAO( int calc_type,
                         int so, int l, double rc, double ep,
                         double *AEW, double AEV[],
                         double *PW,  double *PV )
{
  /******************************************************************
     calc_type
       0, usual
       1, search the upper bound of a scaling factor for norm
       2, making partial wave functions for linear combination
  ******************************************************************/

  static int i,j,po,calc_fail;
  static int loop_max,loop_num;
  static double r,dc,fv[3];
  static double p,p0,p1,v0,v1,v2,Ncore,Norm2;
  static double r2,r3,r4,r5,r6,r7,r8,r9,r10,r11,r12;
  static double c2max,c2min,c0,c2,c4,c6,c8,c10,c12;
  static double DLN,DLN_max,DLN_min,DLN_p,pf1,pf2;
  static double a[ASIZE6][ASIZE6],b[ASIZE6],c[ASIZE6];
  static double a0[ASIZE6][ASIZE6];

  calc_fail = 0;
  loop_max = 10000;

  /****************************************************
        input data from all electron calculation
  ****************************************************/

  Spline_Func(rc,fv,AEW);
  p0 = fv[0];
  p1 = fv[1];
  Spline_Func(rc,fv,AEV);
  v0 = fv[0];
  v1 = fv[1];
  v2 = fv[2];
  Ncore = Norm_Core(AEW,rc);

  b[0] = log(p0/pow(rc,l+1.0));
  b[1] = p1/p0 - (l + 1.0)/rc;
  b[2] = 2.0*v0 - 2.0*ep - 2.0*b[1]*(l + 1.0)/rc - b[1]*b[1];
  b[3] = 2.0*v1 + 2.0*b[1]*(l + 1.0)/rc/rc
        -2.0*b[2]*(l + 1.0)/rc - 2.0*b[1]*b[2];
  b[4] = 2.0*v2 - 4.0*b[1]*(l + 1.0)/rc/rc/rc
         + 4.0*b[2]*(l + 1.0)/rc/rc
         - 2.0*b[3]*(l + 1.0)/rc - 2.0*b[2]*b[2] - 2.0*b[1]*b[3];

  r2  = rc*rc;
  r3  = r2*rc;
  r4  = r2*r2;
  r5  = r4*rc;
  r6  = r4*r2;
  r7  = r6*rc;
  r8  = r4*r4;
  r9  = r8*rc;
  r10 = r8*r2;
  r11 = r10*rc;
  r12 = r10*r2;

  a[0][0] = 1.0;
  a[0][1] = r6;
  a[0][2] = r8; 
  a[0][3] = r10; 
  a[0][4] = r12; 

  a[1][0] = 0.0;
  a[1][1] = 6.0*r5;
  a[1][2] = 8.0*r7; 
  a[1][3] = 10.0*r9; 
  a[1][4] = 12.0*r11;

  a[2][0] = 0.0;
  a[2][1] = 30.0*r4;
  a[2][2] = 56.0*r6; 
  a[2][3] = 90.0*r8; 
  a[2][4] = 132.0*r10;

  a[3][0] = 0.0;
  a[3][1] = 120.0*r3;
  a[3][2] = 336.0*r5; 
  a[3][3] = 720.0*r7; 
  a[3][4] = 1320.0*r9;

  a[4][0] = 0.0;
  a[4][1] = 360.0*r2;
  a[4][2] = 1680.0*r4; 
  a[4][3] = 5040.0*r6; 
  a[4][4] = 11880.0*r8;

  for (i=0; i<=4; i++){
    for (j=0; j<=4; j++){
      a0[i][j] = a[i][j];
    }
  }

  /********************************************************
   find the upper and lower points of the bisection method
  ********************************************************/

  po = 0;
  c2 = -20.0;
  dc = 0.1;
  loop_num = 0;

  DLN_p = Dlog_Norm(l,a,a0,b,c,c2,Ncore,rc,r2,r3,r4);

  do {
    c2 = c2 + dc;
    DLN = Dlog_Norm(l,a,a0,b,c,c2,Ncore,rc,r2,r3,r4);
    if (DLN*DLN_p<0.0){
      po = 1;
      c2max = c2;
      c2min = c2 - dc;
      DLN_max = DLN; 
      DLN_min = DLN_p;
    }
    else{
      DLN_p = DLN;
    }

    /*
    printf("loop_num=%4d DLN=%15.12f\n",loop_num,DLN);
    */

    loop_num++;      

    if (loop_max<loop_num){

      if (calc_type==0 && calc_type==2){  
        printf("********************************************************\n");
        printf("Could not find an appropriate value of the parameter c2\n");
        printf("in the pseudo potential generation by the TM scheme.\n");
        printf("Please check cutoff radii of both pseudo potentials and\n");
        printf("pseudo atomic oribtals\n");
        printf("********************************************************\n");
        exit(0); 
      }
      else if (calc_type==1){
        calc_fail = 1;
        return calc_fail;
      }  
    }

  } while (po==0);

  /****************************************************
                       Bisection
  ****************************************************/

  po = 0;
  loop_num = 0;
  do {
    c2 = 0.50*(c2min + c2max);
    DLN = Dlog_Norm(l,a,a0,b,c,c2,Ncore,rc,r2,r3,r4);
    if (sgn(DLN)==sgn(DLN_min)){
      DLN_min = DLN;
      c2min = c2;
    }
    else{
      DLN_max = DLN;
      c2max = c2;
    }
    if (fabs(DLN)<=1.0e-12) po = 1;

    /*
    printf("loop_num=%4d DLN=%15.12f\n",loop_num,DLN);
    */

    loop_num++;      

    if (loop_max<loop_num){
      printf("********************************************************\n");
      printf("Failure in searching of an appropriate c2\n");
      printf("in the pseudo potential generation by the TM scheme.\n");
      printf("Please check cutoff radii of both pseudo potentials and\n");
      printf("pseudo atomic oribtals\n");
      printf("********************************************************\n");
      exit(0); 
    }

  } while (po==0);

  c0  = c[0];
  c4  = c[5]; 
  c6  = c[1];
  c8  = c[2];
  c10 = c[3];
  c12 = c[4];

  /****************************************************
                       W2 and VPS
  ****************************************************/

  for (i=0; i<Grid_Num; i++){
    r = MRV[i]; 
    if (r<=rc){
      r2  = r*r;
      r3  = r2*r; 
      r4  = r2*r2;
      r5  = r4*r;
      r6  = r4*r2;
      r7  = r6*r;
      r8  = r4*r4;
      r9  = r8*r;
      r10 = r8*r2;
      r11 = r10*r;
      r12 = r10*r2;

      p = c0 + c2*r2 + c4*r4 + c6*r6 + c8*r8 + c10*r10 + c12*r12;
      pf1 = 2.0*c2*r + 4.0*c4*r3 + 6.0*c6*r5 + 8.0*c8*r7
	+ 10.0*c10*r9 + 12.0*c12*r11;
      pf2 = 2.0*c2 + 12.0*c4*r2 + 30.0*c6*r4 + 56.0*c8*r6
	+ 90.0*c10*r8 + 132.0*c12*r10;

      PW[i] = pow(r,(double)l+1.0)*exp(p);

      if (calc_type==0){
	PV[i] = ep + ((double)l + 1.0)*pf1/r + 0.5*(pf2 + pf1*pf1);
      }
      else if (calc_type==2){
	PV[i] =  ((double)l*((double)l+1.0)*pow(r,(double)l-1.0)
                + 2.0*((double)l+1.0)*pow(r,(double)l)*pf1
                + pow(r,(double)l+1.0)*(pf2 + pf1*pf1))*exp(p);
      }
    }
    else{
      PW[i] = AEW[i];
      PV[i] = AEV[i];
    }
  }

  return calc_fail;
}







static void Calc_BoundStates(char *mode, double ***EV, double ****WF,
                             double VAE[], double ****GVPS)
{
  int state_num=0;
  static int i0,i1,i2,imin,p_node,NumL,s,mu,imu;
  static int i,j,k,m,n,l,po,SCF,po_node,po1,po2,CTP_P;
  static int cNode,pNode,Node_min,Node_max,refine_loop;
  static int Enum,so,GL,M,q,MatP,MatP0,MatP1;
  static int nf,fg,Mul,num_k,Smallest_MP;

  static double Sr,Dr,norm_kmax,norm_kmin,norm_k,rmin,rmax,r;
  static double minMapD,tmp0,tmp1,tmp2,cDD,pDD,pDD2,pOEE,Fugou;
  static double DD_max,DD_min,di,alpha;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double Ltmp[ASIZE1],Vpot[ASIZE1];
  static double ratio,Deriv_O,Deriv_I,ep,sum;
  static double kei,DD_old,DD_new,s0,s1,s2,OEE,MaxE,MinE;
  static double Min_ep,Max_ep,Min_D,Max_D,Criterion;
  static double Trial_ep,Trial_D,pTrial_ep,pTrial_D;
  static double p_ep,p_D,cri,scale,sum0,d,dif;
  static double StepE,E_old[ASIZE3][ASIZE2];

  double search_min0,search_max0;
  double search_min, search_max; 
  double search_e[ASIZE11],search_D[ASIZE11];
  int    search_node[ASIZE11],search_invD[ASIZE11];
  double search_e2[ASIZE11],search_D2[ASIZE11];
  int    search_node2[ASIZE11],search_invD2[ASIZE11];
  int    search_MatP[ASIZE11], search_MatP2[ASIZE11];

  int    inode,inode0,ie,ie_min,ie_max,found,really_found,ith,count;
  int    MatP_min,MatP_max;

  search_min0 = search_LowerE;
  search_max0 = search_UpperE;

  /****************************************************
      solve the atomic Kohn-Sham equation with 
      a confinement potential by the Hamming method
  ****************************************************/

  Smallest_MP = Grid_Num/3;
  Criterion = 10e-10;

  CTP_P = 2;
  n = NVPS[0];  
  l = LVPS[0];

  MinE = E[state_num][0][n][l] - 0.50; 
  alpha = 0.1;

  for (so=0; so<SOI_switch; so++){
    for (GL=0; GL<=GVPS_MaxL; GL++){

      printf("<PAO>  Calculating multiple pseudo atomic orbitals for GL=%2d....\n",GL);

      for (inode=0; inode<GVPS_ProNum; inode++) {

        if (strcasecmp(mode,"AE")==0){
          m = GI_VPS[GL][0];
          n = NVPS[m];
          inode0 = inode + n - GL - 1;
          for (i=0; i<Grid_Num; i++) Vpot[i] = VAE[i];
        }
        else if (strcasecmp(mode,"PS")==0){  
          inode0 = inode; 
          for (i=0; i<Grid_Num; i++) Vpot[i] = GVPS[so][GL][0][i];
	}

	search_min = MinE;
	search_max = search_max0;

	search_roughly( Num_Partition, search_e, search_node, search_D, search_MatP,
		        search_min, search_max, 
		        so, GL,CTP_P, Grid_Num, Smallest_MP, alpha,
		        ASIZE1, Mo, Lo, DMo,
		        Mi, Li, DMi, Ltmp, Vpot);

	if (Check_switch==1) {
	  printf("inode=%d\n",inode);
	}

	for (count=0; ; count++) {

	  if (count==0) {
	    for (i=0; i<Num_Partition; i++) {
	      search_e2[i]    = search_e[i];
	      search_node2[i] = search_node[i];
	      search_D2[i]    = search_D[i];
	      search_MatP2[i] = search_MatP[i];
	    }
	  }
	  else {
	    search_roughly(
			   Num_Partition, search_e2, search_node2, search_D2, search_MatP2,
			   search_min,search_max, 
			   so, GL, CTP_P, Grid_Num, Smallest_MP, alpha,
			   ASIZE1, Mo, Lo, DMo,
			   Mi, Li, DMi, Ltmp, Vpot);
	  }

	  ie_min=-100;
	  ie_max=-100;

	  for (ie=0; ie<Num_Partition; ie++) {
	    if (inode0==search_node2[ie]) {
	      ie_min = ie;
	      break;
	    }
	  }

	  for (ie=0; ie<Num_Partition; ie++){
	    if (inode0==search_node2[ie]) {
	      ie_max = ie;
	    }
	  }

	  if (Check_switch==1){
	    printf("node=%d, ie=%d %d\n",inode0,ie_min,ie_max);
	  }

	  if (ie_min==-100 &&
	      ! (search_node2[0]<inode0 && search_node2[Num_Partition-1]>inode0)) {
	    search_min = search_min - (search_max-search_min)*0.5;
	    search_max = search_max - (search_max-search_min)*0.5;
	    continue;
	  }

	  /* argorighm 1, search inversion of sign in search_D */

	  found = 0;

	  for (ie=ie_min; ie<ie_max; ie++){
	    if (search_D2[ie]*search_D2[ie+1]<0.0 && search_node2[ie]==search_node2[ie+1] &&
		search_node2[ie]==inode0 ) {
	      search_invD[found++] = ie;
	    }
	  }

	  if (found) {
	    for (ith=0; ith<found; ith++) {

	      ie=search_invD[ith]; 

	      if (Check_switch==1) {
		printf("inv found ie=%d\n",ie);
	      }

	      really_found=search_bisection( inode0,
                                             search_e2[ie], search_e2[ie+1],
                                             search_D2[ie], search_D2[ie+1],
					     search_MatP2[ie], search_MatP2[ie+1],
					     Criterion,
					     node,Trial_ep, Trial_D,
					     so, GL, CTP_P, Grid_Num, Smallest_MP, alpha,
					     ASIZE1, Mo, Lo, DMo,
					     Mi, Li, DMi, Ltmp, Vpot, 
                                             &EV[so][GL][inode],WF[so][GL][inode]);

	      if (really_found) {
		goto exit_count;
	      }
	    }
	  }

	  /* argorighm 2, search the boundaries */
	  /* upper boundary */
	  if (ie_max<Num_Partition) {
	    if (Check_switch==1) {
	      printf("upper boundary, min, max=%f %f\n",search_e2[ie_max],search_e2[ie_max+1]);
	    }
	    search_min=search_e2[ie_max];
	    search_max=search_e2[ie_max+1];
	    MatP_min=search_MatP[ie_max];
	    MatP_max=search_MatP[ie_max+1];

	    if (search_max-search_min<1.0e-14) {
              printf("can not find wavefunction, GL=%d, node=%d\n",GL,inode);
              printf("This problem might be solved using\n");
              printf("  a SMALLER search.LowerE,\n");
              printf("  a LARGER search.UpperE,\n");
              printf("  a LARGER num.of.partition.\n");
              exit(10);
	    }
	    continue;
	  }

	  /* argorighm 2, search the boundaries */
	  /* lower boundary */
	  if (ie_min>0) {
	    if (Check_switch==0) {
	      printf("lower boundary, min, max=%f %f\n",search_e2[ie_min-1],search_e2[ie_min]);
	    }
	    search_min=search_e2[ie_min-1];
	    search_max=search_e2[ie_min];
	    MatP_min=search_MatP[ie_min-1];
	    MatP_max=search_MatP[ie_min];

	    continue;
	  }

	  /* argorighm 3, search minitely ... */
	  printf("can not find %d\n",inode);
	  break;
      
	} /* count */
      exit_count:  ;
     
      } /* inode */
    } /* GL */
  } /* so */

  printf("\n");

} 




void Simpson(double *PAO)
{
  /****************************************************
                   JRCAT note 97p
  ****************************************************/

  static int i;
  static double s0,s1,s2,sum,sqr_sum;

  sum = 0.0;
  s0 = PAO[0]*PAO[0]*MRV[0]
     + PAO[Grid_Num-1]*PAO[Grid_Num-1]*MRV[Grid_Num-1];

  s1 = 0.0;
  for (i=1; i<=(Grid_Num-2); i=i+2){
    s1 = s1 + PAO[i]*PAO[i]*MRV[i];
  }

  s2 = 0.0; 
  for (i=2; i<=(Grid_Num-3); i=i+2){
    s2 = s2 + PAO[i]*PAO[i]*MRV[i];
  }
  sum = (s0 + 4.0*s1 + 2.0*s2)*dx/3.0;
  sqr_sum = sqrt(sum);

  if (0.0<=PAO[Grid_Num-100]){
    for (i=0; i<=(Grid_Num-1); i++){
      PAO[i] = PAO[i]/sqr_sum;
    }
  }
  else {
    for (i=0; i<=(Grid_Num-1); i++){
      PAO[i] = -PAO[i]/sqr_sum;
    }
  }
}








static void search_roughly(
        int Num_Part,  /* input, array size of search_?[] */
        double search_e[], int search_node[], double search_D[], int search_MatP[],/* output */
        double min, double max,
	int so, int GL, int CTP_P, int Grid_Num, int Smallest_MP,double alpha, /* input */
        int asize1, double Mo[], double Lo[], double DMo[], /* work */
        double Mi[],double Li[], double DMi[],
        double Ltmp[], double Vpot[] )
{
  static int MatP,MatP0,MatP1,i,ie;
  double Trial_ep,tmp0,Trial_D,ratio,kappa;

  /* set kappa for relativistic calculation */

  if (GL==0){
    kappa = -1.0;    /* for j=l+1/2 */
  } 
  else{ 
    if      (so==0)  kappa = (double)(-GL-1);    /* for j=l+1/2 */
    else if (so==1)  kappa = (double)GL;         /* for j=l-1/2 */
  }

  /* rough search */

  for ( ie=0; ie<Num_Part; ie++) {

    Trial_ep = min + (max-min)*ie/(Num_Part-1);

    if (Equation_Type==0){      /* Schrodinger equation */
      Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }
    else if (Equation_Type==1){ /* scalar relativistic equation */
      Hamming_O(2,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      Hamming_I(2,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }
    else if (Equation_Type==2){ /* full relativistic equation */
      Hamming_O(3,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      Hamming_I(3,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }

    if (LL_grid<=UL_grid){
      if (LL_grid<=Grid_Num*MatP_ratio && Grid_Num*MatP_ratio<=UL_grid){
        MatP = Grid_Num*MatP_ratio;         
      }
      else{
        MatP = (LL_grid + UL_grid)/2;
      } 
    }    
    else{
      MatP = LL_grid;
    }

    ratio = Lo[MatP]/Li[MatP];
    Trial_D = Mo[MatP] - ratio*Mi[MatP];

    /* recalculation of the number of nodes */

    for (i=0; i<=MatP; i++){
      Ltmp[i] = Lo[i];
    }
    for (i=(MatP+1); i<Grid_Num; i++){
      Ltmp[i] = ratio*Li[i];
    }
    node = 0;
    for (i=1; i<=(Grid_Num-2); i++){
      if (Ltmp[i]*Ltmp[i+1]<0.0) node++;
    }

    if (Check_switch==1){
      printf("X MPAO GL=%i Trial_ep=%15.12f , ie=%d",GL,Trial_ep,ie);
      printf(" Trial_D=%15.12f  node=%2d  MatP=%2d CTP=%d\n",
               Trial_D,node,MatP,CTP);
    }

    search_e[ie]    = Trial_ep;
    search_node[ie] = node;
    search_D[ie]    = Trial_D;
    search_MatP[ie] = MatP;

  } /* ie: index for energy */

  /* rough search, end */
}



static int search_bisection(
           int NumL, double Min_ep, double Max_ep, double Min_D, double Max_D,
           int Min_MatP, int Max_MatP,
           double Criterion, /* input */
           int retnode, double Trial_ep, double Trial_D, /* output */
	   int so, int GL, int CTP_P, int Grid_Num, int Smallest_MP,double alpha, /* input */
	   int Num_Part, double Mo[], double Lo[], double DMo[], /* work */
	   double Mi[],  double Li[], double DMi[],
	   double Ltmp[], double Vpot[], 
           double *E, double *WF)
{
  int MatP,MatP0,MatP1,i,  po,po2, node,i1; 
  double  pTrial_D, ratio,ep,tmp0,kappa;

  /* set kappa for relativistic calculation */

  if (GL==0){
    kappa = -1.0;    /* for j=l+1/2 */
  } 
  else{ 
    if      (so==0)  kappa = (double)(-GL-1);    /* for j=l+1/2 */
    else if (so==1)  kappa = (double)GL;         /* for j=l-1/2 */
  }

  if (Check_switch==1) {
  printf("bisection start, node=%d,  e=%f %f, D=%f %f MatP=%d %d\n",
     NumL,Min_ep,Max_ep,Min_D,Max_D,Min_MatP,Max_MatP);
  }

  po=0;
  i1=0;
  Trial_D=0;

  do {
    
    Trial_ep = 0.50*(Min_ep + Max_ep);

    if (Equation_Type==0){      /* Schrodinger equation */
      Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      i1++;
      Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }
    else if (Equation_Type==1){ /* scalar relativistic equation */
      Hamming_O(2,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      i1++;
      Hamming_I(2,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }
    else if (Equation_Type==2){ /* full relativistic equation */
      Hamming_O(3,GL,Trial_ep,kappa,Mo,Lo,DMo,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
      i1++;
      Hamming_I(3,GL,Trial_ep,kappa,Mi,Li,DMi,Vpot,Vpot,RNG[NVPS[0]][LVPS[0]]);
    }

    if (LL_grid<=UL_grid){
      if (LL_grid<=Grid_Num*MatP_ratio && Grid_Num*MatP_ratio<=UL_grid){
        MatP = Grid_Num*MatP_ratio;         
      }
      else{
        MatP = (LL_grid + UL_grid)/2;
      } 
    }    
    else{
      MatP = LL_grid;
    }

    ratio = Lo[MatP]/Li[MatP];
    pTrial_D = Trial_D;  
    Trial_D = Mo[MatP] - ratio*Mi[MatP];

    /* recalculation of the number of nodes */

    for (i=0; i<=MatP; i++){
      Ltmp[i] = Lo[i];
    }
    for (i=(MatP+1); i<=(Grid_Num-1); i++){
      Ltmp[i] = ratio*Li[i];
    }

    node = 0;
    for (i=1; i<=(Grid_Num-2); i++){
      if (Ltmp[i]*Ltmp[i+1]<0.0) node++;
    }

    if (Check_switch==1){
    printf("Y MPAO GL=%d NumL=%d Trial_ep=%e Trial_D=%e node=%d MatP=%d CTP=%d\n",
           GL,NumL,Trial_ep,Trial_D,node,MatP,CTP);
    }

    if (fabs(Trial_D)<=Criterion||Max_ep-Min_ep<1.0e-12){

      if (node==NumL){
	po = 1;
	po2 = 0;
            
	ep = Trial_ep;
	*E = ep;

	/**************************************************
          WF is P in JRCAT note 98p and u in AIST note 119p

          ref.  JRCAT note 98p and AIST note 119p
	**************************************************/
           
	for (i=0; i<=MatP; i++){
	  WF[i] = pow(MRV[i],(double)GL+1.0)*Lo[i];
	}
	for (i=(MatP+1); i<=(Grid_Num-1); i++){
	  WF[i] = pow(MRV[i],(double)GL+1.0)*ratio*Li[i];
	}
            
	Simpson(WF);

	printf("<PAO>  GL=%2d  Multiplicity=%2d node=%2d  %15.12f ",
	       GL,NumL,node,Norm_Core0(WF,Grid_Num-1));
	printf("eigenvalue=%15.12f (Hartree)\n",ep);

        if (fabs(Trial_D)>Criterion) {
          printf("warning: fabs(Trial_D)>Criterion, Trial_D=%e ,Criterion= %e at MatP=%d MRV=%lf\n\n",
            Trial_D,Criterion,MatP,MRV[MatP]);
        }

	/*
        printf("MatP=%d MRV=%lf\n",MatP,MRV[MatP]);
	*/

      }

      else{

	/**************************************************
              The number of nodes is not correct.
         In this case, this routine retries to search the
	   true eigenvalue by using a smaller StepE.
	**************************************************/

	po = 3;
      }

    }
    else{
      if (sgn(Trial_D)==sgn(Max_D)){
	Max_ep = Trial_ep;
	Max_D = Trial_D;
      }
      else{
	Min_ep = Trial_ep;
	Min_D = Trial_D;
      }
      if (Check_switch==1){
      printf("Y MPAO new min max=%e %e diff=%e\n",Min_ep,Max_ep,Max_ep-Min_ep);
      }
    } 

  } while(po==0);  


  if (po==1) { return 1; }
  else { return 0; }
}













static double Norm_Core0(double *af, int ic)
{
  int i;
  double sum0,result;

  sum0 = 0.0; 
  for (i=1; i<ic; i++){
    sum0 += af[i]*af[i]*MRV[i];
  }

  result = 0.5*dx*(af[0]*af[0]*MRV[0] + af[ic]*af[ic]*MRV[ic] + 2.0*sum0);

  return result;
}





double Norm_Core(double *af, double rc)
{
  static int i,n,j,l,nf,fg;
  static double r,rmin,rmax,Sr,Dr,sum,dum,fv[3];

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  rmax = rc;
  Sr = rmax + rmin;
  Dr = rmax - rmin;

  sum = 0.0;
  for (i=0; i<=(n-1); i++){
    r = 0.50*(Dr*x[i] + Sr);
    Spline_Func(r,fv,af);
    sum += 0.5*Dr*fv[0]*fv[0]*w[i];
  }

  return sum;
}




double Dlog_Norm(int l,
                 double a[ASIZE6][ASIZE6],
                 double a0[ASIZE6][ASIZE6],
                 double b[ASIZE6], double c[ASIZE6], double c2,
                 double Ncore,
                 double rc, double r2, double r3, double r4)
{
  static int i,j;
  static double c0,c4,c6,c8,c10,c12;
  static double Norm2,DLN;

  for (i=0; i<=4; i++){
    for (j=0; j<=4; j++){
      a[i][j] = a0[i][j];
    }
  }

  c4 = -c2*c2/(2.0*l + 5.0);

  a[0][5] = b[0] - c2*r2 - c4*r4;
  a[1][5] = b[1] - 2.0*c2*rc - 4.0*c4*r3;
  a[2][5] = b[2] - 2.0*c2 - 12.0*c4*r2;
  a[3][5] = b[3] - 24.0*c4*rc;
  a[4][5] = b[4] - 24.0*c4;

  Gauss_LEQ(4,a,c);
  c0  = c[0];
  c6  = c[1];
  c8  = c[2];
  c10 = c[3];
  c12 = c[4];
  c[5] = c4;
  Norm2 = Norm_Core2(l,rc,c0,c2,c4,c6,c8,c10,c12);
  DLN = log(Ncore) - 2.0*c0 - log(Norm2);
  return DLN;
}


double Norm_Core2(int l,double rc,
                  double c0,
                  double c2,
                  double c4,
                  double c6,
                  double c8,
                  double c10,
                  double c12)
{
  static int i,n,j,nf,fg;
  static double r1,r2,r4,r6,r8,r10,r12,p;
  static double rmin,rmax,Sr,Dr,sum,dum,fv[3];

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  rmax = rc;
  Sr = rmax + rmin;
  Dr = rmax - rmin;

  sum = 0.0;
  for (i=0; i<=(n-1); i++){
    r1  = 0.50*(Dr*x[i] + Sr);
    r2  = r1*r1;
    r4  = r2*r2;
    r6  = r4*r2;
    r8  = r4*r4;
    r10 = r8*r2;
    r12 = r10*r2;
    p = c0 + c2*r2 + c4*r4 + c6*r6 + c8*r8 + c10*r10 + c12*r12;
    dum = pow(r1,2.0*((double)l+1.0))*exp(2.0*p - 2.0*c0);
    sum = sum + 0.5*Dr*dum*w[i];
  }
  return sum;
}






void Spline_Func(double R, double fv[3], double *af)
{
  int mp_min,mp_max,m;
  double h1,h2,h3,f1,f2,f3,f4;
  double g1,g2,x1,x2,y1,y2,f,Df,D2f;
  double result;

  mp_min = 4;
  mp_max = Grid_Num - 1;
 
  if (R<MRV[4]){
    m = 5;
  }
  else if (MRV[Grid_Num-1]<R){
    m = Grid_Num - 3;
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
  }

  /****************************************************
                 Spline like interpolation
  ****************************************************/

  h1 = MRV[m-1] - MRV[m-2];
  h2 = MRV[m]   - MRV[m-1];
  h3 = MRV[m+1] - MRV[m];

  f1 = af[m-2];
  f2 = af[m-1];
  f3 = af[m];
  f4 = af[m+1];

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
              calculate the values at R
  ****************************************************/

  g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
  g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

  x1 = R - MRV[m-1];
  x2 = R - MRV[m];
  y1 = x1/h2;
  y2 = x2/h2;

  f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);

  Df = 2.0*y2/h2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + y2*y2*(2.0*f2 + h2*g1)/h2
    + 2.0*y1/h2*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
    - y1*y1*(2.0*f3 - h2*g2)/h2;
    
  D2f =  2.0*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
    + 4.0*y2*(2.0*f2 + h2*g1)
    + 2.0*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1)
    - 4.0*y1*(2.0*f3 - h2*g2);
  D2f = D2f/h2/h2;

  fv[0] = f;
  fv[1] = Df;
  fv[2] = D2f;
}










static double Int_phi0_phi1(double *phi0,double *phi1,int ic)
{
  int i;
  double sum0,result;

  sum0 = 0.0; 
  for (i=1; i<ic; i++){
    sum0 += phi0[i]*phi1[i]*MRV[i];
  }

  result = 0.5*dx*(phi0[0]*phi1[0]*MRV[0] + phi0[ic]*phi1[ic]*MRV[ic] + 2.0*sum0);

  return result;
}


