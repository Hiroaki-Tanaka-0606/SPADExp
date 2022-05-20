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








static double ip_phi0_phi1(double phi0[ASIZE1],double phi1[ASIZE1],double rc);

static double Norm_Core(double af[ASIZE1], double rc);




static double Auxiliary_Func(int order, double r, double AuxF0[11]);

static void Generate_AuxF(double AuxF[11]);

static void Generate_NCWF_PS( int so, int L, int mul,
                              int Num_Node[2][10][20],
                              int Num_AuxF[20], 
                              double pos_node[20],
                              double coes_auxF[20][20],
                              double AuxF0[11],
                              int is[2][10][20][20],
                              int ie[2][10][20][20],
                              int icut[2][10][20],
                              double rs[2][10][20][20],
                              double re[2][10][20][20],
                              double ***EAE, double ****WFAE,
                              double *VAE_Tmp, 
                              double ****WFPS, double ****WFPS0,
                              double ****GVPS);


static void Do_PatchWork(int nloop, int so, int L, int mul,
                         double slope,
                         int Num_Node[2][10][20],
                         int Num_AuxF[20], 
                         double pos_node[20],
                         double coes_auxF[20][20],
                         double AuxF0[11],
                         int is[2][10][20][20],
                         int ie[2][10][20][20],
                         int icut[2][10][20],
                         double rs[2][10][20][20],
                         double re[2][10][20][20],
                         double ***EAE, double ****WFAE,
                         double *VAE_Tmp, 
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
                          double ****GVPS);

static void Generate_TMVPS(double ***EAE, double ****WFAE, double VAE[], 
                           double ***EPS, double ****WFPS, double ****WFPS0,
                           double ****GVPS);

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
  static double max_SF,r,ef;
  static double ep,rc,a2,sum,dnorm,da,dnorm_p;
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

  Generate_TMVPS(EAE, WFAE, V_all_ele, EPS, WFPS, WFPS0, GVPS);

  /*****************************************************
   calculate pseudo-eigenfunctions for the TM potential
  *****************************************************/

  Calc_BoundStates("PS", EPS, WFPS0, V_all_ele, GVPS);

  /*
  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f %15.12f %15.12f\n",
             MRV[i],WFAE[0][0][0][i],WFAE[0][3][0][i],WFPS0[0][0][0][i],WFPS0[0][3][0][i]);

  }
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

  /*****************************************************
        generate a pseudopotential of each state
  *****************************************************/

  printf("A1\n");

  
  Generate_GVPS(EAE, WFAE, V_all_ele, EPS, WFPS, WFPS0, GVPS);

  




  exit(0);

















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
                          double ****GVPS)
{
  int so,L,mul,m,po;
  int i,j,i0,i1,i2;
  int k,k1,k2,k3,k4,k5;
  int k6,k7,k8,k9,k10;
  int knum,knum0,num,ndiv;
  int icut[2][10][20];
  int nstep[20];
  int is[2][10][20][20];
  int ie[2][10][20][20];
  int Num_Node[2][10][20];
  int Num_AuxF[20];
  double rs[2][10][20][20];
  double re[2][10][20][20];
  double adjusted_rcut[2][10][20][20];
  double pos_node[20];
  double coes_auxF[20][20];
  double AuxF0[11];
  double r0,r1,r2,rcut,dr,fv[3];
  double rc,r,r3,r4,r5,r6,r7,r8,r9;
  double start_posnode,end_posnode,dstep_np;
  double start_caf,end_caf,dstep_caf;
  double *VAE_Tmp;

  /* allocation of array */

  VAE_Tmp = (double*)malloc(sizeof(double)*(Grid_Num+10)); 
  for (i=0; i<Grid_Num; i++) VAE_Tmp[i] = VAE[i];

  /*********************************************************
                       set up rs and re
  *********************************************************/

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){
      for (mul=1; mul<GVPS_ProNum; mul++){

        /*********************************************************
          adjust the cutoff position so that the cutoff can not
          be located on a little outside region of a node   
        *********************************************************/

        printf("adjust the cutoff so=%2d L=%2d mul=%2d\n",so,L,mul);


	m = GI_VPS[L][0];
	rc = VPS_Rcut[m];

        /* find i1 for rc */

	po = 0;
	i = 0;

	do {
      
	  if (rc<=MRV[i] && po==0){
	    i1 = i;
	    po = 1;
	  }
      
	  i++;
	} while (po==0 && i<Grid_Num);
       
	if (po==0){
	  printf("The cutoff radius is too large.\n");
	  exit(0);
	}

        /* find i2 for the nearest node with shorter r to i1 */

	po = 0;
	i = i1;

	do {

          if (WFPS0[so][L][mul][i]*WFPS0[so][L][mul][i-1]<0.0){
            i2 = i;
            po = 1;
	  }

          i--;
	} while (po==0 && 1<i);

        /* adjust the cutoff if necessary */ 

        r1 = MRV[i1];
        r2 = MRV[i2];
 
        dr = r1 - r2;

        if (dr<0.3) {
          icut[so][L][mul] = i2 - 5;
          rcut = MRV[icut[so][L][mul]];
        }
        else {
          icut[so][L][mul] = i1;
          rcut = MRV[icut[so][L][mul]];
        }

	/*
        printf("icut=%2d rcut=%15.12f\n",icut[so][L][mul],rcut);
	*/

        /************************************
          find the position of nodes with
          a shorter r compared to rcut
        ************************************/

        printf("find the position of nodes so=%2d L=%2d mul=%2d\n",so,L,mul);

        node = 0;

        is[so][L][mul][node] = 0;
        ie[so][L][mul][node] = 20;

        rs[so][L][mul][node] = MRV[is[so][L][mul][node]];
        re[so][L][mul][node] = MRV[ie[so][L][mul][node]];

        for (i=1; i<=(icut[so][L][mul]-2); i++){

          if (WFPS0[so][L][mul][i]*WFPS0[so][L][mul][i+1]<0.0){

            node++;

            is[so][L][mul][node] = i - 10;
            ie[so][L][mul][node] = i + 10;

            if (is[so][L][mul][node]<0){           
              is[so][L][mul][node] = 0;
	    }
            if ((Grid_Num-1)<ie[so][L][mul][node]){ 
              ie[so][L][mul][node] = Grid_Num - 1;
	    }

            rs[so][L][mul][node] = MRV[is[so][L][mul][node]];
            re[so][L][mul][node] = MRV[ie[so][L][mul][node]];

	  }
	}

        Num_Node[so][L][mul] = node;

        /* set information of rcut */

        is[so][L][mul][node+1] = icut[so][L][mul];
        rs[so][L][mul][node+1] = MRV[icut[so][L][mul]];

      } /* mul */     
    } /* L */
  } /* so */

  /*********************************************************
     construct WFPS

     1.  
     2. 

  *********************************************************/

  /* generate an auxiliary function */

  Generate_AuxF(AuxF0);

  /*
  {
    int nx;
    double xs,xe,stepx,x;
    double f1,f2,f3,f4;
 
    nx = 1000;

    xs = -1.0;
    xe = 1.0;
    stepx = (xe-xs)/(double)(nx-1);  

    for (i=0; i<nx; i++){
      x = xs + (double)i*stepx;

      f1 =  Auxiliary_Func(1,x,AuxF0); 
      f2 =  Auxiliary_Func(2,x,AuxF0); 
      f3 =  Auxiliary_Func(3,x,AuxF0); 
      f4 =  Auxiliary_Func(4,x,AuxF0); 

      printf("%15.12f %15.12f %15.12f %15.12f %15.12f\n",x,f1,f2,f3,f4);
    }
  }
  exit(0);
  */


  /* calculate GVPS */

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){
      for (mul=1; mul<GVPS_ProNum; mul++){

        /* set the number of auxiliary functions */

        printf("set auxiliary functions so=%2d L=%2d mul=%2d\n",so,L,mul);

        for (k=0; k<=Num_Node[so][L][mul]; k++){
          Num_AuxF[k] = 0;
        }
        
        if (mul!=Num_Node[so][L][mul]){

          k = 0;
          m = 0;

          do {

            Num_AuxF[k]++;

            m++; 
            k++; 

            if (Num_Node[so][L][mul]<k) k = 0;
 
          } while ( m<(mul-Num_Node[so][L][mul]) );
	}

        /* one-dimensionalize the multiple loops */



        ndiv = 2;


        knum = 1; 
        for (k=0; k<mul; k++){
          knum *= ndiv;
	}

        /* run the one-dimensionalized loop */

        for (k=0; k<knum; k++){

          /* get each nstep */

          num = 1;
          knum0 = k;
          for (k1=0; k1<mul; k1++){
            nstep[k1] = (knum0%(10*num))/num;
            knum0 -= nstep[k1]*num; 
            num *= ndiv;
          }

          /* set parameters for the position of nodes */

          start_posnode = -0.0;
          end_posnode   =  0.0;


          dstep_np = (end_posnode - start_posnode)/(double)(ndiv-1);

          for (k1=1; k1<=Num_Node[so][L][mul]; k1++){
            pos_node[k1] = start_posnode + nstep[k1-1]*dstep_np;                 
          }

          /* set parameters for coefficients of auxiliary functions */

          start_caf = -0.0;
          end_caf   =  0.0;

          dstep_caf = (end_caf - start_caf)/(double)(ndiv-1); 

          k3 = Num_Node[so][L][mul];
          for (k1=0; k1<=Num_Node[so][L][mul]; k1++){
            for (k2=0; k2<Num_AuxF[k1]; k2++){
              coes_auxF[k1][k2] = start_caf + nstep[k3]*dstep_caf;
              k3++;
            }
	  }

          /*  */

          printf("Generate_NCWF_PS so=%2d L=%2d mul=%2d\n",so,L,mul);
 
          Generate_NCWF_PS( so, L, mul, Num_Node, Num_AuxF, pos_node, coes_auxF, AuxF0,
                            is, ie, icut, rs, re, EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS);





	} /* k */



        
 
        

        





	/*
        if (Num_Node[so][L][mul]==0){
               
	}

        else {
          Do_PatchWork(so, L, mul, Num_Node, is, ie, icut, rs, re,
                       EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS); 
        }
	*/


	/*
          Do_PatchWork(so, L, mul, Num_Node, is, ie, icut, rs, re,
                       EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS); 
	*/


      } /* mul */
    } /* L */
  } /* so */






  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){
      for (mul=1; mul<GVPS_ProNum; mul++){
        for (node=0; node<=Num_Node[so][L][mul]; node++){

          printf("so=%2d L=%2d mul=%2d node=%2d rs=%15.12f re=%15.12f\n",
                  so,L,mul,node,rs[so][L][mul][node],re[so][L][mul][node]);

	}
      }
    }
  }




  for (i=0; i<Grid_Num; i++){
    printf("%15.12f %15.12f %15.12f %15.12f\n",
            MRV[i],WFAE[0][1][1][i],WFPS0[0][1][1][i],WFPS[0][1][1][i]);

  }









  /* freeing of array */

  free(VAE_Tmp);

}


static double Auxiliary_Func(int order, double r, double AuxF0[11])
{
  double p,f,r2,r4,r6,r8,r10;
  
  r2 = r*r;
  r4 = r2*r2;
  r6 = r4*r2;
  r8 = r6*r2;
  r10 = r8*r2; 

  if (1.0<fabs(r)){
    printf("invalid r for Auxiliary_Func\n");
    exit(0);
  }

  if (order<0 || 11<order){
    printf("invalid order for Auxiliary_Func\n");
    exit(0);
  }

  switch (order){
   
  case 0:

    p = 1.0; 
    break;

  case 1:

    p = 1.0;      
    break;

  case 2:

    p = r;      
    break;

  case 3:

    p = 0.5*(3.0*r*r - 1.0);
    break;

  case 4:

    p = 0.5*(5.0*r*r*r - 3.0*r);
    break;

  case 5:

    p = (35.0*r*r*r*r - 30.0*r*r + 3.0)/8.0;
    break;

  case 6:

    p = (63.0*r*r*r*r*r - 70.0*r*r*r + 15.0*r)/8.0;
    break;
  }


  if (order==0){
    f = 0.0;
  }
  else{
    f = (AuxF0[0] + AuxF0[2]*r2 + AuxF0[4]*r4 + AuxF0[6]*r6 + AuxF0[8]*r8 + AuxF0[10]*r10)*p;
  }

  return f;
}


static void Generate_AuxF(double AuxF[11])
{
  int i;
  double a[ASIZE6][ASIZE6];
  double b[ASIZE6];
  
  a[0][0] = 1.0;
  a[0][1] = 1.0;
  a[0][2] = 1.0;
  a[0][3] = 1.0;
  a[0][4] = 1.0;
  a[0][5] =-1.0;

  a[1][0] = 2.0;
  a[1][1] = 4.0;
  a[1][2] = 6.0;
  a[1][3] = 8.0;
  a[1][4] =10.0;
  a[1][5] = 0.0;

  a[2][0] = 2.0;
  a[2][1] =12.0;
  a[2][2] =30.0;
  a[2][3] =56.0;
  a[2][4] =90.0;
  a[2][5] = 0.0;

  a[3][0] = 0.0;
  a[3][1] =24.0;
  a[3][2] =120.0;
  a[3][3] =336.0;
  a[3][4] =720.0;
  a[3][5] = 0.0;

  a[4][0] = 0.0;
  a[4][1] =24.0;
  a[4][2] =360.0;
  a[4][3] =1680.0;
  a[4][4] =5040.0;
  a[4][5] = 0.0;

  Gauss_LEQ(4,a,b);

  for (i=0; i<11; i++) AuxF[i] = 0.0;

  AuxF[0]  = 1.0;
  AuxF[2]  = b[0];
  AuxF[4]  = b[1];
  AuxF[6]  = b[2];
  AuxF[8]  = b[3];
  AuxF[10] = b[4];

}



static void Generate_NCWF_PS( int so, int L, int mul,
                              int Num_Node[2][10][20],
                              int Num_AuxF[20], 
                              double pos_node[20],
                              double coes_auxF[20][20],
                              double AuxF0[11],
                              int is[2][10][20][20],
                              int ie[2][10][20][20],
                              int icut[2][10][20],
                              double rs[2][10][20][20],
                              double re[2][10][20][20],
                              double ***EAE, double ****WFAE,
                              double *VAE_Tmp, 
                              double ****WFPS, double ****WFPS0,
                              double ****GVPS) 
{
  int po,nloop,i;
  double slope_min,slope_max,slope,rc; 
  double NormAE,dNormPS,dNormPS_min,dNormPS_max;

  /* find the upper and lower bound of slope */

  slope_min = 0.7;
  slope_max = 1.3;
   
  i = icut[so][L][mul];
  rc = MRV[i];

  NormAE = Norm_Core0(WFAE[so][L][mul],icut[so][L][mul]);

  po = 0;
  nloop = 0;

  do {

    slope = slope_min;
    Do_PatchWork( nloop, so, L, mul, slope, Num_Node, Num_AuxF, pos_node, coes_auxF, AuxF0,
		  is, ie, icut, rs, re, EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS );

    dNormPS_min = Norm_Core0(WFPS[so][L][mul],icut[so][L][mul]) - NormAE;

    slope = slope_max;
    Do_PatchWork( nloop, so, L, mul, slope, Num_Node, Num_AuxF, pos_node, coes_auxF, AuxF0,
		  is, ie, icut, rs, re, EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS );

    dNormPS_max = Norm_Core0(WFPS[so][L][mul],icut[so][L][mul]) - NormAE;

    if (0.0<dNormPS_min*dNormPS_max){
      slope_min /= 1.1;
      slope_max *= 1.1;

      printf("expand region so=%2d L=%2d mul=%2d\n",so,L,mul);

    }
    else{
      po = 1;
    }

    nloop++;

  } while (po==0 && nloop<100);

  printf("find region so=%2d L=%2d mul=%2d\n",so,L,mul);

  if (po==0){
    printf("0.0<dNormPS_min*dNormPS_max\n");
    exit(0);
  }

  /* refine an optimum slope by a bisection method */

  po = 0;
  nloop = 0;

  do {

    slope = 0.5*(slope_min + slope_max);
    Do_PatchWork( nloop, so, L, mul, slope, Num_Node, Num_AuxF, pos_node, coes_auxF, AuxF0,
                  is, ie, icut, rs, re, EAE, WFAE, VAE_Tmp, WFPS, WFPS0, GVPS );

    dNormPS = Norm_Core0(WFPS[so][L][mul],icut[so][L][mul]) - NormAE;

    printf("L=%2d mul=%2d nloop=%4d slope=%15.12f dNormPS=%18.15f\n",L,mul,nloop,slope,dNormPS);

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








static void Do_PatchWork(int nloop,
                         int so, int L, int mul,
                         double slope,
                         int Num_Node[2][10][20],
                         int Num_AuxF[20], 
                         double pos_node[20],
                         double coes_auxF[20][20],
                         double AuxF0[11],
                         int is[2][10][20][20],
                         int ie[2][10][20][20],
                         int icut[2][10][20],
                         double rs[2][10][20][20],
                         double re[2][10][20][20],
                         double ***EAE, double ****WFAE,
                         double *VAE_Tmp, 
                         double ****WFPS, double ****WFPS0,
                         double ****GVPS)
{
  int i,j,i0,i1,node,k2,po;
  int is_tmp[20];
  int ie_tmp[20];
  double u0,u1,v0,v1,v2,sum;
  double a[ASIZE6][ASIZE6];
  double b[ASIZE6];
  double Q_Coes[20][10];
  double fv[3];
  double r,r0,r1,r2,r3,r4,r5;
  double r6,r7,r8,r9,dr0,dr1;
  double factor,center,rl,rr,width2;

  for (node=0; node<=Num_Node[so][L][mul]; node++){

    i0 = ie[so][L][mul][node  ];
    r0 = re[so][L][mul][node  ];

    i1 = is[so][L][mul][node+1];
    r1 = rs[so][L][mul][node+1];

    if (node==0){
      dr0 = 0.0;
    }
    else {
      dr0 = pos_node[node];
    }

    if (node==Num_Node[so][L][mul]){
      dr1 = 0.0;
    }
    else{
      dr1 = pos_node[node+1];
    }

    /* for r0 */

    r = r0 + dr0;
    r2  = r*r;
    r3  = r2*r;
    r4  = r2*r2;
    r5  = r4*r;
    r6  = r4*r2;
    r7  = r6*r;
    r8  = r4*r4;
    r9  = r8*r;
                      
    a[0][0] = 1.0;
    a[0][1] = r;
    a[0][2] = r2; 
    a[0][3] = r3; 
    a[0][4] = r4; 
    a[0][5] = r5; 
    a[0][6] = r6; 
    a[0][7] = r7; 
    a[0][8] = r8; 
    a[0][9] = r9; 

    a[1][0] = 0.0;
    a[1][1] = 1.0;
    a[1][2] = 2.0*r; 
    a[1][3] = 3.0*r2; 
    a[1][4] = 4.0*r3; 
    a[1][5] = 5.0*r4; 
    a[1][6] = 6.0*r5; 
    a[1][7] = 7.0*r6; 
    a[1][8] = 8.0*r7; 
    a[1][9] = 9.0*r8; 

    a[2][0] = 0.0;
    a[2][1] = 0.0;
    a[2][2] = 2.0; 
    a[2][3] = 6.0*r; 
    a[2][4] = 12.0*r2; 
    a[2][5] = 20.0*r3; 
    a[2][6] = 30.0*r4; 
    a[2][7] = 42.0*r5; 
    a[2][8] = 56.0*r6; 
    a[2][9] = 72.0*r7; 

    a[3][0] = 0.0;
    a[3][1] = 0.0;
    a[3][2] = 0.0; 
    a[3][3] = 6.0; 
    a[3][4] = 24.0*r; 
    a[3][5] = 60.0*r2; 
    a[3][6] = 120.0*r3; 
    a[3][7] = 210.0*r4; 
    a[3][8] = 336.0*r5; 
    a[3][9] = 504.0*r6; 

    a[4][0] = 0.0;
    a[4][1] = 0.0;
    a[4][2] = 0.0; 
    a[4][3] = 0.0; 
    a[4][4] = 24.0; 
    a[4][5] = 120.0*r; 
    a[4][6] = 360.0*r2; 
    a[4][7] = 840.0*r3; 
    a[4][8] = 1680.0*r4; 
    a[4][9] = 3024.0*r5; 
          
    Spline_Func(r0,fv,WFPS0[so][L][mul]);
    u0 = slope*fv[0];
    u1 = slope*fv[1];

    Spline_Func(r0,fv,GVPS[so][L][0]);
    v0 = fv[0];
    v1 = fv[1];
    v2 = fv[2];

    a[0][10] = u0;
    a[1][10] = u1;
    a[2][10] = u0*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r0/r0);
    a[3][10] = u1*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r0/r0)
             + u0*(2.0*v1 - 2.0*(double)L*((double)L+1.0)/r0/r0/r0);
    a[4][10] = a[2][10]*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r0/r0)
             + 4.0*u1*(v1 - (double)L*((double)L+1.0)/r0/r0/r0)
             + 2.0*u0*(v2 + 3.0*(double)L*((double)L+1.0)/r0/r0/r0/r0);

    /* for r1 */

    r = r1 + dr1;
    r2  = r*r;
    r3  = r2*r;
    r4  = r2*r2;
    r5  = r4*r;
    r6  = r4*r2;
    r7  = r6*r;
    r8  = r4*r4;
    r9  = r8*r;
                      
    a[5][0] = 1.0;
    a[5][1] = r;
    a[5][2] = r2; 
    a[5][3] = r3; 
    a[5][4] = r4; 
    a[5][5] = r5; 
    a[5][6] = r6; 
    a[5][7] = r7; 
    a[5][8] = r8; 
    a[5][9] = r9; 

    a[6][0] = 0.0;
    a[6][1] = 1.0;
    a[6][2] = 2.0*r; 
    a[6][3] = 3.0*r2; 
    a[6][4] = 4.0*r3; 
    a[6][5] = 5.0*r4; 
    a[6][6] = 6.0*r5; 
    a[6][7] = 7.0*r6; 
    a[6][8] = 8.0*r7; 
    a[6][9] = 9.0*r8; 

    a[7][0] = 0.0;
    a[7][1] = 0.0;
    a[7][2] = 2.0; 
    a[7][3] = 6.0*r; 
    a[7][4] = 12.0*r2; 
    a[7][5] = 20.0*r3; 
    a[7][6] = 30.0*r4; 
    a[7][7] = 42.0*r5; 
    a[7][8] = 56.0*r6; 
    a[7][9] = 72.0*r7; 

    a[8][0] = 0.0;
    a[8][1] = 0.0;
    a[8][2] = 0.0; 
    a[8][3] = 6.0; 
    a[8][4] = 24.0*r; 
    a[8][5] = 60.0*r2; 
    a[8][6] = 120.0*r3; 
    a[8][7] = 210.0*r4; 
    a[8][8] = 336.0*r5; 
    a[8][9] = 504.0*r6; 

    a[9][0] = 0.0;
    a[9][1] = 0.0;
    a[9][2] = 0.0; 
    a[9][3] = 0.0; 
    a[9][4] = 24.0; 
    a[9][5] = 120.0*r; 
    a[9][6] = 360.0*r2; 
    a[9][7] = 840.0*r3; 
    a[9][8] = 1680.0*r4; 
    a[9][9] = 3024.0*r5; 

    if (node<Num_Node[so][L][mul]){
      Spline_Func(r1,fv,WFPS0[so][L][mul]);
      u0 = slope*fv[0];
      u1 = slope*fv[1];
    }
    else{
      Spline_Func(r1,fv,WFAE[so][L][mul]);
      u0 = fv[0];
      u1 = fv[1];
    }

    if (node<Num_Node[so][L][mul])
      Spline_Func(r1,fv,GVPS[so][L][0]);
    else  
      Spline_Func(r1,fv,VAE_Tmp);

    v0 = fv[0];
    v1 = fv[1];
    v2 = fv[2];

    a[5][10] = u0;
    a[6][10] = u1;
    a[7][10] = u0*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r1/r1);
    a[8][10] = u1*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r1/r1)
             + u0*(2.0*v1 - 2.0*(double)L*((double)L+1.0)/r1/r1/r1);
    a[9][10] = a[7][10]*(2.0*v0 - 2.0*EAE[so][L][mul] + (double)L*((double)L+1.0)/r1/r1)
             + 4.0*u1*(v1 - (double)L*((double)L+1.0)/r1/r1/r1)
             + 2.0*u0*(v2 + 3.0*(double)L*((double)L+1.0)/r1/r1/r1/r1);

    Gauss_LEQ(9,a,b);

    for (i=0; i<=9; i++) Q_Coes[node][i] = b[i];

  } /* node */ 

  /*********************************************************
             construct WFPS by a patchwork method                        
  *********************************************************/

  /* find the position of shifted nodes */

  is_tmp[0] = is[so][L][mul][0];
  ie_tmp[0] = ie[so][L][mul][0];
  is_tmp[Num_Node[so][L][mul]+1] = is[so][L][mul][Num_Node[so][L][mul]+1];

  for (node=1; node<=Num_Node[so][L][mul]; node++){

    r0 = rs[so][L][mul][node] + pos_node[node];
    r1 = re[so][L][mul][node] + pos_node[node];

    po = 0;
    i = 0;
    do {
      if (r0<=MRV[i]){
        po = 1;
        is_tmp[node] = i;
      }
      i++;
    } while (po==0);

    po = 0;
    i = 0;
    do {
      if (r1<=MRV[i]){
        po = 1;
        ie_tmp[node] = i;
      }
      i++;
    } while (po==0);
  }
       
  for (node=0; node<=Num_Node[so][L][mul]; node++){

    if (node==0){
      for (i=is_tmp[node]; i<ie_tmp[node]; i++){
        WFPS[so][L][mul][i] = slope*WFPS0[so][L][mul][i];
      }
    }
    else{
      for (i=is_tmp[node]; i<ie_tmp[node]; i++){
        r = MRV[i] - pos_node[node];
        Spline_Func(r,fv,WFPS0[so][L][mul]);
        WFPS[so][L][mul][i] = slope*fv[0];
      }
    }

    /*
    rl = MRV[ie_tmp[node]];
    rr = MRV[is_tmp[node+1]];
    center = 0.5*(rr + rl);
    width2 = 0.5*(rr - rl);
    */

    for (i=ie_tmp[node]; i<is_tmp[node+1]; i++){

      r = MRV[i];

      /*
      factor = 1.0;
      for (k2=0; k2<Num_AuxF[node]; k2++){
        factor += coes_auxF[node][k2]*Auxiliary_Func(k2+1, (r-center)/width2, AuxF0 );
      }
      factor = 1.0;
      */

      sum = 0.0;
      for (j=0; j<=9; j++) sum += Q_Coes[node][j]*pow(r,(double)j);

      /*
      WFPS[so][L][mul][i] = factor*sum;
      */

      WFPS[so][L][mul][i] = sum;

    }
  }

  for (i=icut[so][L][mul]; i<Grid_Num; i++){
    WFPS[so][L][mul][i] = WFAE[so][L][mul][i];
  }
}







static void Generate_TMVPS(double ***EAE, double ****WFAE, double VAE[], 
                           double ***EPS, double ****WFPS, double ****WFPS0,
                           double ****GVPS)
{
  int so,L,m,po,i;
  double ep,rc;

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<=GVPS_MaxL; L++){

      ep = EAE[so][L][0];
      m = GI_VPS[L][0]; 
      rc = VPS_Rcut[m];

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










double ip_phi0_phi1(double phi0[ASIZE1],double phi1[ASIZE1],double rc)
{
  static int i,n,j,l,nf,fg;
  static double r,rmin,rmax,Sr,Dr,sum,dum;
  static double fv0[3],fv1[3];

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
    Spline_Func(r,fv0,phi0);
    Spline_Func(r,fv1,phi1);
    sum = sum + 0.5*Dr*fv0[0]*fv1[0]*w[i];
  }

  return sum;
}


