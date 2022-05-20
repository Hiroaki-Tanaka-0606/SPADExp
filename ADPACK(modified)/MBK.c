/**********************************************************************
  MBK.c:

     MBK.c is a subroutine to calculate norm-conserving Vanderbilt 
     pseudopotentials using the MBK scheme.

  Log of MBK.c:

     13/Jan/2010  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "adpack.h"

#define debug 0
#define xmin  0.0
#define lmax  4
#define numsb 13

typedef int INTEGER;

static void Calc_new_range(/*input */
                     int *node_c_ene, double *DD_c_ene, int num_nodes,
                    /*output*/
                    double *MinE, double *MaxE , double *dE, int *num_div );

static void Calc_PAO(int state_num,
                     int so, int n, int l, int num_nodes, double Reduce_Num_grid,
                     double Vps[ASIZE1], double pw[ASIZE1]);
static void Calc_PAO_Hamman(int so, int L, double rc, double ep,
                            double PV[ASIZE1], double W[ASIZE1], 
                            double W3[ASIZE1]);
static double Norm_Core(double af[ASIZE1], double rc);
static double Norm_Core2(int l,double rc,
                         double c0,
                         double c2,
                         double c4,
                         double c6,
                         double c8,
                         double c10,
			 double c12);
static double Dlog_Norm(int l,
                        double a[ASIZE6][ASIZE6],
                        double a0[ASIZE6][ASIZE6],
                        double b[ASIZE6], double c[ASIZE6], double c2,
                        double Ncore,
			double rc, double r2, double r3, double r4);
static void Spline_Func(double R, double fv[3], double af[ASIZE1]);
static int Calc_VPS_PAO(int calc_type,
                        int so, int l, double rc, double ep, double scale_factor, 
                        double AEW[ASIZE1], double AEV[ASIZE1],
                        double  PW[ASIZE1], double PEV[ASIZE1]);
static double ip_phi0_phi1(double phi0[ASIZE1],double phi1[ASIZE1],double rc);
static void Hamman(int so,int L,double rc,double ep,double PW[ASIZE1],double PV[ASIZE1]);

static void find_zero_jx();
static void Spherical_Bessel( double x, int Lmax, double *sb, double *dsb );
static void Eigen_lapack_x(double **a, double *ko, int n0, int EVmax);


static double phi0_phi1(double phi0[ASIZE1],
			double phi1[ASIZE1]);

static double phi0_v_phi1(double phi0[ASIZE1],
			  double phi1[ASIZE1],
			  double vpot[ASIZE1]);

static double Int_RadF(double R, double RadF[ASIZE1]);



void dsyevx_(char *JOBZ, char *RANGE, char *UPLO, INTEGER *N, double *A, INTEGER *LDA, 
	     double *VL, double *VU,
	     INTEGER *IL, INTEGER *IU, double *ABSTOL, INTEGER *M, double *W, double *Z, 
	     INTEGER *LDZ, double *WORK,
	     INTEGER *LWORK, INTEGER *IWORK, INTEGER *IFAIL, INTEGER *INFO);


double zero_j[10][20];







void MBK(int state_num)
{ 
  int m,n,l,i,i1,j,k,k1,k2,ip,L,gi;
  int po,po1,numvps,ll,mul,mul1,mul2;
  int m1,n1,loop_num,loop_max,Hamman_flag;
  int num_nodes,so,m2,n2,loop_count;
  double max_SF,r,ef,rc1,diff,ep1,sr,tmp0,tmp,tmp1,tmp2;
  double dep,dep1,ep,ep0,rc,a2,dif,sum,sum2,dnorm,da,dnorm_p,Chi0,Chi1;
  double a_p,lower_a,upper_a,lower_d,upper_d;
  double a11,a12,a13,a21,a22,a23,a31,a32,a33,det;
  double scaling,scaling_min,scaling_max,rc_max;
  double W3[ASIZE1],rj[ASIZE1],beta[ASIZE4][ASIZE1];
  double Chi[ASIZE13][ASIZE4][ASIZE1];
  double phi[8][ASIZE1];
  double coe_sb[numsb+1],cmat[10][numsb+1];
  double c012[10][numsb+1],ia[10][10];
  double a0[ASIZE6][ASIZE6],c[ASIZE6],pe[ASIZE4];
  double a1[ASIZE6][ASIZE6],ia1[ASIZE6][ASIZE6];
  double B[ASIZE4][ASIZE4],iB[ASIZE4][ASIZE4];
  double Lcoes[numsb],Qcoes[numsb];
  double gradF[numsb+1],NormGF,pNormGF,F;
  double norm_all,norm_ps;
  double Vsl[ASIZE1];
  double *V0; /* dummy */
  char *s_vec[20];
  double *sb,*dsb;
  double **BB,*ko;

  /* allocation of arrays */

  sb  = (double*)malloc(sizeof(double)*(lmax+3));
  dsb = (double*)malloc(sizeof(double)*(lmax+3));

  BB = (double**)malloc(sizeof(double*)*ASIZE4);
  for (i=0; i<ASIZE4; i++){
    BB[i] = (double*)malloc(sizeof(double)*ASIZE4);
  } 

  ko = (double*)malloc(sizeof(double)*ASIZE4);

  /****************************************************
   calculate zero points of spherical bessel functions
  ****************************************************/
  
  find_zero_jx();

  /****************************************************
      find the maximum rc among all the cutoff of l
  ****************************************************/

  rc_max = -100.0; 
  for (l=0; l<ASIZE2; l++){

    for (mul=0; mul<NumVPS_L[l]; mul++){
      m = GI_VPS[l][mul];
      rc = VPS_Rcut[m];

      if (rc_max<rc) rc_max = rc; 
    }
  }

  /******************************************************
    find energies used for calculation of unbound states
  ******************************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

	for (mul=0; mul<NumVPS_L[l]; mul++){

	  m = GI_VPS[l][mul];
	  n = NVPS[m];
	  rc = VPS_Rcut[m];
	  ep = E[state_num][so][n][l];

          if ( OcpN[0][0][n][l]<1.0e-12 ){
            E[state_num][so][n][l] = VPS_ene[m];
	  }

	}
      }
    }
  }

  /***********************************************************
     generate nodeless pseudo-wave functions for each chosen 
     state and the corresponding pseudopotentials 
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

	for (mul=0; mul<NumVPS_L[l]; mul++){

	  m = GI_VPS[l][mul];
	  n = NVPS[m];
	  rc = VPS_Rcut[m];
	  ep = E[state_num][so][n][l];

          if (99.99<VPS_ene[m]){
            dep = VPS_ene[m];
	  }
          else {
            dep = 0.0;
	  }

          printf("<VPS>  Generating VPS of SO=%d n=%d l=%d by the MBK scheme....\n",so,n,l);
          printf("<VPS>  Used cutoff radius for VPS=%10.7f (a.u.)\n",rc);

          if ( 1.0e-12<=OcpN[0][0][n][l] ){

            /* calculate pseudopotential and pseudowave function by the TM method */ 

  	    po1 = Calc_VPS_PAO(0,so,l,rc,ep,1.0,PF[so][n][l],V,W2[so][m],V2[so][m]);
            Hamman_flag = 0;
	  }

          else {

            /* calculte a partial all electron wave function by the Hamman method, 
               then calculate pseudopotential and pseudowave function by the TM method */ 

            Hamman(so,l,rc,ep,PF[so][n][l],V_all_ele);
            po1 = Calc_VPS_PAO(0,so,l,rc,ep,1.0,PF[so][n][l],V,W2[so][m],V2[so][m]);
            Hamman_flag = 1;

	  }

          /* calculate a part of a chi function */

 	  for (i=0; i<Grid_Num; i++){

	    r = MRV[i];

            if (Hamman_flag==1){ 
	      if (rc<r)      Chi[so][m][i] = 0.0;
	      else           Chi[so][m][i] = (V2[so][m][i]+dep)*W2[so][m][i];
	    }

            else {
	      if (rc_max<r)  Chi[so][m][i] = 0.0;
	      else           Chi[so][m][i] = (V2[so][m][i]+dep)*W2[so][m][i];
            }

	  }
	}
      }
    }
  }

  /***********************************************************
     construct pseudo-wave functions for each chosen state 
     which satisfies generalized norm-conserving conditions
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

	m = GI_VPS[l][0];
	n = NVPS[m];

	for (mul=1; mul<NumVPS_L[l]; mul++){

	  m1 = GI_VPS[l][mul]; 
	  n1 = NVPS[m1];
	  rc1 = VPS_Rcut[m1];
	  ep1 = E[state_num][so][n1][l];

          if (99.99<VPS_ene[m1]) dep1 = VPS_ene[m1];
          else                   dep1 = 0.0;

          /* calculate excited state of V2[so][m] */ 

          if ( 1.0e-12<=OcpN[0][0][n1][l] ){

  	    Calc_PAO(state_num,so,n1,l,mul,RNG[n1][l],V2[so][m],W3);
            Hamman_flag = 0;
	  }
          else {

            Calc_PAO_Hamman(so,l,rc1,ep1,V2[so][m],W2[so][m1],W3);
            Hamman_flag = 1;
          }

          /* check sign of W3 */ 

          diff = 10000.0;
 	  for (i=0; i<Grid_Num; i++){
            if ( fabs(MRV[i]-rc1)<diff ){
              diff = fabs(MRV[i]-rc1);
              i1 = i;
            }             
	  } 
          i1 = i1 - 2 ;         
          if (i1<0) i1 = 0;

          if (W2[so][m1][i1]*W3[i1]<0.0){
   	    for (i=0; i<Grid_Num; i++){
              W3[i] = -W3[i];
	    }
          }

	  /*
          printf("ABC\n");
 	  for (i=0; i<Grid_Num; i++){
            printf("%15.12f %15.12f %15.12f\n",MRV[i],W3[i],W2[so][m1][i]);
	  }
	  */

          /* calculate difference between W3 and W2 */ 

 	  for (i=0; i<Grid_Num; i++){
            W3[i] = W3[i] - W2[so][m1][i];
	  }          

          /* find initial coefficients by expanding W3 with spherical bessel functions */ 

          for (k=0; k<numsb; k++){

	    for (i=0; i<Grid_Num; i++){
	      r = MRV[i];
	      sr = r/rc1*zero_j[l][k];
	      Spherical_Bessel(sr,lmax,sb,dsb);
	      rj[i] = r*sb[l];
	    }          

	    sr = zero_j[l][k];
	    Spherical_Bessel(sr,lmax,sb,dsb);
            tmp = 0.5*dsb[l]*dsb[l]*rc1*rc1*rc1;

	    coe_sb[k] = ip_phi0_phi1(W3,rj,rc1)/tmp;
	  }

          /* for debugging */ 

          if (debug){

            printf("reconstructed W3 l=%2d mul=%2d\n",l,mul);
	    for (i=0; i<Grid_Num; i++){

	      r = MRV[i];
	      sum = 0.0;

	      for (k=0; k<numsb; k++){

		sr = r/rc1*zero_j[l][k];
		Spherical_Bessel(sr,lmax,sb,dsb);
		sum += r*sb[l]*coe_sb[k]; 
	      }

	      printf("%15.12f %15.12f %15.12f\n",r,W3[i],sum);
	    }
	  }

          /********************************
            construct a constraint matrix
          ********************************/ 

          /* first derivative */

          for (k=0; k<numsb; k++){
            sr = zero_j[l][k]; 
            Spherical_Bessel(sr,lmax,sb,dsb);
            cmat[0][k] = zero_j[l][k]*dsb[l];
	  }          
          cmat[0][numsb] = 0.0;

          /* third derivative */

          for (k=0; k<numsb; k++){
	    sr = zero_j[l][k]; 
	    Spherical_Bessel(sr,lmax,sb,dsb);
	    cmat[1][k] = sr/rc1/rc1*((double)l*((double)l+1.0)-sr*sr)*dsb[l];
	  }
          cmat[1][numsb] = 0.0;

          /* condition of Qij where i!=j */

          for (mul2=0; mul2<mul; mul2++){

  	    m2 = GI_VPS[l][mul2]; 
	    n2 = NVPS[m2];

	    norm_all = ip_phi0_phi1(PF[so][n2][l],PF[so][n1][l],rc1);  
	    norm_ps  = ip_phi0_phi1(W2[so][m2],W2[so][m1],rc1);

	    for (k=0; k<numsb; k++){

	      for (i=0; i<Grid_Num; i++){
		r = MRV[i];
		sr = r/rc1*zero_j[l][k];
		Spherical_Bessel(sr,lmax,sb,dsb);
		rj[i] = r*sb[l];
	      }          

	      cmat[2+mul2][k] = ip_phi0_phi1(W2[so][m2],rj,rc1);  
	    }        
	    cmat[2+mul2][numsb] = norm_ps - norm_all;
	  }

          /***************************************************************
            express c0, c1, c2 by c3, c4, ...,c_(numsb-1) plus some const
          ***************************************************************/ 

          /* calculate the inverse of a0 */

          for (k=0; k<(2+mul); k++){

	    for (i=0; i<(3+mul); i++){
	      for (j=0; j<(3+mul); j++){
		a0[i][j] = 0.0;
	      }
	    }

	    for (i=0; i<(2+mul); i++){
	      for (j=0; j<(2+mul); j++){
		a0[i][j] = cmat[i][j];
	      }
	    }

	    a0[k][2+mul] = 1.0;
	    Gauss_LEQ(1+mul,a0,c);
	    for (i=0; i<(2+mul); i++) ia[i][k] = c[i];
	  }

	  /* construct a matrix showing the linear relation between first few c0,c1,.. and others */

          for (i=0; i<(2+mul); i++){
            for (j=(2+mul); j<=numsb; j++){

              sum = 0.0;              
              for (k=0; k<(2+mul); k++){
                sum -= ia[i][k]*cmat[k][j];
              }

              c012[i][j] = sum;
            }
          }

          /********************************************************
                               minimize Qii^2
          ********************************************************/ 

          /* calculate coefficiencts of the linear and quadratic terms */
          
          for (k=0; k<numsb; k++){

	    for (i=0; i<Grid_Num; i++){
	      r = MRV[i];
	      sr = r/rc1*zero_j[l][k];
	      Spherical_Bessel(sr,lmax,sb,dsb);
	      rj[i] = r*sb[l];
	    }          

	    Lcoes[k] = ip_phi0_phi1(W2[so][m1],rj,rc1);
	  }          

          for (k=0; k<numsb; k++){
	    sr = zero_j[l][k];
	    Spherical_Bessel(sr,lmax,sb,dsb);
            Qcoes[k] = 0.5*dsb[l]*dsb[l]*rc1*rc1*rc1;
	  }

          coe_sb[numsb] = 1.0;

          /* do-loop of optimization */
 
          loop_count = 1; 
          scaling = 0.000001;
          scaling_max = scaling*100.0;
          scaling_min = scaling/1000.0;

          NormGF = 1.0e+10;

          do { 

	    /* correct initial c0, c1, c2, (c3) so that the constraints can be satisfied */  

	    for (k=0; k<(2+mul); k++){

	      sum = 0.0;            
	      for (i=(2+mul); i<=numsb; i++){
		sum += c012[k][i]*coe_sb[i];
	      }
	      coe_sb[k] = sum;
	    }

	    /* calculate gradients of the function F */

	    tmp = 0.0;
	    for (k=0; k<numsb; k++){
	      tmp += 2.0*coe_sb[k]*Lcoes[k] + coe_sb[k]*coe_sb[k]*Qcoes[k];             
	    }          

            F = tmp*tmp;

            pNormGF = NormGF;
            NormGF = 0.0; 

	    for (k=(2+mul); k<numsb; k++){
            
	      sum = Lcoes[k] + coe_sb[k]*Qcoes[k];

	      for (i=0; i<(2+mul); i++){
		sum += c012[i][k]*(Lcoes[i]+coe_sb[i]*Qcoes[i]); 
	      }

	      gradF[k] = 4.0*tmp*sum; 
              NormGF += gradF[k]*gradF[k]; 
	    }          

	    /* update the coefficients */

	    for (k=(2+mul); k<numsb; k++){
              coe_sb[k] += -scaling*gradF[k]; 
	    }

            if (NormGF<pNormGF) scaling *= 2.0;
            else                scaling /= 2.0;

            if ( scaling_max<scaling ) scaling = scaling_max;
            if ( scaling<scaling_min ) scaling = scaling_min;

            if (debug){
              printf("so=%2d l=%2d mul=%2d count=%2d scaling=%15.12f F=%18.15f NormGF=%18.15f\n",
                       so,l,mul,loop_count,scaling,F,NormGF);
	    }

            /* increment loop_count */ 
 
            loop_count++;

	  } while ( 1.0e-12<NormGF );

          /* correct initial c0, c1, c2 so that the constraints can be satisfied */  

	  for (k=0; k<(2+mul); k++){

	    sum = 0.0;            
	    for (i=(2+mul); i<=numsb; i++){
	      sum += c012[k][i]*coe_sb[i];
	    }
	    coe_sb[k] = sum;
	  }

          /************************************************************
                construct a pseudo-wave function which satisfies 
                  the generalized norm conserving conditions 
	  ************************************************************/

          for (i=0; i<Grid_Num; i++){

            r = MRV[i];

            if (rc1<r){
              sum = 0.0;
	    }
            else {

	      sum = 0.0;
              sum2 = 0.0;

	      for (k=0; k<numsb; k++){
 
                tmp2 = zero_j[l][k]/rc1;
		sr = r*tmp2;
		Spherical_Bessel(sr,lmax,sb,dsb);
                tmp = r*sb[l]*coe_sb[k];
		sum += tmp;
                sum2 += tmp2*tmp2*tmp;
	      }
	    }

            W2[so][m1][i] += sum;

            /* add a part of a chi function */

            if (r<=rc1){
              Chi[so][m1][i] += (ep1+dep1)*sum - 0.5*sum2;
	    }
	  }

          /* check whether the conditions are satisfied for debugging */

          if (debug){
	    sum = 0.0;
	    for (k=0; k<numsb; k++){
	      sum += cmat[0][k]*coe_sb[k]; 
	    }
	    printf("The first derivative at the cutoff: %18.15f\n",sum);

	    sum = 0.0;
	    for (k=0; k<numsb; k++){
	      sum += cmat[1][k]*coe_sb[k]; 
	    }
	    printf("The third derivative at the cutoff: %18.15f\n",sum);

            for (mul2=0; mul2<mul; mul2++){

	      m2 = GI_VPS[l][mul2]; 
	      n2 = NVPS[m2];

	      norm_all = ip_phi0_phi1(PF[so][n2][l],PF[so][n1][l],rc1);
	      norm_ps  = ip_phi0_phi1(W2[so][m2],W2[so][m1],rc1);
	      printf("Qij: mul2=%2d all=%18.15f ps=%18.15f\n",mul2,norm_all,norm_ps);
	    }

	    norm_all = ip_phi0_phi1(PF[so][n1][l],PF[so][n1][l],rc1);
	    norm_ps  = ip_phi0_phi1(W2[so][m1],W2[so][m1],rc1);
	    printf("Qii: all=%18.15f ps=%18.15f\n",norm_all,norm_ps);

	    for (i=0; i<Grid_Num; i++){
	      printf("%15.12f %15.12f %15.12f\n",MRV[i],W2[so][m1][i],PF[so][n1][l][i]);
	    }
	  }

        } /* mul */ 
      } /* if (NumVPS_L[l]!=0) */
    } /* l */
  } /* so */

  /****************************************************
              calculate a common local part
  ****************************************************/

  Calc_Vlocal(Vlocal,V0,0);

  /***********************************************************
            calculate the screened local potential   
  ************************************************************/

  /* set OcpN0 */

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){
	for (mul=0; mul<NumVPS_L[l]; mul++){
	  n = NVPS[m];
          if (99.99<VPS_ene[m]){
            OcpN0[state_num][so][n][l] = 0.0;
	  }
	}
      }
    }
  }

  /* calculate hartree and xc parts of the screened local potential */

  Density_V(state_num,OcpN0);
  Hartree(rho_V[state_num],Vh_V);

  if (PCC_switch==0){
    if      (XC_switch==0) XC_CA(rho_V[state_num],Vxc_V,1);
    else if (XC_switch==1) XC4atom_PBE(rho_V[state_num],Vxc_V,1);
  }
  else {

    Density_PCC(rho_PCC,state_num,OcpN0);

    for (i=0; i<Grid_Num; i++){
      rho[1][i] = rho_V[state_num][i] + rho_PCC[i];
    }
    if      (XC_switch==0) XC_CA(rho[1],Vxc_V,1);
    else if (XC_switch==1) XC4atom_PBE(rho[1],Vxc_V,1);
 

  }


  /*
  for (i=0; i<Grid_Num; i++){
    r = MRV[i];
    printf("%15.12f %15.12f %15.12f %15.12f\n",
	   r,0.5*rho_V[state_num][i],0.5*rho_PCC[i],Vxc_V[i]); 
  }  
  exit(0);
  */  

  /* calculate the screened local potential */

  for (i=0; i<Grid_Num; i++){
    Vsl[i] = Vlocal[i] + Vh_V[i] + Vxc_V[i];
  }

  /***********************************************************
    add a remaining contribution of chi, coming from Vsl:
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

	for (mul=0; mul<NumVPS_L[l]; mul++){

	  m = GI_VPS[l][mul];
	  n = NVPS[m];
	  rc = VPS_Rcut[m];
	  ep = E[state_num][so][n][l];

          if ( 1.0e-12<=OcpN[0][0][n][l] ) Hamman_flag = 0;
          else                             Hamman_flag = 1;

          for (i=0; i<Grid_Num; i++){

            r = MRV[i];

            if (Hamman_flag==1){
              if (r<=rc)     Chi[so][m][i] -= Vsl[i]*W2[so][m][i];
	    }
            else{
              if (r<=rc_max) Chi[so][m][i] -= Vsl[i]*W2[so][m][i];
	    }
	  }
	}
      }
    }
  }

  /****************************************************
           enhancement of spin-orbit splitting 
  ****************************************************/

  if (Equation_Type==2){

    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

	for (mul=0; mul<NumVPS_L[l]; mul++){

	  m = GI_VPS[l][mul];
	  n = NVPS[m];

	  for (i=0; i<Grid_Num; i++){

	    Chi0 = Chi[0][m][i];
	    Chi1 = Chi[1][m][i];
	    sum = 0.5*(Chi0 + Chi1);
	    dif = 0.5*(Chi0 - Chi1);
          
	    Chi[0][m][i] = sum + SO_factor[l]*dif;
	    Chi[1][m][i] = sum - SO_factor[l]*dif;
	  }
	}
      }
    }
  }

  if (debug){

    printf("Chi\n");

    l = 0;
    so = 0;
    for (i=0; i<Grid_Num; i++){

      r = MRV[i];
      printf("%15.12f ",r);

      for (mul=0; mul<NumVPS_L[l]; mul++){

	m = GI_VPS[l][mul];

	printf("%15.12f ",Chi[so][m][i]);
      }
      printf("\n");
    }
  }

  /***********************************************************
    calculate the matrix B = <phi|chi> and its inverse matrix
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (l=0; l<ASIZE2; l++){

      if (NumVPS_L[l]!=0){

        /* calculate the matrix B */         

	for (mul1=0; mul1<NumVPS_L[l]; mul1++){

          m1 = GI_VPS[l][mul1];

  	  for (mul2=0; mul2<NumVPS_L[l]; mul2++){

	    m2 = GI_VPS[l][mul2];
            tmp = ip_phi0_phi1(W2[so][m1],Chi[so][m2],rc_max);
            B[mul1][mul2] = tmp;
	  }
	}

        /* symmetrization of B */         

	for (mul1=0; mul1<NumVPS_L[l]; mul1++){
  	  for (mul2=(mul1+1); mul2<NumVPS_L[l]; mul2++){

            tmp1 = B[mul1][mul2];
            tmp2 = B[mul2][mul1];

            tmp = 0.5*(tmp1+tmp2);

            B[mul1][mul2] = tmp;
            B[mul2][mul1] = tmp;
	  }
	}

        /* calculate the inverse matrix of B */

        for (k=0; k<NumVPS_L[l]; k++){

	  for (i=0; i<=NumVPS_L[l]; i++){
	    for (j=0; j<=NumVPS_L[l]; j++){
	      a0[i][j] = 0.0;
	    }
	  }

	  for (i=0; i<NumVPS_L[l]; i++){
	    for (j=0; j<NumVPS_L[l]; j++){
	      a0[i][j] = B[i][j];
	    }
	  }

	  a0[k][NumVPS_L[l]] = 1.0;

	  Gauss_LEQ(NumVPS_L[l]-1,a0,c);
	  for (i=0; i<NumVPS_L[l]; i++) iB[i][k] = c[i];
	}

        if (debug){

	  printf("B l=%2d\n",l);
	  for (i=0; i<NumVPS_L[l]; i++){
	    for (j=0; j<NumVPS_L[l]; j++){
	      printf("%15.12f ",B[i][j]);
	    }
	    printf("\n");
	  }


	  printf("iB l=%2d\n",l);
	  for (i=0; i<NumVPS_L[l]; i++){
	    for (j=0; j<NumVPS_L[l]; j++){
	      printf("%15.12f ",iB[i][j]);
	    }
	    printf("\n");
	  }
	}

        /* calculate |beta> */

	for (mul1=0; mul1<NumVPS_L[l]; mul1++){
          for (i=0; i<Grid_Num; i++){

	    sum = 0.0;
	    for (mul2=0; mul2<NumVPS_L[l]; mul2++){
	      m2 = GI_VPS[l][mul2];
	      sum += iB[mul2][mul1]*Chi[so][m2][i];
	    }

            beta[mul1][i] = sum;
	  }
	}

        /* diagonalize the matrix B */

	for (i=0; i<NumVPS_L[l]; i++){
	  for (j=0; j<NumVPS_L[l]; j++){
	    BB[i][j] = B[i][j];
	  }
	}
       
        Eigen_lapack_x(BB,ko,NumVPS_L[l],NumVPS_L[l]);

        if (debug){

	  for (i=0; i<NumVPS_L[l]; i++){
	    printf("ko l=%2d i=%2d ko=%15.12f\n",l,i,ko[i]);
	  }

	  printf("V l=%2d\n",l);
	  for (i=0; i<NumVPS_L[l]; i++){
	    for (j=0; j<NumVPS_L[l]; j++){
	      printf("%15.12f ",BB[i][j]);
	    }
	    printf("\n");
	  }
	}

        /* calculate non-local projectors */
        
	for (mul1=0; mul1<NumVPS_L[l]; mul1++){

          tmp2 = sqrt(fabs(ko[mul1])); /* for "normalization" */

          for (i=0; i<Grid_Num; i++){

	    sum = 0.0;
	    for (mul2=0; mul2<NumVPS_L[l]; mul2++){
	      sum += beta[mul2][i]*BB[mul2][mul1];
	    }

            VNL_W2[so][l][mul1][i] = tmp2*sum;
	  }

          proj_ene[0][so][l][mul1] = sgn(ko[mul1]); 
	}

      } /* if (NumVPS_L[l]!=0) */ 
    } /* l */
  } /* so */

  if (debug){

    so = 0;
    l = 1;
    mul = 0;

    printf("proj_ene=%15.12f\n",proj_ene[0][so][l][mul]);

    printf("VNL_W2\n");
    for (i=0; i<Grid_Num; i++){
      printf("%15.12f %15.12f\n",MRV[i],VNL_W2[so][l][mul][i]);
    }
  }

  /***********************************************************
    calculate VPS which is used in Log_DeriF to guess 
    the initial Lo and for generation of Blochl projectors 
    for the TM l-pseudopotential in case of 
    NumVPS_L[l]==1 && Blochl_pro_num!=1
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (m=0; m<Number_VPS; m++){
      for (i=0; i<Grid_Num; i++){
        VPS[state_num][so][m][i] = V2[so][m][i] - Vh_V[i] - Vxc_V[i];
      }
    }
  }

  /***********************************************************
    if NumVPS_L[l]==1 && Blochl_pro_num!=1,
    generate Blochl projectors for the TM l-pseudopotential.
  ************************************************************/

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){

      if (NumVPS_L[L]==1 && Blochl_pro_num!=1){

	gi = GI_VPS[L][0];

	/* Set non-local potential and an initial wave function */
	for (i=0; i<Grid_Num; i++){
	  VNL[0][so][L][0][i] = VPS[0][so][gi][i] - Vlocal[i];
	  phi[0][i] = W2[so][gi][i];
	}

	/* phi_m = VNL^m * phi_0 */
	for (m=1; m<Blochl_pro_num; m++){
	  for (i=0; i<Grid_Num; i++){
	    phi[m][i] = pow(VNL[0][so][L][0][i],(double)m)*phi[0][i];
	  }
	}

	/* Normalization */
	for (m=0; m<Blochl_pro_num; m++){

	  tmp0 = 1.0/sqrt(phi0_phi1(phi[m],phi[m]));

	  for (i=0; i<Grid_Num; i++){
	    phi[m][i] = tmp0*phi[m][i];
	  }
	}

	/* Gramm-Schmidt orthogonalization with a norm defined by <f|v|g> */
	for (i=0; i<Grid_Num; i++){
	  VNL_W2[so][L][0][i] = phi[0][i]; 
	}
	pe[0] = 1.0/phi0_v_phi1(VNL_W2[so][L][0],VNL_W2[so][L][0],VNL[0][so][L][0]);

	for (m=1; m<Blochl_pro_num; m++){
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = 0.0;
	  }          
	  for (n=0; n<m; n++){
	    tmp0 = phi0_v_phi1(VNL_W2[so][L][n],phi[m],VNL[0][so][L][0]);  
	    for (i=0; i<Grid_Num; i++){
	      VNL_W2[so][L][m][i] = VNL_W2[so][L][m][i] + VNL_W2[so][L][n][i]*pe[n]*tmp0;  
	    }
	  }
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = phi[m][i] - VNL_W2[so][L][m][i];
	  }
	  pe[m] = 1.0/phi0_v_phi1(VNL_W2[so][L][m],VNL_W2[so][L][m],VNL[0][so][L][0]);
	}

	/* renormalization */
	for (m=0; m<Blochl_pro_num; m++){
	  tmp0 = 1.0/sqrt(phi0_phi1(VNL_W2[so][L][m],VNL_W2[so][L][m]));
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = tmp0*VNL_W2[so][L][m][i];
	  }
	  proj_ene[0][so][L][m] = pe[m]/(tmp0*tmp0);
	}

	/* Calc v*VNL_W2 */
	for (m=0; m<Blochl_pro_num; m++){
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = VNL[0][so][L][0][i]*VNL_W2[so][L][m][i];
	  }
	}

	/* again perform normalization */
	for (m=0; m<Blochl_pro_num; m++){

          tmp2 = sqrt(fabs(proj_ene[0][so][L][m])); 

	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] *= tmp2;
	  }

          proj_ene[0][so][L][m] = sgn(proj_ene[0][so][L][m]);
	}

      } /* if (NumVPS_L[L]==1 && Blochl_pro_num!=1) */
    } /* L */
  } /* so */

  /* for std output */
  printf("\n");

  /* allocation of arrays */

  free(ko);

  for (i=0; i<ASIZE4; i++){
    free(BB[i]);
  } 
  free(BB);

  free(dsb);
  free(sb); 
}


static void Calc_PAO_Hamman(int so, int L, double rc, double ep,
                            double PV[ASIZE1], double W[ASIZE1], 
                            double W3[ASIZE1])
{
  int i,ic,po;
  double kappa,tmp,fugo;
  double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];

  Hamming_O(0,L,ep,kappa,Mo,Lo,DMo,PV,PV,0.95);

  /* find ic corresponding to rc */

  po = 0;
  for (i=0; i<Grid_Num; i++){
    if (rc<MRV[i] && po==0){
      po = 1;
      ic = i;
    }
  }

  /* calculate W3 */

  for (i=0; i<=(ic+10); i++){
    W3[i] = pow(MRV[i],(double)L+1.0)*Lo[i];
  }

  /* normalization of W3 so that the norm is same as that of W. */

  tmp = sqrt(Norm_Core(W,rc))/sqrt(Norm_Core(W3,rc));

  if (0.0<W3[ic]) fugo = 1.0;
  else            fugo =-1.0;  

  for (i=0; i<=(ic+10); i++){
    W3[i] = tmp*fugo*W3[i];
  }

}




static void Hamman(int so, int L, double rc, double ep,
                   double PW[ASIZE1], double PV[ASIZE1])
{
  int i,ic,po;
  double kappa,tmp,fugo;
  double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];

  if (L==0){
    kappa = -1.0;    /* for j=l+1/2 */
  } 
  else{ 
    if      (so==0)  kappa = (double)(-L-1);    /* for j=l+1/2 */
    else if (so==1)  kappa = (double)L;         /* for j=l-1/2 */
  }

  if (Equation_Type==0){      /* Schrodinger equation */
    Hamming_O(0,L,ep,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,0.95);
  }
  else if (Equation_Type==1){ /* scalar relativistic equation */
    Hamming_O(2,L,ep,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,0.95);
  }
  else if (Equation_Type==2){ /* fully relativistic equation */
    Hamming_O(3,L,ep,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,0.95);
  }

  /* find ic corresponding to rc */

  po = 0;
  for (i=0; i<Grid_Num; i++){
    if (rc<MRV[i] && po==0){
      po = 1;
      ic = i;
    }
  }

  /* calculate PW */

  for (i=0; i<=(ic+10); i++){
    PW[i] = pow(MRV[i],(double)L+1.0)*Lo[i];
  }

  /* normalization of PW */

  tmp = 1.0/sqrt(Norm_Core(PW,rc));

  if (0.0<PW[ic]) fugo = 1.0;
  else            fugo =-1.0;  

  for (i=0; i<=(ic+10); i++){
    PW[i] = tmp*fugo*PW[i];
  }

}







double Norm_Core(double af[ASIZE1], double rc)
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

  c4 = -c2*c2/(2.0*(double)l + 5.0);

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


void Spline_Func(double R, double fv[3], double af[ASIZE1])
{

  static int mp_min,mp_max,m;
  static double h1,h2,h3,f1,f2,f3,f4;
  static double g1,g2,x1,x2,y1,y2,f,Df,D2f;
  static double result;

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


int Calc_VPS_PAO(int calc_type,
                 int so, int l, double rc, double ep, double scale_factor,
                 double AEW[ASIZE1], double AEV[ASIZE1],
                 double  PW[ASIZE1], double  PV[ASIZE1])
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
        Input data from all electron calculation
  ****************************************************/

  Spline_Func(rc,fv,AEW);
  p0 = fv[0];
  p1 = fv[1];
  Spline_Func(rc,fv,AEV);
  v0 = fv[0];
  v1 = fv[1];
  v2 = fv[2];
  Ncore = scale_factor*Norm_Core(AEW,rc);

  b[0] = log(p0/pow(rc,(double)l+1.0));
  b[1] = p1/p0 - ((double)l + 1.0)/rc;
  b[2] = 2.0*v0 - 2.0*ep - 2.0*b[1]*((double)l + 1.0)/rc - b[1]*b[1];
  b[3] = 2.0*v1 + 2.0*b[1]*((double)l + 1.0)/rc/rc
        -2.0*b[2]*((double)l + 1.0)/rc - 2.0*b[1]*b[2];
  b[4] = 2.0*v2 - 4.0*b[1]*((double)l + 1.0)/rc/rc/rc
         + 4.0*b[2]*((double)l + 1.0)/rc/rc
         - 2.0*b[3]*((double)l + 1.0)/rc - 2.0*b[2]*b[2] - 2.0*b[1]*b[3];

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

  /****************************************************
           Search the upper and lower points
              of the bisection method
  ****************************************************/

  po = 0;
  c2 = -20.0;
  dc = 0.03;
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

    loop_num++;      

    if (loop_max<loop_num){

      if (calc_type==0 || calc_type==2){  
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
                  find c2 by bisection
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
    if (fabs(DLN)<=0.000000000001) po = 1;

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

      /* calculate the TM pseudofunction */

      PW[i] = pow(r,(double)l+1.0)*exp(p);

      /* make the TM pseudopotential */

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
  for (i=0; i<n; i++){
    r = 0.50*(Dr*x[i] + Sr);
    Spline_Func(r,fv0,phi0);
    Spline_Func(r,fv1,phi1);
    sum = sum + 0.5*Dr*fv0[0]*fv1[0]*w[i];
  }

  return sum;
}


void Calc_new_range(/*input */
		    int *node_c_ene, double *DD_c_ene, int num_nodes, 
		    /*output*/
		    double *MinE, double *MaxE , double *dE, int *num_div )
     /* find min and max of node_p==num_nodes */
{
  int num_min, num_max, found;
  found=0;
  for (num_min=0;num_min<*num_div;num_min++) {
    if ( node_c_ene[num_min]>= num_nodes ) {
      found=1;
      break;
    }
  }
  if (found==0) {
    num_min=0;
  }
  found=0;
  for (num_max=*num_div-1; num_max>=0; num_max--) {
    if ( node_c_ene[num_max] <= num_nodes ) {
      found=1;
      break;
    }
  }
  if (found==0) {
    num_max=*num_div-1;
  }
  /* order num_min < num_max  */
  if (num_min>num_max) {
    int i0=num_max;
    num_max=num_min;
    num_min = i0;
  }
  
  assert(num_max < *num_div);
#if 1
  if (num_max-num_min >= 5 ) {
     int num;
     int found=0;
     int num_smallest;
     double DD_c_val = 1.0e+10;
     /* find the smallest one */
     for (num=num_min;num<= num_max; num++) {
         if (node_c_ene[num]==num_nodes) {
            if (DD_c_val > DD_c_ene[num] ) {
                DD_c_val = DD_c_ene[num]; 
                num_smallest = num;
                found=1;
            }
         }
     }
     /* discrad found, DD_c_val */
     /* use num_smallest */
     if (found) {
        int num_try;
        /* try larger */ 
        num_try=num_smallest+1;
	if (num_try<*num_div)
        if (  node_c_ene[num_smallest-1] == node_c_ene[num_smallest] && 
	      node_c_ene[num_try]!=num_nodes ) {
            num_min = num_smallest-1;
#if 0
           printf("new num_min=%i\n",num_min);
#endif
            goto narrower_range_end;
        } 
        /* try smaller */
        num_try=num_smallest-1;
	if (num_try>=0) 
        if (  node_c_ene[num_smallest] == node_c_ene[num_smallest+1] &&
	      node_c_ene[num_try]!=num_nodes ) {
            num_max = num_smallest+1;
#if 0
            printf("new num_max=%i\n",num_max);
#endif
            goto narrower_range_end;
        }
     }
  narrower_range_end: ;
  }
#endif

  if (num_min>0) { num_min--; }
  if (num_max<*num_div-1) { num_max++;}


  {
    double saveMinE=*MinE;
    *MinE = saveMinE + *dE*num_min;
    *MaxE = saveMinE + *dE*num_max;
  }

  *num_div = 2*(num_max-num_min) + 1;

  *dE = (*MaxE - *MinE)/ ((*num_div) -1 );
/*  *dE ~= (*dE)*0.5; */
#if 0
  printf("num_min,num_max, num_div=%i %i %i\n",num_min,num_max,*num_div);
#endif
} 


void Calc_PAO(int state_num,
              int so, int n, int l, int num_nodes, double Reduce_Num_grid,
              double Vps[ASIZE1], double pw[ASIZE1])
{
  static int num,i,po,p,q,num_div,j,CTP_P;
  static int node_c,node_p,ref_switch,loop_num,loop_max;
  static int find_switch,match_p,match_p0;
  static int floop_num,floop_max;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double MaxE,MinE,trial_ene,DD_c,DD_p,kappa;
  static double Deriv_O,Deriv_I,dE,ratio,tmp0;
  static double DD_MinE,DD_MaxE,Local_Ratio_MatP;
  int num_div_loop;
  double MinE0, MaxE0;
  int *node_c_ene;
  double *DD_c_ene; 
  const double eps=1.0e-10;

  /**************************************************
    search roughly
  ***************************************************/

/*  Local_Ratio_MatP = 0.60; */

  floop_max = 5;

  MinE0 = E[state_num][so][n][l] - 0.4;

  if (Calc_Type!=2){
       MaxE0 = -1.0e-4;
  }
  else{
       MaxE0 = smallest(height_wall,5.0);
  }
  CTP_P=10;

  for (floop_num=0; floop_num< floop_max; floop_num++ ) { /* loop to change Local_Ratio_MatP */

    /* initial values for num_div loop */
    MinE = MinE0;
    MaxE = MaxE0;
    Local_Ratio_MatP = 0.60 + 0.1*floop_num ;

    if (debug) printf("Local_Ratio_MatP=%lf\n",Local_Ratio_MatP);

    num_div =10;
    dE = (MaxE - MinE)/(double)(num_div-1);

    for (num_div_loop=0;num_div_loop<20; num_div_loop++) {
	    /* make dE smaller and smaller */

      DD_c_ene = (double*)malloc(sizeof(double)*(num_div));
      node_c_ene = (int*)malloc(sizeof(int)*(num_div));

#if 0
      printf("min,max,dE=%lf %lf %lf\n",MinE, MaxE, dE);
#endif

      find_switch = 0;
      DD_p = 0.0;
      DD_c = 0.0;

      trial_ene =MinE;
      num=0;
      while ( trial_ene<= MaxE + eps ) { 


	Hamming_O(0,l,trial_ene,kappa,Mo,Lo,DMo,Vps,Vps,Reduce_Num_grid);
	Hamming_I(0,l,trial_ene,kappa,Mi,Li,DMi,Vps,Vps,Reduce_Num_grid);

	/* Matching point */
	if (LL_grid<=UL_grid){
	  if (LL_grid<=(CTP+CTP_P) && (CTP+CTP_P)<=UL_grid){
	    match_p = Grid_Num*Local_Ratio_MatP*Reduce_Num_grid;
	  }
	  else{
	    match_p = (LL_grid + UL_grid)/2;
	  }
	}
	else{
	  printf("Could not find an appropriate matching point in TM.c\n");
	  exit(0);
	}

	/* Deviation between the derivatives for MinE */

	DD_p = DD_c;
	ratio = Lo[match_p]/Li[match_p];
	Deriv_O = Mo[match_p];
	Deriv_I = ratio*Mi[match_p];
	DD_c = Deriv_O - Deriv_I;

	/* recalc # of nodes */

	for (i=0; i<=match_p; i++)           pw[i] = Lo[i];
	for (i=(match_p+1); i<Grid_Num; i++) pw[i] = ratio*Li[i];
	node_p = node_c;
	node_c = 0;
	for (i=1; i<(Grid_Num-3); i++){
	  if (pw[i]*pw[i+1]<0.0) node_c++;
	}

	/* save it */
	assert( 0<=num );
	assert(num<num_div );
	node_c_ene[num] = node_c;
	DD_c_ene[num] = DD_c;

        if (debug){
  	  printf("Rough l=%d trial_ene=%15.12f node_c=%d match_p=%4d DD_c=%15.12f\n",
	         l,trial_ene,node_c,match_p,DD_c);
	}

	if ( (node_p==num_nodes && node_c==num_nodes) && DD_p*DD_c<0.0 ){

	  /*
	    printf("AA find_switch=%3d %10.6f %10.6f %10.6f\n",
	    find_switch,trial_ene,DD_p,DD_c);
	  */

	  find_switch = 1;
	  MaxE = trial_ene ; 
	  MinE = trial_ene - dE;
	  DD_MinE = DD_p;
	  DD_MaxE = DD_c;
	  match_p0 = match_p;
          
          free(DD_c_ene); free(node_c_ene); 

	  goto start_refine;
	}  

	trial_ene += dE;
	num++;
      } /* trial ene */

      /* new trial_ene */
      /* find min and max of node_p==num_nodes */
      Calc_new_range(/*input */   node_c_ene, DD_c_ene, num_nodes,
       /*output*/ &MinE,  &MaxE ,  &dE, &num_div );

      free(DD_c_ene); free(node_c_ene);

    } /* num_div_loop */

   /* failed to to find wf */

  } /* floop_num */

  {
    printf("Error: failed to find node.\n"
	   "What you can do is\n"
	   "1. tune Local_Ratio_MatP.\n"
	   "   increase max. of num_div.\n"
	   "2. report it to authors.\n");
    
    exit(0);
  }    

  /**************************************************
    refine 
  ***************************************************/

 start_refine:

  if (find_switch==1){

    po = 0;
    DD_p = 0.0;
    node_c = 100;
    loop_num = 0; 
    loop_max = 100; 

    do{

      loop_num++;

      trial_ene = 0.5*(MaxE + MinE);
      Hamming_O(0,l,trial_ene,kappa,Mo,Lo,DMo,Vps,Vps,Reduce_Num_grid);
      Hamming_I(0,l,trial_ene,kappa,Mi,Li,DMi,Vps,Vps,Reduce_Num_grid);
      
      /* Deviation between the derivatives for MinE */

      DD_p = DD_c;
      ratio = Lo[match_p0]/Li[match_p0];
      Deriv_O = Mo[match_p0];
      Deriv_I = ratio*Mi[match_p0];
      DD_c = Deriv_O - Deriv_I;

      /* recalc # of nodes */

      for (i=0; i<=match_p0; i++)           pw[i] = Lo[i];
      for (i=(match_p0+1); i<Grid_Num; i++) pw[i] = ratio*Li[i];
      node_p = node_c;
      node_c = 0;
      for (i=0; i<(Grid_Num-1); i++){
        if (pw[i]*pw[i+1]<0.0) node_c++;
      }

      /* bisection */

      if (sgn(DD_c)==sgn(DD_MaxE)){
        MaxE = trial_ene; 
        DD_MaxE = DD_c;
      }
      else {
        MinE = trial_ene; 
        DD_MinE = DD_c;
      }

      /* convergence? */
      if (fabs(DD_c)<1.0e-12) po = 1; 

      /*
      printf("Refine l=%d trial_ene=%15.12f node_c=%d match_p0=%4d DD_c=%15.12f\n",
              l,trial_ene,node_c,match_p0,DD_c);
      */

      if ( po==0 && (loop_num==loop_max || (MaxE-MinE) < 1.0e-15) ) { 
        printf("Warning: not fully convergence in TM.c\n");
        printf("Please compare criterion and fabs(DD_c)\n");
        printf("criterion  = %18.15f\n",1.0e-12);   
        printf("fabs(DD_c) = %18.15f\n",fabs(DD_c));   
        printf("Acceptable or not??\n");
        po = 1;   
      }

    } while(po==0 && loop_num<loop_max);
  }
  else{
    printf("Could not find the eigenstate of VPS\n");
    exit(0);
  }


  if (po==1){
    for (i=0; i<=match_p0; i++){
      pw[i] = pow(MRV[i],(double)l+1.0)*Lo[i];
    }
    for (i=(match_p0+1); i<Grid_Num; i++){
      pw[i] = pow(MRV[i],(double)l+1.0)*ratio*Li[i];
    }
    tmp0 = 1.0/sqrt(ip_phi0_phi1(pw,pw,MRV[Grid_Num-1]));
    for (i=0; i<Grid_Num; i++){
      pw[i] = tmp0*pw[i];

      /*
      printf("%2d %2d %2d %15.12f %15.12f\n",so,n,l,MRV[i],pw[i]);
      */
    }
  }
  else{
    printf("Could not find the second orbital for L=%d in VPS\n",l);
    exit(0);
  }
}




static void Spherical_Bessel( double x, int Lmax, double *sb, double *dsb ) 
{
  int m,n,nmax;
  double *tsb;
  double j0,j1,j0p,j1p,sf,tmp,si,co,ix,ix2;

  if (x<0.0){
    printf("minus x is invalid for Spherical_Bessel\n");
    exit(0);    
  } 

  /* find an appropriate nmax */

  nmax = Lmax + 3*x + 20;
  if (nmax<100) nmax = 100; 

  /* allocate tsb */

  tsb  = (double*)malloc(sizeof(double)*(nmax+1)); 
  
  /* if x is larger than xmin */

  if ( xmin < x ){

    /* initial values*/

    tsb[nmax]   = 0.0;
    tsb[nmax-1] = 1.0e-14;

    /* downward recurrence from nmax-2 to Lmax+2 */

    for ( n=nmax-1; (Lmax+2)<n; n-- ){

      tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      if (1.0e+250<tsb[n-1]){
        tmp = tsb[n-1];        
        tsb[n-1] /= tmp;
        tsb[n  ] /= tmp;
      }

      /*
      printf("n=%3d tsb[n-1]=%18.15f\n",n,tsb[n-1]);
      */
    }

    /* downward recurrence from Lmax+1 to 0 */

    n = Lmax + 3;
    tmp = tsb[n-1];        
    tsb[n-1] /= tmp;
    tsb[n  ] /= tmp;

    for ( n=Lmax+2; 0<n; n-- ){
      tsb[n-1] = (2.0*n + 1.0)/x*tsb[n] - tsb[n+1];

      if (1.0e+250<tsb[n-1]){
        tmp = tsb[n-1];
        for (m=n-1; m<=Lmax+1; m++){
          tsb[m] /= tmp;
        }
      }
    }

    /* normalization */

    si = sin(x);
    co = cos(x);
    ix = 1.0/x;
    ix2 = ix*ix;
    j0 = si*ix;
    j1 = si*ix*ix - co*ix;

    if (fabs(tsb[1])<fabs(tsb[0])) sf = j0/tsb[0];
    else                           sf = j1/tsb[1];

    /* tsb to sb */

    for ( n=0; n<=Lmax+1; n++ ){
      sb[n] = tsb[n]*sf;
    }

    /* derivative of sb */

    dsb[0] = co*ix - si*ix*ix;
    for ( n=1; n<=Lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }

  } 

  /* if x is smaller than xmin */

  else {

    /* sb */

    for ( n=0; n<=Lmax; n++ ){
      sb[n] = 0.0;
    }
    sb[0] = 1.0;

    /* derivative of sb */

    dsb[0] = 0.0;
    for ( n=1; n<=Lmax; n++ ){
      dsb[n] = ( (double)n*sb[n-1] - (double)(n+1.0)*sb[n+1] )/(2.0*(double)n + 1.0);
    }
  }

  /* free tsb */

  free(tsb);
}


static void find_zero_jx()
{
  int num,l,n;
  double *sb,*dsb;
  double x,dx,psb,x0,x1,y0,y1;
  double xt0[numsb],xt1[numsb];

  /* allocation of arrays */

  sb  = (double*)malloc(sizeof(double)*(lmax+3));
  dsb = (double*)malloc(sizeof(double)*(lmax+3));

  /* main calculation */

  for (l=0; l<=lmax; l++){

    /* rough search */

    x = 1.0e-10;
    dx = 1.0e-2;
    num = 0;

    Spherical_Bessel(x,lmax,sb,dsb);

    do {

      x += dx;

      psb = sb[l];
      Spherical_Bessel(x,lmax,sb,dsb);

      if (sb[l]*psb<0.0){
        xt0[num] = x - dx;
        xt1[num] = x;
        num++; 
      }

    } while(num<numsb);

    /* refinement of zero by bisection method */

    for (n=0; n<numsb; n++){

      x0 = xt0[n];
      Spherical_Bessel(x0,lmax,sb,dsb);
      y0 = sb[l]; 
 
      x1 = xt1[n];
      Spherical_Bessel(x1,lmax,sb,dsb);
      y1 = sb[l]; 

      do {

        x = 0.5*(x0+x1);
        Spherical_Bessel(x,lmax,sb,dsb);

        if ( 0.0<=(y0*sb[l]) ){
          x0 = x;
          y0 = sb[l];
        }
        else{
          x1 = x;
          y1 = sb[l];
        }

      } while( 1.0e-15<fabs(sb[l]) ); 

      zero_j[l][n] = x;
    }
  }

  /* freeing of arrays */

  free(sb);
  free(dsb);
}


static void Eigen_lapack_x(double **a, double *ko, int n0, int EVmax)
{

  /*
    F77_NAME(dsyevx,DSYEVX)()
  
    input:  n;
    input:  a[n][n];  matrix A
    output: a[n][n];  eigevectors
    output: ko[n];    eigenvalues 
  */
    
  char *name="Eigen_lapack";

  char  *JOBZ="V";
  char  *RANGE="I";
  char  *UPLO="L";

  INTEGER n=n0;
  INTEGER LDA=n0;
  double VL,VU; /* dummy */
  INTEGER IL,IU; 
  double ABSTOL=1.0e-13;
  INTEGER M;

  double *A,*Z;
  INTEGER LDZ=n;
  INTEGER LWORK;
  double *WORK;
  INTEGER *IWORK;
  INTEGER *IFAIL, INFO;

  int i,j;

  A=(double*)malloc(sizeof(double)*n*n);
  Z=(double*)malloc(sizeof(double)*n*n);

  LWORK=n*8;
  WORK=(double*)malloc(sizeof(double)*LWORK);
  IWORK=(INTEGER*)malloc(sizeof(INTEGER)*n*5);
  IFAIL=(INTEGER*)malloc(sizeof(INTEGER)*n);

  IL = 1;
  IU = EVmax;
 
  for (i=0; i<n; i++) {
    for (j=0; j<n; j++) {
       A[i*n+j] = a[i][j];
    }
  }

#if 0
  printf("A=\n");
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) {
       printf("%f ",A[i*n+j]);
    }
    printf("\n");
  }
  fflush(stdout);
#endif

  dsyevx_( JOBZ, RANGE, UPLO, &n, A, &LDA, &VL, &VU, &IL, &IU,
           &ABSTOL, &M, ko, Z, &LDZ, WORK, &LWORK, IWORK,
           IFAIL, &INFO ); 

  /* store eigenvectors */
  for (i=0;i<EVmax;i++) {
    for (j=0;j<n;j++) {
      /*  a[i][j]= Z[i*n+j]; */
      a[j][i]= Z[i*n+j];
    }
  }

  if (INFO>0) {
    printf("\n%s: error in dsyevx_, info=%d\n\n",name,INFO);
  }
  if (INFO<0) {
     printf("%s: info=%d\n",name,INFO);
     exit(10);
  }
   
  free(IFAIL); free(IWORK); free(WORK); free(Z); free(A);
}


static double phi0_phi1(double phi0[ASIZE1],
			double phi1[ASIZE1])
{
  int i,n,j,l,nf,fg;
  double r,rmin,rmax,Sr,Dr,sum,dum;

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  rmax = MRV[Grid_Num-1];
  Sr = rmax + rmin;
  Dr = rmax - rmin;
  
  sum = 0.0;
  for (i=0; i<=(n-1); i++){
    r = 0.50*(Dr*x[i] + Sr);
    sum += 0.5*Dr*w[i]*Int_RadF(r,phi0)*Int_RadF(r,phi1);
  }

  return sum;
}

static double phi0_v_phi1(double phi0[ASIZE1],
			  double phi1[ASIZE1],
			  double vpot[ASIZE1])
{
  int i,n,j,l,nf,fg;
  double r,rmin,rmax,Sr,Dr,sum,dum;

  n = ASIZE7-2;
  nf = n;
  fg = 1;

  Gauss_Legendre(n,x,w,&nf,&fg);

  rmin = MRV[0];
  rmax = MRV[Grid_Num-1];
  Sr = rmax + rmin;
  Dr = rmax - rmin;
  
  sum = 0.0;
  for (i=0; i<=(n-1); i++){
    r = 0.50*(Dr*x[i] + Sr);
    sum += 0.5*Dr*w[i]*Int_RadF(r,phi0)*Int_RadF(r,vpot)*Int_RadF(r,phi1);
  }

  return sum;
}


static double Int_RadF(double R, double RadF[ASIZE1])
{
  int mp_min,mp_max,m,po;
  double h1,h2,h3,f1,f2,f3,f4;
  double g1,g2,x1,x2,y1,y2,f;
  double rm,y12,y22,df,a,b;
  double result;

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
