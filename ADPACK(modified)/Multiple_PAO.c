/**********************************************************************
  Multiple_PAO.c:

     Multiple_PAO.c is a subroutine to calculate eigenstates
     of atom with confinment pseudo potentials.

  Log of Multiple_PAO.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

static void VP_MPAO(int GL);
static void Simpson(double PAO[ASIZE1]);


static void search_roughly(
        int Num_Part,  /* input, array size of search_?[] */
        double search_e[], int search_node[], double search_D[], int search_MatP[],/* output */
        double min, double max,  int MatP_min,int MatP_max,  /* input */
	int GL, int CTP_P, size_t Grid_Num, int Smallest_MP,double alpha, /* input */
        int asize1, double Mo[], double Lo[], double DMo[], /* work */
        double Mi[],double Li[], double DMi[],
        double Ltmp[] );

static int search_bisection(
           int NumL, double Min_ep, double Max_ep, double Min_D, double Max_D,
           int Min_MatP, int Max_MatP,
           double Criterion, /* input */
           int retnode, double Trial_ep, double Trial_D, /* output */
	   int GL, int CTP_P, size_t Grid_Num, int Smallest_MP,double alpha, /* input */
	   int Num_Part, double Mo[], double Lo[], double DMo[], /* work */
	   double Mi[],  double Li[], double DMi[],
	   double Ltmp[]);





void Multiple_PAO(int state_num)
{
  static int i0,i1,i2,imin,p_node,NumL,s,mu,imu;
  static int i,j,k,n,l,po,SCF,po_node,po1,po2,CTP_P;
  static int cNode,pNode,Node_min,Node_max,refine_loop;
  static int Enum,GL,M,q,MatP,MatP0,MatP1;
  static int nf,fg,Mul,num_k,Smallest_MP;

  static double Sr,Dr,norm_kmax,norm_kmin,norm_k,rmin,rmax,r;
  static double minMapD,tmp0,tmp1,tmp2,cDD,pDD,pDD2,pOEE,Fugou;
  static double DD_max,DD_min,di,alpha;
  static double Mo[ASIZE1],Lo[ASIZE1],DMo[ASIZE1];
  static double Mi[ASIZE1],Li[ASIZE1],DMi[ASIZE1];
  static double Ltmp[ASIZE1];
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

  int    inode,ie,ie_min,ie_max,found,really_found,ith,count;
  int    MatP_min,MatP_max;

  search_min0 = search_LowerE;
  search_max0 = search_UpperE;

  /****************************************************
      solve the atomic Kohn-Sham equation with 
      a confinement potential by the Hamming method
  ****************************************************/

  Smallest_MP = Grid_Num/3;
  Criterion = 1.0e-7;

  CTP_P = 2;
  n = NVPS[0];  
  l = LVPS[0];

  MinE = E[state_num][0][n][l] - 0.50; 
  alpha = 0.1;

  for (GL=0; GL<=MaxL_PAO; GL++){
    
    printf("<PAO>  Calculating multiple pseudo atomic orbitals for GL=%2d....\n",GL);

    VP_MPAO(GL);

    for (inode=0;inode<Num_PAO;inode++) {

      search_min = MinE;
      search_max = search_max0;

      search_roughly(
                   Num_Partition, search_e,search_node,search_D,search_MatP,
                   search_min,search_max, fix_MatP[GL][inode],  fix_MatP[GL][inode], 
                   GL,CTP_P,Grid_Num,Smallest_MP,alpha,
		   ASIZE1, Mo,Lo,DMo,
		   Mi,Li,DMi,
		   Ltmp);

      if (Check_switch==1) {
        printf("inode=%d\n",inode);
      }

      for (count=0; ; count++) {

				if (count==0) {
					for (i=0;i<Num_Partition;i++) {
            search_e2[i]=search_e[i];
            search_node2[i]=search_node[i];
            search_D2[i]=search_D[i];
            search_MatP2[i]=search_MatP[i];
					}

				}
				else {

					search_roughly(
                         Num_Partition, search_e2,search_node2,search_D2, search_MatP2,
                         search_min,search_max, fix_MatP[GL][inode],  fix_MatP[GL][inode],
                         GL,CTP_P,Grid_Num,Smallest_MP,alpha,
												 ASIZE1, Mo,Lo,DMo,
												 Mi,Li,DMi,
												 Ltmp);
				}

        ie_min=-100;
        ie_max=-100;
				for (ie=0;ie<Num_Partition;ie++) {
					if (inode==search_node2[ie]) {
						ie_min=ie;
						break;
					}
				}

				for (ie=0;ie<Num_Partition;ie++) {
					if (inode==search_node2[ie]) {
						ie_max=ie;
					}
				}

        if (Check_switch==1){ printf("node=%d, ie=%d %d\n",inode,ie_min,ie_max); }

        if (ie_min==-100 &&
           ! (search_node2[0]<inode && search_node2[Num_Partition-1]>inode)) {
               search_min=search_min-(search_max-search_min)*0.5;
               search_max=search_max-(search_max-search_min)*0.5;
               continue;
           
        }

	/* argorighm 1, search inversion of sign in search_D */
	found=0;
	for (ie=ie_min;ie<ie_max;ie++) {
	  if (search_D2[ie]*search_D2[ie+1]<0.0 && search_node2[ie]==search_node2[ie+1] &&
	      search_node2[ie]==inode ) {
	    search_invD[found++]=ie;
	  }
	}
	if (found) {
	  for (ith=0;ith<found;ith++) {
	    ie=search_invD[ith]; 
            if (Check_switch==1) {
	    printf("inv found ie=%d\n",ie);
            }
#if 0
            printf("%f %f %f %f\n",search_e2[ie],search_e2[ie+1], search_D2[ie],search_D2[ie+1]);
            really_found=1;
#else
	    really_found=search_bisection(
                                          inode, search_e2[ie],search_e2[ie+1], search_D2[ie],search_D2[ie+1],
                                          search_MatP2[ie],search_MatP2[ie+1],
                                          Criterion,
                                          node,Trial_ep, Trial_D,
                                          GL,CTP_P,Grid_Num,Smallest_MP,alpha,
					  ASIZE1, Mo,Lo,DMo,
					  Mi,Li,DMi,
					  Ltmp);
#endif
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
  printf("\n");

} 



void Simpson(double PAO[ASIZE1])
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

  if (0.0<=PAO[10]){
    for (i=0; i<=(Grid_Num-1); i++){
      PAO[i] = PAO[i]/sqr_sum/MRV[i];
    }
  }
  else {
    for (i=0; i<=(Grid_Num-1); i++){
      PAO[i] = -PAO[i]/sqr_sum/MRV[i];
    }
  }
}





void VP_MPAO(int GL)
{
  static int m,i,j;
  static double sum,Vl,tmp0;

  if (Equation_Type==0 || Equation_Type==1){

    if (0<NumVPS_L[GL]){
      m = GI_VPS[GL][0];
      for (i=0; i<Grid_Num; i++){
	V[i] = VPS[0][0][m][i] + Vh_V[i] + Vxc_V[i];
      }
    }
    else{
      for (i=0; i<Grid_Num; i++){
	sum = 0.0;
	for (j=0; j<Number_VPS; j++){
	  sum = sum + VPS[0][0][j][i];
	}
	tmp0 = (double)Number_VPS;
	Vl = sum/tmp0;
	V[i] = Vl + Vh_V[i] + Vxc_V[i];
      }
    }

  }

  else if (Equation_Type==2){

    if (0<NumVPS_L[GL]){
      m = GI_VPS[GL][0];
      for (i=0; i<Grid_Num; i++){
	V[i] = 0.5*(VPS[0][0][m][i] + VPS[0][1][m][i]) + Vh_V[i] + Vxc_V[i];
      }
    }
    else{
      for (i=0; i<Grid_Num; i++){
	sum = 0.0;
	for (j=0; j<Number_VPS; j++){
	  sum = sum + 0.5*(VPS[0][0][j][i] + VPS[0][1][j][i]);
	}
	tmp0 = (double)Number_VPS;
	Vl = sum/tmp0;
	V[i] = Vl + Vh_V[i] + Vxc_V[i];
      }
    }
  }

}




static void search_roughly(
        int Num_Part,  /* input, array size of search_?[] */
        double search_e[], int search_node[], double search_D[], int search_MatP[],/* output */
        double min, double max,  int MatP_min,int MatP_max,  /* input */
	int GL, int CTP_P, size_t Grid_Num, int Smallest_MP,double alpha, /* input */
        int asize1, double Mo[], double Lo[], double DMo[], /* work */
        double Mi[],double Li[], double DMi[],
        double Ltmp[] )
{
  static int MatP,MatP0,MatP1,i,ie;
  double Trial_ep,tmp0,Trial_D,ratio,kappa;

  /* rough search */
  for ( ie=0;ie<Num_Part; ie++) {

    Trial_ep=min+(max-min)*ie/(Num_Part-1);
    Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,V,V,RNG[NVPS[0]][LVPS[0]]);
    Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,V,V,RNG[NVPS[0]][LVPS[0]]);

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
    for (i=1; i<=(Grid_Num-5); i++){
      if (Ltmp[i]*Ltmp[i+1]<0.0) node++;
    }

    if (Check_switch==1){
      printf("X MPAO GL=%i Trial_ep=%15.12f , ie=%d",GL,Trial_ep,ie);
      printf(" Trial_D=%15.12f  node=%2d  MatP=%2d CTP=%d\n",
               Trial_D,node,MatP,CTP);
    }

    search_e[ie]=Trial_ep;
    search_node[ie]=node;
    search_D[ie]=Trial_D;
    search_MatP[ie]=MatP;
    
  }
  /* rough search, end */
}



static int search_bisection(
           int NumL, double Min_ep, double Max_ep, double Min_D, double Max_D,
           int Min_MatP, int Max_MatP,
           double Criterion, /* input */
           int retnode, double Trial_ep, double Trial_D, /* output */
	   int GL, int CTP_P, size_t Grid_Num, int Smallest_MP,double alpha, /* input */
	   int Num_Part, double Mo[], double Lo[], double DMo[], /* work */
	   double Mi[],  double Li[], double DMi[],
	   double Ltmp[])
{
  int MatP,MatP0,MatP1,i,  po,po2, node,i1; 
  double  pTrial_D, ratio,ep,tmp0,kappa;

  if (Check_switch==1) {
  printf("bisection start, node=%d,  e=%f %f, D=%f %f MatP=%d %d\n",
     NumL,Min_ep,Max_ep,Min_D,Max_D,Min_MatP,Max_MatP);
  }

  po=0;
  i1=0;
  Trial_D=0;
#if 0
  /* check consistency */
    Trial_ep=Min_ep;
    Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,V,V,RNG[NVPS[0]][LVPS[0]]);
    Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,V,V,RNG[NVPS[0]][LVPS[0]]);
    ratio = Lo[MatP]/Li[MatP];
    Trial_D = Mo[MatP] - ratio*Mi[MatP];
    Trial_ep=Max_ep;
    Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,V,V,RNG[NVPS[0]][LVPS[0]]);
    Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,V,V,RNG[NVPS[0]][LVPS[0]]);
    ratio = Lo[MatP]/Li[MatP];
    pTrial_D = Mo[MatP] - ratio*Mi[MatP];
    if (Trial_D*pTrial_D<0.0) {
       if (Check_switch==1) {
       printf("search_bisection: Trial_D= %e %e at MatP=%d\n",Trial_D,pTrial_D,MatP);
       }
    }
    else {
       printf("search_bisection: Trial_D*pTrial_D>0, stop\n");
       exit(10);
    }
#endif

  do {
    
    Trial_ep = 0.50*(Min_ep + Max_ep);
    Hamming_O(0,GL,Trial_ep,kappa,Mo,Lo,DMo,V,V,RNG[NVPS[0]][LVPS[0]]);

    i1++;
    Hamming_I(0,GL,Trial_ep,kappa,Mi,Li,DMi,V,V,RNG[NVPS[0]][LVPS[0]]);

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
    /* Recalculation of the number of nodes */
    for (i=0; i<=MatP; i++){
      Ltmp[i] = Lo[i];
    }
    for (i=(MatP+1); i<=(Grid_Num-1); i++){
      Ltmp[i] = ratio*Li[i];
    }
    node = 0;
    for (i=1; i<=(Grid_Num-5); i++){
      if (Ltmp[i]*Ltmp[i+1]<0.0) node++;
    }

    if (Check_switch==1){
    printf("Y MPAO GL=%d NumL=%d Trial_ep=%e Trial_D=%e node=%d MatP=%d CTP=%d\n",
           GL,NumL,Trial_ep,Trial_D,node,MatP,CTP);
    }

#if 0
    if (fabs(pTrial_D-Trial_D)<10e-14 && Criterion<fabs(Trial_D)){
      /**************************************************
       Failure of the eigenvalue by the bisection method
       In this case, this routine tries to search the
       true eigenvalue for the large eigevalue regime.
      **************************************************/

      po = 2;
    }
    else
#endif
    if (fabs(Trial_D)<=Criterion||Max_ep-Min_ep<1.0e-12){

      if (node==NumL){
	po = 1;
	po2 = 0;
            
	ep = Trial_ep;
	E_PAO[GL][NumL] = ep;
           
	for (i=0; i<=MatP; i++){
	  MPAO[GL][NumL][i] = pow(MRV[i],(double)GL+1.0)*Lo[i];
	}
	for (i=(MatP+1); i<=(Grid_Num-1); i++){
	  MPAO[GL][NumL][i] = pow(MRV[i],(double)GL+1.0)*ratio*Li[i];
	}
            
	Simpson(MPAO[GL][NumL]);

	printf("<PAO>  GL=%2d  Multiplicity=%2d node=%2d  ",GL,NumL,node);
	printf("eigenvalue=%15.12f (Hartree)\n",ep);

        if (fabs(Trial_D)>Criterion) {
          printf("warning: fabs(Trial_D)>Criterion, Trial_D=%e ,Criterion= %e at MatP=%d\n\n",
            Trial_D,Criterion,MatP);
        }
 
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




  
