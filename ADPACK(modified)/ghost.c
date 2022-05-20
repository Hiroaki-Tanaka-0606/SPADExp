/**********************************************************************
  ghost.c:

     ghost.c is a subroutine to check whether ghost states exist or not
     in eigenstates of the separable pseudo potentials (KB and Blochl).

  Log of ghost.c:

     19/Mar/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

static double phi0_phi1(double phi0[ASIZE1], double phi1[ASIZE1], double rmax);
static double Int_RadF(double R, double RadF[ASIZE1]);

void ghost(int state_num)
{  
  static int m,L,loop,num_steps,po,i,nodes;
  static int p,SCF_po,j,N,gi,so;
  static int Grid_Num1,mat_p,num_states;
  static double max_ene,min_ene,dene,ene,mixing,kappa;
  static double tmp0,sum_c,sum_p,r,Reduce_Num_grid;
  static double Deriv_O,Deriv_I,ratio,DD,DD_p,dif_c,dif_p;
  static double Mo[ASIZE1],DMo[ASIZE1],Lo[ASIZE1];
  static double Mi[ASIZE1],DMi[ASIZE1],Li[ASIZE1];
  static double L0[ASIZE1],L0_p[ASIZE1];
  static double Vsl[ASIZE1],U0[ASIZE1],VNLp[ASIZE1];
  static double EigenE[ASIZE2][5];

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){

      if (NumVPS_L[L]!=0){

	printf("\n<VPS>  checking ghost states for so=%i L=%i\n",so,L);

	num_states = 0;         
	num_steps = 50;
	min_ene = -5.0;
	max_ene = -2.0e-2;
	dene = (max_ene - min_ene)/(double)(num_steps-1);
	DD = 0.0;

	gi = GI_VPS[L][0];
	N = NVPS[gi];
	Reduce_Num_grid = RNG[N][L];
	Grid_Num1 = Reduce_Num_grid*Grid_Num;
	mat_p = MatP[so][N][L];

	for (loop=0; loop<num_steps; loop++){
	  ene = min_ene + dene*(double)loop;

	  /****************************************
	        semi-local pseudo potentials
	  ****************************************/

	  m = GI_VPS[L][0];
	  for (i=0; i<Grid_Num1; i++){
	    Vsl[i] = VPS[0][so][m][i] + Vh_V[i] + Vxc_V[i];
	  }
	  Hamming_O(0,L,ene,kappa,Mo,Lo,DMo,Vsl,Vsl,Reduce_Num_grid);
	  Hamming_I(0,L,ene,kappa,Mi,Li,DMi,Vsl,Vsl,Reduce_Num_grid);

	  /* matching point */

	  if (LL_grid<=mat_p && mat_p<=UL_grid){
	    /* it's OK!! */
	  }
	  else{
	    mat_p = (LL_grid + UL_grid)/2;
	    printf("warning!, could not find an appropriate matching point\n");
	  }  

	  /* initial radial function */

	  ratio = Lo[mat_p]/Li[mat_p];
	  for (i=0; i<=mat_p; i++){
	    L0[i] = Lo[i];
	  }
	  for (i=(mat_p+1); i<Grid_Num1; i++){
	    L0[i] = ratio*Li[i];
	  }

	  /* normalization */
	  tmp0 = 1.0/sqrt(phi0_phi1(L0,L0,MRV[Grid_Num1-1]));
	  for (i=0; i<Grid_Num1; i++){
	    L0[i] = tmp0*L0[i];
	  }

	  /****************************************
	      fully separable pseudo potentials
	  ****************************************/

	  DD_p = DD;

	  /* screened Vlocal */
	  for (i=0; i<Grid_Num1; i++){
	    Vsl[i] = Vlocal[i] + Vh_V[i] + Vxc_V[i];
	  }

	  dif_c = 1.0e+5;
	  sum_c = 0.0;
	  mixing = 0.5;
	  SCF_po = 0;
	  j = 0;
	  /* iteration loop */ 
	  do{
	    j++;
	    /* calc Uo */
	    for (i=0; i<Grid_Num1; i++){
	      r = MRV[i];
	      U0[i] = pow(r,(double)L+1.0)*L0[i];
	    }
	    /* sum over projectors */
	    for (i=0; i<Grid_Num1; i++) VNLp[i] = 0.0;
	    sum_p = sum_c;
	    sum_c = 0.0; 
	    for (p=0; p<projector_num[L]; p++){
	      /* <v|phi_i> */

	      tmp0 = phi0_phi1(VNL_W2[so][L][p],U0,Vlocal_maxcut);
	      sum_c = sum_c + tmp0;

	      /* |v>*c*<v|phi_i> */
	      for (i=0; i<Grid_Num1; i++){
		r = MRV[i];
		VNLp[i] = VNLp[i] + VNL_W2[so][L][p][i]/r*proj_ene[0][so][L][p]*tmp0;
	      }
	    }

	    for (i=0; i<Grid_Num1; i++) L0_p[i] = L0[i];
	    Hamming_O(1,L,ene,kappa,Mo,Lo,DMo,Vsl,VNLp,Reduce_Num_grid);
	    Hamming_I(1,L,ene,kappa,Mi,Li,DMi,Vsl,VNLp,Reduce_Num_grid);

	    dif_p = dif_c;
	    dif_c = fabs(sum_p-sum_c);

	    if (dif_c<1.0e-5){
	      SCF_po = 1;
	    }
	    else{

	      /* matching point */

	      if (LL_grid<=mat_p && mat_p<=UL_grid){
		/* it's OK!! */
	      }
	      else{
		mat_p = (LL_grid + UL_grid)/2;
		printf("warning!, could not find an appropriate matching point\n");
	      }  

	      /* deviation between the derivatives */

	      ratio = Lo[mat_p]/Li[mat_p];
	      Deriv_O = Mo[mat_p];
	      Deriv_I = ratio*Mi[mat_p];
	      DD = Deriv_O - Deriv_I;

	      /* initial radial function */

	      for (i=0; i<=mat_p; i++){
		L0[i] = Lo[i];
	      }
	      for (i=(mat_p+1); i<Grid_Num1; i++){
		L0[i] = ratio*Li[i];
	      }

	      /* normalization */
	      tmp0 = 1.0/sqrt(phi0_phi1(L0,L0,MRV[Grid_Num1-1]));
	      for (i=0; i<Grid_Num1; i++){
		L0[i] = tmp0*L0[i];
	      }

	      /* recalc # of nodes */
	      nodes = 0;
	      for (i=1; i<(Grid_Num1-2); i++){
		if (L0[i]*L0[i+1]<0.0) nodes++;
	      }

	      /* simple mixing */

	      if (dif_c<dif_p)   mixing = 1.5*mixing;
	      else               mixing = 0.5*mixing;
 
	      if (mixing<1.0e-4) mixing = 1.0e-4;
	      if (0.9<mixing)    mixing = 0.9;

	      for (i=0; i<Grid_Num1; i++){
		L0[i] = mixing*L0[i] + (1.0 - mixing)*L0_p[i];
	      }
	    }

	    if (300<=j && SCF_po==0){
	      printf("warning!, no convergence (dif=%15.10f) in ghost.c\n",dif_c);
	    } 
	  } while(SCF_po==0 && j<300);

	  if (DD*DD_p<0.0){
	    EigenE[L][num_states] = ene;
	    num_states++;
	  }

	  printf("<VPS>  L=%i,  ene=%10.7f,  dif=%15.10f,  # of nodes=%2d\n",
		 L,ene,DD,nodes);
	}

	/* ghost or not ? */ 
	if (num_states==0){
	  printf("warning!, could not find any bound state for L=%i\n",L);
	}
	else if (num_states==NumVPS_L[L]){
	  printf("<VPS>\n");
	  printf("<VPS>  No ghost state for L=%i\n",L);
	  printf("<VPS>\n");
	}
	else if (num_states==(NumVPS_L[L]+1)){
	  gi = GI_VPS[L][0];
	  N = NVPS[gi];
	  if (EigenE[L][0]<(E[state_num][so][N][L]-0.3)){ 
	    printf("<VPS>\n");
	    printf("<VPS>  For L=%i, a ghost state is found UNDER correct eigenstates.\n",L);
	    printf("<VPS>\n");
	  }
	  else{ 
	    printf("<VPS>\n");
	    printf("<VPS>  For L=%i, a ghost like state is found ABOVE correct eigenstates.\n",L);
	    printf("<VPS>\n");
	  }
	}
	else{
	  gi = GI_VPS[L][0];
	  N = NVPS[gi];
	  if (EigenE[L][0]<(E[state_num][so][N][L]-0.3)){ 
	    printf("<VPS>\n");
	    printf("<VPS>  For L=%i, ghost states are found UNDER correct eigenstates.\n",L);
	    printf("<VPS>\n");
	  }
	  else{ 
	    printf("<VPS>\n");
	    printf("<VPS>  For L=%i, ghost like states are found ABOVE correct eigenstates.\n",L);
	    printf("<VPS>\n");
	  }
	}
      }
    }
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
    sum = sum
           + 
          0.5*Dr*w[i]*Int_RadF(r,phi0)*Int_RadF(r,phi1);
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
