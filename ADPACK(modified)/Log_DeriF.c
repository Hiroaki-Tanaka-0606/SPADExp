/**********************************************************************
  Log_DeriF.c:

     Log_DeriF.c is a subroutine to calculate the logarithmic
     derivatives of radial wave functions which are calculated 
     under the all electron unscreened potential, semi-local
     pseudo potentials, and separable (KB and Blochl) pseudo
     potentials.

  Log of Log_DeriF.c:

     19/Mar/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "adpack.h"

static double phi0_phi1(double phi0[ASIZE1], double phi1[ASIZE1], double rmax);
static double Int_RadF(double R, double RadF[ASIZE1]);

void Log_DeriF()
{  
  int m,L,loop,num_steps,po,i,i_r;
  int p,SCF_po,j,so,renzoku;
  double max_ene,min_ene,dene,ene;
  double Reduce_Num_grid,r,r_i,kappa,geta,sign_Lo;
  double logd0,logd1,logd2,dif_p,dif_c;
  double tmp0,sum_c,sum_p,dif0,dif1,mixing;
  double Mo[ASIZE1],DMo[ASIZE1];
  double Lo[ASIZE1],Lo_p[ASIZE1];
  double Vsl[ASIZE1],Uo[ASIZE1],VNLp[ASIZE1];

  geta = 1.0e-14; 
  num_steps = LogD_num;
  min_ene = LogD_MinE;
  max_ene = LogD_MaxE;
  dene = (max_ene - min_ene)/(double)num_steps;

  /************************************************
    calculate logarithmic derivatives of radial 
    wave functions 
  *************************************************/

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){

      if (NumVPS_L[L]!=0){

	printf("<VPS>  SO=%i L=%i, calculating logarithmic derivatives of wave functions\n",so,L);

        if (so==0) {
 	  printf("<VPS>  Used radius for LogD's=%10.7f (a.u.)\n",LogD_R[L]);
        } 
        else {
          printf("<VPS>  Used radius for LogD's=%10.7f (a.u.)\n",LogD1_R[L]);
        } 

	/************************************************
	 find Reduce_Num_grid which exceeds LogD_R
	*************************************************/

	i = -1;
	po = 0;

        if (so==0) {
  	  do{
	    i++;
	    if (LogD_R[L]<MRV[i]){
	      po = 1;
	      Reduce_Num_grid = (double)(i+5)/(double)Grid_Num;
	      i_r = i;
	    }
	  }while(po==0 && i<Grid_Num);
	}

        else {
          do{
            i++;
            if (LogD1_R[L]<MRV[i]){
              po = 1;
              Reduce_Num_grid = (double)(i+5)/(double)Grid_Num;
              i_r = i;
            }
          }while(po==0 && i<Grid_Num);
        }

	dif0 = 0.0;
	dif1 = 0.0;   
       
	for (loop=0; loop<num_steps; loop++){

	  ene = min_ene + dene*(double)loop;
	  r_i = MRV[i_r];

	  /****************************************
	   RICS note 563p
	   dlog(phi)/dr = (L+Mo/Lo)/r
	  ****************************************/

	  /****************************************
  	    the all electron unscreened potential
	  ****************************************/

          if (L==0){
            kappa = -1.0;    /* for j=l+1/2 */
          } 
          else{ 
            if      (so==0)  kappa = (double)(-L-1);    /* for j=l+1/2 */
            else if (so==1)  kappa = (double)L;         /* for j=l-1/2 */
	  }

	  if (Equation_Type==0){       /* Schrodinger equation */
	    Hamming_O(0,L,ene,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,Reduce_Num_grid);
	  }
	  else if (Equation_Type==1){  /* scalar relativistic equation */
	    Hamming_O(2,L,ene,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,Reduce_Num_grid);
	  }
	  else if (Equation_Type==2){  /* fully relativistic equation */
	    Hamming_O(3,L,ene,kappa,Mo,Lo,DMo,V_all_ele,V_all_ele,Reduce_Num_grid);
	  }

          sign_Lo = sgn(Lo[i_r]);

	  logd0 = ((double)L + Mo[i_r]/(Lo[i_r]+sign_Lo*geta))/r_i;

	  /****************************************
    	         semi-local pseudo potentials
	  ****************************************/

	  m = GI_VPS[L][0];
	  for (i=0; i<Grid_Num; i++){
	    Vsl[i] = VPS[0][so][m][i] + Vh_V[i] + Vxc_V[i];
	  }

	  Hamming_O(0,L,ene,kappa,Mo,Lo,DMo,Vsl,Vsl,Reduce_Num_grid);

	  sign_Lo = sgn(Lo[i_r]); 
	  logd1 = ((double)L + Mo[i_r]/(Lo[i_r]+sign_Lo*geta))/r_i;

	  /****************************************
	      fully separable pseudo potentials
	  ****************************************/

	  /* screened Vlocal */
	  for (i=0; i<Grid_Num; i++){
	    Vsl[i] = Vlocal[i] + Vh_V[i] + Vxc_V[i];
	  }

	  dif_c = 1.0e+5;
	  sum_c = 0.0;
	  mixing = 0.01;
	  SCF_po = 0;
	  j = 0;
          renzoku = 0;

	  /* iteration loop */

	  do {


	    j++;

	    /* calc Uo, its initlal value is calculated by Lo generated
               for the semi-local pseudopotential. */

	    for (i=0; i<Grid_Num; i++) Uo[i] = 0.0;
	    for (i=0; i<i_r; i++){
	      r = MRV[i];
	      Uo[i] = pow(r,(double)L+1.0)*Lo[i];
	    }

	    /* sum over projectors */

	    for (i=0; i<Grid_Num; i++) VNLp[i] = 0.0;

	    sum_p = sum_c;
	    sum_c = 0.0; 

	    for (p=0; p<projector_num[L]; p++){

	      /* <v|phi_i> */
	      tmp0 = phi0_phi1(VNL_W2[so][L][p],Uo,Vlocal_maxcut);
	      sum_c = sum_c + tmp0;

	      /* |v>*c*<v|phi_i> */
	      for (i=0; i<Grid_Num; i++){
		r = MRV[i];
		VNLp[i] = VNLp[i] + VNL_W2[so][L][p][i]/r*proj_ene[0][so][L][p]*tmp0;
	      }
	    }

	    for (i=0; i<i_r; i++) Lo_p[i] = Lo[i];
	    Hamming_O(1,L,ene,kappa,Mo,Lo,DMo,Vsl,VNLp,Reduce_Num_grid);

	    dif_p = dif_c;
	    dif_c = fabs(sum_p-sum_c);

	    /*
	    printf("ene=%10.5f j=%3d mixing=%15.12f dif_c=%15.12f\n",ene,j,mixing,dif_c);
	    */

	    if (dif_c<1.0e-5){
	      SCF_po = 1;
	    }
	    else{

	      /* simple mixing */

	      if (dif_c<dif_p) renzoku++;

	      if (dif_c<=dif_p && 4<renzoku){
                mixing = 1.03*mixing;
	      }
	      else if (dif_p<dif_c){        
                mixing = 0.4*mixing; 
                renzoku = 0;
              }

	      if (mixing<1.0e-4) mixing = 1.0e-4;
	      if (0.4<mixing)    mixing = 0.4;

	      for (i=0; i<i_r; i++){
		Lo[i] = mixing*Lo[i] + (1.0-mixing)*Lo_p[i];
	      }
	    }

	    if (500<=j && SCF_po==0){

	      printf("warning!, no convergence (ene=%15.10f dif=%15.10f) in LogDeriF.c\n",
                      ene,dif_c);

	      /*
	      printf("warning!, no convergence (ene=%15.10f dif=%15.10f mixing=%10.5f) in LogDeriF.c\n",
                      ene,dif_c,mixing);
	      */

	    } 

	  } while(SCF_po==0 && j<500);

          sign_Lo = sgn(Lo[i_r]); 
	  logd2 = ((double)L + Mo[i_r]/(Lo[i_r]+sign_Lo*geta))/r_i;

	  LogDE[loop] = ene;
	  LogDF[so][L][0][loop] = logd0;
	  LogDF[so][L][1][loop] = logd1;
	  LogDF[so][L][2][loop] = logd2;

	  dif0 = dif0 + (logd2 - logd0)*(logd2 - logd0)*dene;
	  dif1 = dif1 + (logd2 - logd1)*(logd2 - logd1)*dene;

	} /* for (loop=0; ... ) */

	LogD_dif[so][L][0] = dif0; 
	LogD_dif[so][L][1] = dif1;

	printf("<VPS>  SO=%i L=%i, I0, dif of log deris for all electrons = %15.10f\n",so,L,dif0);
	printf("<VPS>  SO=%i L=%i, I1, dif of log deris for semi local    = %15.10f\n",so,L,dif1);
      }
    } /* L */
  } /* so */

} 


double phi0_phi1(double phi0[ASIZE1], double phi1[ASIZE1], double rmax)
{
  static int i,n,j,l,nf,fg;
  static double r,rmin,Sr,Dr,sum,dum;

  n = 96;
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
