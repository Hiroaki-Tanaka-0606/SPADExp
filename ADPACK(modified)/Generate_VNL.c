/**********************************************************************
  Generate_VNL.c:

     Generate_VNL.c is a subroutine to generate the non-local parts
     of pseudo potentials, based on

       (1) Kleinman-Bylander form, 
           ref. L.Kleinman and D.M.Bylander,
           Phys.Rev.Lett. 48, 1425 (1982), 
       (2) Blochl form
           ref. P.E.Blochl, Phys.Rev.B 41, 5414 (1990) 
       (3) Vanderbilt form, but one non-local part for each L is used,
           wave functions are calculated using the pseudo potential
           of most core state of each L.
           ref. D.Vanderbilt, Phys.Rev.B 41, 7892 (1990)

  Log of Generate_VNL.c:

     16/Mar/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "adpack.h"

static double phi0_v_phi1(
         double phi0[ASIZE1],
         double phi1[ASIZE1],
         double vpot[ASIZE1]);

static double phi0_phi1(
         double phi0[ASIZE1],
         double phi1[ASIZE1]);

static double Int_RadF(double R, double RadF[ASIZE1]);
static void Modify_VNL(double Vtmp[ASIZE1], int ip);



void Generate_VNL()
{ 
  static int L,i,ip,j,k,gi,po,m,n,so;
  static int N2,p,q,ii,jj,qp,trial,trial_max=1;
  static int po1,po2;
  static double r,tmp0,dL,ef;
  static double sum,tmp,tmp1,tmp2,tmp3;
  static double VNL0,VNL1,dif,*V0;
  static double phi[ASIZE4*ASIZE15][ASIZE1];
  static double phi2[ASIZE4][ASIZE1];
  static double input_phi[ASIZE1];
  static double pe[ASIZE4];
  static double vv[ASIZE1];
  static double A[ASIZE6];
  static double B[ASIZE15][ASIZE4][ASIZE15*ASIZE4];
  static double B2[ASIZE6][ASIZE6];
  static double C2[ASIZE15*ASIZE4];

  /****************************************************
              calculate a common local part
  ****************************************************/

  /* set Vlocal */
  if (Vlocal_switch==1){
    Calc_Vlocal(Vlocal,V0,0);
  }
  else{

    if (Equation_Type==0 || Equation_Type==1){
      for (i=0; i<Grid_Num; i++){
        Vlocal[i] = VPS[0][0][Local_Part_VPS][i];
      }
    }

    else if (Equation_Type==2){
      for (i=0; i<Grid_Num; i++){
        dL = LVPS[Local_Part_VPS];
        Vlocal[i] = ( (dL+1.0)*VPS[0][0][Local_Part_VPS][i] + dL*VPS[0][1][Local_Part_VPS][i])
                    /(2.0*dL+1.0);
      }
    }
  }

  /****************************************************
           enhancement of spin-orbit splitting 
  ****************************************************/

  if (Equation_Type==2){
    for (L=0; L<ASIZE2; L++){

      gi = GI_VPS[L][0];

      for (i=0; i<Grid_Num; i++){
	VNL0 = VPS[0][0][gi][i] - Vlocal[i];
	VNL1 = VPS[0][1][gi][i] - Vlocal[i];
	sum = 0.5*(VNL0 + VNL1);
	dif = 0.5*(VNL0 - VNL1);
          
	VPS[0][0][gi][i] = sum + SO_factor[L]*dif + Vlocal[i];
	VPS[0][1][gi][i] = sum - SO_factor[L]*dif + Vlocal[i];
      }
    }
  }
  
  /****************************************************
                 generate separable form
  ****************************************************/

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){

      /****************************************************
	   Kleinman-Bylander form
	   ref. L.Kleinman and D.M.Bylander,
	   Phys. Rev. Lett. 48, 1425 (1982)
      ****************************************************/

      if (NL_type[L]==1){

	gi = GI_VPS[L][0];
          
	for (i=0; i<Grid_Num; i++){
	  VNL[0][so][L][0][i] = VPS[0][so][gi][i] - Vlocal[i];
	  VNL_W2[so][L][0][i] = VNL[0][so][L][0][i]*W2[so][gi][i];
	}

	/* calculate project energies */
	if (Vlocal_switch==1 || Equation_Type==2){
	  proj_ene[0][so][L][0] = 1.0/phi0_v_phi1(W2[so][gi], W2[so][gi], VNL[0][so][L][0]);
	}
	else if (gi!=Local_Part_VPS){
	  proj_ene[0][so][L][0] = 1.0/phi0_v_phi1(W2[so][gi], W2[so][gi], VNL[0][so][L][0]);
	}
	else{
	  proj_ene[0][so][L][0] = 0.0;
	}
      } 

      /****************************************************
	   Blochl form
	   ref. P.E.Blochl, Phys.Rev.B 41, 5414 (1990) 
      ****************************************************/

      else if (NL_type[L]==2 || NL_type[L]==5){

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
      }

      /****************************************************
         Modified Vanderbilt form for norm-conserving VPS
         ref. D.Vanderbilt, Phys.Rev.B 41, 7892 (1990)
      ****************************************************/

      else if (NL_type[L]==3){

	/* Set non-local potential and an initial wave function */
	for (m=0; m<NumVPS_L[L]; m++){
	  gi = GI_VPS[L][m];
	  for (i=0; i<Grid_Num; i++){
	    VNL[0][so][L][m][i] = VPS[0][so][gi][i] - Vlocal[i];
	    phi[m][i] = W2[so][gi][i];
	  }
	}

	/* Gramm-Schmidt orthogonalization with a norm defined by <f|v|g> */
	for (i=0; i<Grid_Num; i++){
	  VNL_W2[so][L][0][i] = phi[0][i]; 
	}          
	pe[0] = 1.0/phi0_v_phi1(VNL_W2[so][L][0],VNL_W2[so][L][0],VNL[0][so][L][0]);
	for (m=1; m<NumVPS_L[L]; m++){
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = 0.0;
	  }          

	  for (n=0; n<m; n++){
	    /* it's OK for the selection of 0 */ 
	    tmp0 = phi0_v_phi1(VNL_W2[so][L][n],phi[m],VNL[0][so][L][0]); 
	    for (i=0; i<Grid_Num; i++){
	      VNL_W2[so][L][m][i] = VNL_W2[so][L][m][i] + VNL_W2[so][L][n][i]*pe[n]*tmp0;  
	    }
	  }
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = phi[m][i] - VNL_W2[so][L][m][i];
	  }
	  /* it's OK for the selection of 0 */ 
	  pe[m] = 1.0/phi0_v_phi1(VNL_W2[so][L][m],VNL_W2[so][L][m],VNL[0][so][L][0]);
	}

	/* Renormalization */

	for (m=0; m<NumVPS_L[L]; m++){
	  tmp0 = 1.0/sqrt(phi0_phi1(VNL_W2[so][L][m],VNL_W2[so][L][m]));
	  for (i=0; i<Grid_Num; i++){
	    VNL_W2[so][L][m][i] = tmp0*VNL_W2[so][L][m][i];
	  }
	  proj_ene[0][so][L][m] = pe[m]/(tmp0*tmp0);
	}

	/* Calc v*VNL_W2 */
	for (m=0; m<NumVPS_L[L]; m++){
	  for (i=0; i<Grid_Num; i++){
	    /* it's OK for the selection of 0 */ 
	    VNL_W2[so][L][m][i] = VNL[0][so][L][0][i]*VNL_W2[so][L][m][i];
	  }
	}
      }

    } /* L */
  } /* so */

} 


double phi0_phi1(
         double phi0[ASIZE1],
         double phi1[ASIZE1])
{

  static int i,n,j,l,nf,fg;
  static double r,rmin,rmax,Sr,Dr,sum,dum;

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

double phi0_v_phi1(
         double phi0[ASIZE1],
         double phi1[ASIZE1],
         double vpot[ASIZE1])
{

  static int i,n,j,l,nf,fg;
  static double r,rmin,rmax,Sr,Dr,sum,dum;

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


void Modify_VNL(double Vtmp[ASIZE1], int ip)
{
  static int i,j,jp,po,m,n,l;
  static double Vi,dVi,d2Vi,d3Vi,ri,rj,r;
  static double Vj,dVj,d2Vj,d3Vj;
  static double c0,c1,c2,c3,c4,c5,c6,c7;
  static double r2,r3,r4,r5,r6,r7,r8;
  static double dis_r,sum,p,p2,MinD;
  static double a[ASIZE6][ASIZE6],c[ASIZE6];
  static double tmp0,Z,dr;
  static double vm1,vm2,vp1,vp2;

  /****************************************************
    Value, first, second and third derivatives at ri
  ****************************************************/

  dr = 0.001;
  ri = MRV[ip];
  Vi = log(Vtmp[ip]);
  vp1 = log(HokanF(ri+1.0*dr,Vtmp,0));
  vp2 = log(HokanF(ri+2.0*dr,Vtmp,0));
  vm1 = log(HokanF(ri-1.0*dr,Vtmp,0));
  vm2 = log(HokanF(ri-2.0*dr,Vtmp,0));
  dVi  = (vm2 - 8.0*vm1 + 8.0*vp1 - vp2)/(12.0*dr);
  d2Vi = (-vm2 + 16.0*vm1  - 30.0*Vi + 16.0*vp1 - vp2)/(12.0*dr*dr);
  d3Vi = (-vm2 + 2.0*vm1 - 2.0*vp1 + vp2)/(2.0*dr*dr*dr);

  /****************************************************
    Value, first, second and third derivatives at rj
  ****************************************************/

  MinD = 10000.0;
  for (i=0; i<Grid_Num; i++){
    tmp0 = fabs(MRV[i]-(ri+0.3));
    if (tmp0<MinD){
      MinD = tmp0;
      jp = i;
    }
  }

  rj = MRV[jp];
  Vj   = log(1.0e-13);
  dVj  = -1.0e+0;
  d2Vj = -1.0e+1;
  d3Vj = -1.0e+2;

  /****************************************************
                           At ri
  ****************************************************/

  r2 = ri*ri;
  r3 = ri*r2;
  r4 = ri*r3;
  r5 = ri*r4;
  r6 = ri*r5;
  r7 = ri*r6;
 
  a[0][0] = 1.0;
  a[0][1] = ri;
  a[0][2] = r2;
  a[0][3] = r3;
  a[0][4] = r4;
  a[0][5] = r5;
  a[0][6] = r6;
  a[0][7] = r7;

  a[1][0] = 0.0;
  a[1][1] = 1.0;
  a[1][2] = 2.0*ri;
  a[1][3] = 3.0*r2;
  a[1][4] = 4.0*r3;
  a[1][5] = 5.0*r4;
  a[1][6] = 6.0*r5;
  a[1][7] = 7.0*r6;

  a[2][0] = 0.0;
  a[2][1] = 0.0;
  a[2][2] = 2.0;
  a[2][3] = 6.0*ri;
  a[2][4] = 12.0*r2;
  a[2][5] = 20.0*r3;
  a[2][6] = 30.0*r4;
  a[2][7] = 42.0*r5;

  a[3][0] = 0.0;
  a[3][1] = 0.0;
  a[3][2] = 0.0;
  a[3][3] = 6.0;
  a[3][4] = 24.0*ri;
  a[3][5] = 60.0*r2;
  a[3][6] = 120.0*r3;
  a[3][7] = 210.0*r4;

  a[0][8] = Vi;
  a[1][8] = dVi;
  a[2][8] = d2Vi;
  a[3][8] = d3Vi;

  /****************************************************
                         At rj
  ****************************************************/

  r2 = rj*rj;
  r3 = rj*r2;
  r4 = rj*r3;
  r5 = rj*r4;
  r6 = rj*r5;
  r7 = rj*r6;
 
  a[4][0] = 1.0;
  a[4][1] = rj;
  a[4][2] = r2;
  a[4][3] = r3;
  a[4][4] = r4;
  a[4][5] = r5;
  a[4][6] = r6;
  a[4][7] = r7;

  a[5][0] = 0.0;
  a[5][1] = 1.0;
  a[5][2] = 2.0*rj;
  a[5][3] = 3.0*r2;
  a[5][4] = 4.0*r3;
  a[5][5] = 5.0*r4;
  a[5][6] = 6.0*r5;
  a[5][7] = 7.0*r6;

  a[6][0] = 0.0;
  a[6][1] = 0.0;
  a[6][2] = 2.0;
  a[6][3] = 6.0*rj;
  a[6][4] = 12.0*r2;
  a[6][5] = 20.0*r3;
  a[6][6] = 30.0*r4;
  a[6][7] = 42.0*r5;

  a[7][0] = 0.0;
  a[7][1] = 0.0;
  a[7][2] = 0.0;
  a[7][3] = 6.0;
  a[7][4] = 24.0*rj;
  a[7][5] = 60.0*r2;
  a[7][6] = 120.0*r3;
  a[7][7] = 210.0*r4;

  a[4][8] = Vj;
  a[5][8] = dVj;
  a[6][8] = d2Vj;
  a[7][8] = d3Vj;

  /****************************************************
              c0,c1,c2,c3,c4,c5,c6, and c7
  ****************************************************/

  Gauss_LEQ(7,a,c);
  c0 = c[0];
  c1 = c[1];
  c2 = c[2];
  c3 = c[3];
  c4 = c[4];
  c5 = c[5];
  c6 = c[6];
  c7 = c[7];

  /****************************************************
                       modify Vtmp
  ****************************************************/

  for (i=ip; i<=jp; i++){
    r = MRV[i];
    r2 = r*r;
    r3 = r2*r;
    r4 = r2*r2;
    r5 = r4*r;
    r6 = r2*r4;
    r7 = r6*r;
    Vtmp[i] = exp(c0 + c1*r  + c2*r2 + c3*r3 + c4*r4 + c5*r5 + c6*r6 + c7*r7);
  }

  for (i=(jp+1); i<Grid_Num; i++){
    Vtmp[i] = 0.0;
  }

}



