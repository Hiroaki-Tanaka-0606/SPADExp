/**********************************************************************
  Initial_Density.c:

     Initial_Density.c is a subroutine to give the initial density
     based on a hydrogen atom.

  Log of Initial_Density.c:

     31/May/2001  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

double F_iden(double R);

double iden_xv[101];
double iden_rv[101];
double iden[101];
int num;

void Initial_Density()
{
  static int i,j;
  static double dum;
  static double r,tmp0,tmp1,tmp2,Z,b;
  char FN_DEN[ASIZE8];
  FILE *fp;

  if (AtomNum==1)  chcp(FN_DEN, "../idensity/H.aden");
  if (AtomNum==2)  chcp(FN_DEN, "../idensity/He.aden");
  if (AtomNum==3)  chcp(FN_DEN, "../idensity/Li.aden");
  if (AtomNum==4)  chcp(FN_DEN, "../idensity/Be.aden");
  if (AtomNum==5)  chcp(FN_DEN, "../idensity/B.aden");
  if (AtomNum==6)  chcp(FN_DEN, "../idensity/C.aden");
  if (AtomNum==7)  chcp(FN_DEN, "../idensity/N.aden");
  if (AtomNum==8)  chcp(FN_DEN, "../idensity/O.aden");
  if (AtomNum==9)  chcp(FN_DEN, "../idensity/F.aden");
  if (AtomNum==10) chcp(FN_DEN, "../idensity/Ne.aden");
  if (AtomNum==11) chcp(FN_DEN, "../idensity/Na.aden");
  if (AtomNum==12) chcp(FN_DEN, "../idensity/Mg.aden");
  if (AtomNum==13) chcp(FN_DEN, "../idensity/Mg.aden");
  if (AtomNum==14) chcp(FN_DEN, "../idensity/Si.aden");
  if (AtomNum==15) chcp(FN_DEN, "../idensity/P.aden");
  if (AtomNum==16) chcp(FN_DEN, "../idensity/S.aden");
  if (AtomNum==17) chcp(FN_DEN, "../idensity/Cl.aden");
  if (AtomNum==18) chcp(FN_DEN, "../idensity/Ar.aden");
  if (AtomNum==19) chcp(FN_DEN, "../idensity/K.aden");
  if (AtomNum==20) chcp(FN_DEN, "../idensity/Ca.aden");
  if (AtomNum==21) chcp(FN_DEN, "../idensity/Sc.aden");
  if (AtomNum==22) chcp(FN_DEN, "../idensity/Ti.aden");
  if (AtomNum==23) chcp(FN_DEN, "../idensity/V.aden");
  if (AtomNum==24) chcp(FN_DEN, "../idensity/Cr.aden");
  if (AtomNum==25) chcp(FN_DEN, "../idensity/Mn.aden");
  if (AtomNum==26) chcp(FN_DEN, "../idensity/Fe.aden");
  if (AtomNum==27) chcp(FN_DEN, "../idensity/Co.aden");
  if (AtomNum==28) chcp(FN_DEN, "../idensity/Ni.aden");
  if (AtomNum==29) chcp(FN_DEN, "../idensity/Cu.aden");
  if (AtomNum==30) chcp(FN_DEN, "../idensity/Zn.aden");
  if (AtomNum==31) chcp(FN_DEN, "../idensity/Ga.aden");
  if (AtomNum==32) chcp(FN_DEN, "../idensity/Ge.aden");
  if (AtomNum==33) chcp(FN_DEN, "../idensity/As.aden");

  if (AtomNum==79) chcp(FN_DEN, "../idensity/Au.aden");
  if (AtomNum==80) chcp(FN_DEN, "../idensity/Hg.aden");

  fp = fopen(FN_DEN, "r");

  if (fp != NULL){

    num = 100;
    for (i=0; i<num; i++){

      for (j=0; j<=2; j++){
        if (fscanf(fp,"%lf",&dum)==EOF){
          printf("File error in idensity\n");
        }
        else{
          if (j==0)      iden_xv[i] = dum;
          else if (j==1) iden_rv[i] = dum;
          else if (j==2) iden[i] = dum;
        }
      }
    }

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      tmp1 = F_iden(r);
      rho[0][i] = tmp1;
      rho[1][i] = tmp1;
    }

    fclose(fp);
  }
  else{

    Z = (double)AtomNum;
    b = pow(Z,2.0/3.0);
    tmp0 = pow(b,1.5)/sqrt(PI);

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      tmp1 = tmp0*exp(-b*r);
      tmp2 = 0.5*Z*tmp1*tmp1;
      rho[0][i] = tmp2;
      rho[1][i] = tmp2;
    }
  }
}

double F_iden(double R)
{
  static int mp_min,mp_max,m;
  static double h1,h2,h3,f1,f2,f3,f4;
  static double g1,g2,x1,x2,y1,y2,f;
  static double result;

  mp_min = 0;
  mp_max = num - 1;

  if (R<iden_rv[0]){
    f = iden[0];
  }
  else if (iden_rv[num-1]<R){
    f = 0.0;
  }
  else{
    do{
      m = (mp_min + mp_max)/2;
      if (iden_rv[m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    /****************************************************
                   Spline like interpolation
    ****************************************************/


    if (m==1){
      h2 = iden_rv[m]   - iden_rv[m-1];
      h3 = iden_rv[m+1] - iden_rv[m];

      f2 = iden[m-1];
      f3 = iden[m];
      f4 = iden[m+1];

      h1 = -(h2+h3);
      f1 = f4;
    }
    else if (m==(num-1)){
      h1 = iden_rv[m-1] - iden_rv[m-2];
      h2 = iden_rv[m]   - iden_rv[m-1];

      f1 = iden[m-2];
      f2 = iden[m-1];
      f3 = iden[m];

      h3 = -(h1+h2);
      f4 = f1;
    }
    else{
      h1 = iden_rv[m-1] - iden_rv[m-2];
      h2 = iden_rv[m]   - iden_rv[m-1];
      h3 = iden_rv[m+1] - iden_rv[m];

      f1 = iden[m-2];
      f2 = iden[m-1];
      f3 = iden[m];
      f4 = iden[m+1];
    } 

    /****************************************************
                Calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - iden_rv[m-1];
    x2 = R - iden_rv[m];
    y1 = x1/h2;
    y2 = x2/h2;

    f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }
  result = f;
  return result;
}

