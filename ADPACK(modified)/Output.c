/**********************************************************************
  Output.c:

     Output.c is a subroutine to output calculated quantities as files.

  Log of Output.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "adpack.h"

static void Out_Std(FILE *fp, char *filein); 
static void Out_AllLOG(char *filein);
static void Output_Local_Pseudo_Potentials();
static void Output_Nonseparable_Pseudo_Potentials();
static void Output_Separable_Pseudo_Potentials();
static void Output_AllBases();
static void Output_PAOBases();
static void Output_PAOBases2(char *filein);
static void Output_Density();
static void Output_LogD_RadF();
static void Output_TM_Pseudo_Potentials();

void Output(char *filein)
{

  /* All */
  if (Calc_Type==0){
    printf("\nThe following files are generated.\n");
    Out_AllLOG(filein);
    Output_AllBases();
    Output_Density();
  }

  /* VPS */
  else if (Calc_Type==1){
    printf("\nThe following files are generated.\n");
    Output_PAOBases();
    Output_Density();
    Output_Local_Pseudo_Potentials();
    Output_Nonseparable_Pseudo_Potentials();
    Output_Separable_Pseudo_Potentials(filein);
    if (LogD_switch==1) Output_LogD_RadF();

    if (VPP_switch==4){
      Output_TM_Pseudo_Potentials();
    }
  }

  /* PAO */
  else if (Calc_Type==2){
    printf("\nThe following files are generated.\n");
    Output_PAOBases2(filein);
  }

  /* FCORE */
  else if (Calc_Type==3){
    printf("\nThe following files are generated.\n");
    Out_AllLOG(filein);
  }

}



void Out_Std(FILE *fp, char *filein) 
{
  int c,ii,i,j,n,l,SCF;
  char *s_vec[20];
  FILE *fp2;

  fprintf(fp,"***************************************************\n");
  fprintf(fp,"                    Input file\n"                     );
  fprintf(fp,"***************************************************\n\n");

  /* echo the input file */

  fp2 = fopen(filein,"r");
  if (fp2!=NULL){
    for (c=getc(fp2); c!=EOF; c=getc(fp2)){
      putc(c,fp); 
    }
    fclose(fp2); 
  }

  /* fprintf SCF history */

  fprintf(fp,"\n");
  fprintf(fp,"*****************************************************\n");
  fprintf(fp,"    SCF history in all electron calculations      \n");
  fprintf(fp,"*****************************************************\n\n");

  for (SCF=1; SCF<SCF_END; SCF++){
    fprintf(fp," SCF=%4d  Eeigen=%17.13f (Hartree)  NormRD=%17.13f\n",
	    SCF,HisEeigen[SCF],HisNormRD[SCF]);
  }

  /* fprintf eigenvalues */

  fprintf(fp,"\n");
  fprintf(fp,"*****************************************************\n");
  fprintf(fp," Eigenvalues (Hartree) in the all electron calculation\n");
  fprintf(fp,"*****************************************************\n\n");

  if (SOI_switch==1){

    if (Calc_Type==3) j = charge_state_num;
    else              j = 1; 

    for (ii=0; ii<j; ii++){

      if (VPP_switch==3){
	fprintf(fp,"\n charge state = %6.3f\n\n",charge_states[ii]);
      }
      else if (Calc_Type==3){
        if (ii==0)
  	  fprintf(fp," Self consistent\n\n");
        else if (ii==1)
  	  fprintf(fp,"\n Frozen core approximation\n\n");
      }

      for (n=1; n<=max_ocupied_N; n++){
	for (l=0; l<n; l++){
	  if (0.0<OcpN[ii][0][n][l]){
	    fprintf(fp," n=%3d  l=%3d  %22.13f\n",n,l,E[ii][0][n][l]);
	  }
	}
      }
    }
  }
  else if (SOI_switch==2){

    for (ii=0; ii<charge_state_num; ii++){

      if (VPP_switch==3){
	fprintf(fp,"\n charge state = %6.3f\n\n",charge_states[ii]);
      }
      else if (Calc_Type==3){
        if (ii==0)
  	  fprintf(fp," Self consistent\n\n");
        else if (ii==1)
  	  fprintf(fp,"\n Frozen core approximation\n\n");
      }

      fprintf(fp,"                         j=l+1/2                j=l-1/2\n");
      for (n=1; n<=max_ocupied_N; n++){
	for (l=0; l<n; l++){
	  if (0.0<OcpN[ii][0][n][l]){
	    fprintf(fp," n=%3d  l=%3d  %22.13f %22.13f\n",
		    n,l,E[ii][0][n][l],E[ii][1][n][l]);
	  }
	}
      }
    }
  }

  fprintf(fp,"\n\n");
  fprintf(fp,"*****************************************************\n");
  fprintf(fp," Energies (Hartree) in the all electron calculation \n");
  fprintf(fp,"*****************************************************\n\n");

  if (Calc_Type==3){

    for (ii=0; ii<charge_state_num; ii++){

      if (Calc_Type==3){
	if (ii==0)
	  fprintf(fp," Self consistent\n\n");
	else if (ii==1)
	  fprintf(fp,"\n Frozen core approximation\n\n");
      }

      fprintf(fp," Eeigen = %22.13f\n",Eeigen[ii]);
      fprintf(fp," Ekin   = %22.13f\n",Ekin[ii]);
      fprintf(fp," EHart  = %22.13f\n",EHart[ii]);
      fprintf(fp," Exc    = %22.13f\n",Exc[ii]);
      fprintf(fp," Eec    = %22.13f\n",Eec[ii]);
      fprintf(fp," Etot   = Ekin + EHart + Exc + Eec\n");
      fprintf(fp," Etot   = %22.13f\n",Etot[ii]);
    }

  }

  else{
    fprintf(fp," Eeigen = %22.13f\n",Eeigen[0]);
    fprintf(fp," Ekin   = %22.13f\n",Ekin[0]);
    fprintf(fp," EHart  = %22.13f\n",EHart[0]);
    fprintf(fp," Exc    = %22.13f\n",Exc[0]);
    fprintf(fp," Eec    = %22.13f\n",Eec[0]);
    fprintf(fp," Etot   = Ekin + EHart + Exc + Eec\n");
    fprintf(fp," Etot   = %22.13f\n",Etot[0]);
  }

}


static void Out_AllLOG(char *filein)
{
  static int i,n,l,SCF;
  static char file0[ASIZE8] = ".alog";
  char *s_vec[20];
  FILE *fp;

  fnjoint(filepath,filename,file0);
  if ((fp = fopen(file0,"w")) != NULL){
    Out_Std(fp,filein); 
    fclose(fp);
    printf("  %s\n",file0);
  }
  else
    printf("Failure of saving the Eigenvalues.\n");
}










static void Output_Local_Pseudo_Potentials()
{
  static int i,m;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".loc";
  FILE *fp;

  fnjoint(filepath,filename,file0);

  if ((fp = fopen(file0,"w")) != NULL){

    Out_XMIN = Grid_Xmin; 
    Out_XMAX = Grid_Xmax - 0.001;
    dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);

    for (i=0; i<Grid_Num_Output; i++){
      x = Out_XMIN + dx0*(double)i;
      r = exp(x);
      fprintf(fp,"%17.14f %17.14f",x,r);
      for (m=0; m<Number_VPS; m++){
        fprintf(fp,"  %17.14f ",HokanF(r,Vlocal,0));
      }
      fprintf(fp,"\n"); 
    }
    fclose(fp);
    printf("  %s\n",file0);
  }
  else
    printf("Failure of saving the Local_Pseudo_Potentials.\n");
}







static void Output_Nonseparable_Pseudo_Potentials()
{
  static int i,m,so,l,mul;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".nsvps";
  static char fname[ASIZE8];
  FILE *fp;

  /* GVPS */

  if (VPP_switch==4){

    if (Equation_Type==0 || Equation_Type==1){

      sprintf(fname,"%s%s%s",filepath,filename,file0);

      if ((fp = fopen(fname,"w")) != NULL){

	Out_XMIN = Grid_Xmin; 
	Out_XMAX = Grid_Xmax - 0.001;

	dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%17.14f %17.14f",x,r);

          for (l=0; l<=GVPS_MaxL; l++){
            for (mul=0; mul<GVPS_ProNum; mul++){
              fprintf(fp,"  %17.14f ",HokanF(r,GVPS[0][l][mul],0));
	    }
	  }

	  fprintf(fp,"\n"); 
	}
	fclose(fp);
	printf("  %s\n",fname);
      }
      else
	printf("Failure of saving the Nonseparable_Pseudo_Potentials.\n");

    }

    else if (Equation_Type==2){

      for (so=0; so<SOI_switch; so++){

	sprintf(fname,"%s%s%s%i",filepath,filename,file0,so);

	if ((fp = fopen(fname,"w")) != NULL){

	  Out_XMIN = Grid_Xmin; 
	  Out_XMAX = Grid_Xmax - 0.001;

	  dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
	  for (i=0; i<Grid_Num_Output; i++){
	    x = Out_XMIN + dx0*(double)i;
	    r = exp(x);
	    fprintf(fp,"%17.14f %17.14f",x,r);

            for (l=0; l<=GVPS_MaxL; l++){
              for (mul=0; mul<GVPS_ProNum; mul++){
                fprintf(fp,"  %17.14f ",HokanF(r,GVPS[so][l][mul],0));
	      }
	    }

	    fprintf(fp,"\n"); 
	  }
	  fclose(fp);
	  printf("  %s\n",fname);
	}
	else
	  printf("Failure of saving the Nonseparable_Pseudo_Potentials.\n");

      }
    }
  }

  /* non-GVPS */

  else {

    if (Equation_Type==0 || Equation_Type==1){

      sprintf(fname,"%s%s%s",filepath,filename,file0);

      if ((fp = fopen(fname,"w")) != NULL){

	Out_XMIN = Grid_Xmin; 
	Out_XMAX = Grid_Xmax - 0.001;

	dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%17.14f %17.14f",x,r);
	  for (m=0; m<Number_VPS; m++){
	    fprintf(fp,"  %17.14f ",HokanF(r,VPS[0][0][m],0));
	  }
	  fprintf(fp,"\n"); 
	}
	fclose(fp);
	printf("  %s\n",fname);
      }
      else
	printf("Failure of saving the Nonseparable_Pseudo_Potentials.\n");

    }

    else if (Equation_Type==2){

      for (so=0; so<SOI_switch; so++){

	sprintf(fname,"%s%s%s%i",filepath,filename,file0,so);

	if ((fp = fopen(fname,"w")) != NULL){

	  Out_XMIN = Grid_Xmin; 
	  Out_XMAX = Grid_Xmax - 0.001;

	  dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
	  for (i=0; i<Grid_Num_Output; i++){
	    x = Out_XMIN + dx0*(double)i;
	    r = exp(x);
	    fprintf(fp,"%17.14f %17.14f",x,r);
	    for (m=0; m<Number_VPS; m++){
	      fprintf(fp,"  %17.14f ",HokanF(r,VPS[0][so][m],0));
	    }
	    fprintf(fp,"\n"); 
	  }
	  fclose(fp);
	  printf("  %s\n",fname);
	}
	else
	  printf("Failure of saving the Nonseparable_Pseudo_Potentials.\n");

      }
    }
  }

}





static void Output_TM_Pseudo_Potentials()
{
  static int i,m,so,l,mul;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".tmvps";
  static char fname[ASIZE8];
  FILE *fp;

  if (Equation_Type==0 || Equation_Type==1){

    sprintf(fname,"%s%s%s",filepath,filename,file0);

    if ((fp = fopen(fname,"w")) != NULL){

      Out_XMIN = Grid_Xmin; 
      Out_XMAX = Grid_Xmax - 0.001;

      dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
      for (i=0; i<Grid_Num_Output; i++){
	x = Out_XMIN + dx0*(double)i;
	r = exp(x);
	fprintf(fp,"%17.14f %17.14f",x,r);

	for (l=0; l<=GVPS_MaxL; l++){
	  fprintf(fp,"  %17.14f ",HokanF(r,TMVPS[0][l],0));
	}

	fprintf(fp,"\n"); 
      }
      fclose(fp);
      printf("  %s\n",fname);
    }
    else
      printf("Failure of saving the TM_Nonseparable_Pseudo_Potentials.\n");

  }

  else if (Equation_Type==2){

    for (so=0; so<SOI_switch; so++){

      sprintf(fname,"%s%s%s%i",filepath,filename,file0,so);

      if ((fp = fopen(fname,"w")) != NULL){

	Out_XMIN = Grid_Xmin; 
	Out_XMAX = Grid_Xmax - 0.001;

	dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%17.14f %17.14f",x,r);

	  for (l=0; l<=GVPS_MaxL; l++){
	    fprintf(fp,"  %17.14f ",HokanF(r,TMVPS[so][l],0));
	  }

	  fprintf(fp,"\n"); 
	}
	fclose(fp);
	printf("  %s\n",fname);
      }
      else
	printf("Failure of saving the TM_Nonseparable_Pseudo_Potentials.\n");

    }
  }

}



static void Output_Separable_Pseudo_Potentials(char *filein)
{
  static int i,ii,jj,L,m,so,mul,n;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".vps";
  static char fname[ASIZE8];
  FILE *fp;

  /* Schrodinger and scalar relativistic*/
  if (Equation_Type==0 || Equation_Type==1){

    sprintf(fname,"%s%s%s",filepath,filename,file0);

    if ((fp = fopen(fname,"w")) != NULL){

      Out_Std(fp,filein); 

      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"** DATA for factorized norm conserving pseudo potentials **\n");
      fprintf(fp,"***********************************************************\n");

      Out_XMIN = Grid_Xmin;
      Out_XMAX = Grid_Xmax - 0.001;
      dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);

      if (LogD_switch==1){
	fprintf(fp,"\n");
	for (L=0; L<ASIZE2; L++){
	  if (NumVPS_L[L]!=0){
	    fprintf(fp,"<VPS> L=%i, dif of log deris for all electrons = %15.10f\n",
		    L,LogD_dif[0][L][0]);
	    fprintf(fp,"<VPS> L=%i, dif of log deris for semi local    = %15.10f\n",
		    L,LogD_dif[0][L][1]);
	  }
	}
	fprintf(fp,"\n");
      }

      fprintf(fp,"<project.energies\n");
      fprintf(fp,"%2d\n",total_pro_num);
      for (L=0; L<ASIZE2; L++){
	if (projector_num[L]!=0){
	  for (m=0; m<projector_num[L]; m++){
	    fprintf(fp,"%2d  %22.14e\n",L,proj_ene[0][0][L][m]);
	  }
	}
      }
      fprintf(fp,"project.energies>\n\n");

      fprintf(fp,"<Pseudo.Potentials\n");
      for (i=0; i<Grid_Num_Output; i++){
	x = Out_XMIN + dx0*(double)i;
	r = exp(x);
	fprintf(fp,"%17.14e %17.14e",x,r);
	fprintf(fp,"  %17.14e ", HokanF(r,Vlocal,0));
	for (L=0; L<ASIZE2; L++){
	  if (projector_num[L]!=0){
	    for (m=0; m<projector_num[L]; m++){
	      fprintf(fp,"  %17.14e",HokanF(r,VNL_W2[0][L][m],1));
	    }
	  }
	}
	fprintf(fp,"\n"); 
      }
      fprintf(fp,"Pseudo.Potentials>\n");

      /* PCC charge */

      if (PCC_switch==1){
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"***********************************************************\n");
	fprintf(fp,"**            Core electron densities for PCC            **\n");
	fprintf(fp,"***********************************************************\n");

	fprintf(fp,"<density.PCC\n");
	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%17.14e %17.14e",x,r);
	  fprintf(fp,"  %17.14e\n",fabs(HokanF(r,rho_PCC,0)));
	}
	fprintf(fp,"density.PCC>\n");
      }

      /* parameters for LESP */

      /*
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"**                  parameters for LESP                  **\n");
      fprintf(fp,"***********************************************************\n");

      fprintf(fp,"LESP.CoreR    %17.14e\n\n",Core_R[0]);

      fprintf(fp,"<LESP.IntProj\n");
      for (i=0; i<charge_state_num; i++){
	fprintf(fp,"  %17.14e\n",Int_EDPP_Proj[i]);
      }
      fprintf(fp,"LESP.IntProj>\n\n");

      fprintf(fp,"<LESP.CoesEZ\n");
      fprintf(fp," %2d\n",charge_state_num);
      for (i=0; i<charge_state_num; i++){
	fprintf(fp,"  %18.15f\n",CoesEZ[i]);
      }
      fprintf(fp,"LESP.CoesEZ>\n\n");
      */

      /* file close */
      fclose(fp);
      printf("  %s\n",fname);
    }
    else
      printf("Failure of saving the Separable_Pseudo_Potentials.\n");
  }

  /* full relativistic */

  else if (Equation_Type==2){

    sprintf(fname,"%s%s%s",filepath,filename,file0);

    if ((fp = fopen(fname,"w")) != NULL){

      Out_Std(fp,filein); 

      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"** DATA for factorized norm conserving pseudo potentials **\n");
      fprintf(fp,"***********************************************************\n");

      fprintf(fp,"\nj.dependent.pseudo.potentials   on\n");

      Out_XMIN = Grid_Xmin; 
      Out_XMAX = Grid_Xmax - 0.001;
      dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);

      if (LogD_switch==1){
	fprintf(fp,"\n");
	for (L=0; L<ASIZE2; L++){
	  if (NumVPS_L[L]!=0){
	    fprintf(fp,"<VPS> L=%i, dif of log deris for all electrons = %15.10f %15.10f\n",
		    L,LogD_dif[0][L][0],LogD_dif[1][L][0]);
	    fprintf(fp,"<VPS> L=%i, dif of log deris for semi local    = %15.10f %15.10f\n",
		    L,LogD_dif[0][L][1],LogD_dif[1][L][1]);
	  }
	}
	fprintf(fp,"\n");
      }
    
      fprintf(fp,"<project.energies\n");
      fprintf(fp,"%2d\n",total_pro_num);
      for (L=0; L<ASIZE2; L++){
	if (projector_num[L]!=0){
	  for (m=0; m<projector_num[L]; m++){
	    fprintf(fp,"%2d  %22.14e %22.14e\n",L,proj_ene[0][0][L][m],proj_ene[0][1][L][m]);
	  }
	}
      }
      fprintf(fp,"project.energies>\n\n");
     
      fprintf(fp,"<Pseudo.Potentials\n");
      for (i=0; i<Grid_Num_Output; i++){
	x = Out_XMIN + dx0*(double)i;
	r = exp(x);
	fprintf(fp,"%17.14e %17.14e",x,r);
	fprintf(fp,"  %17.14e ", HokanF(r,Vlocal,0));
	for (L=0; L<ASIZE2; L++){
	  if (projector_num[L]!=0){
	    for (m=0; m<projector_num[L]; m++){
	      fprintf(fp,"  %17.14e %17.14e",
                      HokanF(r,VNL_W2[0][L][m],1),
                      HokanF(r,VNL_W2[1][L][m],1));
	    }
	  }
	}
	fprintf(fp,"\n"); 
      }
      fprintf(fp,"Pseudo.Potentials>\n");

      /* PCC charge */

      if (PCC_switch==1){
	fprintf(fp,"\n");
	fprintf(fp,"\n");
	fprintf(fp,"***********************************************************\n");
	fprintf(fp,"**            Core electron densities for PCC            **\n");
	fprintf(fp,"***********************************************************\n");

	fprintf(fp,"<density.PCC\n");
	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%17.14e %17.14e",x,r);
	  fprintf(fp,"  %17.14e\n",fabs(HokanF(r,rho_PCC,0)));
	}
	fprintf(fp,"density.PCC>\n");
      }

      /* potential projectors for EDPP */

      /*
      fprintf(fp,"\n");
      fprintf(fp,"\n");
      fprintf(fp,"***********************************************************\n");
      fprintf(fp,"**                  parameters for LESP                  **\n");
      fprintf(fp,"***********************************************************\n");

      fprintf(fp,"LESP.CoreR    %17.14e\n\n",Core_R[0]);

      fprintf(fp,"<LESP.IntProj\n");
      for (i=0; i<charge_state_num; i++){
	fprintf(fp,"  %17.14e\n",Int_EDPP_Proj[i]);
      }
      fprintf(fp,"LESP.IntProj>\n\n");

      fprintf(fp,"<LESP.CoesEZ\n");
      fprintf(fp," %2d\n",charge_state_num);
      for (i=0; i<charge_state_num; i++){
	fprintf(fp,"  %18.15f\n",CoesEZ[i]);
      }
      fprintf(fp,"LESP.CoesEZ>\n\n");
      */

      /* file close */
      fclose(fp);
      printf("  %s\n",fname);
    }
    else
      printf("Failure of saving the Separable_Pseudo_Potentials.\n");

  }

}


static void Output_AllBases()
{
  static int i,m,n,l,multi,so;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".ao";
  FILE *fp;

  fnjoint(filepath,filename,file0);

  if ((fp = fopen(file0,"w")) != NULL){
		fprintf(fp, "grid.num.output %d\n", Grid_Num_Output);
		fprintf(fp, "maxL.pao %d\n", MaxL_PAO);
		fprintf(fp, "max.N %d\n", max_N);
		
    Out_XMIN = Grid_Xmin; 
    // Out_XMAX = Grid_Xmax;
    // dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);
    Out_XMAX = log(PAO_Rcut+1.0);
    dx0 = (Out_XMAX - Out_XMIN)/Grid_Num_Output;
		/*
    for (so=0; so<SOI_switch; so++){
			// for (n=1; n<=max_ocupied_N; n++){
			// for (l=0; l<n; l++){
				for (n=1; n<=max_N; n++){
					fprintf(fp,"so=%2d n=%2d\n",so,n);
        for (i=0; i<Grid_Num_Output; i++){
          x = Out_XMIN + dx0*(double)i;
          r = exp(x);
          fprintf(fp,"%15.12f %15.12f  ",x,r);
					// for (l=0; l<n; l++){
					for(l=0; (l<n && l<=MaxL_PAO); l++){
						fprintf(fp,"%15.12f ",HokanF(r,PF[so][n][l],1)); 
          }
          fprintf(fp,"\n");
        }
      }
			}*/

		/* beginning of the addition */
		for(l=0; l<=MaxL_PAO; l++){
			fprintf(fp, "<atomic.orbitals.L=%d\n", l);
			for (i=0; i<Grid_Num_Output; i++){
				x = Out_XMIN + dx0*(double)i;
				r = exp(x);
				fprintf(fp,"%15.12f %15.12f ",x,r);
				for(n=l+1; n<=max_N; n++){
					fprintf(fp, "%15.12f ", HokanF(r, PF[0][n][l],1));
				}
				fprintf(fp, "\n");
			}
			fprintf(fp, "atomic.orbitals.L=%d>\n", l);
		}
		/* end of the addition */
    fclose(fp);
    printf("  %s\n",file0);
  }
  else
    printf("Failure of saving the Bases.\n");
}



static void Output_PAOBases()
{
  static int i,m,mul,L,so;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static char file0[ASIZE8] = ".vpao";
  static char fname[ASIZE8];

  FILE *fp;

  if (Equation_Type==0 || Equation_Type==1){

    sprintf(fname,"%s%s%s",filepath,filename,file0);

    if ((fp = fopen(fname,"w")) != NULL){

      Out_XMIN = Grid_Xmin; 
      Out_XMAX = Grid_Xmax;
      dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);

      for (i=0; i<Grid_Num_Output; i++){
	x = Out_XMIN + dx0*(double)i;
	r = exp(x);
	fprintf(fp,"%15.12f %15.12f  ",x,r);

        /* in case of GVPS */
        if (VPP_switch==4){

	  for (L=0; L<=GVPS_MaxL; L++){
	    for (mul=0; mul<GVPS_ProNum; mul++){
	      fprintf(fp,"%15.12f ", HokanF(r,WFPS[0][L][mul],1));
	    }
	  }
        }

        /* else */
        else{
  	  for (m=0; m<Number_VPS; m++){
	    fprintf(fp,"%15.12f ",HokanF(r,W2[0][m],1)); 
	  }
	}

	fprintf(fp,"\n"); 
      }
      fclose(fp);
      printf("  %s\n",fname);
    }
    else
      printf("Failure of saving the Bases.\n");

  }

  else if (Equation_Type==2){

    for (so=0; so<SOI_switch; so++){

      sprintf(fname,"%s%s%s%i",filepath,filename,file0,so);

      if ((fp = fopen(fname,"w")) != NULL){

	Out_XMIN = Grid_Xmin; 
	Out_XMAX = Grid_Xmax;
	dx0 = (Out_XMAX - Out_XMIN)/(double)(Grid_Num_Output-1);

	for (i=0; i<Grid_Num_Output; i++){
	  x = Out_XMIN + dx0*(double)i;
	  r = exp(x);
	  fprintf(fp,"%15.12f %15.12f  ",x,r);

          /* in case of GVPS */
          if (VPP_switch==4){

	    for (L=0; L<=GVPS_MaxL; L++){
	      for (mul=0; mul<GVPS_ProNum; mul++){
	        fprintf(fp,"%15.12f ", HokanF(r,WFPS[so][L][mul],1));
	      }
	    }
	  }

          /* else */
          else{
  	    for (m=0; m<Number_VPS; m++){
	      fprintf(fp,"%15.12f ",HokanF(r,W2[so][m],1)); 
	    }
	  }

	  fprintf(fp,"\n"); 
	}
	fclose(fp);
	printf("  %s\n",fname);
      }
      else
	printf("Failure of saving the Bases.\n");

    }

  }

}

static void Output_PAOBases2(char *filein)
{
  static int i,GL,NumL;
  static double dx0,x,r,Out_XMIN,Out_XMAX;
  static double scaleF;
  static char file0[ASIZE8] = ".pao";
  FILE *fp;

  fnjoint(filepath,filename,file0);

  if ((fp = fopen(file0,"w")) != NULL){

    Out_Std(fp,filein); 

    fprintf(fp,"\n");
    fprintf(fp,"***************************************************\n");
    fprintf(fp,"  Eigen values(Hartree) of pseudo atomic orbitals  \n");
    fprintf(fp,"***************************************************\n\n");

    fprintf(fp,"Eigenvalues\n");
    fprintf(fp,"Lmax=%2i Mul=%2i\n",MaxL_PAO,Num_PAO);
    for (GL=0; GL<=MaxL_PAO; GL++){
      for (NumL=0; NumL<Num_PAO; NumL++){
        fprintf(fp," l mu  %2i %2i   %18.14f\n",GL,NumL,E_PAO[GL][NumL]);
      }
    }

    fprintf(fp,"\n");
    fprintf(fp,"\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"*********** Charge density of valence electrons ***********\n");
    fprintf(fp,"***********************************************************\n");

    if (AtomNum!=0)
      scaleF = ( (double)AtomNum - (total_electron - valence_electron) )/valence_electron;
    else 
      scaleF = 1.0;

    Density_V(0,OcpN0);
    Out_XMIN = Grid_Xmin; 
    Out_XMAX = log(PAO_Rcut+1.0);
    dx0 = (Out_XMAX - Out_XMIN)/Grid_Num_Output;

    fprintf(fp,"<valence.charge.density\n");
    for (i=0; i<Grid_Num_Output; i++){
      x = Out_XMIN + dx0*(double)i;
      r = exp(x);
      fprintf(fp,"%18.15f %18.15f",x,r);
      fprintf(fp,"  %18.15f\n",fabs(Frho_V(0,r))*scaleF);
    }
    fprintf(fp,"valence.charge.density>\n");

    fprintf(fp,"\n");
    fprintf(fp,"\n");
    fprintf(fp,"***********************************************************\n");
    fprintf(fp,"********  DATA for multiple pseudo atomic orbitals  *******\n");
    fprintf(fp,"***********************************************************\n");

    fprintf(fp,"PAO.Lmax  %5d\n",MaxL_PAO);
    fprintf(fp,"PAO.Mul   %5d\n",Num_PAO);

    Out_XMIN = Grid_Xmin; 
    Out_XMAX = log(PAO_Rcut+1.0);
    dx0 = (Out_XMAX - Out_XMIN)/Grid_Num_Output;

    for (GL=0; GL<=MaxL_PAO; GL++){
      fprintf(fp,"<pseudo.atomic.orbitals.L=%d\n",GL);
      for (i=0; i<Grid_Num_Output; i++){
        x = Out_XMIN + dx0*(double)i;
        r = exp(x);
        fprintf(fp,"%18.15f %18.15f ",x,r);
        for (NumL=0; NumL<Num_PAO; NumL++){
          fprintf(fp,"%18.15f ",MPAO_RadialF(GL,NumL,r));
        }
        fprintf(fp,"\n");
      }
      fprintf(fp,"pseudo.atomic.orbitals.L=%d>\n",GL);
    }
    fclose(fp);
    printf("  %s\n",file0);
  }
  else
    printf("Failure of saving the Bases.\n");
}

static void Output_Density()
{
  static int i,j,num,n,l,m,so;
  static double r,sum,p,p2;
  static double rho_V_tmp[ASIZE1];
  static char file0[ASIZE8];
  FILE *fp;

  /**************************************
    read a file for EDPP
  **************************************/
  
  if (VPP_switch==3) {
    i = Restart_load(0);
    if (i==0){
      printf("Could not find a file which stores informations of AE calculation.\n");
      exit(0);
    }  
  }

  /* file extention */
  if (Calc_Type==0) strcpy(file0,".aden");
  if (Calc_Type==1){
    strcpy(file0,".vden");

    for (i=0; i<Grid_Num; i++){
      r = MRV[i];
      sum = 0.0;

      for (so=0; so<SOI_switch; so++){
        for (m=0; m<Number_VPS; m++){
          n = NVPS[m];
          l = LVPS[m];
          p = PF[so][n][l][i];
          p2 = p*p;
          sum = sum + OcpN[0][so][n][l]*p2/4.0/PI/r/r;
        }
      }

      rho_V_tmp[i] = sum;
    }
  }

  num = 100;

  fnjoint(filepath,filename,file0);

  if ((fp = fopen(file0,"w")) != NULL){

    /* write header */
    if (Calc_Type==0) {
       fprintf(fp,"# density of all electrons\n");
       fprintf(fp,"# MXV  MRV  rho rho*4 pi r^2\n");
    }
    if (Calc_Type==1){
       fprintf(fp,"# density of electrons\n");
       if (PCC_switch==0) {
          fprintf(fp,"# MXV  MRV  rho_V  rho  rho_Core (*)* 4 pi r^2\n");
       }
       else {
          fprintf(fp,"# MXV  MRV  rho_V  rho  rho_Core  rho_PCC (*)* 4 pi r^2\n");
       }
    }

    j = 0; 
    for (i=0; i<Grid_Num; i+=Grid_Num/num){
      j++;
      p = 4.0*PI*MRV[i]*MRV[i];
      if (j<=num){
        fprintf(fp,"%17.14f %17.14f",MXV[i],MRV[i]);
        if (Calc_Type==0) fprintf(fp,"  %17.14f %17.14f\n",rho[0][i],rho[0][i]*p);
        if (Calc_Type==1){
          if (PCC_switch==0)
            fprintf(fp,"  %17.14f %17.14f %17.14f %17.14f %17.14f %17.14f\n",
                    fabs(rho_V[0][i]),fabs(rho[0][i]),fabs(rho[0][i]-rho_V_tmp[i]),
                    fabs(rho_V[0][i]*p),fabs(rho[0][i]*p),fabs(rho[0][i]-rho_V_tmp[i])*p);

          else 
            fprintf(fp,"  %17.14f %17.14f %17.14f %17.14f %17.14f %17.14f %17.14f %17.14f\n",
                    fabs(rho_V[0][i]),fabs(rho[0][i]),fabs(rho[0][i]-rho_V_tmp[i]),fabs(rho_PCC[i]),
              fabs(rho_V[0][i]*p),fabs(rho[0][i]*p),fabs(rho[0][i]-rho_V_tmp[i])*p,fabs(rho_PCC[i])*p);

	}
      }
    }
    fclose(fp);
    printf("  %s\n",file0);
  }
  else
    printf("Failure of saving the Density.\n");

}


static void Output_LogD_RadF()
{
  static int i,L,so;
  static char file_exn[ASIZE8] = ".ld";
  static char file0[ASIZE8];
  FILE *fp;

  for (so=0; so<SOI_switch; so++){
    for (L=0; L<ASIZE2; L++){
      if (NumVPS_L[L]!=0){

        if (SOI_switch==1)
  	  sprintf(file0,"%s%i",file_exn,L);
        else if (SOI_switch==2)
  	  sprintf(file0,"%s%i_%i",file_exn,so,L);

	fnjoint(filepath,filename,file0);
	if ((fp = fopen(file0,"w")) != NULL){
	  for (i=0; i<LogD_num; i++){
	    fprintf(fp,"%17.14f %17.14f %17.14f %17.14f\n",
		    LogDE[i],
                    LogDF[so][L][0][i],LogDF[so][L][1][i],LogDF[so][L][2][i]);
	  }
	  fclose(fp);
	  printf("  %s\n",file0);
	}
	else
	  printf("Failure of saving the LogD_RadF.\n");
      }
    }
  }
}






