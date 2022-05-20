#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include "adpack.h"

void Restart_save(int state_num)
{
  FILE *fp;
  int i,j,l,so;
  size_t n;
  double di;
  int i_vec[10];
  double d_vec[10];
  char fname[ASIZE8];

  if ( strlen(restartfile)==0 ) {
    return;
  }
  else if (VPP_switch==3) {
    strcpy(fname,filepath);
    sprintf(fname,"%s%s_%i",fname,restartfile,state_num);
  }
  else {
    /* open a file */
    strcpy(fname,filepath);
    strcat(fname,restartfile);
  }

  if (   (fp=fopen(fname,"wb"))==NULL) {
    printf("can not open %s\n",fname);
    return;
  }

  i_vec[0]=1;  /* version */
  i_vec[1]=0;  /* not used, reserved */
  i_vec[2]=0;  /* not used, reserved */
  i_vec[3]=TwoComp_frag;
  i_vec[4]=Equation_Type;
  i_vec[5]=Calc_Type;
  fwrite(i_vec,sizeof(i_vec[0]), 6, fp);

  /*
    Eeigen =     -134.3607269978033
    Ekin   =      240.4271889616195
    EHart  =      112.7698956006297
    Exc    =      -17.4500373009733
    Eec    =     -577.2762390082630
    Etot   = Ekin + EHart + Exc + Eec
    Etot   =     -241.5291917469871
  */

  d_vec[0] = Eeigen[state_num];
  d_vec[1] = Ekin[state_num];
  d_vec[2] = EHart[state_num];
  d_vec[3] = Exc[state_num];
  d_vec[4] = Eec[state_num];
  d_vec[5] = Etot[state_num];
  d_vec[6] = Vinf;
  fwrite(d_vec,sizeof(d_vec[0]), 7, fp);

  fwrite(rho[0],sizeof(rho[0][0]),Grid_Num, fp);

  for (n=1; n<=max_ocupied_N; n++){
    fwrite(RNG[n],sizeof(RNG[0][0]), n,fp);
  }

  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fwrite(E[state_num][so][n],sizeof(E[0][0][0][0]),n,fp);
    }
  }

#if 0
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n;l++) {
	printf("%lf ", E[state_num][so][n][l]);
      }
      printf("\n");
    }
  }
#endif

  /* MCP */
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fwrite(MCP[so][n],sizeof(MCP[0][0][0]),n,fp);
    }
  }

  /* OcpN */
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fwrite(OcpN[state_num][so][n],sizeof(OcpN[0][0][0][0]),n,fp);
    }
  }

  /* V */
  fwrite(Vcore,sizeof(Vcore[0]),Grid_Num,fp);
  fwrite(Vxc,sizeof(Vxc[0]),Grid_Num,fp);
  fwrite(Vh,sizeof(Vh[0]),Grid_Num,fp);
  fwrite(V,sizeof(V[0]),Grid_Num,fp);
  fwrite(V_all_ele,sizeof(V_all_ele[0]),Grid_Num,fp);

  /* FF */ 
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
	fwrite(FF[so][n][l],sizeof(FF[0][0][0][0]), Grid_Num,fp);
      }
    }
  }

  /* PF */ 
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
	fwrite(PF[so][n][l],sizeof(PF[0][0][0][0]), Grid_Num,fp);
      }
    }
  }

  fclose(fp);
}




int Restart_load(int state_num)
{

  FILE *fp;
  int i,j,l,so;
  size_t n;
  double di;
  int i_vec[10];
  double d_vec[10];

  char fname[ASIZE8];
  /* struct stat fbuf; */

  /* open a file */
  if ( strlen(restartfile)==0 ) {
    return 0;
  }
  else if (VPP_switch==3) {
    strcpy(fname,filepath);
    sprintf(fname,"%s%s_%i",fname,restartfile,state_num);
  }
  else {
    /* open a file */
    strcpy(fname,filepath);
    strcat(fname,restartfile);
  }

  /*
    i=stat(fname,&fbuf);
    printf ("%d\n",i);
  */

  if (   (fp=fopen(fname,"rb"))==NULL) {
    /* printf("can not open %s\n",fname); */
    printf("save results to %s after this calculation\n",fname);
    return 0;
  }

  printf("\nload data from %s\n",fname);

#if 0
  i_vec[0]=1;  /* version */
  i_vec[1]=0;  /* not used, reserved */
  i_vec[2]=0;  /* not used, reserved */
  i_vec[3]=TwoComp_frag;
  i_vec[4]=Equation_Type;
  i_vec[5]=Calc_Type;
#endif
  fread(i_vec,sizeof(i_vec[0]), 6, fp);

  /* check consistency of data */
  if ( i_vec[0]!=1 ) {
    printf("version not supported, version in file in %d\n",i_vec[0]);
    exit(0);
  }

  /*************************
      check Calc_Type
  *************************/

  if ( Calc_Type == 0 || Calc_Type==1 ) {
    /* accept Calc_Type==0 || Calc_Type==1 */
    if (i_vec[5] == 0 ||  i_vec[5] == 1) {
    }
    else {
      printf("Calc_Type is different, Calc_Type(%s) = %d\n",fname,i_vec[5]);
      exit(0);
    }
  }
  if ( Calc_Type == 2 ) {
    /* accept Calc_Type= 2  */
    if ( i_vec[5] ==2 ) {
    }
    else {
      printf("Calc_Type is different, Calc_Type(%s) = %d\n",fname,i_vec[5]);
      exit(0);
    }
  }

  /*************************
     check Equation_Type
  *************************/

  if ( Equation_Type!=i_vec[4] ) {
    printf("Equation_Type is different, Equation_Type(%s) = %d\n",fname,i_vec[4]);
    exit(0);
  }

  /*************************
     check TwoComp_frag
  *************************/

  if ( TwoComp_frag!=i_vec[3] ) {
    printf("TwoComp_frag is different, TwoComp_frag(%s) = %d\n",fname,i_vec[3]);
    exit(0);
  }
  
  fread(d_vec,sizeof(d_vec[0]), 7, fp);
  Eeigen[state_num] = d_vec[0];
  Ekin[state_num]   = d_vec[1];
  EHart[state_num]  = d_vec[2];
  Exc[state_num]    = d_vec[3];
  Eec[state_num]    = d_vec[4];
  Etot[state_num]   = d_vec[5];
  Vinf              = d_vec[6];  

  fread(rho[0],sizeof(rho[0][0]),Grid_Num, fp);

  for (n=1; n<=max_ocupied_N; n++){
    fread(RNG[n],sizeof(RNG[0][0]), n,fp);
  }

  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fread(E[state_num][so][n],sizeof(E[0][0][0][0]),n,fp);
    }
  }
#if 0
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n;l++) {
	printf("%lf ", E[state_num][so][n][l]);
      }
      printf("\n");
    }
  }
#endif

  /* MCP */
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fread(MCP[so][n],sizeof(MCP[0][0][0]),n,fp);
    }
  }

  /* OcpN */
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      fread(OcpN[state_num][so][n],sizeof(OcpN[0][0][0][0]),n,fp);
    }
  }

  /* V */
  fread(Vcore,sizeof(Vcore[0]),Grid_Num,fp);
  fread(Vxc,sizeof(Vxc[0]),Grid_Num,fp);
  fread(Vh,sizeof(Vh[0]),Grid_Num,fp);
  fread(V,sizeof(V[0]),Grid_Num,fp);
  fread(V_all_ele,sizeof(V_all_ele[0]),Grid_Num,fp);

  /* FF */ 
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
	fread(FF[so][n][l],sizeof(FF[0][0][0][0]), Grid_Num,fp);
      }
    }
  }

  /* PF */ 
  for (so=0; so<SOI_switch; so++){
    for (n=1; n<=max_ocupied_N; n++){
      for (l=0; l<n; l++){
	fread(PF[so][n][l],sizeof(PF[0][0][0][0]), Grid_Num,fp);
      }
    }
  }


  fclose(fp);

  return 1;

}
