/**********************************************************************
  readfile.c:

     readfile.c is a subroutine to read an input file.

  Log of readfile.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include "adpack.h"
#include "Inputtools.h"

static void Input_std(char *file);

void readfile(char *argv[])
{ 

  FILE *fp;

  if ((fp = fopen(argv[1],"r")) != NULL){
    Input_std(argv[1]);
    fclose(fp);
  }
  else{
    printf("\n**** Error message in the inputfile ****\n");
    printf("Failure of reading the input file.\n");
    exit(0);
  }
}


void Input_std(char *file)
{
  FILE *fp;
  int i,n,l,j,k,Nmax,Lmax,eq,so;
  int po=0; /* error count */
  int i_vec[20];
  char *s_vec[20];
  char c; 
  double sum0,sum1;
  char buf[500];
  int Lrestartfile; 

  printf("\n*******************************************************\n"); 
  printf("*******************************************************\n"); 
  printf(" Welcome to ADPACK                                     \n"); 
  printf(" Copyright (C), 2002-2007, T.Ozaki                     \n"); 
  printf(" ADPACK comes with ABSOLUTELY NO WARRANTY.             \n"); 
  printf(" This is a free software, and you are welcome to       \n"); 
  printf(" redistribute it under the constitution of the GNU-GPL.\n");
  printf("*******************************************************\n"); 
  printf("*******************************************************\n\n"); 

  /****************************************************
                       open a file
  ****************************************************/

  input_open(file);
  input_string("System.CurrentDirectory",filepath,"./");
  input_string("System.Name",filename,"default");

  input_logical("System.UseRestartfile",&Lrestartfile,1);
  if (Lrestartfile) {
    input_string("System.Restartfile",restartfile,"");
  }
  else {
    strcpy(restartfile,"");
  }

  /****************************************************
                     Calculation type
  ****************************************************/

  s_vec[0]="sch";    s_vec[1]="sdirac1"; s_vec[2]="sdirac2";
  s_vec[3]="dirac1"; s_vec[4]="dirac2";  s_vec[5]="sdirac";
  s_vec[6]="dirac";
  i_vec[0]=0;        i_vec[1]=1;         i_vec[2]=2;
  i_vec[3]=3;        i_vec[4]=4;         i_vec[5]=2;
  i_vec[6]=4;

  input_string2int("eq.type",&eq, 7, s_vec,i_vec);

  switch (eq){
    case 0: { Equation_Type=0; TwoComp_frag=0; } break;
    case 1: { Equation_Type=1; TwoComp_frag=0; } break;
    case 2: { Equation_Type=1; TwoComp_frag=1; } break;
    case 3: { Equation_Type=2; TwoComp_frag=0; } break;
    case 4: { Equation_Type=2; TwoComp_frag=1; } break;
  }

  if (Equation_Type==2)  SOI_switch = 2;
  else                   SOI_switch = 1;

  s_vec[0]="ALL"; s_vec[1]="VPS"; s_vec[2]="PAO"; s_vec[3]="FCORE"; s_vec[4]="FEMLDA"; s_vec[5]="FEMHF"; s_vec[6]="FEMLDA0"; s_vec[7]="ALLFEM";
  i_vec[0]=0    ; i_vec[1]=1    ; i_vec[2]=2    ; i_vec[3]=3      ; i_vec[4]=4;        i_vec[5]=5;       i_vec[6]=6;       i_vec[7]=7; 
  input_string2int("calc.type",&Calc_Type, 8, s_vec,i_vec);
  if (Calc_Type==7) Calc_Type = 4;

  s_vec[0]="LDA"; s_vec[1]="GGA"; s_vec[2]="LDA-VWN";  
  i_vec[0]=0    ; i_vec[1]=1    ; i_vec[2]=2    ;
  input_string2int("xc.type",&XC_switch,3, s_vec,i_vec);

  /****************************************************
                          Atom
  ****************************************************/

  input_double("AtomSpecies",&AtomNum,0.0);

  if (AtomNum<0 || 104<AtomNum){
    printf("\n**** Error message in the inputfile ****\n");
    printf("!! Wrong AtomSpecies !!\n"); 
    po++;
  }

  input_int("max.occupied.N",&max_ocupied_N,0);
	input_int("max.N",&max_N,max_ocupied_N);
  input_double("total.electron",&total_electron,(double)0.0);
  total_electron0 = total_electron;
  input_double("valence.electron",&valence_electron,(double)0.0);

  /* initialize OcpN */
  for (i=0; i<ASIZE15; i++){
    for (so=0; so<2; so++){
      for (n=0; n<=ASIZE3; n++){
        for (l=0; l<ASIZE3; l++){
          OcpN[i][so][n][l] = 0.0;
	}
      }
    }
  }

  if (fp=input_find("<occupied.electrons") ) {
    sum0 = 0.0;
    for (n=1; n<=max_ocupied_N; n++){
      fscanf(fp,"%d",&i);
      if (i!=n){
        printf("\n**** Error message in the inputfile ****\n");
        printf("!! Format error in occupied.electrons !!\n");
        po++;
      } 

      for (l=0; l<n; l++){
        fscanf(fp,"%lf",&OcpN[0][0][n][l]);

        /* find the max L to be occupied  */ 
        if (0.0<OcpN[0][0][n][l] && Occupied_Lmax<l) Occupied_Lmax = l;        

        sum0 = sum0 + OcpN[0][0][n][l];
        if ((double)(4*l+2)<OcpN[0][0][n][l] || OcpN[0][0][n][l]<0.0){
          printf("\n**** Error message in the inputfile ****\n");
          printf("!! Unphysical occupation number !!\n");
          po++;
        }

        if (Equation_Type==2){
          OcpN[0][0][n][l] = 0.50*OcpN[0][0][n][l];
          OcpN[0][1][n][l] = OcpN[0][0][n][l];
	}

      }
    }
    if (!input_last("occupied.electrons>")) {
      /* format error */
      printf("\n**** Error message in the inputfile ****\n");
      printf("!! Format error for occupied.electrons !!\n");
      po++;
    }

    /* set OcpN0 */
    for (i=0; i<ASIZE15; i++){
      for (so=0; so<2; so++){
	for (n=0; n<=ASIZE3; n++){
	  for (l=0; l<ASIZE3; l++){
	    OcpN0[i][so][n][l] = OcpN[i][so][n][l];
	  }
	}
      }
    }

  }

  /* calculation of frozen core effect */
  if (Calc_Type==3){
    if (fp=input_find("<fcore.occupied.electrons") ) {

      sum1 = 0.0;
      for (n=1; n<=max_ocupied_N; n++){
	fscanf(fp,"%d",&i);
	if (i!=n){
	  printf("\n**** Error message in the inputfile ****\n");
	  printf("!! Format error in fcore.occupied.electrons !!\n");
	  po++;
	} 
	for (l=0; l<n; l++){
	  fscanf(fp,"%lf",&OcpN[1][0][n][l]);
	  sum1 = sum1 + OcpN[1][0][n][l];
	  if ((double)(4*l+2)<OcpN[1][0][n][l] || OcpN[1][0][n][l]<0.0){
	    printf("\n**** Error message in the inputfile ****\n");
	    printf("!! Unphysical opupied number !!\n");
	    po++;
	  }

	  if (Equation_Type==2){
	    OcpN[1][0][n][l] = 0.50*OcpN[1][0][n][l];
	    OcpN[1][1][n][l] = OcpN[1][0][n][l];
	  }

	}
      }
      if (!input_last("fcore.occupied.electrons>")) {
	/* format error */
	printf("\n**** Error message in the inputfile ****\n");
	printf("!! Format error for fcore.occupied.electrons !!\n");
	po++;
      }
    }

    total_electron1 = sum1;

    input_int("number.relaxed.states",&Number_Realaxed_States,0);

    if (fp=input_find("<relaxed.NandL") ) {

      Nmax = 0;

      for (i=0; i<Number_Realaxed_States; i++){
	fscanf(fp,"%d %d",&FC_NVPS[i],&FC_LVPS[i]);

	if (Nmax<=FC_NVPS[i]){
	  Nmax = FC_NVPS[i];  
	}
	else{
	  /* format error */
	  printf("\n**** Error message in the inputfile ****\n");
	  printf("!! Please specify 'relaxed.NandL' in ascending\n");
	  printf("!! order for the priciple number n\n");
	  po++;
	}      
      }
      if (!input_last("relaxed.NandL>")){
	/* format error */
	printf("\n**** Error message in the inputfile ****\n");
	printf("!! Format error in 'relaxed.NandL' !!\n");
	po++;
      }
    }

    for (n=1; n<=ASIZE3; n++){
      for (l=0; l<ASIZE3; l++){
        Relaxed_Switch[n][l] = 0;
      }
    }

    for (i=0; i<Number_Realaxed_States; i++){
      n = FC_NVPS[i];
      l = FC_LVPS[i];
      Relaxed_Switch[n][l] = 1;
    }     
  }

  if (10e-13<fabs(total_electron-sum0)){
    printf("\n**** Error message in the inputfile ****\n");
    printf("!! The sum of occupied.electrons is inconsistent with total.electron  !!\n");
    po++;
  } 

  /*****************************************************
   parameters for solving radial differential equations
  *****************************************************/

  input_double("grid.xmin",&Grid_Xmin,(double)-7.0);
  input_double("grid.xmax",&Grid_Xmax,(double)2.5);
  input_int("grid.num",&Grid_Num, 4000);
  if (ASIZE1<SCF_MAX){
    printf("\n**** Error message in the inputfile ****\n");
    printf("ASIZE1<grid.num\n");
    po++;
  } 
  input_int("grid.num.output",&Grid_Num_Output, 2000);

  /****************************************************
                          SCF
  ****************************************************/

  input_int("scf.maxIter",&SCF_MAX, 40);
  if (ASIZE10<SCF_MAX){
    printf("\n**** Error message in the inputfile ****\n");
    printf("ASIZE10<scf.maxIter\n");
    po++;
  }
  s_vec[0]="Simple"; s_vec[1]="GR-Pulay"; s_vec[2]="Pulay"; 
  i_vec[0]=0;        i_vec[1]=1;          i_vec[2]=2; 
  input_string2int("scf.Mixing.Type",&Mixing_switch,3, s_vec,i_vec);
  input_double("scf.Init.Mixing.Weight",&Mixing_weight,(double)0.3);
  Mixing_weight_init = Mixing_weight; 
  input_double("scf.Min.Mixing.Weight",&Min_Mixing_weight,(double)0.001);
  input_double("scf.Max.Mixing.Weight",&Max_Mixing_weight,(double)0.8);
  input_int("scf.Mixing.History",&Num_Mixing_pDM,5);
  input_int("scf.Mixing.StartPulay",&Pulay_SCF,6);
  input_double("scf.criterion",&SCF_criterion,(double)1.0e-9);

  /****************************************************
                     Pseudopotential
  ****************************************************/

  s_vec[0]="BHS"; s_vec[1]="TM"; s_vec[2]="MR"; s_vec[3]="EDPP"; s_vec[4]="GVPS"; s_vec[5]="MBK";
  i_vec[0]=0    ; i_vec[1]=1   ; i_vec[2]=2   ; i_vec[3]=3   ;   i_vec[4]=4   ;   i_vec[5]=5   ; 
  input_string2int("vps.type",&VPP_switch,6, s_vec,i_vec);

  input_int("gvps.maxL",&GVPS_MaxL,3);
  if (ASIZE2<=GVPS_MaxL){
    printf("gvps.maxL should be smaller than ASIZE2.\n"); 
    exit(0);
  }

  input_int("gvps.projector.num",&GVPS_ProNum,3);
  input_double("gvps.cutoffR",&GVPS_Rcut,(double)15.0);

  /* in case of GVPS, correct Grid_Xmax */
  if (VPP_switch==4)  Grid_Xmax = log(GVPS_Rcut);

  input_int("number.vps",&Number_VPS,0);

  if (fp=input_find("<pseudo.NandL") ) {
    Vlocal_maxcut = 0.0; 
    Vlocal_mincut = 1000.0; 

    for (l=0; l<ASIZE2; l++) VPS_RcutL[l] = -1.0;

    Nmax = 0;
    sum0 = 0.0; 

    for (i=0; i<Number_VPS; i++){

      fscanf(fp,"%d %d %d %lf %lf",&j,&NVPS[i],&LVPS[i],&VPS_Rcut[i],&VPS_ene[i]);

      sum0 += OcpN[0][0][NVPS[i]][LVPS[i]]; 

      if (Vlocal_maxcut<VPS_Rcut[i]) Vlocal_maxcut = VPS_Rcut[i];
      if (VPS_Rcut[i]<Vlocal_mincut) Vlocal_mincut = VPS_Rcut[i];

      if (VPS_RcutL[LVPS[i]]<0.0) VPS_RcutL[LVPS[i]] = VPS_Rcut[i];

      if (Nmax<=NVPS[i]){
        Nmax = NVPS[i];  
      }
      else{
        /* format error */
        printf("\n**** Error message in the inputfile ****\n");
        printf("!! Please specify 'pseudo.NandL' in ascending\n");
        printf("!! order for the priciple number n\n");
        po++;
      }      
    }

    if (Equation_Type==2) sum0 = 2.0*sum0;

    if (!input_last("pseudo.NandL>")){
      /* format error */
      printf("\n**** Error message in the inputfile ****\n");
      printf("!! Format error in 'pseudo.NandL' !!\n");
      po++;
    }
  }

  if (1.0e-12<fabs(valence_electron-sum0)){
    printf("\n**** Error message in the inputfile ****\n");
    printf("!! valence.electron is inconsistent. !!\n");
    po++;
  }

  input_int("MR.excited.states.num",&Num_MR_ES,1);
  if (Num_MR_ES<0){
    /* format error */
    printf("\n**** Error message in the inputfile ****\n");
    printf("!! MR.excited.states.num must be equal to or greater than 0 !!\n");
    po++;
  }
  
  input_int("Blochl.projector.num",&Blochl_pro_num,1);
  if (Blochl_pro_num==0){
    /* format error */
    printf("\n**** Error message in the inputfile ****\n");
    printf("!! Blochl.projector.num must be equal to or greater than 1 !!\n");
    po++;
  }
  
  s_vec[0]="Simple"; s_vec[1]="Polynomial";
  i_vec[0]=0    ; i_vec[1]=1   ;
  input_string2int("local.type",&Vlocal_switch, 2, s_vec,i_vec);

  if (Vlocal_switch==0 && VPP_switch==5){
    printf("Only the polynomial local part is available for GVPS.\n");
    exit(0);
  }

  input_int("local.part.vps",&Local_Part_VPS,0);
  if ( Vlocal_switch==0 &&  (Local_Part_VPS<0 || (Number_VPS-1)<Local_Part_VPS) ){
     printf("\n**** Error message in the inputfile ****\n");
     printf("!! Invalid local.part.vps !!\n");
     po++;
  }    
  input_double("local.cutoff",&Vlocal_cutoff,(double)Vlocal_mincut);

  input_double("local.origin.ratio",&Vlocal_origin,(double)3.0);

  /*
  if (AtomNum!=0){
    for (i=0; i<Number_VPS; i++){
      n = NVPS[i];
      l = LVPS[i];

      if (OcpN[0][0][n][l]<10e-12){
	printf("\n**** Error message in the inputfile ****\n");
	printf("!! The ocupancy of the orbital n=%i l=%i which is\n",n,l);
	printf("   specified to make pseudo potentials is zero !!\n");
	printf("   Please provide a finite ocupancy for this orbital\n");
	po++;
      }    
    }
  }
  */
  
  input_logical("log.deri.RadF.calc",&LogD_switch,0);
  input_double("log.deri.MinE",&LogD_MinE,(double)-3.0);
  input_double("log.deri.MaxE",&LogD_MaxE,(double)2.0);
  input_int("log.deri.num",&LogD_num,(int)50);
  input_logical("ghost.check",&ghost_check,0);

  /****************************************************
   For VLESP
  ****************************************************/
 
  /*
  charge_state_num = 3;
  charge_states[0] = 0.0;
  charge_states[1] = -0.5;
  charge_states[2] = 0.5;
  */

  input_int("charge.states.num",&charge_state_num,(int)1);
  if (Calc_Type==3) charge_state_num = 2;

  if (ASIZE15<charge_state_num){
    printf("\n**** Error message in the inputfile ****\n");
    printf("ASIZE15<charge.states.num, please use a large value for ASIZE15.\n");
    po++; 
  }

  if (fp=input_find("<charge.states") ) {
    for (i=0; i<charge_state_num; i++){
      fscanf(fp,"%lf",&charge_states[i]);
    }

    if (!input_last("charge.states>")) {
      /* format error */
      printf("\n**** Error message in the inputfile ****\n");
      printf("!! Format error in charge.states !!\n");
      po++;
    }
  }

  /****************************************************
   Check whether multiple states for each l-component
                   are included or not.
  ****************************************************/
  
  for (i=0; i<ASIZE2; i++) NumVPS_L[i] = 0;
  
  for (i=0; i<Number_VPS; i++){
    l = LVPS[i];
    GI_VPS[l][NumVPS_L[l]] = i;
    NumVPS_L[l]++;
  }

  /****************************************************
                       log.deri.R
  ****************************************************/

  if (  LogD_switch && SOI_switch==1 ) {

    if (fp=input_find("<log.deri.R") ) {
      for (l=0; l<ASIZE2; l++){
        if (NumVPS_L[l]!=0){
          fscanf(fp,"%d %lf",&j,&LogD_R[l]);
          if (j!=l){
            printf("\n**** Error message in the inputfile ****\n");
            printf("!! Format error for log.deri.R !!\n");
            po++;
          }
        }
      }
      if (!input_last("log.deri.R>")){
        /* format error */
        printf("\n**** Error message in the inputfile ****\n");
        printf("!! Format error for log.deri.R !!\n");
        po++;
      }
    }
  }

  else if ( LogD_switch && SOI_switch==2 ) {

    if (fp=input_find("<log.deri.R") ) {
      for (l=0; l<ASIZE2; l++){
        if (NumVPS_L[l]!=0){
          fscanf(fp,"%d %lf %lf",&j,&LogD_R[l],&LogD1_R[l]);
          if (j!=l){
            printf("\n**** Error message in the inputfile ****\n");
            printf("!! Format error for log.deri.R !!\n");
            po++;
          }
        }
      }
      if (!input_last("log.deri.R>")){
        /* format error */
        printf("\n**** Error message in the inputfile ****\n");
        printf("!! Format error for log.deri.R !!\n");
        po++;
      }
    }

  } /* SOI_switch */

  /****************************************************
            factor spin-orbit splitting 
  ****************************************************/

  if (Equation_Type==2){

    if (fp=input_find("<SO.factor") ) {
      for (l=0; l<ASIZE2; l++){
        if (NumVPS_L[l]!=0){
          fscanf(fp,"%d %lf",&j,&SO_factor[l]);
          if (j!=l){
            printf("\n**** Error message in the inputfile ****\n");
            printf("  Format error in SO.factor\n");
            po++;
          }
        }
      }
      if (!input_last("SO.factor>")){
        /* format error */
        printf("\n**** Error message in the inputfile ****\n");
        printf(" Format error in SO.factor\n");
        po++;
      }
    }

    else{
      /* set default */
      for (l=0; l<ASIZE2; l++){
        SO_factor[l] = 1.0;
      }
    }
  }
    
  /****************************************************
                 Pseudo atomic orbitals
  ****************************************************/

  input_int("maxL.pao",&MaxL_PAO,2);
  input_int("num.pao",&Num_PAO,7);
  input_double("radial.cutoff.pao",&PAO_Rcut,(double)5.0);
  input_double("height.of.wall",&height_wall,(double)2000.0);
  input_double("rising.edge",&rising_edge,(double)0.5);
  input_double("search.LowerE",&search_LowerE,(double)-3.000);
  input_double("search.UpperE",&search_UpperE,(double)20.000);
  input_int("num.of.partition",&Num_Partition,(int)300);
  input_double("matching.point.ratio",&MatP_ratio,(double)0.67);

  /* in case of PAO calc., correct Grid_Xmax */
  if (Calc_Type==2 || Calc_Type==0)  Grid_Xmax = log(PAO_Rcut + 2.0);

  /****************************************************
    Core electron density for partial core correction
  ****************************************************/

  input_logical("charge.pcc.calc",&PCC_switch,0);
  input_double("pcc.ratio",&PCC_Ratio,(double)1.0);
  input_double("pcc.ratio.origin",&PCC_Ratio_Origin,(double)6.0);

  /****************************************************
                          Log
  ****************************************************/

  input_logical("Log.print",&Check_switch,0);

  /****************************************************
                    close the file
  ****************************************************/

  input_close();

  /****************************************************
                     Error check
  ****************************************************/

  /*
  if (Calc_Type!=2 && Grid_Num<2000){
    printf("\n**** Error message in the inputfile ****\n");
    printf("You are requested to use over 2000 for grid.num.\n"); 
    printf("Please use a larger value for grid.num.\n"); 
    po++;
  }
  else if (Calc_Type==2 && Grid_Num<6000){
    printf("\n**** Error message in the inputfile ****\n");
    printf("You are requested to use over 6000 for grid.num.\n"); 
    printf("Please use a larger value for grid.num.\n"); 
    po++;
  }

  else if (Calc_Type==2 && 20000<height_wall){
    printf("\n**** Error message in the inputfile ****\n");
    printf("You are requested to use below 20000 for height.of.wall.\n"); 
    printf("Please use a smaller value for height.of.wall.\n"); 
    po++;
  }
  */

  if (po>0 || input_errorCount()>0) {
    printf("\n**** Error message in the inputfile ****\n");
    printf("Invalid values are specified in your inputfile.\n"); 
    printf("Therefore, the execution is terminated.\n"); 
    printf("Please check your input file.\n");
    exit(0);
  } 

}

