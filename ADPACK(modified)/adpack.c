/**************************************************************************
    Copyright (C), 2002-2007, Taisuke Ozaki.
    The program package follows the terms of GNU GPL. 

  adpack.c:

     adpack.c (Atomic Density functional program PACKage) is a program
     package for atomic density functional calculations.

     This program includes
      1) all electron self-consistent calculations for an atom,
      2) generations of norm conserving pseudo potentials by using
         BHS (Bachelet, Hamann, and Schluter) or TM (Troullier and
         Martins) methods,
      3) calculations of charge density for the partial core correction,
      4) calculations of pseudo atomic orbitals for confinement pseudo
         potentials.  

  Log of adpack.c:

     10/Dec/2002  Released by T.Ozaki

**************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "adpack.h"
#include "Inputtools.h"

#ifdef noomp
#include "mimic_omp.h"
#else
#include <omp.h>
#endif

int main(int argc, char *argv[])
{
  int i,k,po,ip;

  /****************************************************
     set the number of threads for openmp
     ./adpack -nt #
  ****************************************************/

  openmp_threads_num = 1; /* default */

  po = 0;
  for (i=0; i<argc; i++){
    if ( strcmp(argv[i],"-nt")==0 ){
      po = 1;
      ip = i;
    }
  }

  if ( (argc-1)<(ip+1) ){
    printf("cannot find the number of threads\n");
    exit(0);
  }

  if ( po==1 ){

    openmp_threads_num = atoi(argv[ip+1]);

    if (openmp_threads_num<=0){
      printf("check the number of threads\n");
      exit(0);
    }
  }

  omp_set_num_threads(openmp_threads_num);

  printf("\nThe number of threads in each node for OpenMP parallelization is %d.\n\n",openmp_threads_num);

  /****************************************************
               Set of input parameters                
  ****************************************************/

  readfile(argv);

  /****************************************************
                     Initialize
  ****************************************************/

  Set_Init();

  /****************************************************
                All electron calculation
  ****************************************************/

  /* calculation of frozen core */
  if (Calc_Type==3){
    All_Electron(0);
    All_Electron(1);
  }

  /* all electron LDA calculation by FEM: original version */
  else if (Calc_Type==4){
    FEM_All_Electron();
  }
  /* all electron HF calculation by FEM */
  else if (Calc_Type==5){
    FEMHF_All_Electron();
  }
  /* all electron LDA calculation by FEM */
  else if (Calc_Type==6){
    FEMLDA_All_Electron();
  }
  else{
    i = Restart_load(0);
    if ( i==0 ) All_Electron(0);
    if ( (i==0 && strlen(restartfile)) )  Restart_save(0);
  }

  /****************************************************
                   Pseudo potentials
  ****************************************************/

  if ( (Calc_Type==1 || Calc_Type==2) && AtomNum!=0){

    Init_VPS();

    if      (VPP_switch==0)  BHS(0);
    else if (VPP_switch==1)  TM(0);
    else if (VPP_switch==2)  MR(0);
    else if (VPP_switch==5)  MBK(0);

    if ( VPP_switch==0 || VPP_switch==1 || VPP_switch==2 ) Generate_VNL();

    if (Calc_Type==1 && LogD_switch==1) Log_DeriF();
    if (Calc_Type==1 && ghost_check==1) ghost(0);
  }
  else if (AtomNum==0){
    Empty_VPS();
  }

  /****************************************************
                   Pseudo atomic orbitals
  ****************************************************/
	 if (Calc_Type==2) Multiple_PAO(0);
	//if (Calc_Type==2 || Calc_Type==0) Multiple_PAO(0);

	/* All-electron atomic orbitals */
	 if(Calc_Type==0) All_Electron_NSCF(0);
		 
  
  /****************************************************
                        Output 
  ****************************************************/

  Output(argv[1]);

  /****************************************************
                     a final message
  ****************************************************/

  printf("\nThe calculation was completed normally.\n");

}



