/**********************************************************************
  Init_VPS.c:

     Init_VPS.c is a subroutine to set parameters for generating
     pseudo potentials.

  Log of Init_VPS.c:

     16/Mar/2003  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "adpack.h"

void Init_VPS()
{ 
  
  static int L,i,gi,po,m;
  
  /****************************************************
     define how to make a non-local pseudo potential
                   for each L-component

       NL_type=0: none
       NL_type=1: KB
       NL_type=2: Blochl
       NL_type=3: modified Vanderbilt
  ****************************************************/
  
  for (i=0; i<ASIZE2; i++){
    NumVPS_L[i] = 0;
    projector_num[i] = 0;
    NL_type[i] = 0;
  }

  if (ASIZE2<(ASIZE14+1)){
    printf("ASIZE2 should be larger than (ASIZE14+1).\n");
    exit(0);    
  }

  for (i=0; i<Number_VPS; i++){
    L = LVPS[i];
    GI_VPS[L][NumVPS_L[L]] = i;
    NumVPS_L[L]++;
  }

  if (VPP_switch==5){

    total_pro_num = 0;
    for (L=0; L<ASIZE2; L++){

      if (NumVPS_L[L]!=1){
        projector_num[L] = NumVPS_L[L];
      }
      else if (NumVPS_L[L]==1){
        projector_num[L] = Blochl_pro_num;
      }

      else{
        projector_num[L] = 0;
      }  

      total_pro_num += projector_num[L];
    }
  }

  else {  

    /* loop for angular momentum */

    for (L=0; L<ASIZE2; L++){

      if (NumVPS_L[L]==0){
	projector_num[L] = 0;
      }

      /* Hamman's scheme */

      else if (NumVPS_L[L]==1 && OcpN[0][0][NVPS[GI_VPS[L][0]]][L]<1.0e-12 ){

	gi = GI_VPS[L][0];
	if (Vlocal_switch==1 || Equation_Type==2){
	  projector_num[L] = Blochl_pro_num;
	  NL_type[L] = 5;
	}
	else if (gi!=Local_Part_VPS){
	  projector_num[L] = Blochl_pro_num;
	  NL_type[L] = 5;
	}
	else{
	  projector_num[L] = 0;
	  NL_type[L] = 0;
	}
      }

      else{

	if (VPP_switch==4){
	  projector_num[L] = GVPS_ProNum;
	  NL_type[L] = 4;
	}

	else {

	  /* KB */

	  if (Blochl_pro_num==1){
	    gi = GI_VPS[L][0];
	    if (Vlocal_switch==1 || Equation_Type==2){
	      projector_num[L] = 1;
	      NL_type[L] = 1;
	    }
	    else if (gi!=Local_Part_VPS){
	      projector_num[L] = 1;
	      NL_type[L] = 1;
	    }
	    else{
	      projector_num[L] = 0;
	      NL_type[L] = 0;
	    }
	  }

	  /* Blochl */
	  else{
	    gi = GI_VPS[L][0];
	    if (Vlocal_switch==1 || Equation_Type==2){
	      projector_num[L] = Blochl_pro_num;
	      NL_type[L] = 2;
	    }
	    else if (gi!=Local_Part_VPS){
	      projector_num[L] = Blochl_pro_num;
	      NL_type[L] = 2;
	    }
	    else{
	      projector_num[L] = 0;
	      NL_type[L] = 0;
	    }
	  }
	}

      } 
    }
  
    total_pro_num = 0;
    for (L=0; L<ASIZE2; L++){
      total_pro_num = total_pro_num + projector_num[L];
      if (NL_type[L]==3 && VPP_switch==0){
	printf("In VPS, you specified multiple states for the same L-component.\n");
	printf("In BHS-type, the multiple states are not supported.\n");
	printf("Please use TM-type for the generation of VPS.\n");
	exit(0);
      }
    }

  }

}







