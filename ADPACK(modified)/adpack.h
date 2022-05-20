/**********************************************************************
  adpack.h:

     adpack.h is a header file of adpack.

  Log of adpack.h:

     22/Nov/2001  Released by T.Ozaki
 
***********************************************************************/
#include <sys/types.h>
 
#define PI      3.1415926535897932384626433832795028841971693993751L
#define ASIZE1   20000   /* max number of grids                           */
#define ASIZE2      20   /* max number of L                               */
#define ASIZE3      20   /* max number of N                               */
#define ASIZE4      15   /* max number of pseudopotentials                */
#define ASIZE5      22   /* max Number of multiple bases                  */
#define ASIZE6     100   /* matrix size in Gauss_LEQ                      */
#define ASIZE7    1000   /* number of grids for Gauss-Legendre            */
#define ASIZE8     100   /* max size of character                         */
#define ASIZE9      20   /* max number of the history densities           */
#define ASIZE10    2000  /* max SCF iterations                            */
#define ASIZE11   15000   /* max number of energy mesh in searching energy */
#define ASIZE12   15000   /* max number of energy mesh in searching energy */
#define ASIZE13      2   /* for kappa of spin-orbit coupling              */
#define ASIZE14      6   /* max L of generated VPS by MR method           */
#define ASIZE15      3   /* max number of reference states in EDPP method */
#define LimitE  1.0e+303 /* limited value in your machine                 */
   
char filepath[ASIZE8],filename[ASIZE8];
char restartfile[ASIZE8];

int Grid_Num_Output,Local_Part_VPS,Number_VPS;
int Grid_Num,openmp_threads_num;;
int MaxL_PAO,Num_PAO,node,CTP,Num_Mixing_pDM,Mixing_switch;
int LL_grid,UL_grid,max_ocupied_N,Occupied_Lmax,PCC_switch,VPP_switch,TwoComp_frag;
int max_N;
int Calc_Type,Equation_Type,SOI_switch,MatP[ASIZE13][ASIZE3][ASIZE2];
int Check_switch,Basis_switch,Multiple_switch,SCF_MAX,XC_switch;
int Number_Realaxed_States;
int Pulay_SCF,MaxO_LF,NVPS[ASIZE4],LVPS[ASIZE4],MCP[ASIZE13][ASIZE3][ASIZE2];
int FC_NVPS[ASIZE4],FC_LVPS[ASIZE4];
int Relaxed_Switch[ASIZE3+1][ASIZE3];
int NumVPS_L[ASIZE2],GI_VPS[ASIZE2][ASIZE4],projector_num[ASIZE2];
int SCF_END,Num_Partition,Vlocal_switch,Blochl_pro_num,total_pro_num;
int NL_type[ASIZE2],LogD_switch,LogD_num,ghost_check;
int Num_MR_ES,charge_state_num,num_EDPP_grid;
int GVPS_MaxL,GVPS_ProNum;
/* manual mode to find PAO, (kino) */
int fix_MatP[ASIZE2][ASIZE5];

double AtomNum;  
double dx,Vinf,Grid_Xmin,Grid_Xmax;
double valence_electron,total_electron;
double total_electron0,total_electron1;
double Ekin[ASIZE15],EHart[ASIZE15],Eeigen[ASIZE15];
double Exc[ASIZE15],Eec[ASIZE15],Etot[ASIZE15];
double Ekin_VPS[ASIZE15],EHart_VPS[ASIZE15],Eeigen_VPS[ASIZE15];
double Exc_VPS[ASIZE15],Eec_VPS[ASIZE15],Etot_VPS[ASIZE15];
double Enl_VPS[ASIZE15],Etot0_VPS[ASIZE15],Eextra[ASIZE15],Eextra2[ASIZE15];
double LogD_MinE,LogD_MaxE,LogD_R[ASIZE2];
double LogD1_R[ASIZE2];
double SO_factor[ASIZE2],VPS_ene[ASIZE4];
double VPS_Rcut[ASIZE4],VPS_RcutL[ASIZE2];
double PAO_Rcut,PCC_Ratio,PCC_Ratio_Origin;
double PF[ASIZE13][ASIZE3][ASIZE2][ASIZE1];
double FF[ASIZE13][ASIZE3][ASIZE2][ASIZE1];
double MRV[ASIZE1],MXV[ASIZE1];
double Mixing_weight,Mixing_weight_init,SCF_criterion;
double Min_Mixing_weight,Max_Mixing_weight;
double search_UpperE,search_LowerE,MatP_ratio;
double Vlocal_cutoff,Vlocal_origin,Vlocal_maxcut,Vlocal_mincut;
double Gaussian_Alp,Av_CoreR,GVPS_Rcut;
double OcpN[ASIZE15][2][ASIZE3+1][ASIZE3];
double OcpN0[ASIZE15][2][ASIZE3+1][ASIZE3];
double Energy_Kin[ASIZE15][2][ASIZE3+1][ASIZE3];
double MPAO[ASIZE2][ASIZE5][ASIZE1];
double E_PAO[ASIZE2][ASIZE5],Vxc[ASIZE1],rho[ASIZE9][ASIZE1];
double Rrho[ASIZE9][ASIZE1],x[ASIZE7],w[ASIZE7];
double V[ASIZE1],V_all_ele[ASIZE1];
double Vcore[ASIZE1],Vh[ASIZE1],History_Uele[ASIZE9],rho_PCC[ASIZE1];
double V1[ASIZE2][ASIZE1],VPS[ASIZE15][ASIZE13][ASIZE2][ASIZE1];
double VNL[ASIZE15][ASIZE13][ASIZE2][ASIZE4][ASIZE1];
double VNL_W2[ASIZE13][ASIZE2][ASIZE4][ASIZE1];
double Vlocal[ASIZE1],E[ASIZE15][ASIZE13][ASIZE3][ASIZE2];
double Prho[ASIZE1],rho_V[ASIZE15][ASIZE1];
double NormRD[ASIZE9],Vxc_V[ASIZE1],Vh_V[ASIZE1],V2[ASIZE13][ASIZE4][ASIZE1];
double W2[ASIZE13][ASIZE4][ASIZE1],HisNormRD[ASIZE10],HisEeigen[ASIZE10];
double ProjectEnl[ASIZE4],height_wall,rising_edge,search_stepE;
double RNG[ASIZE3+1][ASIZE3],proj_ene[ASIZE15][ASIZE13][ASIZE2][ASIZE4];
double LogDE[ASIZE12],LogDF[ASIZE13][ASIZE2][3][ASIZE12];
double LogD_dif[ASIZE13][ASIZE2][2];
double charge_states[ASIZE15],Core_R[ASIZE15];
double Vedpp[ASIZE15][ASIZE1];
double Int_EDPP_Proj[ASIZE15];
double Evps[ASIZE15][ASIZE13][ASIZE4];
double Evps0[ASIZE15][ASIZE13][ASIZE4];
double CoesQ[30];
double CoesEZ[30];
double E_EDPP[ASIZE1];
double Q_EDPP[ASIZE1];
double ****GVPS;
double ****WFPS;
double ***TMVPS;

double PAO_RadialF(int l, double R, double RadF[ASIZE2][ASIZE1]);
double HokanF(double R, double RadF[ASIZE1], int waru_switch);
double MPAO_RadialF(int l, int mul, double R);
double VNLF(int so, int l, int mul, double R);
double Frho_V(int state_num, double R);

void readfile(char *argv[]);
void Set_Init();
void Initial_Density();
void GR_Pulay(int state_num, int SCF_iter);
void Simple_Mixing(int state_num,
                   double rho0[ASIZE1],
                   double rho1[ASIZE1],
                   double rho2[ASIZE1],
                   double Rrho2[ASIZE1]);
void All_Electron(int state_num);
void FEM_All_Electron();
void FEMLDA_All_Electron();
void FEMHF_All_Electron();
void Total_Energy(int state_num);
void Density(int state_num);
void Density_V(int state_num, double OcpN0[ASIZE15][2][ASIZE3+1][ASIZE3]);
void Density_PCC(double rho_PCC[ASIZE1], int state_num, double OcpN0[ASIZE15][2][ASIZE3+1][ASIZE3]);
void Core(int SCF, double Z, double Vcore[ASIZE1]);
void Hartree(double rh[ASIZE1], double h[ASIZE1]);
void XC_Xa(double rh[ASIZE1], double xc[ASIZE1]);
void XC_CA(double rh[ASIZE1], double xc[ASIZE1], int P_switch);
void XC4atom_PBE(double rh[ASIZE1], double xc[ASIZE1], int P_switch);
void XC_VWN(double rh[ASIZE1], double xc[ASIZE1], int P_switch);
void XC_PBE(double dens[2], double GDENS[3][2], double Exc[2],
            double DEXDD[2], double DECDD[2],
            double DEXDGD[3][2], double DECDGD[3][2]);
void XC_PW91C(double dens[2], double Ec[1], double Vc[2]);
void XC_EX(int NSP, double DS0, double DS[2], double EX[1], double VX[2]);
void VP();
void Hamming_O(int eq_type,
               int NL, double ep, double kappa,
               double Mo[ASIZE1],  double Lo[ASIZE1], double DM[ASIZE1], 
               double VPL[ASIZE1], double VNLp[ASIZE1],
               double Reduce_Num_grid);
void Hamming_I(int eq_type, int NL, double ep, double kappa,
               double Mi[ASIZE1],  double Li[ASIZE1], double DM[ASIZE1],
               double VPL[ASIZE1], double VNLp[ASIZE1],
               double Reduce_Num_grid);
void Gauss_LEQ(int n, double a[ASIZE6][ASIZE6], double x[ASIZE6]);
void Gauss_Legendre(int n, double x[ASIZE7], double w[ASIZE7],
                    int *ncof, int *flag);
void BHS(int state_num);
void TM(int state_num);
void MBK(int state_num);
void MBK2(int state_num);
void MR();
void Multiple_PAO();
// void E_NL(int NumVPS);
void Calc_Vlocal(double Vlocal[ASIZE1], double *V0, int local_switch);
void Init_VPS();
void Generate_VNL();
void Log_DeriF();
void ghost(int state_num);
void Output(char *filein);
void Restart_save(int state_num);
int Restart_load(int state_num);
void Make_EDPP();
void Find_LESP();
void Empty_VPS();

void qsort_int(long n, int *a, int *b);
void qsort_double(long n, double *a, double *b);

/* functions defined in addfunc.c */
double isgn(int nu);
double sgn(double nu);
double largest(double a, double b);
double smallest(double a, double b);
void fnjoint(char name1[ASIZE8],char name2[ASIZE8],char name3[ASIZE8]);
void chcp(char name1[ASIZE8],char name2[ASIZE8]);


void All_Electron_NSCF(int state_num);
