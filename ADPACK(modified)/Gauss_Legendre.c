/**********************************************************************
  Gauss_Legendre.c:

     Gauss_Legendre.c is a subrutine to generate radial grids for
     the Gauss_Legendre quadrature. 

  Log of Gauss_Legendre.c:

     10/Dec/2002  Released by T.Ozaki

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "adpack.h"

#define EPS 1.0e-14

static void Gauss_Legendre1(int n, double x[], double w[],
                            int *ncof, int *flag);
static void Gauss_Legendre2( double x1, double x2,  double x[],  double w[], int n );



void Gauss_Legendre(int n, double x[], double w[], int *ncof, int *flag)
{
  int i; 
const int NLISTMAX=20;
static int firstmessage=1;
static int nlistmax=0;
static int nvalue=0;
static int *nlist=NULL;
static double *xlist=NULL, *wlist=NULL;

  /*
  Gauss_Legendre1(n,x,w,&ncof,&flag);
  */

  if ( nlist==NULL ) {
	 xlist = (double*)malloc(sizeof(double)*n*NLISTMAX);
	 wlist = (double*)malloc(sizeof(double)*n*NLISTMAX);
         nlist = (int*)malloc(sizeof(int)*n*NLISTMAX);
         nlistmax=0; 
	 nvalue=n;
  }
  if (n!=nvalue) {
	  printf("Gauss_Legendre: new mesh(%d)!=old mesh(%d), reset them \n",n,nvalue);
	  if (xlist) free(xlist);
          if (wlist) free(wlist);
          if (nlist) free(nlist);
	  xlist = (double*)malloc(sizeof(double)*n*NLISTMAX);
	  wlist = (double*)malloc(sizeof(double)*n*NLISTMAX);
	  nlist = (int*)malloc(sizeof(int)*n*NLISTMAX);
          nlistmax=0;
          nvalue=n;
  }
  for (i=0;i<nlistmax;i++) {
	  if (nlist[i]==n) {
		  memcpy(x,&x[n*i],n);
		  memcpy(w,&w[n*i],n);
		  return;
	  }
  }

  Gauss_Legendre2(-1.0,1.0,x,w,n);

  nlistmax++;
  if (nlistmax>=NLISTMAX) {
	  nlistmax=NLISTMAX-1;
	  if (firstmessage) {
	  printf("better increasing NLISTMAX in Gauss_Legendre, but continue...\n");
	  firstmessage=0;
	  }
     return;
  }
  i=nlistmax-1;
  nlist[i]=n;
  memcpy(&x[n*i],x,n);
  memcpy(&w[n*i],w,n);
}



void Gauss_Legendre1(int n, double x[], double w[],
                     int *ncof, int *flag)
{

  int pointer[23] = { 1,  2,   4,  6,  9 , 12, 16, 20,  25,
		     30, 36,  42,  49, 56, 64, 72, 82,  94,
		     110, 130, 154, 186, 226}; 
                                                                       
  /*
    Abscissae {x } and weights {w } for the interval (-1,+1). 
    i                i
  */
                                                                       
  double gausscof[274][2] = {
    {0.000000000000000 , 0.00000000000000},
    /* 
      n = 2 , start position 1
    */                                                              
    {0.577350269189626 , 1.000000000000000},         
    /* 
      n = 3 , start position 2
    */                                                               
    {0.774596669241483 , 0.555555555555556},         
    {0.000000000000000 , 0.888888888888889},         
    /* 
      n = 4 , start position 4
    */                                                              
    {0.861136311594053 , 0.347854845137454},         
    {0.339981043584856 , 0.652145154862546},         
    /*
      n = 5 , start position 6
    */                                                           
    {0.906179845938664 , 0.236926885056189},         
    {0.538469310105683 , 0.478628670499366},         
    {0.000000000000000 , 0.568888888888889},         
    /*   
      n = 6 , start position 9
    */                                                            
    {0.932469514203152 , 0.171324492379170},         
    {0.661209386466265 , 0.360761573048139},        
    {0.238619186083197 , 0.467913934572691},        
    /*    
      n = 7 , start position 12
    */                                                          
    {0.949107912342759 , 0.129484966168870},        
    {0.741531185599394 , 0.279705391489277},        
    {0.405845151377397 , 0.381830050505119},        
    {0.000000000000000 , 0.417959183673469},        
    /*
      n = 8 , start position 16
    */                                                           
    {0.960289856497536 , 0.101228536290376},        
    {0.796666477413627 , 0.222381034453374},        
    {0.525532409916329 , 0.313706645877887},        
    {0.183434642495650 , 0.362683783378362},        
    /*
      n = 9 , start position 20
    */                                                            
    {0.968160239507626 , 0.081274388361574},        
    {0.836031107326636 , 0.180648160694857},        
    {0.613371432700590 , 0.260610696402935},        
    {0.324253423403809 , 0.312347077040003},        
    {0.000000000000000 , 0.330239355001260},        
    /*
      n = 10 , startposition 25
    */                                                           
    {0.973906528517172 , 0.066671344308688},        
    {0.865063366688985 , 0.149451349150581},        
    {0.679409568299024 , 0.219086362515982},        
    {0.433395394129247 , 0.269266719309996},        
    {0.148874338981631 , 0.295524224714753},        
    /*
      n = 11 , start position 30
    */                                                         
    {0.978228658146057 , 0.055668567116174},          
    {0.887062599768095 , 0.125580369464905},          
    {0.730152005574049 , 0.186290210927734},          
    {0.519096129206812 , 0.233193764591990},          
    {0.269543155952345 , 0.262804544510247},          
    {0.000000000000000 , 0.272925086777901},          
    /*
      n = 12 , start position 36
    */                                                        
    {0.981560634246719 , 0.047175336386512},          
    {0.904117256370475 , 0.106939325995318},          
    {0.769902674194305 , 0.160078328543346},          
    {0.587317954286617 , 0.203167426723066},          
    {0.367831498998180 , 0.233492536538355},          
    {0.125233408511469 , 0.249147045813403},          
    /*
      n = 13 , start position 42
    */                                                            
    {0.984183054718588 , 0.040484004765316},          
    {0.917598399222978 , 0.092121499837728},          
    {0.801578090733310 , 0.138873510219787},          
    {0.642349339440340 , 0.178145980761946},          
    {0.448492751036447 , 0.207816047536889},          
    {0.230458315955135 , 0.226283180262897},          
    {0.000000000000000 , 0.232551553230874},          
    /*
      n = 14 , start position 49
    */                                                           
    {0.986283808696812 , 0.035119460331752},          
    {0.928434883663574 , 0.080158087159760},          
    {0.827201315069765 , 0.121518570687903},          
    {0.687292904811685 , 0.157203167158194},          
    {0.515248636358154 , 0.185538397477938},          
    {0.319112368927890 , 0.205198463721296},          
    {0.108054948707344 , 0.215263853463158},          
    /*
      n = 15 , start position 56
    */                                                            
    {0.987992518020485 , 0.030753241996117},          
    {0.937273392400706 , 0.070366047488108},          
    {0.848206583410427 , 0.107159220467172},          
    {0.724417731360170 , 0.139570677926154},          
    {0.570972172608539 , 0.166269205816994},          
    {0.394151347077563 , 0.186161000015562},          
    {0.201194093997435 , 0.198431485327111},          
    {0.000000000000000 , 0.202578241925561},          
    /*
      n = 16 , start position 64
    */                                                     
    {0.989400934991650 , 0.027152459411754},          
    {0.944575023073233 , 0.062253523938648},          
    {0.865631202387832 , 0.095158511682493},          
    {0.755404408355003 , 0.124628971255534},          
    {0.617876244402644 , 0.149595988816577},          
    {0.458016777657227 , 0.169156519395003},          
    {0.281603550779259 , 0.182603415044924},          
    {0.095012509837637 , 0.189450610455069},          
    /*
      n = 20 , start position 72
    */                                                         
    {0.993128599185094 , 0.017614007139152},
    {0.963971927277913 , 0.040601429800386},          
    {0.912234428251325 , 0.062672048334109},          
    {0.839116971822218 , 0.083276741576704},          
    {0.746331906460150 , 0.101930119817240},          
    {0.636053680726515 , 0.118194531961518},          
    {0.510867001950827 , 0.131688638449176},          
    {0.373706088715419 , 0.142096109318382},          
    {0.227785851141645 , 0.149172986472603},          
    {0.076526521133497 , 0.152753387130725},          
    /*
      n = 24 , start position 82
    */                                                            
    {0.995187219997021 , 0.012341229799987},          
    {0.974728555971309 , 0.028531388628933},          
    {0.938274552002732 , 0.044277438817419},          
    {0.886415527004401 , 0.059298584915436},          
    {0.820001985973902 , 0.073346481411080},          
    {0.740124191578554 , 0.086190161531953},          
    {0.648093651936975 , 0.097618652104113},          
    {0.545421471388839 , 0.107444270115965},          
    {0.433793507626045 , 0.115505668053725},          
    {0.315042679696163 , 0.121670472927803},          
    {0.191118867473616 , 0.125837456346828},          
    {0.064056892862605 , 0.127938195346752},          
    /*
      n = 32 , start position 94
    */                                                            
    {0.997263861849481 , 0.007018610009470},          
    {0.985611511545268 , 0.016274394730905},          
    {0.964762255587506 , 0.025392065309262},          
    {0.934906075937739 , 0.034273862913021},          
    {0.896321155766052 , 0.042835898022226},          
    {0.849367613732569 , 0.050998059262376},          
    {0.794483795967942 , 0.058684093478535},       
    {0.732182118740289 , 0.065822222776361},       
    {0.663044266930215 , 0.072345794108848},        
    {0.587715757240762 , 0.078193895787070},        
    {0.506899908932229 , 0.083311924226946},        
    {0.421351276130635 , 0.087652093004403},        
    {0.331868602282127 , 0.091173878695763},        
    {0.239287362252137 , 0.093844399080804},        
    {0.144471961582796 , 0.095638720079274},         
    {0.048307665687738 , 0.096540088514727},         
    /*
      n = 40 , start position 110
    */                                                            
    {0.998237709710559 , 0.004521277098533},         
    {0.990726238699457 , 0.010498284531152},         
    {0.977259949983774 , 0.016421058381907},         
    {0.957916819213791 , 0.022245849194166},         
    {0.932812808278676 , 0.027937006980023},         
    {0.902098806968874 , 0.033460195282547},         
    {0.865959503212259 , 0.038782167974472},         
    {0.824612230833311 , 0.043870908185673},         
    {0.778305651426519 , 0.048695807635072},         
    {0.727318255189927 , 0.053227846983936},         
    {0.671956684614179 , 0.057439769099391},         
    {0.612553889667980 , 0.061306242492928},         
    {0.549467125095128 , 0.064804013456601},         
    {0.483075801686178 , 0.067912045815233},         
    {0.413779204371605 , 0.070611647391286},         
    {0.341994090825758 , 0.072886582395804},         
    {0.268152185007253 , 0.074723169057968},         
    {0.192697580701371 , 0.076110361900626},         
    {0.116084070675255 , 0.077039818164247},         
    {0.038772417506050 , 0.077505947978424},         
    /*
      n = 48 , start position 130
    */                                            
    {0.998771007252426 , 0.003153346052305},         
    {0.993530172266350 , 0.007327553901276},         
    {0.984124583722826 , 0.011477234579234},         
    {0.970591592546247 , 0.015579315722943},         
    {0.952987703160430 , 0.019616160457355},         
    {0.931386690706554 , 0.023570760839324},         
    {0.905879136715569 , 0.027426509708356},         
    {0.876572020274247 , 0.031167227832798},         
    {0.843588261624393 , 0.034777222564770},         
    {0.807066204029442 , 0.038241351065830},         
    {0.767159032515740 , 0.041545082943464},         
    {0.724034130923814 , 0.044674560856694},         
    {0.677872379632663 , 0.047616658492490},         
    {0.628867396776513 , 0.050359035553854},         
    {0.577224726083972 , 0.052890189485193},         
    {0.523160974722233 , 0.055199503699984},         
    {0.466902904750958 , 0.057277292100403},         
    {0.408686481990716 , 0.059114839698395},         
    {0.348755886292160 , 0.060704439165893},         
    {0.287362487355455 , 0.062039423159892},         
    {0.224763790394689 , 0.063114192286254},         
    {0.161222356068891 , 0.063924238584648},         
    {0.097004699209462 , 0.064466164435950},         
    {0.032380170962869 , 0.064737696812683},         
    /*
      n = 64 , start position 154
    */                                                            
    {0.999305041735772 , 0.001783280721696},         
    {0.996340116771955 , 0.004147033260562},         
    {0.991013371476744 , 0.006504457968978},         
    {0.983336253884625 , 0.008846759826363},         
    {0.973326827789910 , 0.011168139460131},         
    {0.961008799652053 , 0.013463047896718},         
    {0.946411374858402 , 0.015726030476024},         
    {0.929569172131939 , 0.017951715775697},         
    {0.910522137078502 , 0.020134823153530},         
    {0.889315445995114 , 0.022270173808383},         
    {0.865999398154092 , 0.024352702568710},         
    {0.840629296252580 , 0.026377469715054},         
    {0.813265315122797 , 0.028339672614259},         
    {0.783972358943341 , 0.030234657072402},         
    {0.752819907260531 , 0.032057928354851},         
    {0.719881850171610 , 0.033805161837141},         
    {0.685236313054233 , 0.035472213256882},         
    {0.648965471254657 , 0.037055128540240},         
    {0.611155355172393 , 0.038550153178615},         
    {0.571895646202634 , 0.039953741132720},         
    {0.531279464019894 , 0.041262563242623},         
    {0.489403145707052 , 0.042473515123653},         
    {0.446366017253464 , 0.043583724529323},         
    {0.402270157963991 , 0.044590558163756},         
    {0.357220158337668 , 0.045491627927418},         
    {0.311322871990210 , 0.046284796581314},         
    {0.264687162208767 , 0.046968182816210},         
    {0.217423643740007 , 0.047540165714830},         
    {0.169644420423992 , 0.047999388596458},         
    {0.121462819296120 , 0.048344762234802},         
    {0.072993121787799 , 0.048575467441503},         
    {0.024350292663424 , 0.048690957009139},         
    /*
      n = 80 , start position 186
    */                                                              
    {0.999553822651630 , 0.001144950003186},         
    {0.997649864398237 , 0.002663533589512},         
    {0.994227540965688 , 0.004180313124694},         
    {0.989291302499755 , 0.005690922451403},         
    {0.982848572738629 , 0.007192904768117},         
    {0.974909140585727 , 0.008683945269260},         
    {0.965485089043799 , 0.010161766041103},         
    {0.954590766343634 , 0.011624114120797},         
    {0.942242761309872 , 0.013068761592401},         
    {0.928459877172445 , 0.014493508040509},         
    {0.913263102571757 , 0.015896183583725},         
    {0.896675579438770 , 0.017274652056269},         
    {0.878722567678213 , 0.018626814208299},         
    {0.859431406663111 , 0.019950610878141},         
    {0.838831473580255 , 0.021244026115782},         
    {0.816954138681463 , 0.022505090246332},         
    {0.793832717504605 , 0.023731882865930},         
    {0.769502420135041 , 0.024922535764115},         
    {0.744000297583597 , 0.026075235767565},         
    {0.717365185362099 , 0.027188227500486},         
    {0.689637644342027 , 0.028259816057276},         
    {0.660859898986119 , 0.029288369583267},         
    {0.631075773046871 , 0.030272321759557},         
    {0.600330622829751 , 0.031210174188114},         
    {0.568671268122709 , 0.032100498673487},         
    {0.536145920897131 , 0.032941939397645},         
    {0.502804111888784 , 0.033733214984611},         
    {0.468696615170544 , 0.034473120451753},         
    {0.433875370831756 , 0.035160529044747},         
    {0.398393405881969 , 0.035794393953416},         
    {0.362304753499487 , 0.036373749905835},         
    {0.325664370747701 , 0.036897714638276},         
    {0.288528054884511 , 0.037365490238730},         
    {0.250952358392272 , 0.037776364362001},        
    {0.212994502857666 , 0.038129711314477},         
    {0.174712291832646 , 0.038424993006959},         
    {0.136164022809143 , 0.038661759774076},         
    {0.097408398441584 , 0.038839651059051},         
    {0.058504437152420 , 0.038958395962769},         
    {0.019511383256793 , 0.039017813656306},         
    /*
      n = 96 , start position 226
    */                                                          
    {0.999689503883230 , 0.000796792065552},         
    {0.998364375863181 , 0.001853960788946},         
    {0.995981842987209 , 0.002910731817934},         
    {0.992543900323762 , 0.003964554338444},         
    {0.988054126329623 , 0.005014202742927},         
    {0.982517263563014 , 0.006058545504235},         
    {0.975939174585136 , 0.007096470791153},         
    {0.968326828463264 , 0.008126876925698},         
    {0.959688291448742 , 0.009148671230783},         
    {0.950032717784437 , 0.010160770535008},         
    {0.939370339752755 , 0.011162102099838},         
    {0.927712456722308 , 0.012151604671088},         
    {0.915071423120898 , 0.013128229566961},         
    {0.901460635315852 , 0.014090941772314},         
    {0.886894517402420 , 0.015038721026994},         
    {0.871388505909296 , 0.015970562902562},         
    {0.854959033434601 , 0.016885479864245},         
    {0.837623511228187 , 0.017782502316045},         
    {0.819400310737931 , 0.018660679627411},         
    {0.800308744139140 , 0.019519081140145},         
    {0.780369043867433 , 0.020356797154333},         
    {0.759602341176647 , 0.021172939892191},         
    {0.738030643744400 , 0.021966644438744},         
    {0.715676812348967 , 0.022737069658329},         
    {0.692564536642171 , 0.023483399085926},         
    {0.668718310043916 , 0.024204841792364},         
    {0.644163403784967 , 0.024900633222483},         
    {0.618925840125468 , 0.025570036005349},         
    {0.593032364777572 , 0.026212340735672},         
    {0.566510418561397 , 0.026826866725591},         
    {0.539388108324357 , 0.027412962726029},         
    {0.511694177154667 , 0.027970007616848},         
    {0.483457973920596 , 0.028497411065085},         
    {0.454709422167743 , 0.028994614150555},         
    {0.425478988407300 , 0.029461089958167},         
    {0.395797649828908 , 0.029896344136328},         
    {0.365696861472313 , 0.030299915420827},         
    {0.335208522892625 , 0.030671376123669},         
    {0.304364944354496 , 0.031010332586313},         
    {0.273198812591049 , 0.031316425596861},         
    {0.241743156163840 , 0.031589330770727},         
    {0.210031310460567 , 0.031828758894411},         
    {0.178096882367618 , 0.032034456231992},         
    {0.145973714654896 , 0.032206204794030},         
    {0.113695850110665 , 0.032343822568575},         
    {0.081297495464425 , 0.032447163714064},         
    {0.048812985136049 , 0.032516118713868},         
    {0.016276744849602 , 0.032550614492363}         
  };                                                                   
                                                            
  int i, ipoint, index, fraction;

  /****************************************************
     For a given value of n, compute the start
     location for the related abscissae and weights
     in the array gausscof. The start location is
     given by the number stored in the element
     index = pointer[ipoint].
     The number of abscissae (and weights) returned
     is (n+1)/2.
  ****************************************************/

  fraction = 0;

  if (n >= 2 && n <= 16)
    ipoint = n;
  else if (n > 16 && n <= 24) {
    ipoint = 16 + (n-16)/4;
    fraction = (n-16)%4;
  }
  else if (n >= 32 && n <= 48) {
    ipoint = 19 + (n-32)/8;
    fraction = (n-32)%8;
  }
  else if (n >= 64 && n <= 96) {
    ipoint = 22 + (n-64)/16;
    fraction = (n-64)%16;
  }
  else
    fraction = 1;

  /****************************************************
    If value of n is in the set of permitted values
    then transfere the abscissae and weights for the
    related Gauss-Legendre formula to the x[] and
    w[] arrays.
  ****************************************************/

  if (fraction == 0) {
    index = pointer[ipoint -2];
    *ncof = (n + 1)/2;
    for (i = 0; i <= *ncof-1; i++) {
      x[i] = -gausscof[index + i][0];
      w[i] = gausscof[index + i][1];
    } /* End of i-loop. */
    *flag = 1;
  }
  else
    *flag = 0;

  if ((n%2)==0){
    for (i=0; i<=(n/2-1); i++){
      x[n/2+i] = -x[n/2-1-i];
      w[n/2+i] =  w[n/2-1-i];
    }
  }
  else{
    for (i=0; i<=((n-1)/2-1); i++){
      x[(n+1)/2+i] = -x[(n-1)/2-1-i];
      w[(n+1)/2+i] =  w[(n-1)/2-1-i];
    }
  }

}



void Gauss_Legendre2( double x1, double x2,  double x[],  double w[], int n )
{

  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1; 

  m=(n+1)/2;
  xm=0.50*(x2+x1);
  xl=0.50*(x2-x1);

  for (i=1;i<=m;i++) {

    z=cos(PI*(i-0.250)/(n+0.50));

    do {

      p1=1.0;
      p2=0.0;

      for (j=1;j<=n;j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*(double)j-1.0)*z*p2-((double)j-1.0)*p3)/(double)j;
      }

      pp=(double)n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;

    } while (fabs(z-z1) > EPS);

    x[i]=xm-xl*z;     
    x[n+1-i]=xm+xl*z; 
    w[i]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i]=w[i];    
  }

  for (i=1; i<=n; i++) {
    x[i-1] = x[i];    
    w[i-1] = w[i];    
  }  

}

