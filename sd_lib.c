/* library file for sd / panic
   used by sd.c and sd_crunch.c 
   */

int MainSwitch_DEFAULT = 0; 
char *IFN_DEFAULT = "sd.par";
char *OFN_DEFAULT = "sd.dat";
char *OF2N_DEFAULT = "sd.dat2";

/********  global constants for both sd.c and sd_crunch.c **********/  

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
#include <unistd.h>

#include "nrutil.c"
#include "base.c"
#include "readpar.c"


#define SD_LIB_EXIT {_E("sd_lib.c: Exiting to system.\n");exit(-1);}
#define MY_STRLEN 200
float EPSILON = 1.0e-5; 

/* number of walls and wpoints */
#define NW 9   
#define NWP 4  
typedef struct wall { float x1, y1, x2, y2; } wall; 
typedef struct wpoint { float x, y; } wpoint;

/************** constants needed for sd.c only ***************/

#ifdef _FLAG_SD_C_

#include "Xlibext.c" 
#include "Epslib.c"
char *ParticleColorName[]={"yellow", "grey"};
char *SmokeColorName="grey";
XColor BGC_sdef,IC_sdef,PC_sdef[100],SmokeColor_sdef;
int X11_WWi, X11_WHe, BGCCode, ICCode, PaCNum, *PaCCode, SmokeCCode;

#endif /* _FLAG_SD_C_ */

/********** constants needed for sd_crunch.c only ***************/

#ifdef _FLAG_SD_CRUNCH_C_

int ii,ij,ik, RndSeed_main_num, VSmoke_main_num, RndSeed_main_index,
  VSmoke_main_index,*RndSeed_main,RndSeed_main_1st;
float *VSmoke_main,***Tleave,***Tinjured;

#endif /* _FLAG_SD_CRUNCH_C_ */

/********* global parameters -- to be read from parameter file **********/

static int N0, FrictionSwitch, InjurySwitch, ColumnSwitch,
  X11_Margin, X11_InFW, X11_InFH, X11_TLH, X11_GrFH, X11_RightRim, 
  EpsPicMult, Eps_iPF_first, Eps_iPF_max,
  SaveUN, DrawUN, Sleep, Draw, 
  RndSeed, MaxUpdNum, AyS;
static float RoomXSize, RoomYSize, DoorWidth, WallWidth, Dmean,
  deltaD, A, B, A_fire, B_fire, Kappa, Gamma, C_Young, R, R_fire, V0, Tau,
  GaMe, GaTh, GaCM, 
  SmokeStartTime, VSmoke, FCrush_over_1m, 
  ColumnCenterX, ColumnCenterY, ColumnD,
  X11_Magn, 
  EpsXS, EpsYS, EpsMinXMarg, EpsMinYMarg, EpsInFW, EpsInTH, EpsInFH,
  EpsLineWidth, SaveST, DrawST, DrawDMult, 
  JavaXS, JavaYS, JavaMinXMarg, JavaMinYMarg,
  JavaTStep, JavaMaxTime,
  MaxSimTime, Vmax, H, DefaultDeltaT, C_NS, V_ChangeLimit;
static char BackGroundColorName[MY_STRLEN], InfoColorName[MY_STRLEN],
  X11_FontName[MY_STRLEN];  


int *IPar[]={&N0, &FrictionSwitch, &InjurySwitch, &ColumnSwitch,
	     &X11_Margin, &X11_InFW,
	     &X11_InFH, &X11_TLH, &X11_GrFH, &X11_RightRim,
	     &EpsPicMult, &Eps_iPF_first, &Eps_iPF_max, &SaveUN,
	     &DrawUN, &Sleep, &Draw, 
	     &RndSeed, &MaxUpdNum, &AyS}; 
char *IParName[]={"N0", "FrictionSwitch", "InjurySwitch", 
		  "ColumnSwitch", "X11_Margin",
		  "X11_InFW", "X11_InFH", "X11_TLH", "X11_GrFH",
		  "X11_RightRim", "EpsPicMult", "Eps_iPF_first",
		  "Eps_iPF_max", "SaveUN", "DrawUN", "Sleep", 
		  "Draw", "RndSeed", "MaxUpdNum", "AyS"}; 
float *FPar[]={&RoomXSize, &RoomYSize, &DoorWidth, &WallWidth, &Dmean,
	       &deltaD, &A, &B, &A_fire, &B_fire, &Kappa, &Gamma,
	       &C_Young, &R, &R_fire, &V0, &Tau, 
	       &GaMe, &GaTh, &GaCM, &SmokeStartTime, &VSmoke, 
	       &FCrush_over_1m, 
	       &ColumnCenterX, &ColumnCenterY, &ColumnD,
	       &X11_Magn, &EpsXS, &EpsYS,
	       &EpsMinXMarg, &EpsMinYMarg, &EpsInFW, &EpsInTH,
	       &EpsInFH, &EpsLineWidth, &SaveST, &DrawST, 
	       &DrawDMult, 
	       &JavaXS, &JavaYS, &JavaMinXMarg, &JavaMinYMarg,
	       &JavaTStep, &JavaMaxTime,
	       &MaxSimTime, &Vmax, &H, &DefaultDeltaT,
	       &C_NS, &V_ChangeLimit}; 
char *FParName[]={"RoomXSize", "RoomYSize", "DoorWidth", "WallWidth",
	          "Dmean", "deltaD", "A", "B", "A_fire", "B_fire",
		  "Kappa", "Gamma", "C_Young", 
		  "R", "R_fire", "V0", "Tau", "GaMe", "GaTh", "GaCM", 
		  "SmokeStartTime", "VSmoke", 
		  "FCrush_over_1m", 
		  "ColumnCenterX", "ColumnCenterY", "ColumnD",
		  "X11_Magn", "EpsXS", "EpsYS",
		  "EpsMinXMarg", "EpsMinYMarg", "EpsInFW", "EpsInTH",
		  "EpsInFH", "EpsLineWidth", "SaveST",   
		  "DrawST", "DrawDMult", 
		  "JavaXS", "JavaYS", "JavaMinXMarg", "JavaMinYMarg",
		  "JavaTStep", "JavaMaxTime",
		  "MaxSimTime", "Vmax", "H",
		  "DefaultDeltaT", "C_NS", "V_ChangeLimit"}; 
char *SPar[]={BackGroundColorName, InfoColorName, X11_FontName};  
char *SParName[]={"BackGroundColorName", "InfoColorName",
		  "X11_FontName"};   


int IParNum = sizeof(IPar)/sizeof(int*),
    FParNum = sizeof(FPar)/sizeof(float*),
    SParNum = sizeof(SPar)/sizeof(char*);

/******* global variables needed both for sd.c and sd_crunch.c ************/

int MainSwitch, UpdNum, N, GX, GY, G, Mb, Me, *BIndBd, *BInd,
  GaussFlag, NInRoom, *Injured, NInjured, Eps_iPF, JavaXMarg,
  JavaYMarg, JavaImInd, JavaNTStep,
  **JavaXlo,**JavaXhi,**JavaYlo,**JavaYhi,**JavaC,**JavaD;
float *SimTime, XS, YS, *D, *Phi, *X, *Y, *Xprev, *Yprev, *V, *VX,
  *VY, *E, *Vdir, GaussSet1, GaussSet2, *V0of, EpsXMarg, EpsYMarg,
  EpsMagn, JavaMagn; 
char IFN[MY_STRLEN], OFN[MY_STRLEN], OF2N[MY_STRLEN];
FILE *OFP,*OF2P;
wall *W;
wpoint *WP;

struct stat IFStatBuf;
long int IFModTime;

/* MainSwitch (see in sd.c)
   UpdNum     index of present update
   N          number of particles
   IFN, IFP   input file name, input file pointer
   IFStatBuf  input file status buffer
   IFModTime  modification time of input file
   X11_WWi/X11_WHe   
              window width/height
   SimTime[n]    time after the nth update (SimTime[0] = 0.0)

   X[i]       current coordinate of particle i
   Xprev[i]   previous coordinate
   VX[i]      current velocity component of particle i
   V[i]       magnitude of current velocity of particle i
   Phi[i]     preferred direction of particle i
   Vdir[i]    direction of particle i's velocity
   V0of[i]    preferred magnitude of velocity of particle i

   E[n]       is the efficiency at update step no. n
   Mb/Me      mem.alloc. begin/end 
              ie. limits of current time window 
	      (see explanation of AyS in sd.par)
   XS,YS      size of field
   GaussFlag,GaussSet1,GaussSet2 
              are used by the random number generator
   NInRoom    number of particles in room
   Injured[i] =1, if pedestrian i is injured
              =0, if not
   NInjured   number of injured	  
   EpsXMarg   eps x margin
   EpsMagn    magnification of eps image (determined from width of
              margins)
   Eps_iPF    index of eps picture file

   JavaMagn   magnification for java images
   JavaXMarg  x margin
   JavaImInd  current index of image
   JavaNTStep number of time steps
   JavaXlo[n][i] + 256*JavaXhi[n][i] 
              is the x coordinate of particle n at time step i  
   JavaYlo[n][i] + 256*JavaYhi[n][i] 
              is the y coordinate of particle n at time step i  
   JavaC[n][i] the color code of particle n at time step i  
   JavaD[n][i] the diameter of particle n at time step i  
   


BOOK-KEEPING OF PARTICLES:
  GX*GY blocks (all same size)
  the size of a block should be: not less than R, i.e.
  GX = (int)floor(XS/R) i.e. R <= XS/GX
  indexed with 0, 1, ... , G^2 - 1, where G = max(GX,GY)
  block indices of particles stored in BIndBd & BInd
	     
  BIndBd[n] (Block + Index + Board) 
             = 0,1,...,N-1: the number of the nth block's first particle
	     = -1: there's no particle in the nth block
  BInd[i]   = 0,1,...,N-1: the next particle in the block of particle i
            = -1: there are no more particles in this block 
	    */

/*====== global variables needed only for sd.c =========*/

#ifdef _FLAG_SD_C_
float FW_x;
#endif /* _FLAG_SD_C_ */

/* FW_x     x component of force exerted on walls left and right
            from exit  
*/   

/***************** prototypes **********************/

/* functions needed in all cases */

void WallParticleRelation(int iw, int i, float *tr, int *can_see);
void WallPsychForce(int iw, int i, float r, float *fx, float *fy);
void WallYoungForce(int iw, int i, float r, float *fx, float *fy);
void WallTangForce_FS0(int iw, int i, float *fx, float *fy);
void WallTangForce_FS1(int iw, int i, float r, float *fx, float *fy);

void WPointParticleRelation(int iwp, int i, float *r, int *can_see);
void WPointPsychForce(int iwp, int i, float r, float *fx, float *fy);
void WPointYoungForce(int iwp, int i, float r, float *fx, float *fy);
void WPointTangForce_FS0(int iwp, int i, float r, float *fx, float *fy);
void WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy);

void PP_PsychForce(int i1, int i2, float r, float *fx, float *fy);
void PP_YoungForce(int i1, int i2, float r, float *fx, float *fy);
void PP_TangForce_FS0(int i1, int i2, float r, float *fx, float *fy);
void PP_TangForce_FS1(int i1, int i2, float r, float *fx, float *fy);

float DirectionOfExit( int i );
void RemoveParticle( int *n, int i );
void EulTStep( float *tstep, float f ); 


void Start_Bare( int narg, char *argstr[] );
void Init_Bare( char *init_switch );
void Upd();
float GaussRand( float gmean, float gtheta, float gcutmult );
float EMean( char* sw, int cnfreq, float stfreq );



/* functions needed only for sd.c */
#ifdef _FLAG_SD_C_ 

void Start_Demo( int narg, char *argstr[] );
void Init_Demo();
void X11_init();
void Eps_init();
void Java_init();
void Pic();
void X11_Pic();
void EpsDrawParticle( FILE* fp, int i );
void Eps_Pic();
void XDrawParticle( int i, int leftxmargin, int upymargin, float magn);
void Save_Demo();
void Save_Java_in_Loop();
void Save_Java_after_Loop();
void Shutdown_Demo();
void ReInit();

#endif /* _FLAG_SD_C_ */



/* functions needed only for sd_crunch.c */
#ifdef _FLAG_SD_CRUNCH_C_

void Start_Crunch( int narg, char *argstr[] );
void Init_Crunch();
void Save_Crunch_InLoop();
void Save_Crunch_AfterLoop();
void Clean_Crunch_AfterLoop();
void Save_Crunch_AtEnd();

#endif /* _FLAG_SD_CRUNCH_C_ */

/********************* definitions *****************/

void WallParticleRelation(int iw, int i, float *r, int *can_see){
  /* can_see: whether partice i is within the range of wall iw 
     r: distance */


  switch(iw){
  case 0:{ *r = Y[i];           break; }
  case 1:{ *r = WP[0].x-X[i];   break; }
  case 2:{ *r = Y[i]-WP[0].y;   break; }
  case 3:{ *r = X[i]-WP[1].x;   break; }
  case 4:{ *r = RoomYSize-Y[i]; break; }
  case 5:{ *r = X[i]-WP[2].x;   break; }
  case 6:{ *r = WP[2].y-Y[i];   break; }
  case 7:{ *r = WP[3].x-X[i];   break; }
  case 8:{ *r = X[i];           break; }
  }


  switch(iw){
  case 0: {   
      if(Y[i]<=R){ *can_see=1; } 
      else{ *can_see=0; }    
      break;  
  }
  case 1: {
      if((X[i]>=WP[0].x-R)&&(X[i]<=WP[0].x)&&(Y[i]<=WP[0].y)){ *can_see=1; }
      else{ *can_see=0; }
      break;
  }
  case 2: {
      if((X[i]>=WP[0].x)&&(X[i]<=WP[1].x)&&(Y[i]<=WP[0].y+R)){ *can_see=1; }
      else{ *can_see=0; }
      break;
  }
  case 3:{
      if((X[i]>=WP[1].x)&&(X[i]<=WP[1].x+R)&&(Y[i]<=WP[1].y)){ *can_see=1; }
      else{ *can_see=0; }
      break;
  }
  case 4:{
      if(Y[i]>=RoomYSize-R){ *can_see=1; }
      else{ *can_see=0; }
      break;
  }
  case 5:{
      if((X[i]>=WP[2].x)&&(X[i]<=WP[2].x+R)&&(Y[i]>=WP[2].y)){ *can_see=1; }
      else{ *can_see=0; }
      break;
  }
  case 6:{
      if((X[i]>=WP[3].x)&&(X[i]<=WP[2].x)&&(Y[i]>=WP[2].y-R)){ *can_see=1; }
      else{ *can_see=0; }
      break;    
  }
  case 7:{
      if((X[i]<=WP[3].x)&&(X[i]>=WP[3].x-R)&&(Y[i]>=WP[3].y)){ *can_see=1;}
      else{ *can_see=0; }
      break;    
  }
  case 8:{
      if(X[i]<=R){ *can_see=1; }
      else{ *can_see=0; }
      break;    
  }
  }
}

/*------------------------------*/

void WallPsychForce(int iw, int i, float r, float *fx, float *fy){

#define tmp_f (A*exp(-(r-0.5*D[i])/B))

  switch(iw){
  case 0:{ *fx = 0.0;     *fy = tmp_f;   break; }
  case 1:{ *fx = - tmp_f; *fy = 0.0;     break; }
  case 2:{ *fx = 0.0;     *fy = tmp_f;   break; }
  case 3:{ *fx = tmp_f;   *fy = 0.0;     break; }
  case 4:{ *fx = 0.0;     *fy = - tmp_f; break; }
  case 5:{ *fx = tmp_f;   *fy = 0.0;     break; } 
  case 6:{ *fx = 0.0;     *fy = - tmp_f; break; }
  case 7:{ *fx = - tmp_f; *fy = 0.0;     break; }
  case 8:{ *fx = tmp_f;   *fy = 0.0;     break; }
  }

#undef tmp_f
}

/*------------------------------*/

void WallYoungForce(int iw, int i, float r, float *fx, float *fy){

#define tmp_f (2.0*C_Young*(0.5*D[i]-r))

  switch(iw){
  case 0:{ *fx = 0.0;     *fy = tmp_f;   break; }
  case 1:{ *fx = - tmp_f; *fy = 0.0;     break; }
  case 2:{ *fx = 0.0;     *fy = tmp_f;   break; }
  case 3:{ *fx = tmp_f;   *fy = 0.0;     break; }
  case 4:{ *fx = 0.0;     *fy = - tmp_f; break; }
  case 5:{ *fx = tmp_f;   *fy = 0.0;     break; } 
  case 6:{ *fx = 0.0;     *fy = - tmp_f; break; }
  case 7:{ *fx = - tmp_f; *fy = 0.0;     break; }
  case 8:{ *fx = tmp_f;   *fy = 0.0;     break; }
  }

#undef tmp_f
}

/*------------------------------*/

void WallTangForce_FS0(int iw, int i, float *fx, float *fy){

  switch(iw){
  case 0: case 2: case 4: case 6: { 
      *fx = -Gamma*VX[i]; 
      *fy = 0.0;          
      break; 
  }
  case 1: case 3: case 5: case 7: case 8: { 
      *fx = 0.0;          
      *fy = -Gamma*VY[i]; 
      break; 
  }
  }
}

/*------------------------------*/

void WallTangForce_FS1( int iw, int i, float r, float *fx, float *fy ){

#define tmp_delta_r (0.5*D[i]-r)

  /* friction forces */
  switch(iw){
  case 0: case 2: case 4: case 6: { 
      *fx = -Kappa*tmp_delta_r*VX[i]; 
      *fy = 0.0;          
      break; 
  }
  case 1: case 3: case 5: case 7: case 8: { 
      *fx = 0.0;          
      *fy = -Kappa*tmp_delta_r*VY[i]; 
      break; 
  }
  }

#undef tmp_delta_r
}

/*------------------------------------*/

void WPointParticleRelation(int iwp, int i, float *r, int *can_see){
  /* can_see: whether partice i is within the range of wpoint iwp 
     r: distance */
  
  *r = sqrt(SQR(WP[iwp].x-X[i])+SQR(WP[iwp].y-Y[i]));

  switch(iwp){
  case 0:{
    if((X[i]<=WP[0].x)&&(Y[i]>=WP[0].y)){ *can_see=1; }
    else{ *can_see=0; }
    break;
  }
  case 1:{
    if((X[i]>=WP[1].x)&&(Y[i]>=WP[1].y)){ *can_see=1; }
    else{ *can_see=0; }
    break;
  }
  case 2:{
    if((X[i]>=WP[2].x)&&(Y[i]<=WP[2].y)){ *can_see=1; }
    else{ *can_see=0; }
    break;
  }
  case 3:{
    if((X[i]<=WP[3].x)&&(Y[i]<=WP[3].y)){ *can_see=1; }
    else{ *can_see=0; }
    break;
  }
  }
}

/*------------------------------*/

void WPointPsychForce(int iwp, int i, float r, float *fx, float *fy){
  /* exerted by wpoint iwp on particle i */
  
#define tmp_f_over_r (A*exp(-(r-0.5*D[i])/B)/r)

  *fx = (X[i]-WP[iwp].x) * tmp_f_over_r;
  *fy = (Y[i]-WP[iwp].y) * tmp_f_over_r;

#undef tmp_f_over_r
}

/*------------------------------*/

void WPointYoungForce(int iwp, int i, float r, float *fx, float *fy){
  /* exerted by wpoint iwp on particle i */
  
  float rx,ry;

#define tmp_f_over_r ( 2.0*C_Young*(0.5*D[i]-r) / r)

  rx=WP[iwp].x-X[i];
  ry=WP[iwp].y-Y[i];
  *fx = - rx * tmp_f_over_r;
  *fy = - ry * tmp_f_over_r;

#undef tmp_f_over_r
}

/*------------------------------*/

void WPointTangForce_FS0(int iwp, int i, float r, float *fx, float *fy){
  /* exerted by wpoint iwp on particle i */

  float rx,ry,scal_prod_over_rsqr;

  rx = X[i]-WP[iwp].x;
  ry = Y[i]-WP[iwp].y;
  scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(r); 
  *fx = -Gamma * (   ry * scal_prod_over_rsqr );
  *fy = -Gamma * ( - rx * scal_prod_over_rsqr );
}

/*------------------------------*/

void WPointTangForce_FS1(int iwp, int i, float r, float *fx, float *fy){
  /* exerted by wpoint iwp on particle i */

  float rx,ry,scal_prod_over_rsqr;

  rx = X[i]-WP[iwp].x;
  ry = Y[i]-WP[iwp].y;
  scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(r); 
  *fx = -Kappa * (0.5*D[i]-r) * (   ry * scal_prod_over_rsqr );
  *fy = -Kappa * (0.5*D[i]-r) * ( - rx * scal_prod_over_rsqr );
}

/*------------------------------*/

void PP_PsychForce(int i1, int i2, float r, float *fx, float *fy){

  float f_over_r;

  f_over_r = A*exp(-(r-0.5*(D[i1]+D[i2]))/B) / r;
  *fx = (X[i1]-X[i2]) * f_over_r;
  *fy = (Y[i1]-Y[i2]) * f_over_r;
}

/*------------------------------*/

void PP_YoungForce(int i1, int i2, float r, float *fx, float *fy){

  float f_over_r;

  f_over_r = 2.0*C_Young*(0.5*(D[i1]+D[i2])-r) / r;
  *fx = (X[i1]-X[i2]) * f_over_r;
  *fy = (Y[i1]-Y[i2]) * f_over_r;
}

/*------------------------------*/

void PP_TangForce_FS0(int i1, int i2, float r, float *fx, float *fy){
  /* exerted by particle i2 on particle i1 */

  float rx,ry,vx,vy,scal_prod_over_rsqr;

  rx = X[i1]-X[i2];
  ry = Y[i1]-Y[i2];
  vx = VX[i1]-VX[i2];
  vy = VY[i1]-VY[i2];
  scal_prod_over_rsqr = (ry*vx - rx*vy) / SQR(r); 
  *fx = -Gamma * (   ry * scal_prod_over_rsqr );
  *fy = -Gamma * ( - rx * scal_prod_over_rsqr );
}

/*---------------------------*/ 

void PP_TangForce_FS1(int i1, int i2, float r, float *fx, float *fy){
  /* exerted by particle i2 on particle i1 */

  float rx,ry,vx,vy,scal_prod_over_rsqr;

  rx = X[i1]-X[i2];
  ry = Y[i1]-Y[i2];
  vx = VX[i1]-VX[i2];
  vy = VY[i1]-VY[i2];
  scal_prod_over_rsqr = (ry*vx - rx*vy) / SQR(r); 
  *fx = -Kappa * (0.5*(D[i1]+D[i2])-r) * (   ry * scal_prod_over_rsqr );
  *fy = -Kappa * (0.5*(D[i1]+D[i2])-r) * ( - rx * scal_prod_over_rsqr );
}

/*---------------------------*/ 

float DirectionOfExit( int i ){
  /* direction of exit for particle i */

  float dsqr, /* sqr of particle center - door-post distance */
    rsqr; /* sqr of particle's radius */


  /* behind the upper door-post */
  if((Y[i]<=0.5*YS-0.5*DoorWidth+0.5*D[i]+EPSILON)&&(X[i]<=RoomXSize)){
    
          dsqr = SQR(W[1].x2-X[i]) + SQR(W[1].y2-Y[i]);
	  rsqr = SQR(0.5*D[i])+EPSILON;
	  if(dsqr<=rsqr){ 
	          /* very close to the door-post */
	          if(Y[i]<=0.5*YS-0.5*DoorWidth){
		          return( 0.5*PI );
		  }
		  else{
		          return(   0.5*PI 
				  + atan2( W[1].y2-Y[i],W[1].x2-X[i] )
				);  
		  }
	  }
	  else {
	          /* well apart from the door-post */
	          return(   atan2( 1.0, sqrt(dsqr/rsqr-1.0) ) 
                          + atan2( W[1].y2-Y[i],W[1].x2-X[i] )
		        );
	  }
  }
    

  /* behind the lower door-post */
  else if((Y[i]>=0.5*YS+0.5*DoorWidth-0.5*D[i]-EPSILON)&&(X[i]<=RoomXSize)){

          dsqr = SQR(W[6].x2-X[i]) + SQR(W[6].y2-Y[i]);
	  rsqr = SQR(0.5*D[i])+EPSILON;
	  if(dsqr<=rsqr){ 
	          /* very close to the door-post */
	          if(Y[i]>=0.5*YS+0.5*DoorWidth){
		          return( -0.5*PI );
		  }
		  else{
		          return( - 0.5*PI
				  + atan2( W[6].y2-Y[i],W[6].x2-X[i] )
				);
		  }
	  }
	  else {
	          /* well apart from the door-post */
	          return( - atan2( 1.0, sqrt(dsqr/rsqr-1.0) ) 
                          + atan2( W[6].y2-Y[i],W[6].x2-X[i] )
		        );
	  }
  }


  /* in the center or outside */
  else { return 0.0; }
}

/*-----------------------------------------------*/

void RemoveParticle( int *n, int i ){
  /* *n: number of particles now
     i: index of particle to be removed */

  int j;
  

  /* (a) particle i (which is off-board now)
     is removed from the book-keeping
     (block determined by previous coordinates) 
     (b) if i != *n-1 
         (b1) particle *n - 1 is removed from the book-keeping
	      (block determined by previous coordinates)
         (b2) copying all values of particle *n-1 into i's place
         (b3) inserting particle i (that used to be indexed *n-1) into the
              book-keeping, into the block given by the previous
	      coordinates (Xprev[i],Yprev[i]), and not into the block
	      given by (X[i],Y[i])
	      . reason: after this substitution (*n-1 -> i) 
	      particle i will be looked for in the block of
	      (Xprev[i],Yprev[i]) in Upd 
	      because no one tells the main cycle (located in Upd),
	      whether this particle is the result of a substitution 
	      or not 
     (c) decrement particle number  ( *n to *n - 1 ) */


  /* a */
  j = (int)floor(Xprev[i]*GX/XS) + G*(int)floor(Yprev[i]*GY/YS);
  if(BIndBd[j]==i) {
          BIndBd[j] = BInd[i];
  }
  else {
          j = BIndBd[j];
	  while(BInd[j]!=i) {
	          j = BInd[j];
	  }
	  BInd[j] = BInd[i];
  }  



  /* b */
  if(i!=*n-1){

          /* b1 */
    
          j = (int)floor(Xprev[*n-1]*GX/XS) + G*(int)floor(Yprev[*n-1]*GY/YS);
          if(BIndBd[j]==*n-1) {
	          BIndBd[j] = BInd[*n-1];
	  }
	  else {
	          j = BIndBd[j];
	          while(BInd[j]!=*n-1) {
		          j = BInd[j];
		  }
		  BInd[j] = BInd[*n-1];
	  }  



          /* b2 */
          D[i] = D[*n-1];
	  Phi[i] = Phi[*n-1];
	  X[i] = X[*n-1];
	  Y[i] = Y[*n-1];
	  V[i] = V[*n-1];
	  VX[i] = VX[*n-1];
	  VY[i] = VY[*n-1];
	  Xprev[i] = Xprev[*n-1];
	  Yprev[i] = Yprev[*n-1];
	  Vdir[i]=Vdir[*n-1];
	  Injured[i]=Injured[*n-1]; 
	  V0of[i]=V0of[*n-1];

	  

	  /* b3 */
	  j = (int)floor(Xprev[i]*GX/XS) + G*(int)floor(Yprev[i]*GY/YS);
	  if(BIndBd[j]==-1) {
	          BIndBd[j] = i;
		  BInd[i] = -1;
	  }
	  else {
	          j = BIndBd[j];
		  while(BInd[j]!=-1) {
		          j = BInd[j];
		  }
		  BInd[j] = i;
		  BInd[i] = -1;
	  }
  }




  /* c */
  (*n)--;
}

/*-----------------------------------------------------------------------*/

void EulTStep( float *tstep, float f ){
  /* adjusts the time step in a way that the force (fx,fy) doesn't
     change the velocity of particle i by more than V_ChangeLimit */
  
  while( f*(*tstep) >= V_ChangeLimit ){ *tstep *= C_NS; }
}

/********************************/

void Start_Bare( int narg, char *argstr[] ){
  /* reading command line parameters and parameter file */


  /* 1 */

  switch(narg){
  case 1: { 
          MainSwitch = MainSwitch_DEFAULT; 
	  strcpy(IFN,IFN_DEFAULT);
	  strcpy(OFN,OFN_DEFAULT);
	  strcpy(OF2N,OF2N_DEFAULT);
	  break; 
  }
  case 3:{
          MainSwitch = atoi(argstr[1]); 
	  strcpy(IFN,argstr[2]);
	  strcpy(OFN,OFN_DEFAULT);
	  strcpy(OF2N,OF2N_DEFAULT);
	  break; 
  }
#ifdef _FLAG_SD_CRUNCH_C_
  case 6:{
          /* running the program with the loadleveller */
          RndSeed_main_1st = atoi(argstr[1]); 
          MainSwitch = MainSwitch_DEFAULT; 
	  strcpy(IFN,IFN_DEFAULT);
	  strcpy(OFN,OFN_DEFAULT);
	  strcpy(OF2N,OF2N_DEFAULT);
	  break; 
  }
#endif /* _FLAG_SD_CRUNCH_C_ */
  default:{
          _E("Usage:  either \"sd <MainSwitch> <input file name>\",");
	  _E("        or     \"sd\" .\n");
	  SD_LIB_EXIT;
  }}
  fprintf(stderr,"(default values are:\n MainSwitch = %d, input = %s, output = %s, %s)\n", MainSwitch,IFN,OFN,OF2N);
  fflush(stderr);




  /* 2 */

  /* reading parameters */
  readpar ( "start", IFN, IPar, IParName, IParNum, 
	     FPar, FParName, FParNum, SPar, SParName, SParNum ); 


  /* 3 */

}

/*----------------------------------------*/

void Init_Bare( char *init_switch ){
  /* 1 global vars, mem.alloc.
     2 walls, wpoints
     3 particles
   */

  int i,j,ok_flag;



  /* 1 */
  stat(IFN, &IFStatBuf);
  IFModTime = IFStatBuf.st_mtime;
  UpdNum = 0;
  Mb = 0;
  Me = AyS-1;
  SimTime = vector( Mb, Me );
  SimTime[0] = 0.0;
  srand(RndSeed);

  XS = RoomXSize+WallWidth+X11_RightRim+EPSILON;
  YS = RoomYSize;
  N = N0;
  NInRoom = N0;

  GX = (int)MAX(1.0,floor(XS/R));
  GY = (int)MAX(1.0,floor(YS/R));
  G = (int)MAX(GX,GY);

  BIndBd = ivector( 0, SQR(G)-1 );
  BInd = ivector( 0, N0-1 );
  D = vector( 0, N0-1 );
  Phi = vector( 0, N0-1 );
  X = vector( 0, N0-1 );
  Y = vector( 0, N0-1 );
  Xprev = vector( 0, N0-1 );
  Yprev = vector( 0, N0-1 );
  V = vector( 0, N0-1 );
  VX = vector( 0, N0-1 );
  Vdir = vector(0,N0-1);
  VY = vector( 0, N0-1 );
  V0of = vector( 0, N0-1 );
  E = vector( Mb, Me );
  E[0] = 1.0;
  Injured = ivector(0,N0-1);



  /* 2 walls, wpoints */
  /* allocating memory:
   * if there's a column at the door, 
   * the four faces and corners of the column have to be initialized,
   * too  
   */
  W = (wall*)calloc(NW,sizeof(wall));
  WP = (wpoint*)calloc(NWP,sizeof(wpoint));



  /* every wall rotated by PI/2 points towards the inside of the room */

  /* upper part */
  W[0].x1 = 0.0;
  W[0].y1 = 0.0;
  W[0].x2 = XS;
  W[0].y2 = 0.0;

  W[1].x1 = RoomXSize;
  W[1].y1 = 0.0;
  W[1].x2 = W[2].x1 = WP[0].x = RoomXSize;
  W[1].y2 = W[2].y1 = WP[0].y = 0.5*RoomYSize-0.5*DoorWidth;
  W[2].x2 = W[3].x1 = WP[1].x = RoomXSize+WallWidth;
  W[2].y2 = W[3].y1 = WP[1].y = 0.5*RoomYSize-0.5*DoorWidth;
  W[3].x2 = RoomXSize+WallWidth;
  W[3].y2 = 0.0;


  /* lower part */
  W[4].x1 = XS;
  W[4].y1 = RoomYSize;
  W[4].x2 = 0.0;
  W[4].y2 = RoomYSize;
  
  W[5].x1 = RoomXSize+WallWidth;
  W[5].y1 = RoomYSize;
  
  W[5].x2 = W[6].x1 = WP[2].x = RoomXSize+WallWidth;
  W[5].y2 = W[6].y1 = WP[2].y = 0.5*RoomYSize+0.5*DoorWidth;
  W[6].x2 = W[7].x1 = WP[3].x = RoomXSize;
  W[6].y2 = W[7].y1 = WP[3].y = 0.5*RoomYSize+0.5*DoorWidth;
  W[7].x2 = RoomXSize;
  W[7].y2 = RoomYSize;

  
  /* left wall of the room */
  W[8].x1 = 0.0;
  W[8].y1 = RoomYSize;
  W[8].x2 = 0.0;
  W[8].y2 = 0.0;




  /* 3 */

  /* diameters and coordinates */
  for(i=0; i<N; i++) {
      D[i] =   (Dmean + deltaD)
	     - 2.0*deltaD * rand()/(RAND_MAX+1.0);
      X[i] =   0.5*H*D[i]+EPSILON 
	     + (RoomXSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);
      Y[i] =   0.5*H*D[i]+EPSILON 
	     + (RoomYSize-H*D[i]-2.0*EPSILON)*rand()/(RAND_MAX+1.0);


      /* checking whether far enough from the column */
      ok_flag = 1;
      switch(ColumnSwitch){
      default: case 0:{ break; }
      case 1:{ 
	      if(   SQR(X[i]-ColumnCenterX)+SQR(Y[i]-ColumnCenterY)
		 <= SQR(0.5*(D[i]+ColumnD))+EPSILON
		){
		      ok_flag = 0;
		      i--;
	      } 
	      break;
      }
      }


      /* checking distances to already existing particles */
      if(ok_flag==1){
	      for(j=0;j<i;j++) {
		      if(     SQR(X[j] - X[i])
			    + SQR(Y[j] - Y[i])
			 <= SQR( 0.5*H*(D[i]+D[j]) ) + EPSILON
			) {
			      i = i - 1;
			      j = i - 1;
		      }
	      }
      }
  }


  /* book-keeping */
  for(i=0;i<SQR(G);i++){ BIndBd[i] = -1; }
  for(i=0;i<N;i++){ BInd[i] = -1; }

  for(i=0;i<N;i++) {
          j = (int)floor(X[i]*GX/XS) + G * (int)floor(Y[i]*GY/YS);

	  if(BIndBd[j]==-1) {
	          BIndBd[j] = i;
	  }
	  else {
	          j = BIndBd[j];
		  while(BInd[j]!=-1) {
		          j = BInd[j];
                  }
		  BInd[j] = i;
	  }
  }



  /* injuries, velocities and preferred directions */
  NInjured = 0; 
  for(i=0;i<N;i++){ Injured[i] = 0; }

  for(i=0;i<N;i++){
	  Phi[i] = DirectionOfExit( i );
	  Vdir[i] = Phi[i];
	  V[i]=0.0;
	  V0of[i]=V0;
	  VX[i]=0.0;
	  VY[i]=0.0;
  }
}

/*------------------------------*/

void Upd(){

  /* one parallel update step 
     using Euler's method with adaptive stepsize 
     
     0, allocating memory to local arrays
     
     1, forces
        1, walls
	2, wpoints
	3, pairs
	4, column
        5, smoke force & injuries

     2, 
        1. preparing update of the eq. of motion computing vxnew[i], vynew[i]
        2. if falling down is allowed, checking injuries

     3,
        1, current coordinates of particle i stored
        2, updating coordinates of particle i
            (X[i],Y[i]) = (X[i],Y[i]) + (VX[i],VY[i]) * time step  
	  (if any particles have left, storing leaving times)
        3, removing particles that have dropped off (with X[i]>XS) 
	4, updating book-keeping for remainig particles
	   
     4, 1 efficiency in this update
        2 for(i=0..N-1)
          (VX[i],VY[i]) = (VXNew[i],VYNew[i])	
	3 updating preferred directions
	4 time step added to present time value

     5, freeing mem. allocated to local arrays 
  */  


  int allocN,i,j,k,l,mx,my,m,can_see,iwp,iw,j_old,j_new;
  float *fwallx,*fwally,*fwpointx,*fwpointy,*fpairx,*fpairy,
    *fspx,*fspy,*fsumx,*fsumy,*vxnew,*vynew,tstep,tmpr,
    tmp_fpsx,tmp_fpsy,tmp_fyox,tmp_fyoy,tmp_ftax,tmp_ftay,
    tmprsqr,sqrt_fact,ksi,eta,vnew,*ftmagsum,*fsmokex,*fsmokey,
    x_smokefront,tmpf,f_over_r,scal_prod_over_rsqr,rx,ry,
    *fcolx,*fcoly;



  /* 0 */
  allocN=N;
  fwallx=vector(0,allocN-1);
  fwally=vector(0,allocN-1);
  fwpointx=vector(0,allocN-1);
  fwpointy=vector(0,allocN-1);
  fpairx=vector(0,allocN-1);
  fpairy=vector(0,allocN-1);
  fsmokex=vector(0,allocN-1);
  fsmokey=vector(0,allocN-1);

  fspx=vector(0,allocN-1);
  fspy=vector(0,allocN-1);
  fsumx=vector(0,allocN-1);
  fsumy=vector(0,allocN-1);
  vxnew=vector(0,allocN-1);
  vynew=vector(0,allocN-1);
  ftmagsum=vector(0,allocN-1);
  fcolx=vector(0,allocN-1);
  fcoly=vector(0,allocN-1);



  /* 1 */

  /* 1.0 */
  /* default values */
  tstep = DefaultDeltaT;
  for(i=0;i<N;i++){ 
          fwallx[i] = 0.0;
	  fwally[i] = 0.0;
	  fwpointx[i] = 0.0;
	  fwpointy[i] = 0.0;
	  fpairx[i] = 0.0;
	  fpairy[i] = 0.0;
          fsmokex[i] = 0.0;
          fsmokey[i] = 0.0;

	  fspx[i] = 0.0;
	  fspy[i] = 0.0;
	  fsumx[i] = 0.0;
	  fsumy[i] = 0.0;
	  ftmagsum[i] = 0.0;
	  fcolx[i] = 0.0;
	  fcoly[i] = 0.0;
  }
#ifdef _FLAG_SD_C_
  FW_x=0.0;
#endif /* _FLAG_SD_C_ */





  /* 1.1 */
  /* wall force */
  for(i=0;i<N;i++){
      for(iw=0;iw<NW;iw++){

	  WallParticleRelation(iw,i,&tmpr,&can_see);
          if((can_see==1)&&(tmpr<=R)){

	      /* init */
	      tmp_fpsx = tmp_fpsy = 0.0;
	      tmp_fyox = tmp_fyoy = 0.0;
	      tmp_ftax = tmp_ftay = 0.0;

	      /* psychological force */
	      WallPsychForce(iw,i,tmpr,&tmp_fpsx,&tmp_fpsy);
	      /* Young and tangential forces */
	      if(tmpr<=0.5*D[i]){
		  WallYoungForce(iw,i,tmpr,&tmp_fyox,&tmp_fyoy);
		  switch(FrictionSwitch){
		  case 0:{
		      WallTangForce_FS0(iw,i,&tmp_ftax,&tmp_ftay);
		      break;
		  }
		  case 1:{
		      WallTangForce_FS1(iw,i,tmpr,&tmp_ftax,&tmp_ftay);
		      break;
		  }
		  }
	      }
	      /* summing wall forces */
	      if(Injured[i]==0){
		  fwallx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
		  fwally[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
	      }
	      else /* ie. if Injured[i]=1 */ {
		  fwallx[i] += tmp_fyox + tmp_ftax;
		  fwally[i] += tmp_fyoy + tmp_ftay;
	      }

	      /* sum of magnitude of touching forces */
	      if(InjurySwitch==1){
		  ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
	      }

	      /* measuring x component of touching force exerted
		 on walls left and right from exit 
		 -- only in demo mode */ 
#ifdef _FLAG_SD_C_
	      if((iw==1)||(iw==7)){
		  FW_x -= tmp_fyox + tmp_ftax;
	      }
#endif /* _FLAG_SD_C_ */
	  }
      }
  }



  /* 1.2 */
  /* wpoint force */
  for(i=0;i<N;i++){
    for(iwp=0;iwp<NWP;iwp++){

        WPointParticleRelation(iwp,i,&tmpr,&can_see);
	if((can_see==1)&&(tmpr<=R)){

	        /* init */
	        tmp_fpsx = tmp_fpsy = 0.0;
		tmp_fyox = tmp_fyoy = 0.0;
		tmp_ftax = tmp_ftay = 0.0;

	        /* computing forces */  
	        WPointPsychForce(iwp,i,tmpr,&tmp_fpsx,&tmp_fpsy);
		if(tmpr<=0.5*D[i]){
		    WPointYoungForce(iwp,i,tmpr,&tmp_fyox,&tmp_fyoy);
		    switch(FrictionSwitch){
		    case 0:{
			WPointTangForce_FS0(iwp,i,tmpr,&tmp_ftax,&tmp_ftay);
			break;
		    }
		    case 1:{
			WPointTangForce_FS1(iwp,i,tmpr,&tmp_ftax,&tmp_ftay);
			break;
		    }
		    }
		}

		/* summing forces */
		if(Injured[i]==0){
		    fwpointx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
		    fwpointy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
		}
		else /* ie. if Injured[i]=1 */ {
		    fwpointx[i] += tmp_fyox + tmp_ftax;
		    fwpointy[i] += tmp_fyoy + tmp_ftay;
		}

		/* sum of magnitude of touching forces */
		if(InjurySwitch==1){
		    ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
		}

		/* measuring x component of touching force exerted
		   on walls left and right from exit 
		   -- only in demo mode */ 
#ifdef _FLAG_SD_C_
		if((iwp==0)||(iwp==3)){
		    FW_x -= tmp_fyox + tmp_ftax;
		}
#endif /* _FLAG_SD_C_ */
	}
    }
  }



  /* 1.3 */ 
  /* particle-particle forces */
  for(i=0; i<N; i++) {  

    j = (int)floor(X[i]*GX/XS) + G * (int)floor(Y[i]*GY/YS);
    for(k=-1;k<=1;k++){ for(l=-1;l<=1;l++){

	mx = j%G+k; 
	my = j/G+l;
	if((mx>=0)&&(mx<GX)&&(my>=0)&&(my<GY)){

	      m = BIndBd[ (mx+GX)%GX + G * (my%GY) ];
	      /* checking each pair of particles only once */
	      while(m>=i) { m = BInd[m]; }
	      if(m!=-1) {
		do { 

		  tmprsqr = SQR(X[i]-X[m]) + SQR(Y[i]-Y[m]);
		  if( tmprsqr <= SQR(R) ){
		      tmpr = sqrt(tmprsqr);
		      
		      /* init */
		      tmp_fpsx = tmp_fpsy = 0.0;
		      tmp_fyox = tmp_fyoy = 0.0;
		      tmp_ftax = tmp_ftay = 0.0;

		      /* pair forces */
		      /* Force(i,m,...) gives the force exerted by m
			 on i, all forces are symmetric now */
		      PP_PsychForce(i,m,tmpr,&tmp_fpsx,&tmp_fpsy);
		      if(tmpr<=0.5*(D[i]+D[m])){
			  PP_YoungForce(i,m,tmpr,&tmp_fyox,&tmp_fyoy);
			  switch(FrictionSwitch){
			  case 0:{
			      PP_TangForce_FS0(i,m,tmpr,&tmp_ftax,&tmp_ftay);
			      break;
			  }
			  case 1:{
			      PP_TangForce_FS1(i,m,tmpr,&tmp_ftax,&tmp_ftay);
			      break;
			  }
			  }
		      }

		      /* summing forces */
		      if(Injured[i]==0){
			  fpairx[i] += tmp_fpsx + tmp_fyox + tmp_ftax;
			  fpairy[i] += tmp_fpsy + tmp_fyoy + tmp_ftay;
		      }
		      else{ /* ie. if Injured[i]=1 */
			  fpairx[i] += tmp_fyox + tmp_ftax;
			  fpairy[i] += tmp_fyoy + tmp_ftay;
		      }
		      if(Injured[m]==0){
			  fpairx[m] -= tmp_fpsx + tmp_fyox + tmp_ftax;
			  fpairy[m] -= tmp_fpsy + tmp_fyoy + tmp_ftay;
		      }
		      else{ /* ie. if Injured[m]=1 */
			  fpairx[m] -= tmp_fyox + tmp_ftax;
			  fpairy[m] -= tmp_fyoy + tmp_ftay;
		      }

		      /* sum of magnitude of touching forces */
		      if(InjurySwitch==1){
			  ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
			  ftmagsum[m] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
		      }
		  }

		  m = BInd[m]; 
		  while(m>=i) { m = BInd[m]; }
		      
		}while(m!=-1);
	      }
	}
    }}
  }



  /* 1.4 
   * column
   */
  switch(ColumnSwitch){
  default: case 0:{ 
      for(i=0;i<N;i++){
	  fcolx[i] = fcoly[i] = 0.0;
      }
      break; 
  }
  case 1:{
      for(i=0;i<N;i++){
	  tmprsqr = SQR(X[i]-ColumnCenterX)+SQR(Y[i]-ColumnCenterY);
	  if(tmprsqr<=SQR(R)){
	      tmpr=sqrt(tmprsqr);
	      
	      /* init */
	      tmp_fpsx = tmp_fpsy = 0.0;
	      tmp_fyox = tmp_fyoy = 0.0;
	      tmp_ftax = tmp_ftay = 0.0;
	      
	      /* computing forces */
	      /* psychological */
	      f_over_r = A * exp(-(tmpr-0.5*(D[i]+ColumnD))/B) / tmpr;
	      tmp_fpsx = (X[i]-ColumnCenterX) * f_over_r;
	      tmp_fpsy = (Y[i]-ColumnCenterY) * f_over_r;
	      /* touching */
	      if(tmpr<=0.5*(D[i]+ColumnD)){
		  /* Young */
		  f_over_r = 2.0*C_Young*(0.5*(D[i]+ColumnD)-tmpr) / tmpr;
		  tmp_fyox = (X[i]-ColumnCenterX) * f_over_r;
		  tmp_fyoy = (Y[i]-ColumnCenterY) * f_over_r;
		  /* friction */
		  rx = X[i]-ColumnCenterX;
		  ry = Y[i]-ColumnCenterY;
		  scal_prod_over_rsqr = (ry*VX[i] - rx*VY[i]) / SQR(tmpr);
		  switch(FrictionSwitch){
		  case 0:{
		      tmp_ftax = -Gamma * (   ry * scal_prod_over_rsqr );
		      tmp_ftay = -Gamma * ( - rx * scal_prod_over_rsqr );
		      break;
		  }
		  case 1:{
		      tmp_ftax =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr) 
			         * (   ry * scal_prod_over_rsqr );
		      tmp_ftay =   -Kappa * (0.5*(D[i]+ColumnD)-tmpr) 
			         * ( - rx * scal_prod_over_rsqr );
		      break;
		  }
		  }
	      }


	      /* summing forces */
	      if(Injured[i]==0){
		      fcolx[i] = tmp_fpsx + tmp_fyox + tmp_ftax;
		      fcoly[i] = tmp_fpsy + tmp_fyoy + tmp_ftay;
	      }
	      else /* ie. if Injured[i]==1 */ {
		      fcolx[i] = tmp_fyox + tmp_ftax;
		      fcoly[i] = tmp_fyoy + tmp_ftay;
	      }

	      
	      /* sum of magnitude of touching forces */
	      if(InjurySwitch==1){
		      ftmagsum[i] += sqrt(SQR(tmp_fyox)+SQR(tmp_fyoy));
	      }
	  }
      }
      break;
  }
  }
  


  /* 1.5 */
  /* injuries */

  switch(InjurySwitch){
  case 0: { break; }
  case 1:{ 
      
      /* case: people crushed */
      for(i=0;i<N;i++){

	  /* newly injured */
	  if((ftmagsum[i]>FCrush_over_1m*PI*D[i])&&(Injured[i]==0)){
	      Injured[i] = 1;
	      NInjured++;
	      V0of[i] = 0.0;
	  }
      }
      break; 
  }
  case 2: case 3: { 

      /* case: smoke front */
      if(SimTime[UpdNum]>=SmokeStartTime){
	  x_smokefront = (SimTime[UpdNum]-SmokeStartTime)*VSmoke;

	  for(i=0;i<N;i++){
	      /* checking position compared to smoke front */
	      tmpr = X[i] - x_smokefront;

	      /* center of particle behind smoke front: injured */
	      if( tmpr < 0.5*D[i] ){ 
		  if(Injured[i]==0){ 
		      Injured[i] = 1;
		      NInjured++;
		      V0of[i] = 0.0;
		      VX[i] = VY[i] = 0.0;
		  }
	      }
	      /* ahead of front but within its interaction range:
		 trying to escape */ 
	      if( (tmpr>=0.5*D[i])&&(tmpr<=R) ){
		  tmpf = A_fire*exp(-(tmpr-0.5*D[i])/B_fire);
		  fsmokex[i] += cos(Phi[i])*tmpf;
		  fsmokey[i] += sin(Phi[i])*tmpf;
	      }
	  }
      }
      break;
  }
  }



  /* 2 */

  /* 2.1 preparing update of the eq. of motion */

  sqrt_fact = sqrt(tstep/DefaultDeltaT);
  for(i=0;i<N;i++) { 

          /* self-propelling */
          fspx[i] = 1/Tau * (V0of[i]*cos(Phi[i]) - VX[i]);
	  fspy[i] = 1/Tau * (V0of[i]*sin(Phi[i]) - VY[i]);

	  /* noise */
	  if(GaTh!=0.0){ 
	          ksi = GaussRand(GaMe, GaTh, GaCM);
		  eta = 2.0*PI * rand() / (RAND_MAX+1.0);
	  }
	  else{ ksi=0.0; eta=0.0; }


	  /* sum of forces */
	  fsumx[i] =   fspx[i] + fpairx[i] + fwallx[i] + fwpointx[i] 
	             + sqrt_fact * ksi * cos(eta);
	  fsumy[i] =   fspy[i] + fpairy[i] + fwally[i] + fwpointy[i] 
	             + sqrt_fact * ksi * sin(eta);
	  

	  /* adding smoke force */
	  if((InjurySwitch==2)||(InjurySwitch==3)){
	      fsumx[i] += fsmokex[i];
	      fsumy[i] += fsmokey[i];
	  }
	  /* adding force of column */
	  switch(ColumnSwitch){
	  default: case 0:{ break; }
	  case 1:{
	      fsumx[i] += fcolx[i];
	      fsumy[i] += fcoly[i];
	      break;
	  }
	  }


	  /* time step adjustment for velocity change */
	  EulTStep( &tstep, sqrt(SQR(fsumx[i])+SQR(fsumy[i])) );


	  /* new velocity */
	  if(  (Injured[i]==1)
	     &&((InjurySwitch==1)||(InjurySwitch==3))
	    ){
	          vxnew[i] = 0.0;
		  vynew[i] = 0.0;
	  }
	  else{
	          vxnew[i] = VX[i] + fsumx[i] * tstep;
		  vynew[i] = VY[i] + fsumy[i] * tstep;
	  }

	    
	  /* checking new velocity */
	  vnew = sqrt( SQR(vxnew[i]) + SQR(vynew[i]) );
	  if(vnew > Vmax) {
	          vxnew[i] = vxnew[i]/vnew * Vmax;
		  vynew[i] = vynew[i]/vnew * Vmax;
	  }
  }




  /* 3 */

  for(i=0; i<N; i++) { 

          /* .1 */
	  Xprev[i] = X[i];
	  Yprev[i] = Y[i];

          /* .2 */
	  X[i] += VX[i] * tstep;
	  Y[i] += VY[i] * tstep;

	  if((Xprev[i]>RoomXSize)&&(X[i]<=RoomXSize)){ NInRoom++; }
	  if((Xprev[i]<=RoomXSize)&&(X[i]>RoomXSize)){
	      NInRoom--;
#ifdef _FLAG_SD_CRUNCH_C_
	      Tleave[RndSeed_main_index][VSmoke_main_index][NInRoom] = 
		  SimTime[UpdNum];
#endif /* _FLAG_SD_CRUNCH_C_ */
	  }
  }

  
  /* .3 and .4 */
  for(i=0;i<N;i++){

	  /* (a) if the particle is on the board, its book-keeping
	     arrays are modified only if its block has changed during
	     the last update
	     (b) if the particle is off-board, it will be removed */


          /* a */
	  if(X[i]<XS){
	          j_old =   (int)floor(Xprev[i]*GX/XS) 
		          + G*(int)floor(Yprev[i]*GY/YS);
	          j_new = (int)floor(X[i]*GX/XS) + G*(int)floor(Y[i]*GY/YS);
		  if( j_new != j_old ) {

		      /* deleting particle i from its old block */
		      j = j_old;
		      if(BIndBd[j]==i) {
		          BIndBd[j] = BInd[i];
		      }
		      else {
		          j = BIndBd[j];
			  while(BInd[j]!=i) {
			          j = BInd[j];
			  }
			  BInd[j] = BInd[i];
		      }
		  

		      /* inserting particle i into its new block */
		      j = j_new;
		      if(BIndBd[j]==-1) {
		          BIndBd[j] = i;
			  BInd[i] = -1;
		      }
		      else {
		          j = BIndBd[j];
			  while(BInd[j]!=-1) {
			          j = BInd[j];
			  }
			  BInd[j] = i;
			  BInd[i] = -1;
		      }
		  }	  
	  }
	  else{ 
	          RemoveParticle( &N, i ); 
		  i--;
	  }
  }
  

  
  /* 4 */

  /* 4.1 */
  E[UpdNum+1] = 0.0;
  for(i=0;i<N;i++) {
          E[UpdNum+1] += VX[i] * cos(Phi[i]) + VY[i] * sin(Phi[i]); 
  }
  if(N>0){ E[UpdNum+1] /= N; }


  /* 4.2 */
  for(i=0;i<N;i++){
          VX[i] = vxnew[i];
          VY[i] = vynew[i];
	  V[i] = sqrt(SQR(VX[i])+SQR(VY[i]));
	  Vdir[i] = atan2(VY[i],VX[i]);
          Phi[i] = DirectionOfExit( i );
  }


  /* 4.3 */
  SimTime[UpdNum+1] = SimTime[UpdNum] + tstep;
  UpdNum++;


  /*
if(NInjured>0){
  fprintf(stdout,"t[%d]=%g\n",UpdNum,SimTime[UpdNum]);
  fflush(stdout);
}
*/


  /* 5 */
  free_vector(fwallx,0,allocN-1);
  free_vector(fwally,0,allocN-1);
  free_vector(fwpointx,0,allocN-1);
  free_vector(fwpointy,0,allocN-1);
  free_vector(fpairx,0,allocN-1);
  free_vector(fpairy,0,allocN-1);
  free_vector(fsmokex,0,allocN-1);
  free_vector(fsmokey,0,allocN-1);

  free_vector(fspx,0,allocN-1);
  free_vector(fspy,0,allocN-1);
  free_vector(fsumx,0,allocN-1);
  free_vector(fsumy,0,allocN-1);
  free_vector(vxnew,0,allocN-1);
  free_vector(vynew,0,allocN-1);
  free_vector(ftmagsum,0,allocN-1);
  free_vector(fcolx,0,allocN-1);
  free_vector(fcoly,0,allocN-1);


  /* 6 */
}

/*------------------------------*/

float GaussRand( float gmean, float gtheta, float gcutmult ){
  /* generates a random number (x) with
     P(x) = exp[- (x-gmean)^2 / (2*gtheta)], if x is in 
            [gmean - gcutmult*sqrt(gtheta), gmean + gcutmult*sqrt(gtheta)]
          = 0                              , if not */
     
  if( (GaussFlag==1) && (fabs(GaussSet2-gmean) <= gcutmult*sqrt(gtheta)) ) {
          GaussFlag = 0;
	  return GaussSet2;
  }
  else {
          float v1,v2,rsq,fac;

	  GaussFlag = 0;
	  do {
	          do {
		          v1 = 1.0 - 2.0*(rand()/(RAND_MAX+1.0));
			  v2 = 1.0 - 2.0*(rand()/(RAND_MAX+1.0));
		  } while((rsq=v1*v1+v2*v2) >= 1.0);
		  fac = sqrt(-2.0*gtheta*log(rsq)/rsq);
		  GaussSet1 = v1*fac;
		  GaussSet2 = v2*fac;
	  } while(    (fabs(GaussSet1-gmean) > gcutmult*sqrt(gtheta)) 
		   && (fabs(GaussSet2-gmean) > gcutmult*sqrt(gtheta)) );

	  if(fabs(GaussSet1-gmean) <= gcutmult*sqrt(gtheta)) {
	          GaussFlag = 1;
		  return GaussSet1;
	  }
	  else {
	          GaussFlag = 0;
		  return GaussSet2;
	  }
  } 
} 

/*------------------------------*/

float EMean( char* sw, int unfreq, float stfreq ) {
  /* calculates the mean value of the efficiency of the system for the last
     few update steps -- NOTE: use this function only when UpdNum > 0 

     if unfreq != 0, the average will be calculated for the last unfreq
     updates (the present one included)
     if unfreq == 0, the average will be calculated for the shortest
     possible time interval exceeding stfreq */

  int i, start;
  float e_mean, f;


  if(strcmp(sw,"un")==0) { start = UpdNum - unfreq; }
  else /* i.e. if(strcmp(sw,"st")==0) */ { 
          start = Mb; /* start from beginning of present time window */
	  f = floor( SimTime[UpdNum] / stfreq );
	  while( f - floor( SimTime[start] / stfreq ) > 1.0 ) { start++; }
	  if( start==UpdNum ) { start--; }
  }
  e_mean = 0.0;
  for(i=start+1; i<=UpdNum; i++) {
          e_mean += E[i] * ( SimTime[i] - SimTime[i-1] );
  }	  
  e_mean /= SimTime[UpdNum] - SimTime[start];


  e_mean /= V0; 
  return e_mean;
}


/*==============================*/

/*******************************/
#ifdef _FLAG_SD_C_
/*******************************/

void Start_Demo( int narg, char *argstr[] ){

  /* standard start */
  Start_Bare( narg, argstr );

  /* additional start */
  /* ... */
}

/*--------------------------------------------------*/

void Init_Demo(){
  /* 1 general
     2 special 
     */


  _E("Initializing, please wait... \n"); 


  /* 1 */
  Init_Bare("demo");



  /* 2 */

  /* opening files */
  if(!(OFP=fopen(OFN,"w"))){
          fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",OFN);
	  SD_LIB_EXIT;
  }
  fprintf(OFP,"UpdNum, SimTime, N, <E>\n");
  fflush(OFP);

  if(!(OF2P=fopen(OF2N,"w"))){
          fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",OF2N);
	  SD_LIB_EXIT;
  }
  fprintf(OF2P,"\n");
  fflush(OF2P);



  /* init visual or data output */
  switch(MainSwitch){
  case 0: default: {
      X11_init();
      break;
  }
  case 1:{
      Eps_init();
      break;
  }
  case 2:{
      Java_init();
      break;
  }
  }


  _E("... finished.\n");
}

/*------------------------------*/

void X11_init(){
  int ii,last_ok;
  XColor sdef,edef;
  

  /* general */
  X11_WWi = X11_InFW + (int)(X11_Magn*XS) + X11_Margin;
  X11_WHe = (int)MAX( X11_InFH, X11_Magn*YS+X11_GrFH + 3*X11_Margin );
  g_win( "open", " self-driven", "sd", 0, 0, X11_WWi, X11_WHe, 4);
  g_font( "open", X11_FontName );


  /* colors */
  if( !XAllocNamedColor(display,cmap,BackGroundColorName,&edef,&sdef) ) {
      fprintf(stderr,"Error: couldn't allocate color: %s\n",
	      BackGroundColorName);
      SD_LIB_EXIT;
  }
  BGCCode = sdef.pixel;

  if( !XAllocNamedColor(display,cmap,InfoColorName,&edef,&sdef) ) {
      fprintf(stderr,"Error: couldn't allocate color: %s\n",InfoColorName);
      SD_LIB_EXIT;
  }
  ICCode = sdef.pixel;
	  

  PaCNum = sizeof(ParticleColorName)/sizeof(char*);
  PaCCode = ivector(0,PaCNum-1);
  for(ii=0,last_ok=0;ii<PaCNum;ii++){
      if( !XAllocNamedColor(display,cmap,ParticleColorName[ii],&edef,&sdef) ) {
	  fprintf(stderr,"WARNING: couldn't allocate color: %s\n",
		  ParticleColorName[ii]);
	  fprintf(stderr,"Using %s instead\n",ParticleColorName[last_ok]);
	  PaCCode[ii]=PaCCode[last_ok];
      }
      else{
	  PaCCode[ii] = sdef.pixel;
	  last_ok=ii;
      }
  }
}

/*------------------------------------------*/

void Eps_init(){
  XColor edef;
  int ii;


  EpsMagn = MIN( (EpsXS-2.0*EpsMinXMarg-EpsInFW)/XS,
 	         (EpsYS-2.0*EpsMinYMarg)/YS
	       );
  EpsXMarg = 0.5*(EpsXS-EpsMagn*XS-EpsInFW);
  EpsYMarg = 0.5*(EpsYS-EpsMagn*YS);
  Eps_iPF = Eps_iPF_first;


  /* using X11 for a short time here (just for the color codes) */
  g_win( "open", " self-driven", "sd", 0, 0, 50, 50, 4);
  g_font( "open", X11_FontName );

  PaCNum = sizeof(ParticleColorName)/sizeof(char*);
  XLookupColor(display,cmap,BackGroundColorName,&edef,&BGC_sdef);
  XLookupColor(display,cmap,InfoColorName,&edef,&IC_sdef);
  for(ii=0;ii<PaCNum;ii++){
      XLookupColor(display,cmap,ParticleColorName[ii],&edef,&(PC_sdef[ii]));
  }
  XLookupColor(display,cmap,SmokeColorName,&edef,&SmokeColor_sdef);


  /* closing X11 */
  /* don't... the program stops, if you do this */
  /*
  g_win( "close", " self-driven", "sd", 0, 0, 50, 50, 4);
  g_font( "close", X11_FontName );
  */
}

/*------------------------------------------*/

void Java_init(){

  JavaXMarg=10;
  JavaYMarg=10;
  JavaMagn = MIN((JavaXS-20.0)/XS,(JavaYS-20.0)/YS);
  JavaImInd = 1;
  JavaNTStep = 1 + (int)ceil(JavaMaxTime/JavaTStep);
fprintf(stderr,"XS,YS,M,XM,YM=%g,%g,%g,%d,%d\n",XS,YS,JavaMagn,JavaXMarg,JavaYMarg);

  JavaXlo=imatrix(1,N0,1,JavaNTStep);
  JavaXhi=imatrix(1,N0,1,JavaNTStep);
  JavaYlo=imatrix(1,N0,1,JavaNTStep);
  JavaYhi=imatrix(1,N0,1,JavaNTStep);
  JavaC=imatrix(1,N0,1,JavaNTStep);
  JavaD=imatrix(1,N0,1,JavaNTStep);
}

/*------------------------------------------*/

void Pic(){

    if(  (UpdNum==0)
       ||(  (UpdNum>0) 	 
          &&(  (  (DrawUN != 0)
	        &&(UpdNum % DrawUN == 0)  
	       )
             ||(  (DrawUN == 0)
	        &&(   floor( SimTime[UpdNum] / DrawST )
		    > floor( SimTime[UpdNum-1] / DrawST )
	          )     
	       )
	    ) 
	 )
      ) {
            switch(MainSwitch){
	    case 0: default: {
	            X11_Pic();
		    break;
	    }
	    case 1: {
	            Eps_Pic();
		    break;
	    }
	    }
    }
}

/*------------------------------*/

void X11_Pic(){
  /* 1 cleaning the whole window
     2 drawing particles (smoke front, column)
     3 walls
     4 cleaning the info surface, drawing info 
     5 showing it, time delay
     */

  int i,disp_height;
  char disp_str[MY_STRLEN];
  /*  float pmean;*/
  float x;


  /* 1 */
  XSetForeground( display, gc, BGCCode );
  XFillRectangle( display, pix1, gc, 0, 0, X11_WWi, X11_WHe );


  /* 2 */
  for(i=0;i<N;i++){
	  XDrawParticle( i, X11_InFW, X11_Margin, X11_Magn); 
  }
  


  /* 2.B */
  /* smoke front, if needed */
  if(  ((InjurySwitch==2)||(InjurySwitch==3))
     &&(SimTime[UpdNum]>=SmokeStartTime)){
      XSetForeground( display, gc, ICCode );
      x =   X11_InFW + X11_Magn*(SimTime[UpdNum]-SmokeStartTime)*VSmoke;
      for(i=0;i<=X11_Magn*YS/6.0;i++) { 
	  XDrawLine(display, pix1, gc, 
		    x, X11_Margin + 6*i,
		    x, X11_Margin + 6*i+3
		   ); 
      }
  }

  /* 2.C */
  /* column */
  switch(ColumnSwitch){
  default: case 0:{ break; }
  case 1:{
          XSetForeground( display, gc, ICCode );          
	  XDrawArc(display, pix1, gc,
		   (int)floor(X11_InFW+X11_Magn*(ColumnCenterX-0.5*ColumnD)),
		   (int)floor(X11_Margin+X11_Magn*(ColumnCenterY-0.5*ColumnD)),
		   (int)floor(X11_Magn*ColumnD),
		   (int)floor(X11_Magn*ColumnD),
		   0, 23040		   
		  );
	  break;
  }
  }



  /* 3 */
  XSetForeground( display, gc, ICCode );
  for(i=0;i<NW;i++){
          XDrawLine(display,pix1,gc,
		    (int)floor(X11_InFW+X11_Magn*W[i].x1),
		    (int)floor(X11_Margin+X11_Magn*W[i].y1),
		    (int)floor(X11_InFW+X11_Magn*W[i].x2),
		    (int)floor(X11_Margin+X11_Magn*W[i].y2)
		   );
  }





  /* 4 */
  XSetForeground( display, gc, BGCCode );

  /* cleaning the x=XS end of the field to allow particles 
     leave the screen gradually */
  XFillRectangle( display, pix1, gc, 
		  (int)floor(X11_InFW+X11_Magn*XS), 0, 
		  (int)floor(X11_WWi-X11_InFW-X11_Magn*XS), X11_WHe );

  /* writing info */
  XSetForeground( display, gc, ICCode );
  disp_height = X11_Margin + X11_TLH;

  disp_height += X11_TLH;
  sprintf( disp_str, "t [%6d] = %.1f", UpdNum, SimTime[UpdNum] );
  XDrawString( display, pix1, gc, X11_Margin, disp_height, 
	       disp_str, (signed int)strlen(disp_str) );

  disp_height += X11_TLH;
  sprintf( disp_str, "N (in room) = %d", NInRoom );
  XDrawString( display, pix1, gc, X11_Margin, disp_height, 
	       disp_str, (signed int)strlen(disp_str) );

  disp_height += X11_TLH;
  sprintf( disp_str, "N_injured = %d", NInjured );
  XDrawString( display, pix1, gc, X11_Margin, disp_height, 
	       disp_str, (signed int)strlen(disp_str) );

  disp_height += X11_TLH;
  sprintf( disp_str, "V0 = %g", V0 );
  XDrawString( display, pix1, gc, X11_Margin, disp_height, 
	       disp_str, (signed int)strlen(disp_str) );

  disp_height += X11_TLH;
  sprintf( disp_str, "FWall_x = %.2f", FW_x );
  XDrawString( display, pix1, gc, X11_Margin, disp_height, 
	       disp_str, (signed int)strlen(disp_str) );



  /* 5 */
  h_show(X11_WWi,X11_WHe);
  sleep(Sleep);
}

/*------------------------------*/

void EpsDrawParticle( FILE* fp, int i ) {
  /* putting particle i onto the eps image 
    pm: picture magnification
    */
  
  float theta,v,x,y,vx,vy,d,magn,lxm,uym;


  theta = atan2(vy,vx); 
  v = sqrt(SQR(vx)+SQR(vy));
  x = X[i];
  y=Y[i];
  vx=VX[i];
  vy=VY[i];
  d=D[i]; 
  lxm = EpsXMarg + EpsInFW;
  uym = EpsYMarg;
  magn = EpsMagn;



  /* drawing particle */
  switch( Draw ) {
  default: case 0: {

          d = D[i];
	  x = X[i];
	  y = Y[i];

          EpsFillCircle( fp,
			 lxm + magn * x,          
			 uym + magn * y, 
			 0.5*magn*d
		       );
	  break;
  }
  case 1: {

          d = D[i];
	  x = X[i];
	  y = Y[i];

          EpsFillCircle( fp,
			 lxm + magn * x,          
			 uym + magn * y, 
			 0.5*magn*d
		       );
	  EpsSetRgb( fp, 
		     (double)(IC_sdef.red), 
		     (double)(IC_sdef.green), 
		     (double)(IC_sdef.blue) 
		   );	      
	  if(EpsLineWidth>0.0){
	      EpsSetLinewidth( fp, EpsLineWidth );
	      EpsDrawCircle( fp,
			     lxm + magn * x,          
			     uym + magn * y, 
			     0.5*magn*d
			   );
	  }
	  break;
  }
  case 2: {

          d = DrawDMult * D[i];
	  x = X[i];
	  y = Y[i];

          EpsFillCircle( fp,
			 lxm + magn * x,          
			 uym + magn * y, 
			 0.5*magn*d
		       );
	  EpsSetRgb( fp, 
		     (double)(IC_sdef.red), 
		     (double)(IC_sdef.green), 
		     (double)(IC_sdef.blue) 
		   );	      
	  if(EpsLineWidth>0.0){
	      EpsSetLinewidth( fp, EpsLineWidth );
	      EpsDrawCircle( fp,
			     lxm + magn * x,          
			     uym + magn * y, 
			     0.5*magn*d
			   );
	  }
	  break;
  }
  }
}

/*--------------------------------------*/

void Eps_Pic(){
  /* writing the current configuration in eps format into EpsPicMult
     files (identical copies) */

  int ii,i;
  float lw,disp_height,x_smokefront;
  char fn[MY_STRLEN],sh_com[MY_STRLEN],disp_str[MY_STRLEN];
  FILE *fp;



  /* 1 opening file */
  sprintf(fn,"sd.%d.eps",Eps_iPF);
  if(!(fp=fopen(fn,"w"))){
      fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",fn);
      SD_LIB_EXIT;
  }
  fprintf(stdout,"Writing %s...\n",fn);
  

  /* 2 init and background */
  EpsInit( fp, 0, 0, EpsXS-1, EpsYS-1 );
  EpsSetFont( fp, "Times-Bold", EpsInTH );
  EpsSetRgb( fp, BGC_sdef.red, BGC_sdef.green,
	     BGC_sdef.blue );
  EpsFillRectangle( fp, 0, 0, EpsXS, EpsYS );

  

  /* 2b smoke, if needed */
  if((InjurySwitch==2)||(InjurySwitch==3)){
      if(SimTime[UpdNum]>=SmokeStartTime){
	  x_smokefront = (SimTime[UpdNum]-SmokeStartTime)*VSmoke;
      
	  EpsSetRgb( fp, 
		     SmokeColor_sdef.red,
		     SmokeColor_sdef.green,
		     SmokeColor_sdef.blue
		   );
	  EpsFillRectangle( fp, 
			    EpsXMarg+EpsInFW, 
			    EpsYMarg, 
			    EpsXMarg+EpsInFW+EpsMagn*x_smokefront, 
			    EpsYMarg+EpsMagn*RoomYSize
			  );
      }
  }




  /* 3 particles */

  /* 3a: not injured */
  EpsSetRgb( fp,  
	     (double)(PC_sdef[0].red), 
	     (double)(PC_sdef[0].green), 
	     (double)(PC_sdef[0].blue)
	   );	      
  for(i=0;i<N;i++){ 
      if(Injured[i]==0){
	  EpsDrawParticle( fp, i ); 
      }
  }

  /* 3b: injured */
  EpsSetRgb( fp,  
	     (double)(PC_sdef[1].red), 
	     (double)(PC_sdef[1].green), 
	     (double)(PC_sdef[1].blue)
	   );	      
  for(i=0;i<N;i++){ 
      if(Injured[i]==1){
	  EpsDrawParticle( fp, i ); 
      }
  }



  /* 4 drawing walls */
  EpsSetRgb( fp, IC_sdef.red, IC_sdef.green, IC_sdef.blue );
  EpsSetLinewidth( fp, 0.0 );
  lw = EpsMagn*0.5*WallWidth;
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW-lw,
		    EpsYMarg-lw,
		    EpsXMarg+EpsInFW,
		    EpsYMarg+EpsMagn*YS+lw
		  );
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW-lw,
		    EpsYMarg-lw,
		    EpsXMarg+EpsInFW+EpsMagn*XS,
		    EpsYMarg
		  );
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW-lw,
		    EpsYMarg+EpsMagn*YS,
		    EpsXMarg+EpsInFW+EpsMagn*XS,
		    EpsYMarg+EpsMagn*YS+lw
		  );
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW+EpsMagn*RoomXSize,
		    EpsYMarg-lw,
		    EpsXMarg+EpsInFW+EpsMagn*(RoomXSize+WallWidth),
		    EpsYMarg+EpsMagn*0.5*(RoomYSize-DoorWidth)
		  );
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW+EpsMagn*RoomXSize,
		    EpsYMarg+EpsMagn*0.5*(RoomYSize+DoorWidth),
		    EpsXMarg+EpsInFW+EpsMagn*(RoomXSize+WallWidth),
		    EpsYMarg+EpsMagn*RoomYSize+lw
		  );

  /* drawing column, if needed */
  if(ColumnSwitch==1){
          EpsSetRgb( fp, IC_sdef.red, IC_sdef.green, IC_sdef.blue );
          EpsFillCircle( fp,
			 EpsXMarg+EpsInFW+EpsMagn*ColumnCenterX,
			 EpsYMarg+EpsMagn*ColumnCenterY,
			 0.5*EpsMagn*ColumnD
		       );
  }





  /* 5 removing rubbish from outer rim */
  EpsSetRgb( fp, BGC_sdef.red, BGC_sdef.green, BGC_sdef.blue ); 
  EpsFillRectangle( fp, 
		    EpsXMarg+EpsInFW+EpsMagn*XS,
		    EpsYMarg,
		    EpsXS,
		    EpsYMarg+EpsMagn*YS
		  );


  /* 6 text */
  EpsSetRgb( fp, IC_sdef.red, IC_sdef.green, IC_sdef.blue ); 
  sprintf( disp_str,"t = %d", (int)floor(SimTime[UpdNum]) );
  disp_height = EpsYS - (EpsYMarg + 1.2*EpsInTH);
  EpsDrawString( fp, 0.0, EpsXMarg, disp_height, disp_str );
  
  disp_height -= 1.2*EpsInTH;
  sprintf( disp_str, "N = %d", NInRoom );
  EpsDrawString( fp, 0.0, EpsXMarg, disp_height, disp_str );

  disp_height -= 1.2*EpsInTH;
  sprintf( disp_str, "V0 = %g", V0 );
  EpsDrawString( fp, 0.0, EpsXMarg, disp_height, disp_str );

  /* number of injured */
  if((InjurySwitch==1)||(InjurySwitch==2)||(InjurySwitch==3)){
      disp_height -= 1.2*EpsInTH;
      sprintf( disp_str, "Inj.: %d", NInjured );
      EpsDrawString( fp, 0.0, EpsXMarg, disp_height, disp_str );
  }      



  /* 10 closing file */
  fflush(fp);
  fclose(fp);
  fprintf(stdout,"finished.\n");
  Eps_iPF++;

  /* 11 creating copies of eps file */
  if(EpsPicMult>1){
    for(ii=1;ii<=EpsPicMult-1;ii++){
        sprintf(sh_com,"ln -s sd.%d.eps sd.%d.eps",Eps_iPF-1,Eps_iPF-1+ii);
	fprintf(stdout,"%s\n",sh_com);
	system(sh_com);
	_O("ok.\n");
    }	 
    Eps_iPF+=EpsPicMult-1;
  }
}

/*------------------------------*/

void XDrawParticle( int i, int leftxmargin, int upymargin, float magn){ 

  /* - drawing the particle  */

  int lxm = leftxmargin, uym = upymargin;
  float d,x,y;



  /* particle color */
  switch(Injured[i]){
  case 0: { XSetForeground( display, gc, PaCCode[0] ); break; }
  case 1: { XSetForeground( display, gc, PaCCode[1] ); break; }
  }




  /* drawing particle */
  switch( Draw ) {
  default: case 0: {

          d = D[i];
	  x = X[i];
	  y = Y[i];

          XFillArc(display, pix1, gc, 
		   (int)floor(lxm + magn * (x - d/2)), 
		   (int)floor(uym + magn * (y - d/2)), 
		   (int)floor(magn * d), 
		   (int)floor(magn * d), 
		   0, 23040
		   );

	  break;
  }
  case 1: {

          d = D[i];
	  x = X[i];
	  y = Y[i];

          XFillArc(display, pix1, gc, 
		   (int)floor(lxm + magn * (x - d/2)), 
		   (int)floor(uym + magn * (y - d/2)), 
		   (int)floor(magn * d), 
		   (int)floor(magn * d), 
		   0, 23040
		   );

	  XSetForeground( display, gc, ICCode );
          XDrawArc(display, pix1, gc, 
		   (int)floor(lxm + magn * (x - d/2)), 
		   (int)floor(uym + magn * (y - d/2)), 
		   (int)floor(magn * d), 
		   (int)floor(magn * d), 
		   0, 23040
		   );

	  break;
  }
  case 2: {

          d = DrawDMult * D[i];
	  x = X[i];
	  y = Y[i];

          XFillArc(display, pix1, gc, 
		   (int)floor(lxm + magn * (x - d/2)), 
		   (int)floor(uym + magn * (y - d/2)), 
		   (int)floor(magn * d), 
		   (int)floor(magn * d), 
		   0, 23040
		   );

	  XSetForeground( display, gc, ICCode );
          XDrawArc(display, pix1, gc, 
		   (int)floor(lxm + magn * (x - d/2)), 
		   (int)floor(uym + magn * (y - d/2)), 
		   (int)floor(magn * d), 
		   (int)floor(magn * d), 
		   0, 23040
		   );

	  break;
  }
  case 3: {

          d = D[i];
	  x = X[i];
	  y = Y[i];

	  if(sqrt(SQR(VX[i])+SQR(VY[i]))>=0.5){
	          XFillArc(display, pix1, gc, 
			   (int)floor(lxm + magn * (x - d/2)), 
			   (int)floor(uym + magn * (y - d/2)), 
			   (int)floor(magn * d), 
			   (int)floor(magn * d), 
			   0, 23040
			   );
	  }

	  break;
  }
  }

}

/*-----------------------------------------------*/

void Save_Demo(){

  char sw[MY_STRLEN];
  float simtime_now,simtime_now_minus_1,e_now,e_now_minus_1,
    e_mean;

  if(  (UpdNum==0)
     ||(  (UpdNum>0) 	 
        &&(  (  (SaveUN != 0)
	      &&(UpdNum % SaveUN == 0)  
	     )
           ||(  (SaveUN == 0)
	      &&(   floor( SimTime[UpdNum] / SaveST )
		  > floor( SimTime[UpdNum-1] / SaveST )
	        )     
	      )
	   )
       )
    ) { 


          if(UpdNum>0){
	          if(SaveUN!=0){strcpy(sw,"un");}
		  else /* ie. if(SaveUN==0) */ { strcpy(sw,"st"); }
		  e_mean=EMean(sw,SaveUN,SaveST);
	  }
	  else /* ie. if(UpdNum==0) */ { e_mean=1.0; }
	  fprintf(OFP,"%d\t%g\t%d\t%g\n",
		  UpdNum,SimTime[UpdNum],N,e_mean);
	  fflush(OFP); 



	  /* closing present time window, opening new time window */
	  if(UpdNum>0){

	          simtime_now = SimTime[UpdNum];
		  simtime_now_minus_1 = SimTime[UpdNum-1];
		  e_now = E[UpdNum];
		  e_now_minus_1 = E[UpdNum-1];

		  free_vector(SimTime,Mb,Me);
		  free_vector(E,Mb,Me);

		  Mb = UpdNum-1;
		  Me = UpdNum-1 + AyS-1;
		  SimTime = vector(Mb,Me);
		  SimTime[UpdNum-1] = simtime_now_minus_1;
		  SimTime[UpdNum] = simtime_now;
		  E = vector(Mb,Me);
		  E[UpdNum-1] = e_now_minus_1;
		  E[UpdNum] = e_now;
	  }
  }
}

/*------------------------------*/

void Save_Java_in_Loop(){

  int i,intx,inty,intd,i_lastgood;
  float simtime_now,simtime_now_minus_1,e_now,e_now_minus_1,
    e_mean;

  if(  (JavaImInd<=JavaNTStep)
     &&(  (UpdNum==0)
        ||(  (UpdNum>0) 	 
           &&(   floor( SimTime[UpdNum] / JavaTStep )
	       > floor( SimTime[UpdNum-1] / JavaTStep )
	     )
          )
       )
    ) { 

          /* saving data */

          /* meaningful values (ie. particles still on the board) */
          for(i=0;i<N;i++){
	      intx=(int)floor(JavaXMarg+JavaMagn*X[i]);
	      JavaXlo[i+1][JavaImInd]=intx%256;
	      JavaXhi[i+1][JavaImInd]=intx/256;

	      inty=(int)floor(JavaYMarg+JavaMagn*Y[i]);
	      JavaYlo[i+1][JavaImInd]=inty%256;
	      JavaYhi[i+1][JavaImInd]=inty/256;

	      JavaD[i+1][JavaImInd]=(int)floor(JavaMagn*D[i]);

	      JavaC[i+1][JavaImInd]=Injured[i];
	  }
	  /* dummy values (for particles that have already left) */
	  i_lastgood=i-1;
	  for(;i<N0;i++){
	      JavaXlo[i+1][JavaImInd]=JavaXlo[i_lastgood+1][JavaImInd];
	      JavaXhi[i+1][JavaImInd]=JavaXhi[i_lastgood+1][JavaImInd];
	      JavaYlo[i+1][JavaImInd]=JavaYlo[i_lastgood+1][JavaImInd];
	      JavaYhi[i+1][JavaImInd]=JavaYhi[i_lastgood+1][JavaImInd];

	      JavaD[i+1][JavaImInd]=JavaD[i_lastgood+1][JavaImInd];
	      JavaC[i+1][JavaImInd]=JavaC[i_lastgood+1][JavaImInd];
	  }
	  JavaImInd++;



	  /* closing present time window, opening new time window */
	  if(UpdNum>0){

	          simtime_now = SimTime[UpdNum];
		  simtime_now_minus_1 = SimTime[UpdNum-1];
		  e_now = E[UpdNum];
		  e_now_minus_1 = E[UpdNum-1];

		  free_vector(SimTime,Mb,Me);
		  free_vector(E,Mb,Me);

		  Mb = UpdNum-1;
		  Me = UpdNum-1 + AyS-1;
		  SimTime = vector(Mb,Me);
		  SimTime[UpdNum-1] = simtime_now_minus_1;
		  SimTime[UpdNum] = simtime_now;
		  E = vector(Mb,Me);
		  E[UpdNum-1] = e_now_minus_1;
		  E[UpdNum] = e_now;
	  }
  }
}

/*------------------------------*/

void Save_Java_after_Loop(){
  char c,fn[MY_STRLEN];
  FILE *fp;
  int i,t,intx1,inty1,intx2,inty2;


  /****************************************************
   * DATA FILE FORMAT:
   *
   * byte order in a word: lo,hi
   *
   *
   * item                       symbol       size [bytes]  (2bytes=1word)
   * ----------------------------------------------
   * initial number of particles  
   *                            N0           2
   * number of time steps       NTSTEP       2
   * time step [msec]           TSTEP        2
   * number of lines            NL           2
   * number of boxes            NB           2
   * ---------------------------------
   * particle coordinates,diameters,colors
   * X[1][t=1],Y[1][t=1],D[1][t=1],C[1][t=1]...
   * ...X[1][t=NTSTEP],Y[1][NTSTEP],D[1][t=NTSTEP].C[1][NTSTEP]...
   * ...X[N][t=NTSTEP],Y[N][NTSTEP],D[N][t=NTSTEP],C[N][NTSTEP]
   *                                         (2+2+1+1)*N*NTSTEP
   * line coordinates (walls in this case)
   * X1(wall1),Y1(wall1),X2(wall1),Y2(wall1)...Y2(wallNW)
   *                                         8*NW
   * coordinates of boxes (blanks in this case)
   * X1(blank1),Y1(blank1),X2(blank1),Y2(blank1)...Y2(blankNB)
   *                                         8*NB
   * --------------------------------------------------
   * extra switch               ES           1
   *
   * switch(ES){
   * case 0: default: {
   *      no further bytes 
   * }
   * case 1: (smoke) {
   *      v_smoke [milli pixel/timestep]      
   *                            VSMOKE       2
   *      t_smoke_start [timestep]  
   *                            TSMOKESTART  2
   * }
   * case 4: there is a cylindrical column in the room {
   *   COLUMN_CENTER_XLO,COLUMN_CENTER_XHI,  
   *   COLUMN_CENTER_YLO,COLUMN_CENTER_YHI,
   *   COLUMN_CENTER_DLO,COLUMN_CENTER_DHI   6 
   * }
   * }
   ****************************************/



  /* init */
  sprintf(fn,"sd.javadat");
  fp=fopen(fn,"w");
  if(!fp){ 
      fprintf(stderr,"ERROR: file %s could not be opened for writing.\n",fn);
      fprintf(stderr,"       exiting to system.\n");
      exit(-1);
  }

  /* parameters */
  fprintf(fp,"%c%c",(char)(N0%256),(char)(N0/256));
  fprintf(fp,"%c%c",(char)(JavaNTStep%256),(char)(JavaNTStep/256));
  fprintf(fp,"%c%c",
	  (char)(((int)(floor(JavaTStep/0.001)))%256),
	  (char)(((int)(floor(JavaTStep/0.001)))/256)
	 );
  fprintf(fp,"%c%c",(char)(NW%256),(char)(NW/256));
  fprintf(fp,"%c%c",(char)(5),(char)(0));

  /* particle data */
  for(i=1;i<=N0;i++){
      for(t=1;t<=JavaNTStep;t++){
	  fprintf(fp,"%c%c%c%c%c%c",
		  (char)(JavaXlo[i][t]),(char)(JavaXhi[i][t]),
		  (char)(JavaYlo[i][t]),(char)(JavaYhi[i][t]),
		  (char)(JavaD[i][t]),(char)(JavaC[i][t])
		  );
      }
  }
  /* walls */
  for(i=0;i<NW;i++){
      intx1=(int)floor(JavaXMarg+JavaMagn*W[i].x1);
      inty1=(int)floor(JavaYMarg+JavaMagn*W[i].y1);
      intx2=(int)floor(JavaXMarg+JavaMagn*W[i].x2);
      inty2=(int)floor(JavaYMarg+JavaMagn*W[i].y2);
      fprintf(fp,"%c%c%c%c%c%c%c%c",
	      (char)(intx1%256),(char)(intx1/256),
	      (char)(inty1%256),(char)(inty1/256),
	      (char)(intx2%256),(char)(intx2/256),
	      (char)(inty2%256),(char)(inty2/256)
	     );
  }
  /* 1st blank */
  intx1=(int)floor(MAX(0.0,JavaXMarg+JavaMagn*(W[0].x1-0.6*(Dmean+deltaD))));
  inty1=(int)floor(MAX(0.0,JavaYMarg+JavaMagn*(W[0].y2-0.6*(Dmean+deltaD))));
  intx2=(int)floor(JavaXMarg+JavaMagn*W[0].x2);
  inty2=(int)floor(JavaYMarg+JavaMagn*W[0].y1);
  fprintf(fp,"%c%c%c%c%c%c%c%c",
	  (char)(intx1%256),(char)(intx1/256),
	  (char)(inty1%256),(char)(inty1/256),
	  (char)(intx2%256),(char)(intx2/256),
	  (char)(inty2%256),(char)(inty2/256)
	 );
  /* 2nd blank */
  intx1=(int)floor(JavaXMarg+JavaMagn*W[1].x1);
  inty1=(int)floor(JavaYMarg+JavaMagn*W[1].y1);
  intx2=(int)floor(JavaXMarg+JavaMagn*W[3].x1);
  inty2=(int)floor(JavaYMarg+JavaMagn*W[3].y1);
  fprintf(fp,"%c%c%c%c%c%c%c%c",
	  (char)(intx1%256),(char)(intx1/256),
	  (char)(inty1%256),(char)(inty1/256),
	  (char)(intx2%256),(char)(intx2/256),
	  (char)(inty2%256),(char)(inty2/256)
	 );
  /* 3rd blank */
  intx1=(int)floor(MAX(0.0,JavaXMarg+JavaMagn*(W[4].x2-0.6*(Dmean+deltaD))));
  inty1=(int)floor(JavaYMarg+JavaMagn*W[4].y2);
  intx2=(int)floor(JavaXMarg+JavaMagn*W[4].x1);
  inty2=(int)floor(MIN(JavaYS,JavaYMarg+JavaMagn*(W[4].y1+0.6*(Dmean+deltaD))));
  fprintf(fp,"%c%c%c%c%c%c%c%c",
	  (char)(intx1%256),(char)(intx1/256),
	  (char)(inty1%256),(char)(inty1/256),
	  (char)(intx2%256),(char)(intx2/256),
	  (char)(inty2%256),(char)(inty2/256)
	 );
  /* 4th blank */
  intx1=(int)floor(JavaXMarg+JavaMagn*W[7].x1);
  inty1=(int)floor(JavaYMarg+JavaMagn*W[7].y1);
  intx2=(int)floor(JavaXMarg+JavaMagn*W[5].x1);
  inty2=(int)floor(JavaYMarg+JavaMagn*W[5].y1);
  fprintf(fp,"%c%c%c%c%c%c%c%c",
	  (char)(intx1%256),(char)(intx1/256),
	  (char)(inty1%256),(char)(inty1/256),
	  (char)(intx2%256),(char)(intx2/256),
	  (char)(inty2%256),(char)(inty2/256)
	 );
  /* 5th blank */
  intx1=(int)floor(MAX(0.0,JavaXMarg+JavaMagn*(W[8].x2-0.6*(Dmean+deltaD))));
  inty1=(int)floor(JavaYMarg+JavaMagn*W[8].y2);
  intx2=(int)floor(JavaXMarg+JavaMagn*W[8].x1);
  inty2=(int)floor(JavaYMarg+JavaMagn*W[8].y1);
  fprintf(fp,"%c%c%c%c%c%c%c%c",
	  (char)(intx1%256),(char)(intx1/256),
	  (char)(inty1%256),(char)(inty1/256),
	  (char)(intx2%256),(char)(intx2/256),
	  (char)(inty2%256),(char)(inty2/256)
	 );





  /* extra switch */
  switch(InjurySwitch){
  case 0: case 1: default: {
      switch(ColumnSwitch){
      default: case 0: {
	      fprintf(fp,"%c",(char)(0)); 
	      break;
      }
      case 1:{
	      fprintf(fp,"%c",(char)(4)); 
	      fprintf(fp,"%c%c%c%c%c%c",
		    (char)(((int)floor(JavaXMarg+JavaMagn*ColumnCenterX))%256),
		    (char)(((int)floor(JavaXMarg+JavaMagn*ColumnCenterX))/256),
		    (char)(((int)floor(JavaYMarg+JavaMagn*ColumnCenterY))%256),
		    (char)(((int)floor(JavaYMarg+JavaMagn*ColumnCenterY))/256),
		    (char)(((int)floor(JavaMagn*ColumnD))%256),
		    (char)(((int)floor(JavaMagn*ColumnD))/256)
		      );
	      break;
      }
      break;
      }
  }
  case 2: case 3:{
      fprintf(fp,"%c",(char)(1)); 
      fprintf(fp,"%c%c",
	      (char)(((int)floor(1000.0*JavaMagn*VSmoke*JavaTStep))%256),
	      (char)(((int)floor(1000.0*JavaMagn*VSmoke*JavaTStep))/256)
	     );
      fprintf(fp,"%c%c",
	      (char)(((int)floor(SmokeStartTime/JavaTStep))%256),
	      (char)(((int)floor(SmokeStartTime/JavaTStep))/256)
	     );
      break;
  }
  }
}

/*------------------------------*/

void Shutdown_Demo(){}

/*------------------------------*/

void ReInit(){
  /* if parameter file has changed, re-read  */

  stat(IFN, &IFStatBuf);
  if( IFStatBuf.st_mtime != IFModTime ) { 
          IFModTime = IFStatBuf.st_mtime;
	  readpar ( "re", IFN, IPar, IParName, IParNum, 
		     FPar, FParName, FParNum, 
		     SPar, SParName, SParNum ); 
  }
}

/*******************************/
#endif /* _FLAG_SD_C_ */
/*******************************/

/*========================================*/

/*******************************/
#ifdef _FLAG_SD_CRUNCH_C_
/*******************************/

void Start_Crunch( int narg, char *argstr[] ){

  RndSeed_main_1st = RndSeed;
  /* ! this may change in 'Start_Bare' */


  /* standard start */
  Start_Bare( narg, argstr );
}

/*--------------------------------------------------*/

void Init_Crunch(){

  _E("Initializing, please wait... \n"); 
  Init_Bare("crunch");

  _E("... finished.\n");
}

/*------------------------------*/

void Save_Crunch_InLoop(){

  char sw[MY_STRLEN];
  float simtime_now,simtime_now_minus_1,
    e_now,e_now_minus_1,e_mean;

  if(  (UpdNum==0)
     ||(  (UpdNum>0) 	 
        &&(  (  (SaveUN != 0)
	      &&(UpdNum % SaveUN == 0)  
	     )
           ||(  (SaveUN == 0)
	      &&(   floor( SimTime[UpdNum] / SaveST )
		  > floor( SimTime[UpdNum-1] / SaveST )
	        )     
	      )
	   )
       )
    ) { 

          if(UpdNum>0){
	          if(SaveUN!=0){strcpy(sw,"un");}
		  else{ strcpy(sw,"st"); }
		  e_mean=EMean(sw,SaveUN,SaveST);
	  }
	  else{ e_mean=1.0; }

	  /*
	    fprintf(stdout,"%d\t%g\t%d\t%g\n",
		  UpdNum,SimTime[UpdNum],N,e_mean);
	  fflush(stdout); 
	  */


	  /* closing present time window, opening new time window */
	  if(UpdNum>0){

	          simtime_now = SimTime[UpdNum];
		  simtime_now_minus_1 = SimTime[UpdNum-1];
		  e_now = E[UpdNum];
		  e_now_minus_1 = E[UpdNum-1];

		  free_vector(SimTime,Mb,Me);
		  free_vector(E,Mb,Me);

		  Mb = UpdNum-1;
		  Me = UpdNum-1 + AyS-1;
		  SimTime = vector(Mb,Me);
		  SimTime[UpdNum-1] = simtime_now_minus_1;
		  SimTime[UpdNum] = simtime_now;
		  E = vector(Mb,Me);
		  E[UpdNum-1] = e_now_minus_1;
		  E[UpdNum] = e_now;
	  }
  }
}

/*------------------------------*/

void Save_Crunch_AfterLoop(){

  char fn[MY_STRLEN];
  FILE *fp;
  int i;


  /* 1 */
  /* ALL CASES: saving leaving times */
  sprintf(fn,"%s_RndSeed=%d_VSmoke=%g_N0=%d_",OFN,RndSeed,VSmoke,N0);
  if(!(fp=fopen(fn,"w"))){
      fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",fn);
      SD_LIB_EXIT;
  }
  fprintf(fp,"#index\tleaving time\n\n");
  for(i=0;i<N0;i++){
      fprintf(fp,"%d\t%g\n",
	      i,Tleave[RndSeed_main_index][VSmoke_main_index][N0-1-i]);
  }
  fflush(fp);
  fclose(fp);


  switch(InjurySwitch){
  case 0: { break; }
  case 1: case 2: case 3: {

      /* case 1: people crushed */
      /* case 2, case 3: smoke */

      /* saving injury data */
      sprintf(fn,"%s_RndSeed=%d_VSmoke=%g_N0=%d_",OF2N,RndSeed,VSmoke,N0);
      if(!(fp=fopen(fn,"w"))){
	  fprintf(stderr,"sd_lib.c: Couldn't open %s for writing.\n",fn);
	  SD_LIB_EXIT;
      }
      fprintf(fp,"# N pedestrians remained in the room until time t:\n");
      fprintf(fp,"# %d\t%g\n",N,SimTime[UpdNum]);
      fprintf(fp,"\n");
      fprintf(fp,"# index\tinjured?\tx\ty\td\n");
      fprintf(fp,"\n");
      for(i=0;i<N;i++){
	  fprintf(fp,"%d\t%d\t%.4f\t%.4f\t%.4f\n",
		  i,Injured[i],X[i],Y[i],D[i]
		 );
      }
      fflush(fp);
      fclose(fp);

      break;
  }
  }
}

/*------------------------------*/

void Clean_Crunch_AfterLoop(){

  free_vector(SimTime,Mb,Me);
  free_ivector(BIndBd, 0, SQR(G)-1 );
  free_ivector(BInd, 0, N0-1 );
  free_vector(D, 0, N0-1 );
  free_vector(Phi, 0, N0-1 );
  free_vector(X, 0, N0-1 );
  free_vector(Y, 0, N0-1 );
  free_vector(Xprev, 0, N0-1 );
  free_vector(Yprev, 0, N0-1 );
  free_vector(V, 0, N0-1 );
  free_vector(VX, 0, N0-1 );
  free_vector(VY, 0, N0-1 );
  free_vector(E, Mb, Me );
  free_vector(Vdir,0,N0-1);
  free_vector(V0of,0,N0-1);
  free_ivector(Injured,0,N0-1);  
}

/*------------------------------*/

void Save_Crunch_AtEnd(){}

/*******************************/
#endif /* _FLAG_SD_CRUNCH_C_ */
/*******************************/
