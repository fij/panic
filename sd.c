/* self-driven / panic / main file for demo version 
   MainSwitch options:
   0 simple demo
   1 eps images instead of X11
   2 creating data file for java
   */

/*********************/

#define _FLAG_SD_C_
#include "sd_lib.c"

int main( int NArg, char * ArgStr[] ) {

  Start_Demo( NArg, ArgStr );
  Init_Demo();

  switch(MainSwitch){
  case 0: {
    
          do { 
	          Pic();
		  Save_Demo();
		  Upd(); 
		  ReInit();

	  } while( /* (N>0)&& */
		  ( UpdNum < MaxUpdNum )
	 	   &&( SimTime[UpdNum] < MaxSimTime )
		 );
	  
	  Shutdown_Demo();
	  break;
  } 
  case 1: {
    
          do { 
	          Pic();
                  Save_Demo();
		  Upd(); 

	  } while( /* (N>0)&& */
		  ( UpdNum < MaxUpdNum )
	 	   &&( SimTime[UpdNum] < MaxSimTime )
		  &&( Eps_iPF <= Eps_iPF_max )
		 );
	  
	  Shutdown_Demo();
	  break;
  } 
  case 2: {
    

          do { 
		  Save_Java_in_Loop();
		  Upd(); 
	  } while( JavaImInd < JavaNTStep );
	  
	  Save_Java_after_Loop();
	  Shutdown_Demo();
	  break;
  } 
  default: {
         
          fprintf(stderr,"Error: no such MainSwitch: %d\n",MainSwitch);
	  _E("Exiting to system.\n");
	  exit(-1);
  }}

  return 0;
}
