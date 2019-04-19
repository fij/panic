/* self-driven / panic / main file for number-crunching version */

/*********************/

#define _FLAG_SD_CRUNCH_C_
#include "sd_lib.c"

int main( int NArg, char * ArgStr[] ) {

  Start_Crunch( NArg, ArgStr );


  /* rnd seeds */
  RndSeed_main_num=1;
  RndSeed_main=ivector(1,RndSeed_main_num);
  for(ii=1;ii<=RndSeed_main_num;ii++){
	RndSeed_main[ii]=RndSeed_main_1st + ii-1;
  }

  /* VSmoke values */
  VSmoke_main_num=3;
  VSmoke_main=vector(1,VSmoke_main_num);
  VSmoke_main[1]=0.05;
  VSmoke_main[2]=0.25;
  VSmoke_main[3]=0.5;



  /* */ 
  Tleave=f3tensor(1,RndSeed_main_num,1,VSmoke_main_num,0,N0-1);
  for(ii=1;ii<=RndSeed_main_num;ii++){
      for(ij=1;ij<=VSmoke_main_num;ij++){
	  for(ik=0;ik<=N0-1;ik++){
	      Tleave[ii][ij][ik]=0.0;
	  }
      }
  }



  /* points of parameter space to be checked */
  for(RndSeed_main_index=1;
      RndSeed_main_index<=RndSeed_main_num;
      RndSeed_main_index++){

    for(VSmoke_main_index=1;
	VSmoke_main_index<=VSmoke_main_num;
	VSmoke_main_index++){


          VSmoke = VSmoke_main[VSmoke_main_index];
	  RndSeed = RndSeed_main[RndSeed_main_index];
          Init_Crunch();

	  /* loop */
          do { 
		  Save_Crunch_InLoop();
		  Upd(); 

	  } while(   ( N > NInjured )
		  && ( UpdNum < MaxUpdNum )
	 	  &&( SimTime[UpdNum] < MaxSimTime )
		 );
	  
	  Save_Crunch_AfterLoop();
	  Clean_Crunch_AfterLoop();
    }
  }

  Save_Crunch_AtEnd();
_E("END\n");

	return 0;
}
