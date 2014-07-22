
#include "../paul.h"

int numPlanets( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){
   thePlanets[0].M     = 1.0; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = 0.0; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = 0.0;
}

void movePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

