
#include "../paul.h"

int numPlanets( void ){
   return(2);
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   thePlanets[0].M     = 1.0; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = 0.0; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = 1e-3; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = 1.0; 
   thePlanets[1].r     = 1.0; 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = 0.025;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[1].phi += thePlanets[1].omega*dt;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

