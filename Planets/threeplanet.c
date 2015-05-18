
#include "../paul.h"

static double q_planet = 1.0;
static double Mach = 1.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 4; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;

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

   double r1 = pow( 2. , -1./3. );
   double r2 = pow( 2. , 1./3.  );
   double r3 = .5;

   thePlanets[1].M     = q_planet; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = pow( r1 , -1.5 ); 
   thePlanets[1].r     = r1; 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = 0.5*r1/Mach;

   thePlanets[2].M     = q_planet; 
   thePlanets[2].vr    = 0.0; 
   thePlanets[2].omega = pow( r2 , -1.5 ); 
   thePlanets[2].r     = r2;
   thePlanets[2].phi   = 0.0; 
   thePlanets[2].eps   = 0.5*r2/Mach;

   thePlanets[3].M     = q_planet; 
   thePlanets[3].vr    = 0.0; 
   thePlanets[3].omega = pow( r3 , -1.5 ); 
   thePlanets[3].r     = r3;
   thePlanets[3].phi   = 0.0; 
   thePlanets[3].eps   = 0.5*r3/Mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[1].phi += thePlanets[1].omega*dt;
   thePlanets[2].phi += thePlanets[2].omega*dt;
   thePlanets[3].phi += thePlanets[3].omega*dt;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

