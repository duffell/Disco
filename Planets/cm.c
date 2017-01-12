
#include "../paul.h"

static double q_planet = 1.0;
static double smooth = 0.0;
static double Mach = 1.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   smooth = 0.5/Mach;

}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
 
   double q = q_planet;
   double mu = q/(1.+q);

   double om = pow( a , -1.5 );

   thePlanets[0].M     = 1.-mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = om; 
   thePlanets[0].r     = a*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].eps   = smooth;// + 1.0;//0.025; 

   thePlanets[1].M     = mu;  
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om;  
   thePlanets[1].r     = a*(1.-mu); 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = smooth;// + 1.0;//0.025;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){
   thePlanets[0].phi += thePlanets[0].omega*dt;
   thePlanets[1].phi += thePlanets[1].omega*dt;
   double eps = smooth + exp(-t/30.);
   thePlanets[0].eps = eps;
   thePlanets[1].eps = eps;
}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

