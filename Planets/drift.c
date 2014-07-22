
#include "../paul.h"

static double rate = 0.0;
static double dexp = 0.0;
static double q = 0.0;
static double t_max = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2;

   rate  = theDomain->theParList.Drift_Rate;
   dexp  = theDomain->theParList.Drift_Exp;
   q     = theDomain->theParList.Mass_Ratio;
   t_max = theDomain->theParList.t_max;

}

int planet_motion_analytic( void ){
   return(1);
}

double drift_pos( double R , double t ){
   return( pow( 1. + R*( t - t_max )/dexp , dexp ) );
}

void initializePlanets( struct planet * thePlanets ){

   double r = drift_pos( rate , 0.0 );

   thePlanets[0].M     = 1.0; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = 0.0; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = q; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = pow(r,-1.5); 
   thePlanets[1].r     = r; 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = 0.025;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){

   double r = drift_pos( rate , t+dt );

   thePlanets[1].r     = r;
   thePlanets[1].omega = pow(r,-1.5);
   thePlanets[1].phi  += thePlanets[1].omega*dt;

}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

