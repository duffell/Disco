
#include "../paul.h"

static double q_planet = 1.0; 
static double Mach = 1.0;
static double e_planet = 0.0;

void setPlanetParams( struct domain * theDomain ){

   theDomain->Npl = 2; 
   q_planet = theDomain->theParList.Mass_Ratio;
   Mach = theDomain->theParList.Disk_Mach;
   e_planet = theDomain->theParList.Eccentricity;

}

double root0( double E , double e ){ 
   return( E  - e*sin(E) );
}

double root1( double E , double e ){ 
   return( 1. - e*cos(E) );
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
   double e  = e_planet;
   double R = a*(1.-e);
   double om = pow( a , -1.5 )*sqrt(1.-e*e)/(1.-e)/(1.-e);
 
   double q = q_planet;
   double mu = q/(1.+q);

   thePlanets[0].M     = 1.0 - mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = R*mu; 
   thePlanets[0].phi   = M_PI; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = mu; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om; 
   thePlanets[1].r     = R*(1.0-mu); 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = 0.5/Mach;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){

   double TOL = 1e-8;

   double r0   = thePlanets[0].r + thePlanets[1].r; 
   double phi0 = thePlanets[1].phi;

   double vr = thePlanets[0].vr + thePlanets[1].vr; 
   double omega = thePlanets[1].omega; 

   double l = r0*r0*omega;
   double en = 0.5*vr*vr - 1./r0 + 0.5*l*l/r0/r0;

   double a = 1./2./fabs(en);
   double b = l/sqrt(2.*fabs(en));
   double f = sqrt(fabs(a*a-b*b));
   double e = f/a; 

   double x0 = r0*cos(phi0);
   double y0 = r0*sin(phi0);

   double E0 = atan2( y0/b , (x0+f)/a );
   double M0 = E0 - e*sin(E0);
   double M = M0 + l*dt/a/b;

//Newton-Rapheson to solve M = E - e*sin(E)
      double E = M;  //Guess value for E is M.
      double ff = root0( E , e ) - M; 
      while( fabs(ff) > TOL ){
         double dfdE = root1( E , e ); 
         double dE = -ff/dfdE;
         E += dE;
         ff = root0( E , e ) - M; 
      } 
      double x = a*cos(E)-f;
      double y = b*sin(E);
      double R   = sqrt(x*x+y*y);
      double phi = atan2(y,x);

   vr = sqrt( fabs( 2.*en + 2./R - l*l/R/R ) );
   if( y<0.0 ) vr *= -1.;

   double mu = q_planet/(1.+q_planet);

   thePlanets[1].r   = R*(1.-mu);
   thePlanets[1].phi = phi;
   thePlanets[1].omega = l/R/R; 
   thePlanets[1].vr = vr*(1.-mu);

   thePlanets[0].r   = R*mu;
   thePlanets[0].phi = phi+M_PI;
   thePlanets[0].omega = l/R/R;
   thePlanets[0].vr  = vr*mu;

}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

