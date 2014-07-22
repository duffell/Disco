
#include "../paul.h"

double root0( double E , double e ){ 
   return( E  - e*sin(E) );
}

double root1( double E , double e ){ 
   return( 1. - e*cos(E) );
}

int numPlanets( void ){
   return(2);
}

int planet_motion_analytic( void ){
   return(1);
}

void initializePlanets( struct planet * thePlanets ){

   double a  = 1.0;
   double e  = 0.1;
   double r  = a*(1.-e);
   double om = pow( a , -1.5 )*sqrt(1.-e*e)/(1.-e)/(1.-e);
 
   double q = 3e-6;//1.25e-4;//3e-6;
   double mu = q/(1.+q);

   thePlanets[0].M     = 1.0 - mu; 
   thePlanets[0].vr    = 0.0; 
   thePlanets[0].omega = 0.0; 
   thePlanets[0].r     = 0.0; 
   thePlanets[0].phi   = 0.0; 
   thePlanets[0].eps   = 0.0;

   thePlanets[1].M     = mu; 
   thePlanets[1].vr    = 0.0; 
   thePlanets[1].omega = om; 
   thePlanets[1].r     = r; 
   thePlanets[1].phi   = 0.0; 
   thePlanets[1].eps   = 0.025;

}

void movePlanets( struct planet * thePlanets , double t , double dt ){

   double TOL = 1e-8;

   double r0   = thePlanets[1].r; 
   double phi0 = thePlanets[1].phi;

   double vr = thePlanets[1].vr; 
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
      double r   = sqrt(x*x+y*y);
      double phi = atan2(y,x);

   vr = sqrt( fabs( 2.*en + 2./r - l*l/r/r ) );
   if( y<0.0 ) vr *= -1.;

   thePlanets[1].r   = r;
   thePlanets[1].phi = phi;
   thePlanets[1].omega = l/r/r; 
   thePlanets[1].vr = vr;

}

void forcePlanets( struct planet * thePlanets , double dt ){
   //Silence is golden.
}

