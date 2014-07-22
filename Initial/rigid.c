
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r = x[0];
   double phi = x[1];

   double rho = 1.0;
   double omega = 1.0;

   double xc = r*cos(phi) - 0.75;
   double yc = r*sin(phi);
   double r2 = xc*xc + yc*yc;

   double Pp  = 0.01 + .5*rho*omega*omega*r*r;
   if( r2 < 0.01 ) Pp += 1.0;

   double X = 0.0; 
   if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;
}
