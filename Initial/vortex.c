
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double rho = 1.0;
   double P0  = 1.0;
   double R   = 1.0;
   double Om  = 1.0;

   double r   = x[0];
   double phi = x[1];

   double X = 0.0;
   if( cos(phi) > 0.0 ) X = 1.0;

   double Pp,omega;
   omega = Om*exp(-.5*r*r/R/R);
   Pp = P0 - .5*rho*Om*Om*R*R*exp(-r*r/R/R);

   prim[RHO] = 1.0;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
