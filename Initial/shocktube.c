
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r = x[0];

   double rho = 1.0;
   double Pp  = 1.0;

   if( r>0.5 ){
      rho = 0.125;
      Pp  = 0.1;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
