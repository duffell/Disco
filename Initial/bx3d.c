
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];
   double z = x[2];
   r = sqrt(r*r+z*z);

   double rho  = 1.0;
   double Pp   = 0.1;

   double b = 2.0;

   if( r<0.1 ){
      Pp  = 10.0;
   }

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = b;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
