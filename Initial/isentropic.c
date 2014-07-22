
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   double r   = x[0];
   double phi = x[1];
   double z   = x[2];
   double xx = r*cos(phi);
   double yy = r*sin(phi);//-.5;
   double R2 = xx*xx + yy*yy + z*z;
  // double xx = r*cos(phi);
  // double yy = r*sin(phi)-.5;
   prim[RHO] = 1.0 + 3.0*exp(-80.*R2);
   prim[PPP] = pow(prim[RHO],5./3.);
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
