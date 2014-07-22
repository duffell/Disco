
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){
   prim[RHO] = 1.0;
   prim[PPP] = 1e-6;
   prim[URR] = -1.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
