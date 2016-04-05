
#include "../paul.h"
 
static int ranky=0;

void setICparams( struct domain * theDomain ){
   ranky = theDomain->rank;
}

void initial( double * prim , double * x ){
   prim[RHO] = 1.0+(double)ranky;
   prim[PPP] = 1.0;
   prim[URR] = 0.0;
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;
}
