
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double Bx   = 1.0;
   double By   = 0.0;
   double P0   = 1.0;
   double rho0 = 1.0;

   prim[RHO] = rho0; 
   prim[PPP] = P0;
   prim[URR] = 0.0; 
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0; 

   prim[BRR] =  Bx*cos(phi) + By*sin(phi);
   prim[BPP] = -Bx*sin(phi) + By*cos(phi);
   prim[BZZ] = 0.0; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
