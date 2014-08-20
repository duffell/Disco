
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double r0 = 0.1; 
   double r1 = 0.115;
   double Bx = 5.0/sqrt(4.*M_PI);
   double v0 = 2.0; 
   double Pp = 1.0; 

   double f = 0.0; 
   if( r<r0 ) f = 1.0; else if( r<r1 ) f = (r1-r)/(r1-r0);
   double rho  = 1.0 + 9.*f;

   double omega = v0/r0*f;
   if( r>r0 ) omega = v0/r*f;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;

   prim[BRR] =  Bx*cos(phi);
   prim[BPP] = -Bx*sin(phi);
   prim[BZZ] = 0.0;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
