
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double v0 = 0.1;
   double P0 = 10.0;
   double R  = 1.0;
   double n  = 6.0;

   double X; 
   double omega;
   double Pp;
   double rho;

   if( r<R ){
      X = 0.0;
      omega = 2.0;
      Pp = P0 + 2.*r*r;
      rho = 1.0;
   }else{
      X = 1.0;
      omega = 1.0;
      Pp = P0 + r*r + R*R;
      rho = 2.0;
   }

   double vr = v0*cos(n*phi)*exp(-80.*pow(r-R,2.));

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
