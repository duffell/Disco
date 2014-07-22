
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double r1 = 0.75;
   double r2 = 2.0;
   double rc = 1.25;

   double P0  = 0.9;
   double v0  = 0.5;

   double X; 
   double omega;
   double Pp;

   if( r<r1 ){
      X = -1.0;
      omega = pow(1./r,1.5);
      Pp = -1./r;
   }else if( r<r2 ){
      X = 0.0;
      omega = r1*r1*pow(1./r,3.5);
      Pp = -(1./5.)*pow(r1,4.)/pow(r,5.) - (4./5.)/r1;
   }else{
      X = 1.0;
      omega = r1*r1/r2/r2*pow(1./r,1.5);
      Pp = -pow(r1/r2,4.)*(1./r-1./r2) - (1./5.)*pow(r1,4.)/pow(r2,5.) - (4./5.)/r1;
   }

   Pp += P0 + 1./r;

   double vr = v0*cos(2.*phi)*exp(-10.*pow(r-rc,2.));

   //if( cos(phi) > 0.0 ) X = 1.0; 

   prim[RHO] = 1.0;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
