
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r = x[0];
   double phi = x[1];

   double rho = 1.0;
   double Pp  = 1./r;

   double omega = 0.0;
   double v_rad = 0.0;

   double dsinp = sin(phi)+1-2./r;

   if( fabs(dsinp) < .2 && cos(phi)>0.0 ){
      rho += 1e4*exp(-.5*pow(dsinp/.1,2.) );
      omega = sqrt(2.)/r/r;
      v_rad = -.5*r*r*cos(phi)*omega;
   }
   

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = v_rad;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
