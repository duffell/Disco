
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double s   = x[0];
   double phi = x[1];
   double z   = x[2];
 
   double r = sqrt(s*s+z*z);

   double r0 = 0.1; 
   double rho = 1.0;
   double Om0 = 1.5;
   double omega = Om0*exp(-.5*r*r/r0/r0);

   double Pp = (1.1)*rho*Om0*Om0*r0*r0/(2.*exp(1.));
   double Bx = 0.0;//0.01*sqrt(Pp);//2e-3;//0.001;

   double v = omega*s;
   double Bp = sqrt(rho)*v;
   double dP = -.5*Bp*Bp;

   double vr = 0.0;

   prim[RHO] = rho;
   prim[PPP] = Pp+dP;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;

   prim[BRR] =  Bx*cos(phi);
   prim[BPP] = -Bx*sin(phi) + Bp;
   prim[BZZ] = 0.0;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

