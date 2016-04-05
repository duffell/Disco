
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double Rin = 0.5;
   double Rout = 1.5; 

   double rho = 1.0;
   double R = 1.0;
   double P0 = 1.0;
   double P1 = 0.5;
   double k = 15.0;
   double Bz = 0.001;
   if( r < 0.75 || r > 1.25 ) Bz = 0.0;
   double Pp = P0+P1*tanh(k*(r-R)/R) - .5*Bz*Bz;
//   double Pp  = P0-P1*cos(M_PI*(r-Rin)/R);
//   double omega = sqrt( P1/rho*M_PI/R/r*sin(M_PI*(r-Rin)/R) );
   double om2r = P1/rho*k/R/pow( cosh(k*(r-R)/R) , 2. );
   double omega = sqrt( om2r/r );


   double nz = 1.0;
   double np = 6.0;
   double vr = 0.01*exp(-.5*pow(8.*(r-R),2.))*sin(2.*M_PI*nz*z)*sin(np*phi)*omega*r;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

