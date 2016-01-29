
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double Rl = 0.3;
   double B0 = 0.01;
   double Om = 1.0;
   double P0 = 1.0;

   double x0 = 0.25;


   double omega = Om;
   double Pp = P0 + .5*Om*Om*r*r;

   double xl = r*cos(phi)-x0;
   double yl = r*sin(phi);//-x0;

   double rl = sqrt(xl*xl+yl*yl);
   double xx = .5*M_PI*rl/Rl;

   double Bz = B0*pow(cos(xx),4.);
   if( rl > Rl ) Bz = 0.0; 

   double dP = -.5*Bz*Bz;

   prim[RHO] = 1.0; 
   prim[PPP] = Pp+dP;
   prim[URR] = 0.0; 
   prim[UPP] = omega;
   prim[UZZ] = 0.0; 

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
