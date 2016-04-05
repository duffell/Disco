
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double Rl = 0.15;
   double B0 = 1e-4;
   double Om = 1.0;
   double P0 = 0.01;

   double x0 = 0.25;

   double omega = Om;
   double Pp = P0 + .5*Om*Om*r*r;

   double xl = r*cos(phi)-x0;
   double yl = r*sin(phi);

   double rl = sqrt(xl*xl+yl*yl);
   double xx = M_PI*rl/Rl;

   double Bp = B0*pow(sin(xx),2.)*sqrt(2.*rl/Rl);
   if( rl > Rl ) Bp = 0.0; 

   double dP = B0*B0*( - (rl/Rl)*pow(sin(xx),4.) - (1./16./M_PI)*( 12.*xx - 8.*sin(2.*xx) + sin(4.*xx) ) ); 
   if( rl > Rl ) dP = 0.0;

   double Bx = -Bp*yl/rl;
   double By =  Bp*xl/rl;

   prim[RHO] = 1.0; 
   prim[PPP] = Pp+dP;
   prim[URR] = 0.0; 
   prim[UPP] = omega;
   prim[UZZ] = 0.0; 

   prim[BRR] =  Bx*cos(phi) + By*sin(phi);
   prim[BPP] = -Bx*sin(phi) + By*cos(phi);
   prim[BZZ] = 0.0; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
