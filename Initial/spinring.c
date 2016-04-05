
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double Rl = 0.45;
   double B0 = 1e-10;//1e-6;//0.0001;
   double Om = 1.0;
   double P0 = 1.0;

   double omega = Om;
   double Pp = P0 + .5*Om*Om*r*r;

   double xl = r*cos(phi);
   double yl = r*sin(phi);
   double zl = z;

   double rl = sqrt(xl*xl+zl*zl);
   double xx = M_PI*rl/Rl;

   double Bp = B0*pow(sin(xx),2.)*pow(cos(M_PI*yl/Rl),2.)*sqrt(2.*rl/Rl);
   if( rl > Rl || fabs(yl) > .5*Rl ) Bp = 0.0; 


   double Bx = -Bp*zl/rl;
   double Bz =  Bp*xl/rl;

   prim[RHO] = 1.0; 
   prim[PPP] = Pp;
   prim[URR] = 0.0; 
   prim[UPP] = omega;
   prim[UZZ] = 0.0; 

   prim[BRR] = Bx*cos(phi);
   prim[BPP] =-Bx*sin(phi);
   prim[BZZ] = Bz; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
