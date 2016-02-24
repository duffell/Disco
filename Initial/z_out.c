
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double r   = x[0];
   double phi = x[1];

   double gam = 5./3.;
   double s    = 0.01;
   double e    = 1.0;
   double mdot = 1.0;

   double u = 1.0;
   
   u = sqrt( 2.*e - 2.*gam/(gam-1.)*s*pow(mdot/u/r,gam-1.) );
   u = sqrt( 2.*e - 2.*gam/(gam-1.)*s*pow(mdot/u/r,gam-1.) );
   u = sqrt( 2.*e - 2.*gam/(gam-1.)*s*pow(mdot/u/r,gam-1.) );
   u = sqrt( 2.*e - 2.*gam/(gam-1.)*s*pow(mdot/u/r,gam-1.) );

   double vr  = u;
   double rho = mdot/(u*r);
   double Pp  = s*pow(rho,gam);

   double Rl = 0.3;
   double B0 = 0.01;

   double x0 = 0.4;


   double xl = r*cos(phi)-x0;
   double yl = r*sin(phi);//-x0;

   double rl = sqrt(xl*xl+yl*yl);
   double xx = .5*M_PI*rl/Rl;

   double Bz = B0*pow(cos(xx),4.);
   if( rl > Rl ) Bz = 0.0; 

   double dP = -.5*Bz*Bz;

   prim[RHO] = rho; 
   prim[PPP] = Pp+dP;
   prim[URR] = vr; 
   prim[UPP] = 0.0;
   prim[UZZ] = 0.0; 

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz; 

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
