
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

void initial( double * prim , double * x ){

   double r = x[0];

   double rho = 1.0/pow(r,1.5);
   double Pp  = rho/Mach/Mach/gam;
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = 1.5/Mach/Mach/r/r;

   double omega = sqrt( omega02 - omegaP2 );

   double X = 0.0; 
   if( r > 1.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/Mach/Mach*sqrt(r);
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
