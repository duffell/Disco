
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

   double rho = 1.0;
   double Pp  = rho/Mach/Mach/gam;
   double n = 8.0;
   //double omega = ( pow( r , n-1.5 ) + 1. )/( pow( r , n ) + 1. );
   double omega = 1./pow( pow( r , 1.5*n ) + 1. , 1./n );

   double X = 0.0; 
   if( r > 1.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;//-1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
