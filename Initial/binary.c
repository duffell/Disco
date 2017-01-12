
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

double get_cs2( double );

void initial( double * prim , double * x ){

   double r = x[0];
   double z = x[2];

   double sint = z/sqrt(r*r+z*z);
   double cs2 = get_cs2( r );

   double rho = 1.0*exp(-sint*sint*Mach*Mach);
   double Pp  = rho*cs2/gam;
   //double n = 1.5;
   //double omega = ( pow( r , n-1.5 ) + 1. )/( pow( r , n ) + 1. );
   //double omega = 1./pow( pow( r , 1.5*n ) + 0.3 , 1./n );
   double omega2  = pow( r , -3. );
   double omega2P = 0.0*cs2/r/r;
   double omega = sqrt( omega2 - omega2P );

   if( omega > 1. ) omega = 1.;

   double X = 0.0; 
   if( r > 1.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;//-1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
