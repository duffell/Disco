
#include "../paul.h"

static double gam   = 0.0; 
static double alpha = 0.0; 
static double Mach  = 0.0; 

void setICparams( struct domain * theDomain ){
   gam   = theDomain->theParList.Adiabatic_Index;
   alpha = theDomain->theParList.viscosity;
   Mach  = theDomain->theParList.Disk_Mach;
}

double get_cs2( double * );

void initial( double * prim , double * x ){ 

   double r = x[0];
   double cs2 = get_cs2( x ); 
 
   double kappa = 0.5; 
   double nu_exp = 0.5; 
   double pexp = kappa + 2.*nu_exp;

   double rho = 1.0*pow( r , -kappa );
   double Pp = rho*cs2;
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = pexp*cs2/r/r;

   double omega = sqrt( omega02 - omegaP2 );
   double nu = alpha*cs2/omega;

   double X = 0.0; 
   if( r*cos(x[1]) > 0.0 ) X = 1.0; 

   prim[RHO] = rho; 
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0; 
   if( NUM_N>0 ) prim[NUM_C] = X; 

}
