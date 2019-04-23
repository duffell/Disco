
#include "../paul.h"

static double gam   = 0.0; 
static double alpha = 0.0; 
static double Mach  = 0.0; 
static int alpha_flag = 0;

void setICparams( struct domain * theDomain ){
   gam   = theDomain->theParList.Adiabatic_Index;
   alpha = theDomain->theParList.viscosity;
   Mach  = theDomain->theParList.Disk_Mach;
   alpha_flag = theDomain->theParList.alpha_flag;
}

double get_cs2( double * );

void initial( double * prim , double * x ){ 

   double r = x[0];
   double cs2 = get_cs2( x ); 
 
   double kappa = 0.5; 
   kappa = 0.0;
   double nu_exp = 0.5; 
   double pexp = kappa + 2.*nu_exp;

   double rho = 1.0*pow( r , -kappa );
   if( r < 1. ) rho = 0.3;
   double Pp = rho*cs2;
   double omega02 = 1.0/pow(r,3.);
   if( r < 1. ) omega02 = 1.0;
   double omegaP2 = pexp*cs2/r/r;
   if( r < 1. ) omegaP2 = pexp*cs2;

   double omega = sqrt( omega02 - omegaP2 );
   double nu = alpha*cs2*pow(r,1.5);
   if( !alpha_flag ) nu = alpha;
   if( r < 1. ) nu = 0.0;

   double X = 0.0; 
//   if( r*cos(x[1]) > 0.0 ) X = 1.0; 
   if( r>1. ) X = 1.0;

//   if( r < .1 ) omega *= pow(r/.1,3.);

   prim[RHO] = rho; 
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0; 
   if( NUM_N>0 ) prim[NUM_C] = X; 

}
