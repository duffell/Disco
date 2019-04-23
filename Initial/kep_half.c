
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

double get_cs2( double * );

void initial( double * prim , double * x ){

   double r = x[0];

   double xx = r*cos(x[1]);
   double yy = r*sin(x[1]);
   double sc_r2 = (xx-.5)*(xx-.5) + yy*yy;

   double k = 0.5;

   double rho = 1.0/pow(r,k) ;//+ 1.*exp(-20.*sc_r2);
   double cs2 = get_cs2( x );
   double Pp = rho*cs2/gam;
   double omega02 = 0.5/pow(r,3.) ;//+.5-.5/pow(1.-r,3.);
   double omegaP2 = omega02*( (1.+k)/Mach/Mach );

   double omega = sqrt( omega02 - omegaP2 );

   double X = 0.0; 
   if( r > 0.45 ) X = 1.0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = -1.5*nu/r;
   prim[UPP] = omega;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
