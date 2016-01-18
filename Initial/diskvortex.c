
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

void initial( double * prim , double * x_c ){

   double r = x_c[0];
   double phi = x_c[1];

   double rho = 1.0;

   double x = r*cos(phi);
   double y = r*sin(phi);
   double script_r = sqrt( (x-1)*(x-1) + y*y );
   double R = 0.1;
   double domega = .1*exp(-.5*script_r*script_r/R/R);
   double dP = -.5*rho*(domega*domega*R*R);

   double Pp = rho/Mach/Mach/gam;
   double omega02 = 1.0/pow(r,3.);
   double omegaP2 = 0.0;

   double omega = sqrt( omega02 - omegaP2 );
   double omega2 = domega*((x-1)*x + y*y)/script_r/r;
   double v2 = -domega*y/r;

   double X = 0.0; 
   if( r > 1.0 ) X = 1.0; 

   prim[RHO] = rho;
   prim[PPP] = Pp + dP;
   prim[URR] = -1.5*nu/r+v2;
   prim[UPP] = omega+omega2;;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = X;

}
