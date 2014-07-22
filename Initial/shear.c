
#include "../paul.h"

static double nu = 0.0;
static double t = 0.0;

void setICparams( struct domain * theDomain ){
   nu = theDomain->theParList.viscosity;
   t  = theDomain->theParList.t_min;
}

void initial( double * prim , double * x ){

   double v0 = 0.001;

   double r   = x[0];
   double phi = x[1];
   double xx = r*cos(phi)-1.;

   double vy = v0*exp(-xx*xx/(4.*nu*t))/sqrt(4.*M_PI*nu*t);

   double vr = vy*sin(phi);
   double om = vy*cos(phi)/r;

   prim[RHO] = 1.0;
   prim[PPP] = 1.0;
   prim[URR] = vr;
   prim[UPP] = om;
   prim[UZZ] = 0.0;
   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}
