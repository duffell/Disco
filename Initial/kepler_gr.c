
#include "../paul.h"

static double gam  = 0.0;
static double nu   = 0.0;
static double Mach = 0.0;
static double M = 1.0;

void setICparams( struct domain * theDomain ){
   gam  = theDomain->theParList.Adiabatic_Index;
   nu   = theDomain->theParList.viscosity;
   Mach = theDomain->theParList.Disk_Mach;
}

void initial(double *prim, double *x)
{
    double r = x[0];
    double z = x[2];
    double R = sqrt(r*r+z+z);

    double g00 = -1 + 2*M/r;
    double g0r = 2*M/r;
    double grr = 1 + 2*M/r;
    double gpp = r*r;

    double u0, ur, up;

    if(r > 6*M)
    {
        u0 = sqrt(r/(r-3*M));
        ur = 0.0;
        up = sqrt(M/(r*r*r - 3*M*r*r));
    }
    else
    {
        double x = 6*M/r - 1.0;
        u0 = 2.0 * (sqrt(2.0)*r - M*sqrt(x*x*x)) / (3.0*(r-2.0*M));
        ur = -sqrt(x*x*x) / 3.0;
        up = 2.0*sqrt(3.0)*M/(r*r);
    }

    double lr = g0r*u0 + grr*ur;
    double lp = gpp*up;

    double uIsco2 = 0.5;
    double u_cs2 = uIsco2 / (Mach*Mach);
    double cs2 = u_cs2 / (1.0 + u_cs2);

    double rho = 1.0;
    double Pp = rho / (1.0/cs2 - gam/(gam-1));

    prim[RHO] = rho;
    prim[PPP] = Pp;
    prim[URR] = lr;
    prim[UPP] = lp;
    prim[UZZ] = 0.0;
    
    if( NUM_N>0 ) 
        prim[NUM_C] = r>10.0 ? 1.0 : 0.0;
}
