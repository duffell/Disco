//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../paul.h"

static double Mach = 0.0;

void setICparams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
}

double I(double r){
   double v=0.0;
   if( r>1.2 && r<3.8 ) v = 1.0;
   return(v);
}

double floor( double x ){
   return( (double)(int)x );
}

double Rs( double x ){
   double R0 = 1.0;
   double H0 = 0.2;
   return( 1.2 + (floor((x-R0)/H0)-.5)*H0 );
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double rho = 100.0;
   double omega = pow( r , -1.5 );
   double Pp = rho/Mach/Mach;

   double dvp = 2e-2*rnd1 - 1e-2;
   double vz = 2e-2*rnd2 - 1e-2;

   double A0 = 0.37;
   double H0 = 0.2;
   double R0 = 1.0;
   double Bz = A0*sin(2.*M_PI*(r-R0)/H0)/r*I(r)/sqrt(Rs(r));

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = 0.0;
   prim[UPP] = omega*(1.+dvp);
   prim[UZZ] = vz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

