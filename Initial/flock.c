//http://www.aanda.org/articles/aa/pdf/2010/08/aa12443-09.pdf
#include "../paul.h"

void setICparams( struct domain * theDomain ){
}

void initial( double * prim , double * x ){

   double rnd1 = (double)rand()/(double)(RAND_MAX);
   double rnd2 = (double)rand()/(double)(RAND_MAX);

   double r   = x[0];
   double phi = x[1];
   double z   = x[2];

   double n = 4.0;

   double rho = 1.0;
   double omega = pow( r , -1.5 );
   double cs = 0.1*omega*r;
   double Pp = cs*cs*rho/(5./3.);

   double vr = 1e-3*rnd1 - 5e-4;
   double vz = 1e-3*rnd2 - 5e-4;

   double Bz = 0.05513/n;
   if( r < 2. || r > 3. ) Bz = 0.0;

   prim[RHO] = rho;
   prim[PPP] = Pp;
   prim[URR] = vr;
   prim[UPP] = omega;
   prim[UZZ] = vz;

   prim[BRR] = 0.0;
   prim[BPP] = 0.0;
   prim[BZZ] = Bz;

   if( NUM_N>0 ) prim[NUM_C] = 0.0;

}

