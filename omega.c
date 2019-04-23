
#include "paul.h"

static double Mach = 0.0;
static int Npl = 0;
static struct planet * thePlanets = NULL;

double get_pot( struct planet * , int , double * );

void setDiskParams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
   Npl = theDomain->Npl;
   thePlanets = theDomain->thePlanets;
}

double mesh_om( double r ){
   double n = 4.0;
   double omega = 1./pow( pow( r , 1.5*n ) + 1. + 0.0*pow(.25/r,2.*n) , 1./n );
//   double omega = pow(r,-1.5);
//   double omega = 1.0;
   return( omega );
}

double get_om( double r ){
   return(0.0);
   return( 1.0/pow(r,1.5) );
//   return( 10.*exp(-.5*r*r) );
}
  
double get_om1( double r ){
   return(0.0);
   return( -1.5/pow(r,2.5) );
//   return( -10.*r*exp(-.5*r*r) );
}

double get_cs2( double * x ){
//   double r = x[0];
   double pot = get_pot( thePlanets , Npl , x );

//   double nu = .5;
//if( r<1. ) r=1.;
//   return( 1./Mach/Mach/pow(r,2.*nu) );
   return( pot / Mach / Mach );

//   return( 1./Mach/Mach );
//   return(1.0);
}

