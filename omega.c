
#include "paul.h"

static double Mach = 0.0;

void setDiskParams( struct domain * theDomain ){
   Mach = theDomain->theParList.Disk_Mach;
}

double get_om( double r ){
   return(0.0);
//   return( 1.0/pow(r,1.5) );
//   return(0.0);
//   return( 10.*exp(-.5*r*r) );
}
  
double get_om1( double r ){
   return(0.0);
//   return( -1.5/pow(r,2.5) );
//   return( -10.*r*exp(-.5*r*r) );
}

double get_cs2( double r ){
   return( 1./Mach/Mach );
}

