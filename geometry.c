
#include "paul.h"

static double phi_max = 0.0;

void setGeometryParams( struct domain * theDomain ){
   phi_max = theDomain->phi_max;
}

double get_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp < 0.0 ) dp += phi_max;
   while( dp > phi_max) dp -= phi_max;
   return(dp);
}

double get_signed_dp( double phip , double phim ){
   double dp = phip-phim;
   while( dp <-.5*phi_max ) dp += phi_max;
   while( dp > .5*phi_max ) dp -= phi_max;
   return(dp);
}

double get_moment_arm( double * xp , double * xm ){
   double rp = xp[0];
   double rm = xm[0];
   double r2 = .5*(rp*rp+rm*rm);
   return( sqrt(r2) );
}

double get_dL( double * xp , double * xm , int dim ){
   double r   = .5*(xp[0]+xm[0]);
   double dphi = xp[1]-xm[1];
   while( dphi < 0.0 ) dphi += phi_max;
   while( dphi > phi_max ) dphi -= phi_max;
   if( dim==0 ) return( r*dphi );
   else if( dim==1 ) return( xp[0]-xm[0] );
   else return( xp[2]-xm[2] );
}

double get_dA( double * xp , double * xm , int dim ){
   double r  = .5*(xp[0]+xm[0]);
   double dr   = xp[0]-xm[0];
   double dphi = xp[1]-xm[1];
   while( dphi < 0.0 ) dphi += phi_max;
   while( dphi > phi_max ) dphi -= phi_max;
   double dz   = xp[2]-xm[2];
   if( dim==0 ) return( dr*dz );
   else if( dim==1 ) return( r*dphi*dz );
   else return( r*dr*dphi );
}

double get_dV( double * xp , double * xm ){
   double r  = .5*(xp[0]+xm[0]);
   double dr   = xp[0]-xm[0];
   double dphi = get_dp(xp[1],xm[1]);
   double dz   = xp[2]-xm[2];

   return( r*dr*dphi*dz );
}
