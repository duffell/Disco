
#include "paul.h"

int num_diagnostics(void){
   return(2);
}

void get_diagnostics( double * x , double * prim , double * Q ){
   double r = x[0];
   Q[0] = prim[RHO];
   Q[1] = 2.*M_PI*r*prim[RHO]*prim[URR];
}

void zero_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);

   int i,q;
   for( i=0 ; i<Nr ; ++i ){
      for( q=0 ; q<Nq ; ++q ){
         int iq = i*Nq + q;
         theTools->Qr[iq] = 0.0;
      }
   }
   theTools->t_avg = 0.0;

}

void avg_diagnostics( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   double dt = theTools->t_avg;

   int i,q; 
   for( i=0 ; i<Nr ; ++i ){
      for( q=0 ; q<Nq ; ++q ){
         int iq = i*Nq + q; 
         theTools->Qr[iq] /= dt; 
      }    
   }
   theTools->t_avg = 0.0; 

}

double get_dV( double * , double * );

void add_diagnostics( struct domain * theDomain , double dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Nq = theDomain->num_tools;
   struct diagnostic_avg * theTools = &(theDomain->theTools);
   struct cell ** theCells = theDomain->theCells;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;

   double temp_sum[Nr*Nq];
   memset( temp_sum , 0 , Nr*Nq*sizeof(double) );
   int i,j,k,q;
   int kmin = 0;
   int kmax = Nz;

   for( j=0 ; j<Nr ; ++j ){
      double dV_tot = 0.0;
      for( k=kmin ; k<kmax ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = theCells[jk]+i;
            double phip = c->piph;
            double phim = phip - c->dphi;
            double xp[3] = {r_jph[j  ] , phip , z_kph[k  ]};  
            double xm[3] = {r_jph[j-1] , phim , z_kph[k-1]};
            double xc[3];
            for( q=0 ; q<3 ; ++q ) xc[q] = .5*(xp[q]+xm[q]);
            double dV = get_dV(xp,xm);
            double Q[Nq];
            get_diagnostics( xc , c->prim , Q );
            for( q=0 ; q<Nq ; ++q ) temp_sum[ j*Nq + q ] += Q[q]*dV;
            dV_tot += dV;
         }
      }
      for( q=0 ; q<Nq ; ++q ){
         theTools->Qr[ j*Nq + q ] += temp_sum[ j*Nq + q ]*dt/dV_tot;
      }
   }

   theTools->t_avg += dt;

}

