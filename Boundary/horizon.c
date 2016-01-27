
#include "../paul.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double * , double );
void prim2cons( double * , double * , double * , double );

void boundary_trans( struct domain * theDomain , int dim ){

   struct cell ** theCells = theDomain->theCells;

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Ng = theDomain->Ng;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   double M = 1.0; //TODO: get from parfile

   if( dim==1 && dim_rank[0] == dim_size[0]-1 ){
      int j;
      for( j=Nr-1 ; j>Nr-1-Ng ; --j ){
         int i,k;
         for( k=0 ; k<Nz ; ++k ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x );
            }
         }
      }
   }

   if( dim==2 && dim_rank[1] == 0 ){
      int i,j,k;
      for( k=0 ; k<Ng ; ++k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x ); 
            }    
         }    
      } 
   }
   if( dim==2 && dim_rank[1] == dim_size[1]-1 ){ 
      int i,j,k;
      for( k=Nz-1 ; k>Nz-1-Ng ; --k ){
         for( j=0 ; j<Nr ; ++j ){
            int jk = j+Nr*k;
            for( i=0 ; i<Np[jk] ; ++i ){
               struct cell * c = &(theCells[jk][i]);
               double phi = c->piph - .5*c->dphi;
               double x[3] = { .5*(r_jph[j]+r_jph[j-1]) , phi , .5*(z_kph[k]+z_kph[k-1]) };
               initial( c->prim , x );
            }
         }
      }
   }

   //Horizon jazz
   double rh = 2*M;
   double cut = 0.75*rh;

   //Check if the cut-out region is even on this process.
   if(r_jph[-1] < cut && z_kph[-1] < cut && z_kph[Nz-1] > -cut)
   {
       int i,j,k;
       for(k=0; k<Nz; k++)
       {
           double zp = z_kph[k];
           double zm = z_kph[k-1];
           double zc = fabs(zp) < fabs(zm) ? zp : zm;

           if(fabs(zc) > cut)
               continue;

           for(j=0; j<Nr; j++)
           {
               double rp = r_jph[j];
               double rm = r_jph[j-1];
               double rc = fabs(rp) < fabs(rm) ? rp : rm;

               if(rc*rc + zc*zc > cut)
                   break;

               int jk = j+Nr*k;
               for(i=0; i<Np[jk]; i++)
               {
                  struct cell * c = &(theCells[jk][i]);
                  double phi = c->piph - .5*c->dphi;
                  double x[3] = {0.5*(rm+rp), phi, 0.5*(zm+zp)};
                  initial(c->prim, x);
               }
           }
       }
   }
}

