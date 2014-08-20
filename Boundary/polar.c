
#include "../paul.h"
#include <string.h>

void initial( double * , double * );
double get_dV( double * , double * );
void cons2prim( double * , double * , double , double );
void prim2cons( double * , double * , double , double );

void boundary_trans( struct domain * theDomain , struct face * theFaces , int * nn , int dim ){

   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Ng = theDomain->Ng;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

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

   if( NUM_Q >= 8 ){
      if( dim==1 && dim_rank[0] == 0 ){
         int j;
         for( j=0 ; j<Ng ; ++j ){
            int i,k;
            for( k=0 ; k<Nz ; ++k ){
               int jk = j+Nr*k;
               for( i=0 ; i<Np[jk] ; ++i ){
                  struct cell * c = &(theCells[jk][i]);
                  c->prim[BRR] = 0.0;
                  c->prim[BPP] = 0.0;
                  c->cons[BRR] = 0.0;
                  c->cons[BPP] = 0.0;
                  c->RKcons[BRR] = 0.0;
                  c->RKcons[BPP] = 0.0;
               }
            }
         }
      }
   }
/*
   if( dim==1 && (dim_rank[0]==0 || dim_rank[0]==dim_size[0]-1) ){
      int k,LR;
      for( LR=0 ; LR<2 ; ++LR ){
         if( ( dim_rank[0]==0 && LR==0 ) || ( dim_rank[0]==dim_size[0]-1 && LR==1 ) ){
            int j0 = 0;
            int j1 = 1;
            if( LR==1 ){
               j0 = Nt-1;
               j1 = Nt-2;
            }
            for( k=0 ; k<Np ; ++k ){
               int jk0 = j0+Nt*k;
               int jk1 = j1+Nt*k;
               Nr[jk0] = Nr[jk1];
               theCells[jk0] = realloc( theCells[jk0] , Nr[jk0]*sizeof(struct cell) );
               memcpy( &(theCells[jk0][0]) , &(theCells[jk1][0]) , Nr[jk0]*sizeof(struct cell) );
               int i;
               for( i=0 ; i<Nr[jk0] ; ++i ){
                  struct cell * c = &( theCells[jk0][i] );
                  c->prim[UU2] *= .5;
                  double rp = c->riph;
                  double rm = 0.0;
                  if( i>0 ) rm = theCells[jk0][i-1].riph;
                  //double r = (2./3.)*(rp*rp*rp-rm*rm*rm)/(rp*rp-rm*rm);
                  double xp[3] = {rp,t_jph[j0]  ,p_kph[k]  };
                  double xm[3] = {rm,t_jph[j0-1],p_kph[k-1]};
                  double dV = get_dV( xp , xm );
                  xp[1] = t_jph[j1  ];
                  xm[1] = t_jph[j1-1];
                  double dV2 = get_dV( xp , xm );
                  //prim2cons( c->prim , c->cons , r , dV );
                  int q;
                  for( q=0 ; q<NUM_Q ; ++q ){
                     c->cons[q]   *= dV/dV2;
                     c->RKcons[q] *= dV/dV2;
                  }
               }
            }  
         } 
      }
   }
*/
}

