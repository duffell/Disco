enum{_HLL_,_HLLC_,_HLLD_};

#include "paul.h"

int mesh_motion = 0;
int riemann_solver = 0;
int visc_flag = 0;

void setRiemannParams( struct domain * theDomain ){
   mesh_motion = theDomain->theParList.Mesh_Motion;
   riemann_solver = theDomain->theParList.Riemann_Solver;
   visc_flag = theDomain->theParList.visc_flag;
}

void prim2cons( double * , double * , double , double );
void flux( double * , double * , double , double * );
void getUstar( double * , double * , double , double , double , double * , double * );
void get_Ustar_HLLD( double , double * , double * , double * , double * , double , double * );
void vel( double * , double * , double * , double * , double * , double * , double , double * );
double get_signed_dp( double , double );
void visc_flux( double * , double * , double * , double , double * );

void solve_riemann( double * , double * , double *, double * , double * , double * , double , double * , double , double , int );

void riemann_phi( struct cell * cL , struct cell * cR, double r , double dAdt ){

   double primL[NUM_Q];
   double primR[NUM_Q];

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5*cL->gradp[q]*cL->dphi;
      primR[q] = cR->prim[q] - .5*cR->gradp[q]*cR->dphi;
   }

   double n[3] = {0.0,1.0,0.0};

   solve_riemann( primL , primR , cL->cons , cR->cons , cL->gradp , cR->gradp , r , n , r*cL->wiph , dAdt , 0 );

}

void riemann_trans( struct face * F , double dt , int dim ){

   struct cell * cL = F->L;
   struct cell * cR = F->R;
   double dAdt      = F->dA*dt;
   double dxL       = F->dxL;
   double dxR       = F->dxR;
   double r         = F->cm[0];
   double phi       = F->cm[1];

   double primL[NUM_Q];
   double primR[NUM_Q];

   double phiL = cL->piph - .5*cL->dphi;
   double phiR = cR->piph - .5*cR->dphi;
   double dpL = get_signed_dp(phi,phiL);
   double dpR = get_signed_dp(phiR,phi);

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + cL->grad[q]*dxL + cL->gradp[q]*dpL;
      primR[q] = cR->prim[q] - cR->grad[q]*dxR - cR->gradp[q]*dpR;
   }
   
   double n[3] = {0.0,0.0,0.0};
   if( dim==1 ) n[0] = 1.0; else n[2] = 1.0;

   solve_riemann( primL , primR , cL->cons , cR->cons , cL->grad , cR->grad , r , n , 0.0 , dAdt , dim );

}


void solve_riemann( double * primL , double * primR , double * consL , double * consR , double * gradL , double * gradR , double r , double * n , double w , double dAdt , int dim ){

   int q;

   double Fl[NUM_Q];
   double Fr[NUM_Q];
   double Ul[NUM_Q];
   double Ur[NUM_Q];
   
   double Flux[NUM_Q];

   if( riemann_solver == _HLL_ || riemann_solver == _HLLC_ ){
      double Sl,Sr,Ss;
      double Bpack[5];
      vel( primL , primR , &Sl , &Sr , &Ss , n , r , Bpack );

      if( w < Sl ){
         flux( primL , Fl , r , n );
         prim2cons( primL , Ul , r , 1.0 );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fl[q] - w*Ul[q];
         }
      }else if( w > Sr ){
         flux( primR , Fr , r , n );
         prim2cons( primR , Ur , r , 1.0 );

         for( q=0 ; q<NUM_Q ; ++q ){
            Flux[q] = Fr[q] - w*Ur[q];
         }
      }else{
         if( riemann_solver == _HLL_ ){
            double Fstar;
            double Ustar;
            double aL =  Sr;
            double aR = -Sl;
 
            prim2cons( primL , Ul , r , 1.0 );
            prim2cons( primR , Ur , r , 1.0 );
            flux( primL , Fl , r , n );
            flux( primR , Fr , r , n );

            for( q=0 ; q<NUM_Q ; ++q ){
               Fstar = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
               Ustar = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );
  
               Flux[q] = Fstar - w*Ustar;
            }
         }else{
            double Ustar[NUM_Q];
            double Uk[NUM_Q];
            double Fk[NUM_Q];
            if( w < Ss ){
               prim2cons( primL , Uk , r , 1.0 );
               getUstar( primL , Ustar , r , Sl , Ss , n , Bpack ); 
               flux( primL , Fk , r , n ); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sl*( Ustar[q] - Uk[q] ) - w*Ustar[q];
               }    
            }else{
               prim2cons( primR , Uk , r , 1.0 );
               getUstar( primR , Ustar , r , Sr , Ss , n , Bpack ); 
               flux( primR , Fk , r , n ); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sr*( Ustar[q] - Uk[q] ) - w*Ustar[q];
               } 
            } 
         }
      }
   }else{
      double Ustar[NUM_Q];
      get_Ustar_HLLD( w , primL , primR , Flux , Ustar , r , n );
      for( q=0 ; q<NUM_Q ; ++q ) Flux[q] -= w*Ustar[q];
   }

   if( visc_flag ){
      double vFlux[NUM_Q];
      double prim[NUM_Q];
      double gprim[NUM_Q];
      for( q=0 ; q<NUM_Q ; ++q ){
         prim[q] = .5*(primL[q]+primR[q]);
         gprim[q] = .5*(gradL[q]+gradR[q]);
         if( dim==0 ) gprim[q] /= r;
         vFlux[q] = 0.0;
      }
      visc_flux( prim , gprim , vFlux , r , n );
      for( q=0 ; q<NUM_Q ; ++q ) Flux[q] += vFlux[q];
   }

   for( q=0 ; q<NUM_Q ; ++q ){
      consL[q] -= Flux[q]*dAdt;
      consR[q] += Flux[q]*dAdt;
   }

}




