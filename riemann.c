enum{_HLL_,_HLLC_,_HLLD_};

#include "paul.h"

static int mesh_motion = 0;
static int riemann_solver = 0;
static int visc_flag = 0;
static int use_B_fields = 1;

int set_B_flag( void );

void setRiemannParams( struct domain * theDomain ){
   mesh_motion = theDomain->theParList.Mesh_Motion;
   riemann_solver = theDomain->theParList.Riemann_Solver;
   visc_flag = theDomain->theParList.visc_flag;
   use_B_fields = set_B_flag();
}

void prim2cons( double * , double * , double , double );
void flux( double * , double * , double , double * );
void getUstar( double * , double * , double , double , double , double * , double * );
void get_Ustar_HLLD( double , double * , double * , double * , double * , double , double * );
void vel( double * , double * , double * , double * , double * , double * , double , double * );
double get_signed_dp( double , double );
void visc_flux( double * , double * , double * , double , double * );

void solve_riemann( double * , double * , double *, double * , double * , double * , double , double * , double , double , int , double * , double * );

void riemann_phi( struct cell * cL , struct cell * cR, double r , double dAdt ){

   double primL[NUM_Q];
   double primR[NUM_Q];

   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      primL[q] = cL->prim[q] + .5*cL->gradp[q]*cL->dphi;
      primR[q] = cR->prim[q] - .5*cR->gradp[q]*cR->dphi;
   }

   double n[3] = {0.0,1.0,0.0};

   if( use_B_fields && NUM_Q > BPP ){
      double Bp = .5*(primL[BPP]+primR[BPP]);
      primL[BPP] = Bp;
      primR[BPP] = Bp;
   }

   double E,B;
   solve_riemann( primL , primR , cL->cons , cR->cons , cL->gradp , cR->gradp , r , n , r*cL->wiph , dAdt , 0 , &E , &B );

   if( NUM_EDGES == 4 ){
      cL->E[0] = .5*E;
      cL->E[1] = .5*E;

      cR->E[2] = .5*E;
      cR->E[3] = .5*E;

      cL->B[0] = .5*B;
      cL->B[1] = .5*B;

      cR->B[2] = .5*B;
      cR->B[3] = .5*B;
   }
   if( NUM_EDGES == 12 ){

   }
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

   if( use_B_fields ){
      int BTRANS;
      if( dim==1 ) BTRANS = BRR; else BTRANS = BZZ;
      double Bavg = .5*(primL[BTRANS]+primR[BTRANS]);
      primL[BTRANS] = Bavg;
      primR[BTRANS] = Bavg;
   }

   double E,B;
   solve_riemann( primL , primR , cL->cons , cR->cons , cL->grad , cR->grad , r , n , 0.0 , dAdt , dim , &E , &B );

//  GET ACTUAL AREAS HERE ///////
   double dA  = F->dA;
   double dAL = r*cL->dphi;
   double dAR = r*cR->dphi;

   if( NUM_EDGES == 4 ){ 
      cL->E[1] += .5*E*dA/dAL;
      cL->E[3] += .5*E*dA/dAL;

      cR->E[0] += .5*E*dA/dAR;
      cR->E[2] += .5*E*dA/dAR;

      cL->B[1] += .5*B*dA/dAL;
      cL->B[3] += .5*B*dA/dAL;

      cR->B[0] += .5*B*dA/dAR;
      cR->B[2] += .5*B*dA/dAR;
   }

}


void solve_riemann( double * primL , double * primR , double * consL , double * consR , double * gradL , double * gradR , double r , double * n , double w , double dAdt , int dim , double * E_riemann , double * B_riemann ){

   int q;

   double Flux[NUM_Q];
   double Ustr[NUM_Q];
   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = 0.0;
      Ustr[q] = 0.0;
   }

   if( riemann_solver == _HLL_ || riemann_solver == _HLLC_ ){
      double Sl,Sr,Ss;
      double Bpack[5];
      vel( primL , primR , &Sl , &Sr , &Ss , n , r , Bpack );

      if( w < Sl ){
         flux( primL , Flux , r , n );
         prim2cons( primL , Ustr , r , 1.0 );

      }else if( w > Sr ){
         flux( primR , Flux , r , n );
         prim2cons( primR , Ustr , r , 1.0 );

      }else{
         if( riemann_solver == _HLL_ ){
            double Fl[NUM_Q];
            double Fr[NUM_Q];
            double Ul[NUM_Q];
            double Ur[NUM_Q];
   
            double aL =  Sr;
            double aR = -Sl;
 
            prim2cons( primL , Ul , r , 1.0 );
            prim2cons( primR , Ur , r , 1.0 );
            flux( primL , Fl , r , n );
            flux( primR , Fr , r , n );

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
               Ustr[q] = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );
            }
         }else{
            double Uk[NUM_Q];
            double Fk[NUM_Q];
            if( w < Ss ){
               prim2cons( primL , Uk , r , 1.0 );
               getUstar( primL , Ustr , r , Sl , Ss , n , Bpack ); 
               flux( primL , Fk , r , n ); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sl*( Ustr[q] - Uk[q] );
               }    
            }else{
               prim2cons( primR , Uk , r , 1.0 );
               getUstar( primR , Ustr , r , Sr , Ss , n , Bpack ); 
               flux( primR , Fk , r , n ); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sr*( Ustr[q] - Uk[q] );
               } 
            } 
         }
      }
   }else{
      get_Ustar_HLLD( w , primL , primR , Flux , Ustr , r , n );
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
      consL[q] -= (Flux[q] - w*Ustr[q])*dAdt;
      consR[q] += (Flux[q] - w*Ustr[q])*dAdt;
   }

   if( use_B_fields && NUM_Q > BPP ){
      if( dim==0 ){
         *E_riemann = Flux[BRR]*r;
      }else if( dim==1 ){
         *E_riemann = -Flux[BPP]*r;
      }

      *B_riemann = Ustr[BRR]*r*r;
   }

}




