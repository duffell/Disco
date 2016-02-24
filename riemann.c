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

   if( !use_B_fields && riemann_solver == _HLLD_ ){
      printf("Ya dun goofed.\nRiemann Solver = HLLD,\nHydro does not include magnetic fields.\n");
      exit(1);
   }
   if( !use_B_fields && (NUM_C > 5 || NUM_EDGES > 0 || NUM_FACES > 0 || NUM_AZ_EDGES > 0 ) ){
      printf("Warning:  You are not solving MHD equations but possibly storing more variables than you need.\nNum Conserved Vars = %d\nFaces = %d\nEdges = %d\nAzimuthal Edges = %d\nCode might still work fine, this is just a warning.\n",NUM_C,NUM_FACES,NUM_EDGES,NUM_AZ_EDGES);
   }

}

void prim2cons( double * , double * , double * , double );
void flux( double * , double * , double * , double * );
void getUstar( double * , double * , double * , double , double , double * , double * );
void get_Ustar_HLLD( double , double * , double * , double * , double * , double , double * );
void vel( double * , double * , double * , double * , double * , double * , double * , double * );
double get_signed_dp( double , double );
void visc_flux( double * , double * , double * , double * , double * );

void solve_riemann( double * , double * , double *, double * , double * , double * , double * , double * , double , double , int , double * , double * , double * , double * );

void riemann_phi( struct cell * cL , struct cell * cR, double * x , double dAdt ){

   double r = x[0];
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

   double Er,Ez,Br,Bz;
   solve_riemann( primL , primR , cL->cons , cR->cons , cL->gradp , cR->gradp , x , n , r*cL->wiph , dAdt , 0 , &Ez , &Br , &Er , &Bz );

   if( NUM_EDGES == 4 ){
      cL->E[0] = .5*Ez;
      cL->E[1] = .5*Ez;

      cR->E[2] = .5*Ez;
      cR->E[3] = .5*Ez;

      cL->B[0] = .5*Br;
      cL->B[1] = .5*Br;

      cR->B[2] = .5*Br;
      cR->B[3] = .5*Br;
   }
   if( NUM_EDGES == 8 ){

      cL->E[0] = .5*Ez;
      cL->E[1] = .5*Ez;

      cR->E[2] = .5*Ez;
      cR->E[3] = .5*Ez;

      cL->B[0] = .5*Br;
      cL->B[1] = .5*Br;

      cR->B[2] = .5*Br;
      cR->B[3] = .5*Br;

      cL->E[4] = .5*Er;
      cL->E[5] = .5*Er;

      cR->E[6] = .5*Er;
      cR->E[7] = .5*Er;

      cL->B[4] = .5*Bz;
      cL->B[5] = .5*Bz;

      cR->B[6] = .5*Bz;
      cR->B[7] = .5*Bz;
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

   double Erz,Brz,Ephi,buffer;
   solve_riemann( primL , primR , cL->cons , cR->cons , cL->grad , cR->grad , F->cm , n , 0.0 , dAdt , dim , &Erz , &Brz , &Ephi , &buffer );
 
   double fracL = F->dphi / cL->dphi;
   double fracR = F->dphi / cR->dphi;

   if( NUM_EDGES >= 4 && dim==1 ){ 
      cL->E[1] += .5*Erz*fracL;
      cL->E[3] += .5*Erz*fracL;

      cR->E[0] += .5*Erz*fracR;
      cR->E[2] += .5*Erz*fracR;

      cL->B[1] += .5*Brz*fracL;
      cL->B[3] += .5*Brz*fracL;

      cR->B[0] += .5*Brz*fracR;
      cR->B[2] += .5*Brz*fracR;
   }
   if( NUM_AZ_EDGES == 4 && dim==1 ){
      cL->E_phi[1] = Ephi;
      cR->E_phi[0] = Ephi;
   }
   if( NUM_EDGES == 8 && dim==2){
      cL->E[5] += .5*Erz*fracL;
      cL->E[7] += .5*Erz*fracL;

      cR->E[4] += .5*Erz*fracR;
      cR->E[6] += .5*Erz*fracR;

      cL->B[5] += .5*Brz*fracL;
      cL->B[7] += .5*Brz*fracL;

      cR->B[4] += .5*Brz*fracR;
      cR->B[6] += .5*Brz*fracR;
   }
   if( NUM_AZ_EDGES == 4 && dim==2 ){
      cL->E_phi[3] = Ephi;
      cR->E_phi[2] = Ephi;
   }
}


void solve_riemann( double * primL , double * primR , double * consL , double * consR , double * gradL , double * gradR , double * x , double * n , double w , double dAdt , int dim , double * E1_riemann , double * B1_riemann , double * E2_riemann , double * B2_riemann ){

   int q;
   double r = x[0];

   double Flux[NUM_Q];
   double Ustr[NUM_Q];
   for( q=0 ; q<NUM_Q ; ++q ){
      Flux[q] = 0.0;
      Ustr[q] = 0.0;
   }

   if( riemann_solver == _HLL_ || riemann_solver == _HLLC_ ){
      double Sl,Sr,Ss;
      double Bpack[5];
      vel( primL , primR , &Sl , &Sr , &Ss , n , x , Bpack );

      if( w < Sl ){
         flux( primL , Flux , x , n );
         prim2cons( primL , Ustr , x , 1.0 );

      }else if( w > Sr ){
         flux( primR , Flux , x , n );
         prim2cons( primR , Ustr , x , 1.0 );

      }else{
         if( riemann_solver == _HLL_ ){
            double Fl[NUM_Q];
            double Fr[NUM_Q];
            double Ul[NUM_Q];
            double Ur[NUM_Q];
   
            double aL =  Sr;
            double aR = -Sl;
 
            prim2cons( primL , Ul , x , 1.0 );
            prim2cons( primR , Ur , x , 1.0 );
            flux( primL , Fl , x , n );
            flux( primR , Fr , x , n );

            for( q=0 ; q<NUM_Q ; ++q ){
               Flux[q] = ( aL*Fl[q] + aR*Fr[q] + aL*aR*( Ul[q] - Ur[q] ) )/( aL + aR );
               Ustr[q] = ( aR*Ul[q] + aL*Ur[q] + Fl[q] - Fr[q] )/( aL + aR );
            }
         }else{
            double Uk[NUM_Q];
            double Fk[NUM_Q];
            if( w < Ss ){
               prim2cons( primL , Uk , x , 1.0 );
               getUstar( primL , Ustr , x , Sl , Ss , n , Bpack ); 
               flux( primL , Fk , x , n ); 

               for( q=0 ; q<NUM_Q ; ++q ){
                  Flux[q] = Fk[q] + Sl*( Ustr[q] - Uk[q] );
               }    
            }else{
               prim2cons( primR , Uk , x , 1.0 );
               getUstar( primR , Ustr , x , Sr , Ss , n , Bpack ); 
               flux( primR , Fk , x , n ); 

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
      visc_flux( prim , gprim , vFlux , x , n );
      for( q=0 ; q<NUM_Q ; ++q ) Flux[q] += vFlux[q];
   }

   for( q=0 ; q<NUM_Q ; ++q ){
      consL[q] -= (Flux[q] - w*Ustr[q])*dAdt;
      consR[q] += (Flux[q] - w*Ustr[q])*dAdt;
   }

   if( use_B_fields && NUM_Q > BZZ ){
      if( dim==0 ){
         *E1_riemann = Flux[BRR]*r;   //Ez
         *B1_riemann = Ustr[BRR]*r*r; // r*Br
         *E2_riemann = Flux[BZZ];    //Er
         *B2_riemann = Ustr[BZZ]*r;  //-r*Bz
      }else if( dim==1 ){
         *E1_riemann = -Flux[BPP]*r;  //Ez
         *B1_riemann = Ustr[BRR]*r*r; // r*Br
         *E2_riemann = 1.0*Flux[BZZ];     //Ephi
      }else{
         *E1_riemann = -Flux[BPP]*r;   //Er
         *B1_riemann = Ustr[BZZ]*r;  //-r*Bz
         *E2_riemann = 1.0*-Flux[BRR]*r;  //Ephi
      }
   }

}




