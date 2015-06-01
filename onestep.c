
#include "paul.h"

void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clean_pi( struct domain * );
void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );
void adjust_RK_planets( struct domain * , double );
void move_cells( struct domain * , double );
void calc_dp( struct domain * );
void calc_prim( struct domain * );
void B_faces_to_cells( struct domain * , int , int );

void setup_faces( struct domain * , struct face ** , int * , int );
void phi_flux( struct domain * , double dt );
void trans_flux( struct domain * , struct face * , int , double dt , int );
void add_source( struct domain * , double dt );
void avg_Efields( struct domain * , int );
void update_B_fluxes( struct domain * , int * , int , double );
void subtract_advective_B_fluxes( struct domain * );

void movePlanets( struct planet * , double , double );
int planet_motion_analytic(void);

void boundary_r( struct domain * );
void boundary_trans( struct domain * , struct face * , int * , int );
void exchangeData( struct domain * , int );

int get_num_rzFaces( int , int , int );

void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step , double global_dt ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
 
   if( first_step ) set_wcell( theDomain );
   adjust_RK_cons( theDomain , RK );

   //struct face * theFaces_1 = NULL;
   //struct face * theFaces_2 = NULL;

   int NRZ1 = get_num_rzFaces( Nr , Nz , 1 );
   int NRZ2 = get_num_rzFaces( Nr , Nz , 2 );
   int * nfr = (int *) malloc( (NRZ1+1)*sizeof(int) );
   int * nfz = (int *) malloc( (NRZ2+1)*sizeof(int) );;

   phi_flux( theDomain , dt );

   setup_faces( theDomain , &(theDomain->theFaces_1) , nfr , 1 );
   trans_flux( theDomain , theDomain->theFaces_1 , nfr[NRZ1] , dt , 1 );

   if( Nz > 1 ){
      setup_faces( theDomain , &(theDomain->theFaces_2) , nfz , 2 );
      trans_flux( theDomain , theDomain->theFaces_2 , nfz[NRZ2] , dt , 2 );
   }

   avg_Efields( theDomain , nfr[NRZ1] );
   subtract_advective_B_fluxes( theDomain );
   update_B_fluxes( theDomain , nfr , NRZ1 , dt );

   add_source( theDomain , dt );
   if( first_step ){
      move_cells( theDomain , dt );
//      check_flipped( theDomain );
//      flip_fluxes( theDomain );
   }

   if( !planet_motion_analytic() || !last_step ){
      adjust_RK_planets( theDomain , RK );
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
   }
   clean_pi( theDomain );
   calc_dp( theDomain );
//   if( theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , nfr[NRZ1] , 1 );
//   }

   calc_prim( theDomain ); //ORDERING??? AFTER?

   if( last_step ){
      AMR( theDomain );
   }

   boundary_trans( theDomain , theDomain->theFaces_1 , nfr , 1 );
   exchangeData( theDomain , 0 );
   if( Nz > 1 ){
      int Periodic = theDomain->theParList.Z_Periodic;
      if( !Periodic ) boundary_trans( theDomain , theDomain->theFaces_2 , nfz , 2 );
      exchangeData( theDomain , 1 );
   }

   free( nfr );
   free( nfz );

   if( theDomain->theFaces_1 ) free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );

}

