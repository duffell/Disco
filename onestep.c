
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
void B_faces_to_cells( struct domain * , int );

void setup_faces( struct domain * , int );
void phi_flux( struct domain * , double dt );
void trans_flux( struct domain * , double dt , int );
void add_source( struct domain * , double dt );

void avg_Efields( struct domain * );
void update_B_fluxes( struct domain * , double );
void subtract_advective_B_fluxes( struct domain * );
void check_flipped( struct domain * , int );
void flip_fluxes( struct domain * , int );

void movePlanets( struct planet * , double , double );
int planet_motion_analytic(void);

void boundary_r( struct domain * );
void boundary_trans( struct domain * , int );
void exchangeData( struct domain * , int );

//int get_num_rzFaces( int , int , int );
int set_B_flag( void );

void onestep( struct domain * theDomain , double RK , double dt , int first_step , int last_step , double global_dt ){

   int Nz = theDomain->Nz;
   int bflag = set_B_flag();
 
   if( first_step ) set_wcell( theDomain );
   adjust_RK_cons( theDomain , RK );

   phi_flux( theDomain , dt );

   setup_faces( theDomain , 1 );
   trans_flux( theDomain , dt , 1 );

   if( Nz > 1 ){
      setup_faces( theDomain , 2 );
      trans_flux( theDomain , dt , 2 );
   }

   if( bflag && NUM_EDGES >= 4 ){
      avg_Efields( theDomain );
      subtract_advective_B_fluxes( theDomain );
      update_B_fluxes( theDomain , dt );
   }

   add_source( theDomain , dt );

   if( first_step ){
      move_cells( theDomain , dt );
      if( bflag ){
         check_flipped( theDomain , 0 );
         flip_fluxes( theDomain , 0 );
         if( Nz>1 ){
            check_flipped( theDomain , 1 );
            flip_fluxes( theDomain , 1 );
         }
      }
   }

   if( !planet_motion_analytic() || first_step ){
      adjust_RK_planets( theDomain , RK );
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
   }
   clean_pi( theDomain );
   calc_dp( theDomain );

   if( bflag && theDomain->theParList.CT ){
      B_faces_to_cells( theDomain , 1 );
   }

   calc_prim( theDomain ); //ORDERING??? AFTER?

   if( last_step ){
      AMR( theDomain );
   }

   boundary_trans( theDomain , 1 );
   exchangeData( theDomain , 0 );
   if( Nz > 1 ){
      int Periodic = theDomain->theParList.Z_Periodic;
      if( !Periodic ) boundary_trans( theDomain , 2 );
      exchangeData( theDomain , 1 );
   }

   if( theDomain->theFaces_1 ) free( theDomain->theFaces_1 );
   if( theDomain->theFaces_2 ) free( theDomain->theFaces_2 );

}

