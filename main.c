
#include "paul.h"

int mpiSetup( struct domain * , int , char *[] );
void setupGrid( struct domain * );
void timestep( struct domain * , double );
void setupCells( struct domain * );
void regrid( struct domain * );
void exchangeData( struct domain * , int );
double getmindt( struct domain * );
void calc_prim( struct domain * );

void read_par_file( struct domain * );
 
int  set_B_flag( void );
void set_B_fields( struct domain * );

void setupDomain( struct domain * );
void freeDomain( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );

void start_clock( struct domain * );
void generate_log( struct domain * );

int main( int argc , char * argv[] ){
 
   MPI_Init(&argc,&argv);
   struct domain theDomain = {0};
   start_clock( &theDomain ); 
   read_par_file( &theDomain );
   
   int error = mpiSetup(&theDomain,argc,argv);
   if( error==1 ) return(0);

   if(theDomain.rank==0) remove("abort");

   setupGrid( &theDomain );   
   setupDomain( &theDomain );
 
   setupCells( &theDomain );
/*
   if( theDomain.theParList.Initial_Regrid && !(theDomain.theParList.restart_flag) ) regrid( &theDomain );
*/
   if( theDomain.Nr > 1 ) exchangeData( &theDomain , 0 );
   if( theDomain.Nz > 1 ) exchangeData( &theDomain , 1 );

   int restart_flag = theDomain.theParList.restart_flag;
   if( set_B_flag() && NUM_FACES >= 3 && !restart_flag) 
       set_B_fields( &theDomain );

   if( theDomain.rank==0 && !(theDomain.theParList.restart_flag) ){
      FILE * rFile = fopen("report.dat","w");
      fclose(rFile);
   }

   while( !(theDomain.final_step) ){

      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );

   }

   possiblyOutput( &theDomain , 1 );
   generate_log( &theDomain );
   MPI_Barrier(theDomain.theComm);
   freeDomain( &theDomain );
   MPI_Finalize();

   return(0);

}

