
#include "paul.h"

double get_moment_arm( double * , double * );
double get_dV( double * , double * );

int num_diagnostics( void );
void initializePlanets( struct planet * );

void setICparams( struct domain * );
void setHydroParams( struct domain * );
void setGeometryParams( struct domain * );
void setRiemannParams( struct domain * );
void setPlanetParams( struct domain * );

void setupDomain( struct domain * theDomain ){

   srand(theDomain->rank);
   rand();

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   theDomain->theCells = (struct cell **) malloc( Nr*Nz*sizeof(struct cell *) );
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      theDomain->theCells[jk] = (struct cell *) malloc( Np[jk]*sizeof(struct cell) );
   }

   setPlanetParams( theDomain );
   int Npl = theDomain->Npl;
   theDomain->thePlanets = (struct planet *) malloc( Npl*sizeof(struct planet) );
   initializePlanets( theDomain->thePlanets );

   double num_tools = num_diagnostics();
   theDomain->num_tools = num_tools;
   theDomain->theTools.t_avg = 0.0;
   theDomain->theTools.Qr = (double *) calloc( Nr*num_tools , sizeof(double) );

   int i;
   double Pmax = theDomain->theParList.phimax;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      double p0 = Pmax*(double)rand()/(double)RAND_MAX;
      double dp = Pmax/(double)Np[jk];
      for( i=0 ; i<Np[jk] ; ++i ){
         double phi = p0+dp*(double)i;
         if( phi > Pmax ) phi -= Pmax;
         theDomain->theCells[jk][i].piph = phi;
         theDomain->theCells[jk][i].dphi = dp;
      }
   }

   theDomain->t       = theDomain->theParList.t_min;
   theDomain->t_init  = theDomain->theParList.t_min;
   theDomain->t_fin   = theDomain->theParList.t_max;
   theDomain->phi_max = theDomain->theParList.phimax;

   theDomain->N_rpt = theDomain->theParList.NumRepts;
   theDomain->N_snp = theDomain->theParList.NumSnaps;
   theDomain->N_chk = theDomain->theParList.NumChecks;

   theDomain->count_steps = 0;
   theDomain->final_step = 0;

   theDomain->nrpt=-1;
   theDomain->nsnp=-1;
   theDomain->nchk=-1;

   setICparams( theDomain );
   setHydroParams( theDomain );
   setGeometryParams( theDomain );
   setRiemannParams( theDomain );

}

void initial( double * , double * ); 
void prim2cons( double * , double * , double , double );
void cons2prim( double * , double * , double , double );
void restart( struct domain * ); 
void calc_dp( struct domain * );
void set_wcell( struct domain * );
void adjust_gas( struct planet * , double * , double * , double );
void set_B_fields( struct domain * );

void setupCells( struct domain * theDomain ){

   int restart_flag = theDomain->theParList.restart_flag;
   if( restart_flag ) restart( theDomain );

   calc_dp( theDomain );

   int i,j,k;
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;
   double * r_jph = theDomain->r_jph;
   double * z_kph = theDomain->z_kph;
   int atmos = theDomain->theParList.include_atmos;

   for( j=0 ; j<Nr ; ++j ){
      for( k=0 ; k<Nz ; ++k ){
         int jk = j+Nr*k;
         for( i=0 ; i<Np[jk] ; ++i ){
            struct cell * c = &(theCells[jk][i]);
            double phip = c->piph;
            double phim = phip-c->dphi;
            c->wiph = 0.0; 
            double xp[3] = {r_jph[j  ],phip,z_kph[k  ]};
            double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
            double r = get_moment_arm( xp , xm );
            double dV = get_dV( xp , xm );
            double phi = c->piph-.5*c->dphi;
            double x[3] = { r , phi , .5*(z_kph[k]+z_kph[k-1])};
            if( !restart_flag ) initial( c->prim , x ); 
            if( atmos ){
               int p;
               for( p=0 ; p<Npl ; ++p ){
                  double gam = theDomain->theParList.Adiabatic_Index;
                  adjust_gas( theDomain->thePlanets+p , x , c->prim , gam );
               }
            }
            prim2cons( c->prim , c->cons , r , dV );
            cons2prim( c->cons , c->prim , r , dV );
         }    
      }    
   }

   set_wcell( theDomain );

}


/*
void clear_cell( struct cell * c ){
   int q;
   for( q=0 ; q<NUM_Q ; ++q ){
      c->prim[q]   = 0.0;
      c->cons[q]   = 0.0;
      c->RKcons[q] = 0.0;
      c->grad[q]   = 0.0;
      c->gradr[q]  = 0.0;
   }
   c->riph = 0.0;
   c->RKriph = 0.0;
   c->dr = 0.0;
   c->wiph = 0.0;
}
*/

void freeDomain( struct domain * theDomain ){

   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int jk;
   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      free( theDomain->theCells[jk] );
   }
   free( theDomain->theCells );
   free( theDomain->Np );
   theDomain->r_jph--;
   free( theDomain->r_jph );
   theDomain->z_kph--;
   free( theDomain->z_kph );
   free( theDomain->thePlanets );
   free( theDomain->theTools.Qr );

}

void check_dt( struct domain * theDomain , double * dt ){

   double t = theDomain->t;
   double tmax = theDomain->t_fin;
   int final=0;
   if( t + *dt > tmax ){
      *dt = tmax-t;
      final=1;
   }

   if( theDomain->rank==0 ){
      FILE * abort = NULL;
      abort = fopen("abort","r");
      if( abort ){ final = 1; fclose(abort); }
   }

   MPI_Allreduce( MPI_IN_PLACE , &final , 1 , MPI_INT , MPI_SUM , theDomain->theComm );
   if( final ) theDomain->final_step = 1;

}

void report( struct domain * );
void snapshot( struct domain * , char * );
void output( struct domain * , char * );

void possiblyOutput( struct domain * theDomain , int override ){

   double t = theDomain->t;
   double t_min = theDomain->t_init;
   double t_fin = theDomain->t_fin;
   double Nrpt = theDomain->N_rpt;
   double Nsnp = theDomain->N_snp;
   double Nchk = theDomain->N_chk;
   int LogOut = theDomain->theParList.Out_LogTime;
   int n0;

   n0 = (int)( t*Nrpt/t_fin );
   if( LogOut ) n0 = (int)( Nrpt*log(t/t_min)/log(t_fin/t_min) );
   if( theDomain->nrpt < n0 || override ){
      theDomain->nrpt = n0;
      //longandshort( &theDomain , &L , &S , &iL , &iS , theDomain.theCells[0] , 0 , 0 );
      report( theDomain );
      if( theDomain->rank==0 ) printf("t = %.3e\n",t);
   }
   n0 = (int)( t*Nchk/t_fin );
   if( LogOut ) n0 = (int)( Nchk*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nchk < n0 && Nchk>0) || override ){
      theDomain->nchk = n0;
      char filename[256];
      if( !override ){
         if(theDomain->rank==0) printf("Creating Checkpoint #%04d...\n",n0);
         sprintf(filename,"checkpoint_%04d",n0);
         output( theDomain , filename );
      }else{
         if(theDomain->rank==0) printf("Creating Final Checkpoint...\n");
         output( theDomain , "output" );
      }
   }
   n0 = (int)( t*Nsnp/t_fin );
   if( LogOut ) n0 = (int)( Nsnp*log(t/t_min)/log(t_fin/t_min) );
   if( (theDomain->nsnp < n0 && Nsnp>0) || override ){
      theDomain->nsnp = n0;
      char filename[256];
      if(!override) sprintf( filename , "snapshot_%04d" , n0 );
      else sprintf( filename , "snapshot" );
      //snapshot( theDomain , filename );
   }

}


