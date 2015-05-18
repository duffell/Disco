
#include "paul.h"

void start_clock( struct domain * theDomain ){
   theDomain->Wallt_init = time(NULL);
}
 
int count_cells( struct domain * theDomain ){

   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;

   int Nc=0;
   int jmin = Ng;
   int jmax = Nr-Ng;
   if( dim_rank[0] == 0 ) jmin = 0; 
   if( dim_rank[0] == dim_size[0]-1 ) jmax = Nr;
   int kmin = Ng;
   int kmax = Nz-Ng;
   if( dim_rank[1] == 0 ) kmin = 0; 
   if( dim_rank[1] == dim_size[1]-1 ) kmax = Nz;

   int j,k;
   for( j=jmin ; j<jmax ; ++j ){
      for( k=kmin ; k<kmax ; ++k ){
         Nc += Np[j+Nr*k];
      }
   }
   MPI_Allreduce( MPI_IN_PLACE , &Nc ,  1 , MPI_INT    , MPI_SUM , theDomain->theComm );
   return(Nc);
}

void generate_log( struct domain * theDomain ){
   time_t endtime = time(NULL);
   int seconds = (int) (endtime - theDomain->Wallt_init);
   
   int Nc = count_cells( theDomain );
   int Nt = theDomain->count_steps;

   double avgdt = (double)seconds/2./(double)Nc/(double)Nt;

   int size = theDomain->size;

   if( theDomain->rank==0 ){
      FILE * logfile = fopen("times.log","w");
      fprintf(logfile,"Run using %d MPI process",size);
      if( theDomain->size > 1 ) fprintf(logfile,"es");
      fprintf(logfile,".\n");
      fprintf(logfile,"Total time = %d sec\n",seconds);
      fprintf(logfile,"Number of cells = %d\n",Nc);
      fprintf(logfile,"Number of timesteps = %d (x%d)\n",Nt,2);
      fprintf(logfile,"Megazones per second = %.2e\n",1./(avgdt*1e6));
      fprintf(logfile,"Megazones per CPU second = %.2e\n",1./(avgdt*1e6*size));
      fprintf(logfile,"Time/zone/step = %.2e microseconds\n",(avgdt*1e6));
      fclose(logfile);
   }
}

/*
void profiler_start( clock_t * prevtime , clock_t * currtime ){
   *prevtime = clock();
   *currtime = clock();
   printf("\n***\nProfiler Running...\n");
}

void profiler_report( char * event , clock_t * prevtime , clock_t * currtime ){
   *currtime = clock();
   printf("%s:\t\t%d ticks\n",event,(int)(*currtime-*prevtime) );
   *prevtime = *currtime;
}

void profiler_end(void){
   printf("***\n\n");
}
*/
